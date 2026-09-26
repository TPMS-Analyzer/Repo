"""Numerical behavior checks independent of MATLAB availability."""
import unittest
import numpy as np
from TPMS_Analyzer import (TPMS_TYPES, analyze, calculate_tpms_function,
    calculate_traversable_pore_diameter, calculate_wetted_surface_area,
    check_tpms_connection, create_mask, find_level_set_constant,
    calculate_porosity, resolve_display_digits, format_significant, geometry_mesh)


class NumericalTests(unittest.TestCase):
    def test_field_landmarks(self):
        # Analytically known values at (0,0,0).
        expected = [0, 0, 1, 3, 3, 13, 1, 2]
        for name, value in zip(TPMS_TYPES, expected):
            self.assertAlmostEqual(float(calculate_tpms_function(name, 0., 0., 0., 1.)), value)

    def test_mask_threshold_inequalities(self):
        f = np.array([-1., 0., 1.])
        np.testing.assert_array_equal(create_mask(f, 0, 'Solid'), [False, False, True])
        np.testing.assert_array_equal(create_mask(f, 1, 'Sheet'), [True, True, True])

    def test_six_connectivity_and_axis_order(self):
        for direction, axis in enumerate((1, 0, 2)):
            region = np.zeros((7, 7, 7), bool)
            sl = [3, 3, 3]
            sl[axis] = slice(None)
            region[tuple(sl)] = True
            for other in range(3):
                self.assertEqual(check_tpms_connection(region, other), direction == other)
        diagonal = np.zeros((7, 7, 7), bool)
        diagonal[np.arange(7), np.arange(7), np.arange(7)] = True
        self.assertFalse(check_tpms_connection(diagonal, 0))

    def test_internal_area_excludes_box_faces(self):
        solid = np.zeros((5, 5, 5), bool)
        solid[2, 2, 2] = True
        self.assertEqual(calculate_wetted_surface_area(solid, .2), 6*.2**2)
        self.assertEqual(calculate_wetted_surface_area(np.ones_like(solid), .2), 0.)

    def test_straight_channel_known_clearance(self):
        # A one-voxel-wide X channel has center clearance 1-0.5 voxels.
        solid = np.ones((7, 7, 7), bool)
        solid[3, :, 3] = False
        pore = calculate_traversable_pore_diameter(solid, 1.)
        self.assertLessEqual(pore.radius_xyz[0], .5)
        self.assertLess(.5-pore.radius_xyz[0], .1)
        np.testing.assert_array_equal(pore.radius_xyz[1:], [0., 0.])
        self.assertEqual(pore.diameter, 0.)
        self.assertEqual(pore.limiting_direction, 'Y')

    def test_no_void_and_no_solid(self):
        self.assertEqual(calculate_traversable_pore_diameter(np.ones((3, 3, 3), bool), 1.).diameter, 0.)
        with self.assertRaisesRegex(ValueError, 'No solid'):
            calculate_traversable_pore_diameter(np.zeros((3, 3, 3), bool), 1.)

    def test_display_precision(self):
        for text, digits in [('2.54', 4), ('2.540', 5), ('0.002540', 5), ('2.540e-3', 5), ('100', 4)]:
            self.assertEqual(resolve_display_digits(text, 'auto'), digits)
        self.assertEqual(format_significant(2.5, 4), '2.500')
        self.assertEqual(format_significant(0., 4), '0')
        for setting in ['1', '16', '2.5', 'nan']:
            with self.assertRaises(ValueError):
                resolve_display_digits('2.54', setting)

    def test_threshold_and_all_geometry_combinations(self):
        for network in ['Solid', 'Sheet']:
            for tpms in TPMS_TYPES:
                with self.subTest(network=network, tpms=tpms):
                    r = analyze(network, tpms, grid_size=24)
                    self.assertAlmostEqual(r.actual_porosity+r.solid_fraction, 100.)
                    # Small symmetric grids can have tied levels exceeding .01 pp.
                    self.assertLess(abs(r.actual_porosity-70), 1.)
                    self.assertEqual(len(r.rows), 17)
                    self.assertGreaterEqual(r.pore.diameter, 0.)
                    self.assertEqual(r.pore.diameter, min(r.pore.diameter_xyz))
                    vertices, faces = geometry_mesh(r)
                    self.assertGreater(len(faces), 0)
                    self.assertTrue(np.all(vertices >= 0))
                    self.assertTrue(np.all(vertices <= r.alpha+1e-6))
                    if r.pore.radius > 0:
                        self.assertTrue(np.all(np.isfinite(r.pore.center_index)))
                        self.assertFalse(create_mask(r.field, r.c, network)[tuple(r.pore.center_index.astype(int))])

    def test_physical_scaling(self):
        a = analyze(tpms_type='Gyroid', alpha=1., grid_size=24)
        b = analyze(tpms_type='Gyroid', alpha=2., grid_size=24)
        self.assertEqual(a.actual_porosity, b.actual_porosity)
        self.assertEqual(b.wetted_area, 4*a.wetted_area)
        np.testing.assert_allclose(b.pore.diameter_xyz, 2*a.pore.diameter_xyz)

    def test_preview_mesh_is_lighter_without_changing_analysis(self):
        r = analyze(tpms_type='Gyroid', grid_size=60)
        original_field = r.field.copy()
        full_vertices, full_faces = geometry_mesh(r, max_display_points=1000)
        vertices, faces = geometry_mesh(r)
        self.assertLess(len(faces), len(full_faces)/2)
        self.assertGreater(len(vertices), 0)
        self.assertTrue(np.all(vertices >= 0))
        self.assertTrue(np.all(vertices <= r.alpha+1e-6))
        np.testing.assert_array_equal(r.field, original_field)

    def test_grid_rounding_and_validation(self):
        self.assertEqual(analyze(grid_size=20.5).grid_size, 21)
        for kwargs in [dict(target_porosity=0), dict(target_porosity=100),
                       dict(alpha=0), dict(grid_size=19), dict(grid_size=np.nan)]:
            with self.assertRaises(ValueError):
                analyze(**kwargs)

    def test_plot_render_and_transparency(self):
        import io
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from TPMS_Analyzer import plot_geometry, set_preview_transparency
        r = analyze(grid_size=24)
        fig = Figure(figsize=(5, 5))
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111, projection='3d')
        lattice = plot_geometry(ax, r, geometry_mesh(r))
        self.assertEqual(len(ax.collections), 2)  # combined surfaces and center dot
        self.assertGreater(lattice._tpms_face_count, lattice._tpms_lattice_count)
        self.assertAlmostEqual(lattice._facecolor3d[0, 3], .22)
        self.assertAlmostEqual(lattice._facecolor3d[-1, 3], 1.)
        shaded_rgb = lattice._facecolor3d[:, :3].copy()
        self.assertGreater(np.ptp(shaded_rgb[:lattice._tpms_lattice_count, 0]), .1)
        set_preview_transparency(lattice, False)
        self.assertTrue(np.all(lattice._facecolor3d[:, 3] == 1.))
        set_preview_transparency(lattice, True)
        self.assertAlmostEqual(lattice._facecolor3d[0, 3], .22)
        self.assertAlmostEqual(lattice._facecolor3d[-1, 3], 1.)
        np.testing.assert_array_equal(lattice._facecolor3d[:, :3], shaded_rgb)
        self.assertGreater(lattice._edgecolor3d[0, 3], 0.)
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        self.assertGreater(buf.tell(), 1000)


if __name__ == '__main__':
    unittest.main()
