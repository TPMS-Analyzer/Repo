"""Analytical geometry checks for uncapped isosurface area."""
import unittest
from unittest.mock import patch

import numpy as np

from TPMS_Analyzer import analyze, calculate_wetted_surface_area as area


class WettedAreaTests(unittest.TestCase):
    def setUp(self):
        self.n = 25
        self.dx = 1/(self.n-1)
        q = np.linspace(0., 1., self.n)
        self.y, self.x, self.z = np.meshgrid(q, q, q, indexing='ij')

    def test_plane_area_excludes_unit_cell_caps(self):
        self.assertAlmostEqual(area(self.x, .43, ' Solid ', self.dx), 1., places=6)

    def test_oblique_plane_avoids_staircase_area(self):
        field = self.x+self.y
        measured = area(field, 1., 'solid', self.dx)
        self.assertAlmostEqual(measured, np.sqrt(2), places=6)
        mask = field > 1.
        old = sum(np.count_nonzero(np.diff(mask, axis=a)) for a in range(3))*self.dx**2
        self.assertGreater(old, measured*1.3)

    def test_sheet_two_separate_interfaces(self):
        field = self.x-.5
        self.assertAlmostEqual(area(field, .2, 'Sheet', self.dx), 2., places=6)
        self.assertEqual(area(field, 0., 'sheet', self.dx), 0.)
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            area(field, -.2, 'sheet', self.dx)

    def test_levels_at_or_outside_extrema(self):
        for c in (-1., 0., 1., 2.):
            self.assertEqual(area(self.x, c, 'solid', self.dx), 0.)
        self.assertEqual(area(np.ones_like(self.x), 1., 'solid', self.dx), 0.)
        # -c is outside the field but +c intersects it; count only one surface.
        self.assertAlmostEqual(area(self.x, .3, 'sheet', self.dx), 1., places=6)

    def test_physical_area_scaling(self):
        base = area(self.x+self.y, 1., 'solid', self.dx)
        self.assertAlmostEqual(area(self.x+self.y, 1., 'solid', 3*self.dx), 9*base)

    def test_sphere_convergence(self):
        errors = []
        for n in (20, 48):
            q = np.linspace(-1., 1., n)
            y, x, z = np.meshgrid(q, q, q, indexing='ij')
            measured = area(x*x+y*y+z*z, .6**2, 'solid', 2/(n-1))
            errors.append(abs(measured-4*np.pi*.6**2))
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[1]/(4*np.pi*.6**2), .01)

    def test_invalid_inputs_and_nonfinite_triangle_area(self):
        with self.assertRaisesRegex(ValueError, 'Unknown network'):
            area(self.x, .5, 'invalid', self.dx)
        with self.assertRaises(ValueError):
            area(self.x, np.nan, 'solid', self.dx)
        with self.assertRaises(ValueError):
            area(self.x, .5, 'solid', 0.)
        broken = self.x.copy()
        broken[0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            area(broken, .5, 'solid', self.dx)
        mesh = (np.array([[0., 0., 0.], [np.nan, 0., 0.], [0., 1., 0.]]),
                np.array([[0, 1, 2]]), None, None)
        with patch('skimage.measure.marching_cubes', return_value=mesh):
            with self.assertRaisesRegex(ValueError, 'Nonfinite triangle area'):
                area(self.x, .5, 'solid', self.dx)

    def test_degenerate_triangles_contribute_zero(self):
        mesh = (np.array([[0., 0., 0.], [0., 0., 0.], [0., 1., 0.]]),
                np.array([[0, 1, 2]]), None, None)
        with patch('skimage.measure.marching_cubes', return_value=mesh):
            self.assertEqual(area(self.x, .5, 'solid', self.dx), 0.)

    def test_analysis_and_output_use_new_area(self):
        result = analyze(tpms_type='Gyroid', grid_size=24)
        expected = area(result.field, result.c, result.network_type, result.alpha/23)
        self.assertEqual(result.wetted_area, expected)
        rows = {name: value for name, value, _ in result.rows}
        self.assertIn('Wetted area (isosurface)', rows)
        self.assertAlmostEqual(float(rows['Surface area / volume']),
                               expected/result.alpha**3, delta=.005)


if __name__ == '__main__':
    unittest.main()
