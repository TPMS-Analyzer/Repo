"""Compare the optimized renderer with Matplotlib's original projection."""
import types
import unittest

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from TPMS_Analyzer import analyze, geometry_mesh, plot_geometry, set_preview_transparency


class PreviewRenderingTests(unittest.TestCase):
    def test_render_matches_original_during_navigation_and_toggles(self):
        result = analyze(grid_size=24)
        fig = Figure(figsize=(4, 4), dpi=80)
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_subplot(projection='3d')
        preview = plot_geometry(ax, result, geometry_mesh(result))
        for transparent in (True, False, True):
            preview = set_preview_transparency(preview, transparent)
            for elev, azim, low, high in ((30, 52.5, 0., 2.54),
                                         (15, 120, .3, 2.84),
                                         (65, -40, .5, 2.0)):
                with self.subTest(transparent=transparent, elev=elev, azim=azim):
                    ax.view_init(elev, azim)
                    ax.set(xlim=(low, high), ylim=(low, high), zlim=(low, high))
                    fast = preview.do_3d_projection
                    preview.do_3d_projection = types.MethodType(
                        Poly3DCollection.do_3d_projection, preview)
                    canvas.draw()
                    reference = np.asarray(canvas.buffer_rgba()).copy()
                    preview.do_3d_projection = fast
                    canvas.draw()
                    np.testing.assert_array_equal(canvas.buffer_rgba(), reference)

    def test_explicit_depth_order_uses_compatible_fallback(self):
        result = analyze(grid_size=24)
        fig = Figure(figsize=(3, 3), dpi=60)
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_subplot(projection='3d')
        preview = plot_geometry(ax, result, geometry_mesh(result))
        preview.set_sort_zpos(.5)
        canvas.draw()
        reference = Poly3DCollection.do_3d_projection(preview)
        self.assertEqual(preview.do_3d_projection(), reference)


if __name__ == '__main__':
    unittest.main()
