"""Faster projection of the fixed triangular TPMS preview.

Keep Matplotlib's shading and painter ordering, but batch the projection and
sorting work. Unsupported collection features use Matplotlib's own renderer.
"""
import numpy as np
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class TrianglePreviewCollection(Poly3DCollection):
    def do_3d_projection(self):
        # These internal buffers differ across Matplotlib releases. Limit the
        # optimized path to the fixed triangle/color layout we understand.
        triangles = getattr(self, '_tpms_triangles', None)
        count = len(triangles) if triangles is not None else 0
        vec = getattr(self, '_vec', None)
        faces = getattr(self, '_facecolor3d', None)
        edges = getattr(self, '_edgecolor3d', None)
        if (not count or vec is None or vec.shape != (4, count*3)
                or faces is None or faces.shape != (count, 4)
                or edges is None or edges.shape != (count, 4)
                or getattr(self, '_A', None) is not None
                or getattr(self, '_codes3d', None) is not None
                or getattr(self, '_axlim_clip', False)
                or getattr(self, '_sort_zpos', None) is not None
                or not getattr(self, '_closed', False)):
            return super().do_3d_projection()

        projected = self.axes.M @ vec
        with np.errstate(divide='ignore', invalid='ignore'):
            projected = projected[:3] / projected[3]
        if not np.all(np.isfinite(projected)):
            return super().do_3d_projection()
        xyz = projected.reshape(3, count, 3)
        depths = self._zsortfunc(xyz[2], axis=1)
        # Python's original descending sort is stable for equal depths.
        order = np.argsort(-depths, kind='stable')
        xy = xyz[:2].transpose(1, 2, 0)[order]
        self._facecolors2d = faces[order]
        self._edgecolors2d = edges[order]
        # An ndarray activates PolyCollection's fast path for closed polygons.
        PolyCollection.set_verts(self, xy, self._closed)
        return float(np.min(xyz[2]))
