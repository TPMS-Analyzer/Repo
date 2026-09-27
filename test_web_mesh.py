import unittest
import numpy as np

from TPMS_Analyzer import analyze, geometry_mesh
from web_mesh import detailed_mesh, split_surface_parts


class WebMeshTests(unittest.TestCase):
    def test_detail_and_planar_caps(self):
        result = analyze(grid_size=60)
        original = result.field.copy()
        vertices, faces = detailed_mesh(result)
        _, coarse = geometry_mesh(result)
        self.assertGreater(len(faces), 2*len(coarse))
        parts = split_surface_parts(vertices, faces, result.alpha)
        self.assertEqual(sum(len(f) for _, _, f, _ in parts), len(faces))
        self.assertEqual(sum(flat for _, _, _, flat in parts), 6)
        for name, v, f, flat in parts:
            self.assertTrue(np.all((v >= 0) & (v <= result.alpha)))
            normals = np.cross(v[f[:, 1]]-v[f[:, 0]], v[f[:, 2]]-v[f[:, 0]])
            self.assertTrue(np.all(np.linalg.norm(normals, axis=1) > 0))
            if flat:
                axis = 'XYZ'.index(name[0])
                boundary = 0 if 'min' in name else result.alpha
                self.assertTrue(np.all(v[:, axis] == boundary))
                self.assertTrue(np.all(normals[:, axis]*(-1 if boundary == 0 else 1) > 0))
        np.testing.assert_array_equal(result.field, original)

    def test_sheet_caps_preserve_face_count(self):
        result = analyze(network_type='Sheet', tpms_type='Gyroid', grid_size=40)
        v, f = detailed_mesh(result)
        parts = split_surface_parts(v, f, result.alpha)
        self.assertEqual(sum(len(group) for _, _, group, _ in parts), len(f))


if __name__ == '__main__':
    unittest.main()
