"""Web-only display mesh preparation; numerical and desktop code stay separate."""
import numpy as np

from TPMS_Analyzer import geometry_mesh


def detailed_mesh(result):
    # The desktop preview uses at most 35 samples per axis. WebGL can display
    # the full default 150-point grid without that coarse surface sampling.
    vertices, faces = geometry_mesh(result, max_display_points=150)
    vertices = np.clip(vertices.astype(np.float64), 0., result.alpha)
    tolerance = result.alpha*1e-6
    vertices[np.abs(vertices) < tolerance] = 0.
    vertices[np.abs(vertices-result.alpha) < tolerance] = result.alpha
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1]-triangles[:, 0],
                       triangles[:, 2]-triangles[:, 0])
    # Padding/clamping can collapse boundary triangles. Do not send their
    # nearly zero normals into the lighting shader.
    keep = np.linalg.norm(normals, axis=1) > (result.alpha/(result.grid_size-1))**2*1e-8
    return vertices, faces[keep]


def split_surface_parts(vertices, faces, alpha):
    """Separate flat cut faces from curved surfaces for correct sharp edges."""
    triangles = vertices[faces]
    cap_masks = []
    for axis in range(3):
        for boundary in (0., alpha):
            cap_masks.append(np.all(triangles[:, :, axis] == boundary, axis=1))
    is_cap = np.logical_or.reduce(cap_masks)
    groups = [('TPMS surface', faces[~is_cap], False)]
    for index, mask in enumerate(cap_masks):
        group = faces[mask].copy()
        if len(group):
            axis = index//2
            normals = np.cross(vertices[group[:, 1]]-vertices[group[:, 0]],
                               vertices[group[:, 2]]-vertices[group[:, 0]])
            desired = -1 if index % 2 == 0 else 1
            reverse = normals[:, axis]*desired < 0
            group[reverse] = group[reverse, ::-1]
        groups.append((f"{'XYZ'[index//2]} {'min' if index % 2 == 0 else 'max'} cut face", group, True))
    parts = []
    for name, group, flat in groups:
        if not len(group):
            continue
        used, remapped = np.unique(group, return_inverse=True)
        # Independent vertex buffers prevent smooth normals from bleeding
        # across the sharp curve-to-cap and cube-edge boundaries.
        parts.append((name, vertices[used], remapped.reshape(-1, 3), flat))
    return parts
