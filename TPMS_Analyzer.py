"""Python translation of TPMS_Analyzer.m. Run this file to open the GUI.

Arrays follow MATLAB meshgrid order (Y, X, Z). Analysis uses the original
single-cell, nonperiodic voxel approximation, not an analytical pore size.
"""
from __future__ import annotations

import math
import re
import queue
import threading
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

TPMS_TYPES = ('Koch', 'Gyroid', 'Diamond', 'Primitive', 'IWP', 'Neovius', 'FRD', 'PMY')
CONNECTIVITY = ndimage.generate_binary_structure(3, 1)


def validate_inputs(porosity, alpha, grid_size):
    if not np.isfinite(porosity) or not 0 < porosity < 100:
        raise ValueError('Porosity must be between 0 and 100%.')
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError('Unit cell size must be greater than zero.')
    if not np.isfinite(grid_size) or grid_size < 20:
        raise ValueError('Grid size must be at least 20.')


def resolve_display_digits(alpha_text, setting):
    if setting.strip().lower() == 'auto':
        match = re.match(r'^[+-]?(\d+\.?\d*|\.\d+)', alpha_text.strip())
        digits = re.sub(r'[^0-9]', '', match.group(0) if match else '').lstrip('0')
        if not digits:
            raise ValueError('Enter a positive unit-cell size.')
        return min(15, len(digits) + 1)
    try:
        value = float(setting)
        if not math.isfinite(value) or value != math.floor(value) or not 2 <= value <= 15:
            raise ValueError
        return int(value)
    except ValueError:
        raise ValueError('Display sig. digits must be auto or an integer from 2 to 15.') from None


def format_significant(value, digits):
    if not np.isfinite(value):
        return 'NaN' if np.isnan(value) else ('Inf' if value > 0 else '-Inf')
    return '0' if value == 0 else format(float(value), f'#.{digits}g')


def calculate_tpms_function(tpms_type, x, y, z, alpha):
    X, Y, Z = (2 * np.pi * a / alpha for a in (x, y, z))
    cx, cy, cz = np.cos(X), np.cos(Y), np.cos(Z)
    sx, sy, sz = np.sin(X), np.sin(Y), np.sin(Z)
    kind = tpms_type.strip().lower()
    if kind == 'koch':
        return np.cos(2*X)*sy*cz + cx*np.cos(2*Y)*sz + sx*cy*np.cos(2*Z)
    if kind == 'gyroid':
        return sx*cy + sy*cz + sz*cx
    if kind == 'diamond':
        return cx*cy*cz - sx*sy*sz
    if kind == 'primitive':
        return cx + cy + cz
    if kind == 'iwp':
        return 2*(cx*cy + cy*cz + cz*cx) - (np.cos(2*X)+np.cos(2*Y)+np.cos(2*Z))
    if kind == 'neovius':
        return 3*(cx+cy+cz) + 4*cx*cy*cz
    if kind == 'frd':
        a, b, c = np.cos(2*X), np.cos(2*Y), np.cos(2*Z)
        return 4*cx*cy*cz - (a*b + b*c + c*a)
    if kind == 'pmy':
        return 2*cx*cy*cz + np.sin(2*X)*sy + sx*np.sin(2*Z) + np.sin(2*Y)*sz
    raise ValueError(f'Unknown TPMS type: {tpms_type}')


def create_mask(field, c, network_type):
    if network_type.lower() == 'solid':
        return field > c
    if network_type.lower() == 'sheet':
        return np.abs(field) <= c
    raise ValueError('Unknown network type.')


def calculate_porosity(field, c, network_type):
    mask = create_mask(field, c, network_type)
    return 100 * (mask.size - np.count_nonzero(mask)) / mask.size


def find_level_set_constant(field, network_type, target_porosity):
    if network_type.lower() == 'solid':
        margin = 16 * np.spacing(max(1., np.max(np.abs(field))))
        low, high, increasing = field.min()-margin, field.max()+margin, True
    elif network_type.lower() == 'sheet':
        low, high, increasing = 0., np.max(np.abs(field)), False
    else:
        raise ValueError('Unknown network type.')
    p_low = calculate_porosity(field, low, network_type)
    p_high = calculate_porosity(field, high, network_type)
    if not min(p_low, p_high) <= target_porosity <= max(p_low, p_high):
        raise ValueError('Target porosity cannot be generated on this grid.')
    c, best_error = low, abs(p_low-target_porosity)
    if abs(p_high-target_porosity) < best_error:
        c, best_error = high, abs(p_high-target_porosity)
    for _ in range(60):
        mid = (low+high)/2
        if mid == low or mid == high:
            break
        p_mid = calculate_porosity(field, mid, network_type)
        error = abs(p_mid-target_porosity)
        if error < best_error:
            c, best_error = mid, error
        if best_error <= .01:
            return float(c)
        if (increasing and p_mid < target_porosity) or (not increasing and p_mid > target_porosity):
            low = mid
        else:
            high = mid
    return float(c)


def calculate_wetted_surface_area(mask, dx):
    return sum(np.count_nonzero(np.diff(mask, axis=a)) for a in range(3)) * dx**2


def spanning_labels(labels, direction_index):
    if direction_index not in (0, 1, 2):
        raise ValueError('Invalid traversal direction.')
    axis = (1, 0, 2)[direction_index]
    common = np.intersect1d(np.unique(np.take(labels, 0, axis)),
                            np.unique(np.take(labels, -1, axis)))
    return common[common != 0]


def check_tpms_connection(region, direction_index):
    """direction_index is zero-based: X=0, Y=1, Z=2."""
    labels, _ = ndimage.label(region, CONNECTIVITY)
    return bool(spanning_labels(labels, direction_index).size)


@dataclass
class PoreResult:
    diameter: float
    radius: float
    diameter_xyz: np.ndarray
    radius_xyz: np.ndarray
    limiting_direction: str
    center_index: np.ndarray  # zero-based Y, X, Z


def calculate_traversable_pore_diameter(mask, dx):
    mask = np.asarray(mask, dtype=bool)
    void = ~mask
    radii = np.zeros(3)
    empty = PoreResult(0., 0., radii.copy(), radii, 'None', np.full(3, np.nan))
    if not void.any():
        return empty
    if not mask.any():
        raise ValueError('No solid voxels exist; pore diameter is undefined.')
    # MATLAB bwdist returns single precision. Preserve that rounding here.
    clearance = ndimage.distance_transform_edt(void).astype(np.float32)
    clearance = np.maximum(clearance-np.float32(.5), 0) * np.float32(dx)
    clearance[mask] = 0
    maximum_radius = float(clearance.max())
    if maximum_radius <= 0:
        return empty
    for direction in range(3):
        if not check_tpms_connection(void, direction):
            continue
        low, high = 0., maximum_radius
        for _ in range(50):
            trial = (low+high)/2
            if check_tpms_connection(void & (clearance >= trial), direction):
                low = trial
            else:
                high = trial
            if high-low < dx/10:
                break
        radii[direction] = low
    limiting = int(np.argmin(radii))
    radius = float(radii[limiting])
    center = np.full(3, np.nan)
    if radius > 0:
        labels, _ = ndimage.label(void & (clearance >= radius), CONNECTIVITY)
        components = spanning_labels(labels, limiting)
        # MATLAB enumerates components and voxels in column-major order.
        flat = labels.ravel(order='F')
        members = [(np.flatnonzero(flat == label), label) for label in components]
        members.sort(key=lambda item: item[0][0])
        best = (False, -np.inf, -np.inf)
        shape = np.array(mask.shape)
        for indices, _ in members:
            coords = np.column_stack(np.unravel_index(indices, mask.shape, order='F'))
            boundary = dx*np.minimum(coords, shape-1-coords).min(axis=1)
            central = np.sum((coords-(shape-1)/2)**2, axis=1)
            inside = boundary >= radius
            has_inside = bool(inside.any())
            score = np.inf if has_inside else boundary.max()
            candidates = np.flatnonzero(inside if has_inside else boundary == score)
            chosen = candidates[np.argmin(central[candidates])]
            rank = (has_inside, score, -central[chosen])
            if rank > best:
                best, center = rank, coords[chosen]
    return PoreResult(2*radius, radius, 2*radii, radii, 'XYZ'[limiting], center)


@dataclass
class AnalysisResult:
    network_type: str
    tpms_type: str
    target_porosity: float
    alpha: float
    grid_size: int
    digits: int
    field: np.ndarray
    c: float
    actual_porosity: float
    solid_fraction: float
    wetted_area: float
    pore: PoreResult

    @property
    def rows(self):
        p = self.pore
        mean = np.mean(p.diameter_xyz)
        error = np.ptp(p.diameter_xyz)/mean*100 if mean > 0 else np.nan
        values = [
            ('Network type', self.network_type, ''), ('TPMS type', self.tpms_type, ''),
            ('Target porosity', self.target_porosity, '%'),
            ('Calculated porosity', self.actual_porosity, '%'),
            ('Solid volume fraction', self.solid_fraction, '%'),
            ('Level-set constant, c', self.c, '-'),
            ('Wetted area (voxel)', self.wetted_area, 'mm^2'),
            ('Unit cell volume', self.alpha**3, 'mm^3'),
            ('Surface area / volume', self.wetted_area/self.alpha**3, '1/mm'),
            *[(f'Pore diameter {d}', v, 'mm') for d, v in zip('XYZ', p.diameter_xyz)],
            ('Directional error', error, '%'), ('Limiting direction', p.limiting_direction, ''),
            ('Pore radius', p.radius, 'mm'), ('Pore diameter, Dp', p.diameter, 'mm'),
            ('Dp / cell size', p.diameter/self.alpha, '-')]
        return [(name, value if isinstance(value, str) else format_significant(value, self.digits), unit)
                for name, value, unit in values]


def analyze(network_type='Solid', tpms_type='Koch', target_porosity=70.,
            alpha=2.54, grid_size=150, digits=4):
    if not np.isfinite(grid_size):
        raise ValueError('Grid size must be at least 20.')
    grid_size = math.floor(grid_size+.5)  # MATLAB rounds positive halves up.
    validate_inputs(target_porosity, alpha, grid_size)
    coordinate = np.linspace(0., alpha, grid_size)
    x, y, z = np.meshgrid(coordinate, coordinate, coordinate, indexing='xy', sparse=True)
    field = calculate_tpms_function(tpms_type, x, y, z, alpha)
    c = find_level_set_constant(field, network_type, target_porosity)
    mask = create_mask(field, c, network_type)
    solid_fraction = 100*np.count_nonzero(mask)/mask.size
    porosity = 100*np.count_nonzero(~mask)/mask.size
    dx = alpha/(grid_size-1)
    area = calculate_wetted_surface_area(mask, dx)
    pore = calculate_traversable_pore_diameter(mask, dx)
    return AnalysisResult(network_type, tpms_type, target_porosity, alpha, grid_size,
                          digits, field, c, porosity, solid_fraction, area, pore)


def geometry_mesh(result, max_display_points=35):
    """Build a capped display mesh without changing the analysis grid.

    Matplotlib sorts every triangle on each 3D interaction. Sampling the
    field before marching cubes keeps rotation and zoom responsive. The
    sampled coordinates include both box boundaries, and are mapped back to
    physical coordinates after interpolation.
    """
    from skimage.measure import marching_cubes
    field = result.field
    if max_display_points < 2:
        raise ValueError('Display resolution must be at least 2.')
    step = max(1, math.ceil((result.grid_size-1)/(max_display_points-1)))
    indices = np.r_[np.arange(0, result.grid_size-1, step), result.grid_size-1]
    field = field[np.ix_(indices, indices, indices)]
    volume = field-result.c if result.network_type.lower() == 'solid' else result.c-np.abs(field)
    outside = -max(1., float(np.max(np.abs(volume))))
    padded = np.pad(volume.astype(np.float32), 1, constant_values=outside)
    if padded.max() <= 0:
        return np.empty((0, 3)), np.empty((0, 3), dtype=int)
    vertices, faces, _, _ = marching_cubes(padded, 0, allow_degenerate=False)
    vertices = np.clip(vertices-1, 0, len(indices)-1)
    for axis in range(3):
        vertices[:, axis] = np.interp(vertices[:, axis], np.arange(len(indices)), indices)
    vertices *= result.alpha/(result.grid_size-1)
    vertices = vertices[:, [1, 0, 2]]  # YXZ to physical XYZ
    triangles = vertices[faces]
    valid = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                    triangles[:, 2]-triangles[:, 0]), axis=1) > 0
    return vertices, faces[valid]


def plot_geometry(ax, result, mesh, transparent=True):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    ax.clear()
    vertices, faces = mesh
    surfaces = [vertices[faces]]
    lattice_count = len(faces)
    p = result.pore
    if p.radius > 0 and np.all(np.isfinite(p.center_index)):
        center = p.center_index[[1, 0, 2]] * result.alpha/(result.grid_size-1)
        u = np.linspace(0, 2*np.pi, 49)
        v = np.linspace(-np.pi/2, np.pi/2, 49)
        sx = np.outer(np.cos(u), np.cos(v))
        sy = np.outer(np.sin(u), np.cos(v))
        sz = np.outer(np.ones_like(u), np.sin(v))
        sphere = center + p.radius*np.stack((sx, sy, sz), axis=-1)
        a, b, c, d = sphere[:-1, :-1], sphere[1:, :-1], sphere[1:, 1:], sphere[:-1, 1:]
        surfaces.append(np.concatenate((np.stack((a, b, c), axis=2).reshape(-1, 3, 3),
                                        np.stack((a, c, d), axis=2).reshape(-1, 3, 3))))
        ax.scatter(*center, color='black', s=15)
        ax.text(center[0], center[1], center[2]+p.radius,
                 f'  D_p = {format_significant(p.diameter, result.digits)} mm\n  Representative position',
                 color=(.75, .05, .05), weight='bold', fontsize=11)
    # A single collection sorts lattice and sphere triangles together on every
    # redraw; separate collections can incorrectly swap their depth ordering.
    preview = None
    if any(len(surface) for surface in surfaces):
        triangles = np.concatenate(surfaces)
        colors = np.empty((len(triangles), 4))
        colors[:lattice_count] = (.55, .47, .28, .22 if transparent else 1.)
        colors[lattice_count:] = (.90, .15, .10, 1.)
        preview = Poly3DCollection(triangles, facecolors=colors,
                                   linewidths=0, shade=True)
        preview.set_edgecolor('none')
        preview._tpms_lattice_count = lattice_count
        preview._tpms_face_count = len(triangles)
        ax.add_collection3d(preview)
    alpha = result.alpha
    ax.set(xlim=(0, alpha), ylim=(0, alpha), zlim=(0, alpha))
    for name in 'xyz':
        getattr(ax, f'set_{name}label')(f'{name.upper()} (mm)', fontweight='bold', fontfamily='Times New Roman')
    ax.set_title(f'{result.network_type} - {result.tpms_type}, Porosity = '
                 f'{format_significant(result.actual_porosity, result.digits)}%', weight='bold')
    ax.set_box_aspect((1, 1, 1))
    # MATLAB and Matplotlib measure azimuth from different starting axes.
    ax.view_init(elev=30, azim=52.5)
    ax.grid(False)
    return preview


def set_preview_transparency(preview, transparent):
    if preview is None:
        return
    colors = np.empty((preview._tpms_face_count, 4))
    count = preview._tpms_lattice_count
    colors[:count] = (.55, .47, .28, .22 if transparent else 1.)
    colors[count:] = (.90, .15, .10, 1.)
    preview.set_facecolor(colors)


class TPMSAnalyzer:
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        self.root = root
        self.messages = queue.Queue()
        self.transparent = True
        self.lattice = None
        self.result = None
        self.busy = False
        root.title('TPMS Geometry Analyzer')
        w = min(1440, max(640, root.winfo_screenwidth()-80))
        h = min(960, max(520, root.winfo_screenheight()-110))
        root.geometry(f'{w}x{h}+40+45')
        root.minsize(640, 520)
        root.configure(bg='white')
        bg, fg, accent = '#f5f7fa', '#1f2938', '#215287'
        style = ttk.Style(root)
        style.theme_use('clam')
        style.configure('TPMS.Treeview', font=('Arial', 11), rowheight=27, foreground=fg)
        style.configure('TPMS.Treeview.Heading', font=('Arial', 11, 'bold'))
        self.input_panel = tk.LabelFrame(root, text='TPMS Input', font=('Arial', 14, 'bold'), bg=bg, fg=fg)
        self.variables, self.labels, self.controls = [], [], []
        captions = ['Network type', 'TPMS type', 'Target porosity (%)', 'Unit cell size (mm)',
                    'Grid points per axis', 'Display sig. digits']
        for i, (caption, value) in enumerate(zip(captions, ['Solid', 'Koch', '70', '2.54', '150', 'auto'])):
            var = tk.StringVar(value=value)
            self.variables.append(var)
            self.labels.append(tk.Label(self.input_panel, text=caption, font=('Arial', 12), bg=bg, fg=fg, anchor='w'))
            control = (ttk.Combobox(self.input_panel, textvariable=var, state='readonly',
                                   values=('Solid', 'Sheet') if i == 0 else TPMS_TYPES, font=('Arial', 12))
                       if i < 2 else tk.Entry(self.input_panel, textvariable=var, font=('Arial', 12), fg=fg))
            self.controls.append(control)
        self.add_tooltip(self.controls[4], 'Number of grid points along each axis (minimum 20).')
        self.add_tooltip(self.controls[5], 'auto = unit-cell significant digits + 1 (maximum 15);\n'
                         'or enter 2 to 15. Press CALCULATE to apply.')
        self.calculate_button = tk.Button(root, text='CALCULATE', command=self.calculate,
                                           bg=accent, fg='white', font=('Arial', 13, 'bold'))
        self.transparency_button = tk.Button(root, text='Transparency: ON', command=self.toggle_transparency,
                                              font=('Arial', 11))
        self.add_tooltip(self.transparency_button,
                         'Switch the lattice between transparent and opaque without recalculating.')
        self.result_panel = tk.LabelFrame(root, text='Results', font=('Arial', 14, 'bold'), bg=bg, fg=fg)
        self.table = ttk.Treeview(self.result_panel, columns=('Parameter', 'Value', 'Unit'),
                                  show='headings', style='TPMS.Treeview', selectmode='browse')
        for col in self.table['columns']:
            self.table.heading(col, text=col)
            self.table.column(col, anchor='w', stretch=False)
        self.table.tag_configure('alternate', background='#f0f5fa')
        self.scrollbar = ttk.Scrollbar(self.result_panel, orient='vertical', command=self.table.yview)
        self.table.configure(yscrollcommand=self.scrollbar.set)
        self.status = tk.Label(self.result_panel, text='Press CALCULATE to begin.', anchor='w', justify='left',
                               font=('Arial', 11), bg=bg, fg=accent)
        self.plot_frame = tk.Frame(root, bg='white')
        self.figure = Figure(figsize=(7, 7), dpi=100, facecolor='white')
        self.axes = self.figure.add_subplot(111, projection='3d')
        self.axes.set(xlabel='X (mm)', ylabel='Y (mm)', zlabel='Z (mm)')
        self.axes.set_box_aspect((1, 1, 1))
        self.axes.view_init(30, -37.5)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame, pack_toolbar=False)
        self.toolbar.pack(side='top', fill='x')
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        root.bind('<Configure>', self.layout)
        self._poll_id = root.after(100, self.poll)
        root.protocol('WM_DELETE_WINDOW', self.close)

    def add_tooltip(self, widget, text):
        import tkinter as tk
        tip = [None]

        def hide(event=None):
            if tip[0] is not None:
                tip[0].destroy()
                tip[0] = None

        def show(event=None):
            hide()
            tip[0] = tk.Toplevel(self.root)
            tip[0].wm_overrideredirect(True)
            tip[0].geometry(f'+{widget.winfo_rootx()}+{widget.winfo_rooty()+widget.winfo_height()+4}')
            tk.Label(tip[0], text=text, justify='left', background='#fffbdc',
                     relief='solid', borderwidth=1, padx=6, pady=4).pack()

        widget.bind('<Enter>', show, add='+')
        widget.bind('<Leave>', hide, add='+')
        widget.bind('<ButtonPress>', hide, add='+')

    def layout(self, event=None):
        if event is not None and event.widget is not self.root:
            return
        w, h = self.root.winfo_width(), self.root.winfo_height()
        left = min(max(430, min(580, math.floor(.36*w+.5))), max(200, w-36))
        self.input_panel.place(x=18, y=18, width=left, height=325)
        button_y = 355
        bw = math.floor(.52*(left-12))
        self.calculate_button.place(x=18, y=button_y, width=bw, height=42)
        self.transparency_button.place(x=18+bw+12, y=button_y, width=left-bw-12, height=42)
        rh = max(60, h-427)
        self.result_panel.place(x=18, y=409, width=left, height=rh)
        lw = round(.51*(left-36))
        for i, (label, control) in enumerate(zip(self.labels, self.controls)):
            label.place(x=16, y=22+i*45, width=lw, height=25)
            control.place(x=18+lw+8, y=19+i*45, width=max(50, left-lw-44), height=32)
        tw = left-28
        self.table.place(x=12, y=12, width=tw-16, height=max(25, rh-96))
        self.scrollbar.place(x=12+tw-16, y=12, width=16, height=max(25, rh-96))
        for col, width in zip(self.table['columns'], [max(90, tw-181), 100, 65]):
            self.table.column(col, width=width)
        self.status.place(x=12, y=max(38, rh-80), width=tw, height=48)
        if w-(18+left+76)-48 < 100:
            self.plot_frame.place_forget()
        else:
            self.plot_frame.place(x=left+42, y=18, width=w-left-60, height=h-36)

    def calculate(self):
        if self.busy:
            return
        self.table.delete(*self.table.get_children())
        try:
            network, tpms, porosity, alpha, grid, digits = [v.get() for v in self.variables]
            porosity_value, alpha_value, grid_value = float(porosity), float(alpha), float(grid)
            if not np.isfinite(grid_value):
                raise ValueError('Grid size must be at least 20.')
            validate_inputs(porosity_value, alpha_value, math.floor(grid_value+.5))
            digits_value = resolve_display_digits(alpha, digits)
        except ValueError as exc:
            self.status.configure(text=f'ERROR: {exc}')
            return
        self.busy = True
        self.calculate_button.configure(state='disabled')
        self.status.configure(text='Calculating...')

        def work():
            try:
                result = analyze(network, tpms, porosity_value, alpha_value, grid_value, digits_value)
                mesh = geometry_mesh(result)
                self.messages.put((result, mesh, None))
            except Exception as exc:
                self.messages.put((None, None, str(exc)))
        threading.Thread(target=work, daemon=True).start()

    def poll(self):
        try:
            result, mesh, error = self.messages.get_nowait()
        except queue.Empty:
            pass
        else:
            try:
                if error:
                    raise RuntimeError(error)
                self.result = result
                for i, row in enumerate(result.rows):
                    self.table.insert('', 'end', values=row, tags=('alternate',) if i % 2 else ())
                self.lattice = plot_geometry(self.axes, result, mesh, self.transparent)
                self.toolbar.update()
                self.canvas.draw_idle()
                fmt = lambda value: format_significant(value, result.digits)
                self.status.configure(text=f'Dp = {fmt(result.pore.diameter)} mm  |  Porosity = {fmt(result.actual_porosity)} %\n'
                                           f'{result.network_type} / {result.tpms_type} | {result.digits} significant digits')
                print(f'\nTPMS ANALYSIS RESULT ({result.digits} significant digits)', flush=True)
                for name, value, unit in result.rows:
                    print(f'{name:26s} : {value} {unit}')
                print('Radius X / Y / Z           : ' + ' / '.join(fmt(v) for v in result.pore.radius_xyz) + ' mm')
                if np.all(np.isfinite(result.pore.center_index)):
                    center = result.pore.center_index[[1, 0, 2]]*result.alpha/(result.grid_size-1)
                    print('Representative pore center [X Y Z]: ' + ' '.join(fmt(v) for v in center) + ' mm')
            except Exception as exc:
                self.status.configure(text=f'ERROR: {exc}')
                print(f'TPMS ANALYZER ERROR: {exc}')
            finally:
                self.busy = False
                self.calculate_button.configure(state='normal')
        self._poll_id = self.root.after(100, self.poll)

    def toggle_transparency(self):
        self.transparent = not self.transparent
        self.transparency_button.configure(text=f'Transparency: {"ON" if self.transparent else "OFF"}')
        set_preview_transparency(self.lattice, self.transparent)
        self.canvas.draw_idle()

    def close(self):
        self.root.after_cancel(self._poll_id)
        self.root.destroy()


def main():
    import tkinter as tk
    root = tk.Tk()
    TPMSAnalyzer(root)
    root.mainloop()


if __name__ == '__main__':
    main()
