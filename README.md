# TPMS Geometry Analyzer — desktop and Streamlit

## Run in a web browser

This branch adds a Streamlit interface using the existing numerical engine from `main`. The desktop application is still available.

```text
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

Open the local URL printed by Streamlit (normally `http://localhost:8501`). Choose inputs in the sidebar and click **Calculate**. The page includes all 17 result rows, a CSV download, and a Plotly 3D preview with lighting, transparency, and the representative red pore sphere. The browser view does not require Tkinter.

Calculation results and sampled display geometry are kept in each user's session. Transparency changes reuse those results. The dense analysis field is released after calculation. The web interface allows 20–250 analysis grid points per axis; memory and processing time increase approximately with the cube of this setting. The desktop engine's calculations are unchanged.

### Deploy with Streamlit Community Cloud

In [Streamlit Community Cloud](https://share.streamlit.io), create an app using:

- Repository: `TPMS-Analyzer/Repo`
- Branch: `streamlit-web-app`
- Main file path: `streamlit_app.py`
- Python: `3.12` (the version used for validation)

Authorize access to this private repository if requested. Dependencies are declared in the root `requirements.txt`, and the theme is in `.streamlit/config.toml`. The branch is ready to deploy; creating this branch does not publish a website.

See the [official deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy).

### Web validation

```text
python -m unittest -v test_tpms_analyzer.py test_streamlit_app.py
```

Tests cover numerical agreement with the desktop engine, sphere coordinates through transparency changes, calculation submission, retained results, invalid input handling, and CSV export. Streamlit's AppTest checks application behavior; browser interaction and WebGL rendering still depend on the client browser.

## Desktop application


Install Python 3.10 or later (with Tcl/Tk), open a terminal in this folder, and run:

```text
python -m pip install -r requirements.txt
python TPMS_Analyzer.py
```

On Windows, `py` may be used instead of `python`. Standard python.org Windows installers include Tkinter. On Linux, your distribution may require its `python3-tk` package.

The GUI keeps the original title, input labels, defaults, eight TPMS choices, Solid/Sheet choices, calculate button, transparency toggle, results table, and 3D preview with a red representative pore sphere. Drag the plot to rotate it; the plot toolbar provides navigation and image saving. Results also print to the terminal. Narrow windows hide the preview and preserve the controls and scrollable table. Calculations run in a background thread so the controls remain responsive.

## Preserved calculation details

- Exact eight level-set expressions, including the source's two-term Diamond expression.
- Inclusive 0-to-cell-size grid; MATLAB Y/X/Z array convention.
- Solid `F > c`, Sheet `abs(F) <= c`.
- Voxel-based threshold bisection, 60 iterations, 0.01 percentage-point stopping tolerance, and closest sampled threshold fallback.
- Internal solid/void face counting for wetted area; no external-face or periodic-wrap contributions.
- Single-cell Euclidean distance to solid, minus half a voxel; six-neighbor connectivity; directional radius bisection with 50 iterations and `dx/10` tolerance.
- `Dp = min(Dx, Dy, Dz)` and the original representative-center selection, including column-major tie order.
- All 17 result rows, significant-digit formatting, and original validation rules. Positive grid values ending in .5 round upward as in MATLAB.

## Rendering and precision

The interface uses Tkinter and Matplotlib, so native widget appearance, toolbar, lighting, and transparency rendering are not pixel-identical to MATLAB. The original layout and interaction are retained. The lattice has the original gold color and opacity (0.22 or 1); the red sphere stays opaque. Matplotlib transparent surfaces may show depth-order artifacts.

The display uses scikit-image marching cubes, with the boundary solid portions capped by a padded-and-clamped display mesh. Pore openings remain open. Mesh triangulation differs from MATLAB `isosurface`/`isocaps`; display mesh padding is never used for analysis. Wetted area remains the source's voxel face count, not triangle area.

SciPy replaces MATLAB Image Processing Toolbox operations. Distance results are cast to single precision to follow `bwdist`. Small floating-point differences, especially at exact level-set ties, can affect voxel counts and thresholds. No MATLAB installation was available for a direct numerical side-by-side comparison; do not assume bit-for-bit identity.

The original pore measurement is a voxel-based, nonperiodic approximation. The plotted sphere is a representative position, not an exact throat. Large grids require substantial memory and rendering time, with storage growing approximately as the cube of the grid size.

## Programmatic use

```python
from TPMS_Analyzer import analyze
result = analyze(network_type='Sheet', tpms_type='Gyroid',
                 target_porosity=70, alpha=2.54, grid_size=150)
print(result.rows)
print(result.pore.diameter_xyz)  # X, Y, Z in mm
```

Run the included checks with `python -m unittest -v test_tpms_analyzer.py`.
