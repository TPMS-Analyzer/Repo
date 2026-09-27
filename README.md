# TPMS Geometry Analyzer

Explore triply periodic minimal surface (TPMS) structures in 3D and calculate their porosity, wetted surface area, and representative pore dimensions.

**[Open the TPMS Analyzer website — recommended](https://tpmsanalyzer.streamlit.app/)**

The hosted website is the easiest way to use the app: no Python installation is needed.

## Choose how to run the app

| Option | How it runs | Best for |
| --- | --- | --- |
| **Hosted website (recommended)** | Open [tpmsanalyzer.streamlit.app](https://tpmsanalyzer.streamlit.app/) in your browser. | Getting started immediately without installing software. |
| **Local Streamlit app** | Run Python on your computer and use the app in your browser. | The web interface with calculations performed on your own computer. |
| **Local desktop app** | Run Python to open a standalone application window. | Using the desktop interface without a browser. |

Both local options use the same underlying numerical calculations. The Streamlit interface includes a detailed interactive preview, surface finishes, contrasting color schemes, and CSV results downloads.

## Using the app

1. Choose a **network type**: Solid or Sheet.
2. Choose a **TPMS type**: Koch, Gyroid, Diamond, Primitive, IWP, Neovius, FRD, or PMY.
3. Enter the **target porosity (%)** and **unit cell size (mm)**.
4. Set the **grid points per axis**. The default is 150; larger grids require more memory and calculation time.
5. Leave **display significant digits** at `auto`, or enter a value from 2 to 15.
6. Press **Calculate** to display the geometry and results.

In the web interface, drag to rotate, scroll to zoom, and use the chart toolbar to pan or reset the view. Change transparency, choose a Satin, Glossy, or Matte finish, and select a surface color. For two-tone options, the first color applies to the curved surface and the second to the flat cut faces. The red sphere shows the representative pore location and diameter.

The results include calculated porosity, solid volume fraction, wetted area, unit-cell volume, and directional pore diameters. Use **Download results (CSV)** in the web interface to save the results.

## Local installation

**Anaconda is recommended**, using a separate environment with **Python 3.12**.

### 1. Download the application

Download the `main` branch from [this repository](https://github.com/TPMS-Analyzer/Repo) using **Code → Download ZIP**, then extract it. Keep the application files together in the extracted folder.

### 2. Create a Python environment

Install [Anaconda Distribution](https://www.anaconda.com/download). Open **Anaconda Prompt** on Windows, or a terminal with conda available on macOS/Linux, and run:

```bash
conda create -n tpms-analyzer python=3.12 tk pip
conda activate tpms-analyzer
```

### 3. Install the required packages

In the same terminal, change to the extracted application folder. Replace the example path below with the folder on your computer:

```bash
cd "path/to/extracted/application"
python -m pip install -r requirements.txt
```

The requirements file installs the packages needed for both local interfaces:

| Package | Purpose |
| --- | --- |
| NumPy | Numerical arrays and geometry calculations |
| SciPy | Distance and connectivity calculations |
| scikit-image | Building the 3D surface mesh |
| Matplotlib | Desktop 3D preview |
| Streamlit | Local browser interface |
| Plotly | Interactive web 3D preview |

The desktop interface also needs **Tkinter/Tcl-Tk**, provided by the `tk` package in the conda command above. It is not installed through pip. Package version requirements are listed in [requirements.txt](requirements.txt).

## Run locally through Streamlit

From the application folder, with the environment activated, run:

```bash
python -m streamlit run streamlit_app.py
```

Open the local address printed in the terminal, normally **http://localhost:8501**. Keep the terminal running while using the app. Press **Ctrl+C** in that terminal to stop it.

This runs the application on your computer; it does not publish a website.

## Run locally as a desktop app

From the same folder and environment, run:

```bash
python TPMS_Analyzer.py
```

A standalone window opens with model inputs, calculation results, a transparency toggle, and a 3D preview. Use the plot toolbar for navigation and image saving. Results also print in the terminal.

For either local option, activate the environment again when opening a new terminal:

```bash
conda activate tpms-analyzer
```

## Troubleshooting

- **A package is missing:** activate `tpms-analyzer`, then run `python -m pip install -r requirements.txt` from the application folder.
- **The desktop window does not open:** run `python -m tkinter` to check Tkinter. If it is missing, install it in the active conda environment with `conda install tk`. The desktop app requires a graphical desktop session.
- **The browser does not open automatically:** copy the local URL printed by Streamlit into your browser.
- **Calculations take too long or use too much memory:** reduce the grid points per axis. The web interface accepts values from 20 to 250.
- **A transparent surface is difficult to inspect:** turn transparency off to examine the curved surface and flat cut faces.

## Understanding the results

Results are numerical estimates from a sampled unit-cell grid, so grid resolution affects accuracy. The displayed surface may use fewer samples than the calculation grid. Wetted area is calculated from internal solid/void interfaces in the grid, rather than from the displayed triangles.

The red sphere represents the calculated pore size at a representative location; it is not an exact reconstruction of a physical pore throat. Pore calculations use a single unit cell without periodic wrapping.
