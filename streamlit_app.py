"""Run with: python -m streamlit run streamlit_app.py."""
import csv
import io

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from TPMS_Analyzer import TPMS_TYPES, analyze, resolve_display_digits
from web_mesh import detailed_mesh, split_surface_parts


def calculate_model(network, tpms, porosity, cell_text, grid, precision):
    cell_size = float(cell_text)
    digits = resolve_display_digits(cell_text, precision)
    result = analyze(network, tpms, porosity, cell_size, grid, digits)
    vertices, faces = detailed_mesh(result)
    center = result.pore.center_index[[1, 0, 2]] * cell_size/(result.grid_size-1)
    # Retain only results and display geometry in this user's session; release
    # the large analysis field after each calculation.
    return dict(vertices=vertices, faces=faces, center=center,
                parts=split_surface_parts(vertices, faces, cell_size),
                radius=result.pore.radius, diameter=result.pore.diameter,
                rows=result.rows, porosity=result.actual_porosity,
                area=result.wetted_area, alpha=cell_size, digits=digits,
                network=network, tpms=tpms, grid=result.grid_size,
                revision=repr((network, tpms, porosity, cell_size, grid, digits)))


SURFACE_COLORS = {
    'Gold (default)': '#b19a57',
    'Silver': '#b8c2cc',
    'Copper': '#b87333',
    'Blue': '#397dcc',
    'Teal': '#279b92',
    'Green': '#57964b',
    'Purple': '#8963bc',
    'Graphite': '#505862',
}


def preview_figure(model, transparent, finish='Satin', surface_color='Gold (default)'):
    materials = {'Matte': (.08, .85, .05), 'Satin': (.6, .4, .18),
                 'Glossy': (1.2, .22, .35)}
    specular, roughness, fresnel = materials[finish]
    lighting = dict(ambient=.42, diffuse=.78, specular=specular,
                    roughness=roughness, fresnel=fresnel,
                    facenormalsepsilon=1e-15, vertexnormalsepsilon=1e-15)
    light = dict(x=4*model['alpha'], y=6*model['alpha'], z=8*model['alpha'])
    fig = go.Figure()
    for name, vertices, faces, flat in model['parts']:
        fig.add_trace(go.Mesh3d(
            x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            color=SURFACE_COLORS[surface_color], opacity=.28 if transparent else 1.,
            flatshading=flat, name=name, hoverinfo='skip', showlegend=False,
            lighting=lighting, lightposition=light))
    if model['radius'] > 0 and np.all(np.isfinite(model['center'])):
        u = np.linspace(0, 2*np.pi, 97)
        v = np.linspace(-np.pi/2, np.pi/2, 97)
        radius, center = model['radius'], model['center']
        x = center[0] + radius*np.outer(np.cos(u), np.cos(v))
        y = center[1] + radius*np.outer(np.sin(u), np.cos(v))
        z = center[2] + radius*np.outer(np.ones_like(u), np.sin(v))
        fig.add_trace(go.Surface(
            x=x, y=y, z=z, surfacecolor=np.zeros_like(x),
            colorscale=[[0, '#e33726'], [1, '#e33726']],
            showscale=False, opacity=1., name='Representative pore',
            lighting=dict(ambient=.4, diffuse=.8, specular=.55, roughness=.25),
            lightposition=light,
            hovertemplate=f"Representative diameter: {model['diameter']:.6g} mm<extra></extra>"))
    # Do not put the WebGL clipping planes exactly on the cut faces. Floating
    # point clipping there causes speckles even on perfectly planar triangles.
    margin = .03*model['alpha']
    axis = dict(range=[-margin, model['alpha']+margin],
                tickvals=np.linspace(0., model['alpha'], 6),
                showbackground=True, showgrid=False,
                backgroundcolor='#f6f8fb', gridcolor='#dce3ec')
    fig.update_layout(
        height=620, margin=dict(l=0, r=0, t=20, b=0),
        paper_bgcolor='white', uirevision=model['revision'],
        scene=dict(xaxis=dict(axis, title='X (mm)'),
                   yaxis=dict(axis, title='Y (mm)'),
                   zaxis=dict(axis, title='Z (mm)'), aspectmode='cube',
                   dragmode='orbit', camera=dict(eye=dict(x=1.5, y=1.5, z=1.2))))
    return fig


def results_csv(rows):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(['Parameter', 'Value', 'Unit'])
    writer.writerows(rows)
    return output.getvalue()


def main():
    st.set_page_config(page_title='TPMS Geometry Analyzer', page_icon='🧊', layout='wide')
    st.title('TPMS Geometry Analyzer')
    st.caption('Explore periodic surfaces, porosity, and representative pore dimensions.')
    with st.sidebar:
        st.header('Model inputs')
        with st.form('model_inputs'):
            network = st.selectbox('Network type', ['Solid', 'Sheet'])
            tpms = st.selectbox('TPMS type', TPMS_TYPES)
            porosity = st.number_input('Target porosity (%)', min_value=.1,
                                       max_value=99.9, value=70., step=1.)
            cell_text = st.text_input('Unit cell size (mm)', '2.54')
            grid = st.number_input('Grid points per axis', min_value=20,
                                   max_value=250, value=150, step=10,
                                   help='Higher values use more time and memory. The web interface limits this to 250.')
            precision = st.text_input('Display significant digits', 'auto',
                                      help='auto, or a number from 2 to 15.')
            submitted = st.form_submit_button('Calculate', type='primary', width='stretch')
        st.caption('Calculations run when you press Calculate.')
    if submitted:
        try:
            with st.spinner('Calculating geometry and pore dimensions…'):
                model = calculate_model(network, tpms, porosity, cell_text, grid, precision)
            st.session_state['analysis'] = model
        except (ValueError, RuntimeError, MemoryError) as exc:
            st.error(f'Calculation failed: {exc}')
    model = st.session_state.get('analysis')
    if model is None:
        st.info('Choose model inputs in the sidebar, then press Calculate to display the 3D geometry.')
        return
    st.subheader(f"{model['network']} · {model['tpms']}")
    st.caption(f"Displayed result: cell size {model['alpha']:g} mm · {model['grid']} grid points per axis")
    metrics = st.columns(3)
    digits = model['digits']
    metrics[0].metric('Calculated porosity', f"{model['porosity']:.{digits}g}%")
    metrics[1].metric('Representative diameter', f"{model['diameter']:.{digits}g} mm")
    metrics[2].metric('Wetted area', f"{model['area']:.{digits}g} mm²")
    view, results = st.columns([1.6, 1])
    with view:
        transparent = st.toggle('Transparent lattice', value=True)
        finish = st.selectbox('Surface finish', ['Satin', 'Glossy', 'Matte'])
        surface_color = st.selectbox('Surface color', list(SURFACE_COLORS))
        st.plotly_chart(preview_figure(model, transparent, finish, surface_color), width='stretch',
                        theme=None, key='tpms_preview',
                        config=dict(scrollZoom=True, displaylogo=False))
        st.caption('Drag to rotate. Scroll to zoom. Use the chart toolbar to pan or reset the camera.')
    with results:
        st.subheader('Results')
        st.dataframe([dict(Parameter=n, Value=v, Unit=u) for n, v, u in model['rows']],
                     hide_index=True, width='stretch', height=620)
        st.download_button('Download results (CSV)', results_csv(model['rows']),
                           file_name=f"TPMS_{model['network']}_{model['tpms']}.csv", mime='text/csv')
    st.caption('The red sphere marks the representative pore location from the numerical analysis. '
               'The detailed preview uses up to 150 grid samples per axis; results use the full analysis grid.')


if __name__ == '__main__':
    main()
