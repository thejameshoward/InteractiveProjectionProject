import magic
import base64

from io import BytesIO, StringIO
from PIL import PngImagePlugin
from pathlib import Path

import pandas as pd
import numpy as np

from dash import Dash, dcc, html, Input, Output, State, callback, ctx, no_update
from plotly import graph_objects as go

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from components.axes_selector import get_axes_selection_div, failure_div

from styles import *
from styles import button_style, SCATTERPLOT_TEXT_FONT_COLOR, upload_button_style
#
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)

#app = Dash(__name__)

df = None

selected_molecules = []

RDKIT_IMAGE_WIDTH = 250
RDKIT_IMAGE_HEIGHT = 250

app.layout = html.Div([ # Big block

                # This is the row div for the selection and example file button
                html.Div([

                    # Selector div
                    html.Div([
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div(
                                            ['Drag and drop or ',
                                            html.A('Select Files')]
                                            ),
                            style=upload_button_style,
                            multiple=True),
                            ],

                        style={'width': '45%', 'display': 'inline-block'}), # Style for row

                    html.Div([html.P("or...")], style={'width': '10%', 'display': 'inline-block'}),

                    html.Div([html.Button('Use example file', id='example-button', n_clicks=1, style=button_style)], style={'width': '45%','display': 'inline-block'}),
                ], style={'margin': 'auto', 'justify':'center', 'align':'center', 'textAlign':'center'}),


                # This is the display div
                html.Div(id='output-data-upload'),


            ], style={'backgroundColor':BG_COLOR})


def parse_contents(contents, filename, date):
    print('Called parse_contents()')
    # Content string is an iobytes string
    content_type, content_string = contents.split(',')

    # Decode the base string
    decoded = base64.b64decode(content_string)

    global df

    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    if 'SMILES' not in df.columns:
        return failure_div('(no SMILES column)')

    if 'NAME' not in df.columns:
        return failure_div('(no NAME column)')

    #print(f'Found smiles column')
    df['IMAGE'] = None
    df['IMAGE'] = df['SMILES'].apply(get_smiles_image)
    #print('Done ingesting dataframe')

    return get_axes_selection_div(filename, df, contents)

def get_smiles_image(smiles: str, **kwargs):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    drawer = rdMolDraw2D.MolDraw2DCairo(width=RDKIT_IMAGE_WIDTH, height=RDKIT_IMAGE_HEIGHT)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    im = BytesIO(drawer.GetDrawingText())
    im = PngImagePlugin.PngImageFile(im)
    buff = BytesIO()
    im.save(buff, format='png')
    encoded = base64.b64encode(buff.getvalue()).decode("utf-8")

    return 'data:image/png;base64, ' + encoded

@callback(Output('output-data-upload', 'children'),
          Input('upload-data', 'contents'),
          Input('example-button', 'n_clicks'),
          State('upload-data', 'filename'),
          State('upload-data', 'last_modified'),
          prevent_initial_call=True)
def update_output(list_of_contents, n_clicks, list_of_names, list_of_dates):

    #if isinstance(list_of_names, int):
    #    list_of_names = [list_of_names]

    # If file selection is used contents becomes the bytes object, names becomes the list of
    # file name strings, dates becomes
    trigger_id = ctx.triggered_id
    print(f'type(list_of_contents): {type(list_of_contents)}')
    print(f'list_of_names: {list_of_names}')
    print(f'dates: {list_of_dates}')
    print(f'n_clicks: {n_clicks}')
    print(f'trigger: {trigger_id}')
    print('\n')

    if trigger_id == 'example-button':
        p = Path(f'{Path(__file__).parent.name}/data/EXAMPLE_UMAP.csv')

        with open(p, 'rb') as infile:
            data = infile.read()
        content_bytes = base64.b64encode(data)
        content_string = content_bytes.decode("utf-8")

        mime = magic.Magic(mime=True)
        mime_type = mime.from_file(str(p.absolute()))
        content_type = "".join(["data:", mime_type, ";base64"])

        contents = "".join([content_type, ",", content_string])

        children = [parse_contents(contents, 'EXAMPLE_UMAP.csv', 0)]
        return children

    elif trigger_id == 'upload-data':
        if list_of_contents is not None:
            children = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_names, list_of_dates)]
            return children


    return children

@callback(Output('graph-div', 'children'),
          Input('x_selected', 'value'),
          Input('y_selected', 'value'),
          Input('z_selected', 'value'),
          Input('color_selected', 'value'),
          prevent_initial_call=True)
def update_graph(xcol, ycol, zcol, color_selected):
    print('In plotting')
    cols = [xcol, ycol, zcol]
    for col in cols:
        if col is None:
            print(f'Found None columns')
            return None

    global df

    # Get the figure object
    fig = get_scatter_trace(df, 3, xcol, ycol, zcol, color_col=color_selected)

    # Return a list since the output target is a div's children
    return [dcc.Graph(figure=fig, id='3dplot', style={'backgroundColor':BG_COLOR}, clear_on_unhover=True), dcc.Tooltip(id="graph-tooltip", style={'width':'150px', 'backgroundColor':BG_COLOR})]


def hover_template():
    return ("image %{customdata[0]}<br>")

def get_scatter_trace(df: pd.DataFrame,
                      ptsize: int,
                      xcol: str,
                      ycol: str,
                      zcol: str,
                      color_col: str = None) -> go.Scatter3d:
    x = df[xcol]
    y = df[ycol]
    z = df[zcol]

    if color_col is None:
        color = 'red'
    else:
        color = df[color_col]

    margin=dict(l=0, r=0, t=0, b=0)

    scene = dict(
        xaxis=dict(nticks=10, color=SCATTERPLOT_TEXT_FONT_COLOR, title=xcol),
        yaxis=dict(nticks=10, color=SCATTERPLOT_TEXT_FONT_COLOR, title=ycol),
        zaxis=dict(nticks=10, color=SCATTERPLOT_TEXT_FONT_COLOR, title=zcol),
        aspectmode='cube',
        dragmode='orbit',  # use orbital (not turntable) 3d rotation
        bgcolor=BG_COLOR,
        #annotations=player_unames
    )

    #fig = go.Figure()
    #fig.add_trace(main_trace)

    try:
        customdata = list(zip(df['IMAGE'].to_list(), df['NAME'].to_list()))
    except KeyError as e:
        customdata = list(zip(x, y, z))

    fig = go.Figure(go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=ptsize,
            #color=ptcolors[show_inds],
            #opacity=styles.SCATTERPLOT_PTS_OPACITY,
            color=color,
            colorscale='viridis',
            #line=dict(width=0)  # no lines between markers
        ),
        #hovertemplate=hover_template(),
        #customdata=[d for i, d in enumerate(hoverdata)]
        hoverinfo='none',
        hovertemplate=None,
        customdata=customdata,
    ))

    fig.update_layout(
        uirevision='constant',  # don't reset axes when updating plot
        showlegend=False,
        scene=scene,
        margin=margin
    )

    return fig

@app.callback(Output("graph-tooltip", "show"),
              Output("graph-tooltip", "bbox"),
              Output("graph-tooltip", "children"),
              Input("3dplot", "hoverData"),
              prevent_initial_call=True)
def display_hover(hoverData):

    if hoverData is None:
        return False, no_update, no_update

    pt = hoverData["points"][0]
    bbox = pt["bbox"]

    #if isinstance(hoverData['points'][0]['customdata'][0], float):
    #    children = [
    #    html.Div([
    #        html.P(hoverData['points'][0]['customdata'][0], style={"width": "100%", 'display': 'block', 'margin': '0 auto', 'backgroundColor':BG_COLOR}),
    #    ])
    #]
    #else:
    children = [
        html.Div([
            html.Img(src=hoverData['points'][0]['customdata'][0], style={"width": "100%", 'display': 'block', 'margin': '0 auto'}),
            html.P(hoverData['points'][0]['customdata'][1], style={"width": "100%", 'display': 'block', 'margin': '0 auto', 'backgroundColor':BG_COLOR}),
        ])
    ]


    return True, bbox, children

if __name__ == '__main__':
    app.run(debug=True)

