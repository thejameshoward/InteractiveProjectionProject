from dash import dcc, html, dash_table, Input, Output, State, callback, ctx

from styles import BG_COLOR

def get_axes_selection_div(filename, df, contents) -> html.Div:
    # Parent Div
    _ = html.Div([
            html.Div([
                html.H5(filename),
                #html.H6(datetime.datetime.fromtimestamp(date)),

                dcc.Dropdown([x for x in df.columns], id='x_selected', placeholder='Select X', style={'color':'white'}),

                dcc.Dropdown([x for x in df.columns], id='y_selected', placeholder='Select Y', style={'color':'white'}),
                #html.P("Select Z:"),
                dcc.Dropdown([x for x in df.columns], id='z_selected', placeholder='Select Z', style={'color':'white'}),
                #html.P("Color by:"),
                dcc.Dropdown([x for x in df.columns], id='color_selected', placeholder='Select color', style={'color':'white'}),

                # For debugging, display the raw contents provided by the web browser
                #html.Div('Raw Content'),
                #html.Pre(contents[0:200] + '...', style={
                #    'whiteSpace': 'pre-wrap',
                #    'wordBreak': 'break-all'
                #})
            ], style={'width':'50%'}),

            # horizontal line
            html.Hr(),

            html.Div([], id='graph-div', style={'backgroundColor':'blue'}),

    ],)

    return _


def failure_div(reason: str) -> html.Div:
    reason = f'Invalid CSV file {reason}'

    return html.Div(
        [
            html.H1(reason, style={'textAlign': 'center'}),
            dcc.Textarea(
                id='failure-text-area',
                value='Your CSV file is improperly formatted. Please ensure that you have\n \
                    the following columns. \n\n \
                        NAME: name or ID of your compounds\n \
                        SMILES: the smiles string that describes your compounds',
                style={'width': '100%', 'height': 300, 'textAlign':'center', 'backgroundColor':'#333333', 'color':'white'}
            )
         ]
         )