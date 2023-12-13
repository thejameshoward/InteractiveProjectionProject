""" Main scatterplot figure. """

from typing import Dict, Any, List

import numpy as np
from plotly import graph_objects as go
from dash import Output, Input, no_update, dcc

from src.app import styles

def scatterplot():
    return dcc.Graph(
        id='scatterplot',
        className='scatterplot-graph',
        clear_on_unhover=True,
        style={'height': styles.SCATTERPLOT_HEIGHT}
    )