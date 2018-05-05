### CRISPR-SURF Dash Application
import dash
from dash.dependencies import Input, Output, State, Event
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import plotly.graph_objs as go
from plotly import tools
from flask import Flask, request, redirect, url_for, render_template, jsonify, send_from_directory, send_file, session
import requests
import re
import subprocess as sb
import numpy as np
import pandas as pd
import sys
import os
import glob
import base64
import urllib
import ast
import random
import gzip
import uuid
import json
import cPickle as cp
import csv
import time
import io

# Import CRISPR-SURF Functions
from CRISPR_SURF_Core_WebApp import gaussian_pattern, crispr_surf_deconvolution, crispr_surf_find_gamma, crispr_surf_deconvolved_signal, crispr_surf_statistical_significance, crispr_surf_sgRNA_summary_table_update, complete_beta_profile, crispr_surf_significant_regions, crispr_surf_IGV, str2bool, reverse_complement, total_count_normalization, median_normalization, normalize_sgRNA_counts

_ROOT = os.path.abspath(os.path.dirname(__file__))

class Dash_responsive(dash.Dash):
    def __init__(self, *args, **kwargs):
        super(Dash_responsive,self).__init__(*args, **kwargs)

    def index(self, *args, **kwargs):  # pylint: disable=unused-argument
        scripts = self._generate_scripts_html()
        css = self._generate_css_dist_html()
        config = self._generate_config_html()
        title = getattr(self, 'title', 'CRISPR-SURF')
        return '''
        <!DOCTYPE html>
        <html>
            <head>
                <meta charset="UTF-8">
                <title>{}</title>
                {}
                <link rel="stylesheet" href="/static/CRISPR-SURF.css">
                <link rel="stylesheet" href="/static/Loading-State.css">
                <script src="/static/jquery-3.3.1.min.js"></script>
            </head>
            <body>
                <div id="react-entry-point">
                    <div class="_dash-loading">
                        Loading...
                    </div>
                </div>
                <footer>
                    {}
                    {}
                </footer>
            </body>
        </html>
        '''.format(title, css, config, scripts)

def get_data(path):
        return os.path.join(_ROOT, path)

# Logos
crisprsurf_logo = get_data('crisprsurf_logo.png')
crisprsurf_logo_image = base64.b64encode(open(crisprsurf_logo, 'rb').read())

mgh_logo = get_data('mgh.png')
mgh_logo_image = base64.b64encode(open(mgh_logo, 'rb').read())

mitbe_logo = get_data('mitbe.png')
mitbe_logo_image = base64.b64encode(open(mitbe_logo, 'rb').read())

hms_logo = get_data('hms.png')
hms_logo_image = base64.b64encode(open(hms_logo, 'rb').read())

# Generate ID to initialize CRISPR-SURF instance
server = Flask(__name__)
server.secret_key = '\x1f\x1a*\x88\xbe\xb3"\xc9\xa5n\xf3\x9c)\xd5\xec\xa6C\xf2b:\xf1\xa7\x8e\xa8'
app = Dash_responsive('crisprsurf-compute', server = server, url_base_pathname = '/compute/', csrf_protect=False)

app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

app2 = Dash_responsive(name = 'crisprsurf-app-precomputed', server = server, url_base_pathname = '/precomputed/', csrf_protect=False)

app.server.config['UPLOADS_FOLDER']='/tmp/UPLOADS_FOLDER'
app.server.config['RESULTS_FOLDER']='/tmp/RESULTS_FOLDER'
app.server.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024 #200 Mb max

# app.config['suppress_callback_exceptions'] = True

@app.server.route('/static/<path:path>')
def static_file(path):
    static_folder = os.path.join(os.getcwd(), 'static')
    return send_from_directory(static_folder, path)

@server.route('/help')
def help():
	return render_template('help.html')

@server.route('/')
def index():
    newpath = 'CRISPR-SURF_' + str(uuid.uuid4())

    UPLOADS_FOLDER = os.path.join(app.server.config['UPLOADS_FOLDER'], newpath)
    RESULTS_FOLDER = os.path.join(app.server.config['RESULTS_FOLDER'], newpath)

    if not os.path.exists(UPLOADS_FOLDER):
        os.makedirs(UPLOADS_FOLDER)

    if not os.path.exists(RESULTS_FOLDER):
        os.makedirs(RESULTS_FOLDER)

    param_dict = {'example_data_clicks':0,'deconvolution-clicks':0, 'significance-clicks':0, 'checkbutton1':1, 'checkbutton2':1,}

    data_dict = {'gammas2betas':None, 'guideindices2bin':None, 'computed':False, 'time_per_gamma': 0}

    data_dict2 = {'gammas2betas':None, 'computed':False, 'fdr':0.05}

    with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
        json_string = json.dumps(param_dict)
        f.write(json_string + '\n')

    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
        json_string = json.dumps(data_dict)
        f.write(json_string + '\n')

    with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
        json_string = json.dumps(data_dict2)
        f.write(json_string + '\n')

    return render_template('index.html',newpath=newpath)

#Import some other useful functions
def generate_table(dataframe, max_rows = 100):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))]
    )

# file upload function
def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)

    if filename.endswith('.csv'):
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        return df
    else:
        return None

### Compute Layout ###
app.layout = html.Div([

	dcc.Location(id='url', refresh=False),

	dcc.Interval(id='common-interval-1', interval=1000000),
	dcc.Interval(id='common-interval-2', interval=1000000),
	# dcc.Interval(id='common-interval-3', interval=1000000),
	# dcc.Interval(id='common-interval-4', interval=1000000),

	html.Div(id = 'custom-loading-states-1',
		children = [

		html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

		], style = {'display':'none'}),

	html.Div(id = 'custom-loading-states-2',
		children = [

		html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

		], style = {'display':'none'}),

	html.Div(id = 'custom-loading-states-3',
		children = [

		html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

		], style = {'display':'none'}),

	html.Div(id = 'custom-loading-states-4',
		children = [

		html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

		], style = {'display':'none'}),

	html.Img(src='data:image/png;base64,{}'.format(crisprsurf_logo_image), width = '50%'),
	html.H2('CRISPR Screening Uncharacterized Region Function'),

    dcc.Tabs(
        tabs=[
            {'label': 'Step 1: Upload sgRNA Table', 'value': 'upload'},
            {'label': 'Step 2: Analysis Parameters', 'value': 'parameters'},
            {'label': 'Step 3: Run CRISPR-SURF', 'value': 'deconvolution'},
            {'label': 'Step 4: Find Regions', 'value': 'significance'},
            {'label': 'Step 6: Download Results', 'value': 'download'}
        ],
        value='upload',
        id='tabs',
        style = {'font-weight':'bold'}
    ),

    html.Div(id = 'upload-container-total', children = [

        # Drag and Drop upload component
        html.H3('Upload sgRNA Summary Table'),

        html.Button(id = 'example-data-button', children = 'Load Example Data', n_clicks = 0),

        html.Label(id = 'upload-log', children = 'Upload Status: Incomplete', style = {'font-weight':'bold'}),

        html.Label('Supported File Formats: .csv and .xls', style = {'font-weight':'bold'}),

        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ], style = {'font-weight':'bold', 'font-size': 20}),
            style={
                'width': '100%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '2px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
            },
            multiple=False),

        dt.DataTable(id='sgRNA-table-view', filterable=True, sortable=True, column_widths = [250]*100000, rows=[{}]),

        ], className = 'row', style = {'display':'none'}),

    html.Div(id = 'parameters-notification-container-total', children = [html.H2('Please Complete Step 1')], style = {'display':'none'}),

    html.Div(id = 'parameters-container-total', children = [

        html.Div([

            html.H3('Basic Parameters'),

            html.Div([

                html.Div([

                    html.Div([

                        html.Label('CRISPR Perturbation Type', style = {'font-weight':'bold'}),
                        dcc.Dropdown(
                            id = 'pert',
                            options=[
                                {'label': 'CRISPR-Cas9', 'value': 'cas9'},
                                {'label': 'CRISPR-Cpf1', 'value': 'cpf1'},
                                {'label': 'CRISPR Interference', 'value': 'crispri'},
                                {'label': 'CRISPR Activation', 'value': 'crispra'},
                                {'label': 'Base Editor', 'value': 'be'}
                            ],
                            value = 'cas9',
                        ),

                        ], className = 'four columns'),

                    html.Div([

                        html.Label('Genome', style = {'font-weight':'bold'}),
                        dcc.Dropdown(
                            id = 'genome',
                            options=[
                                {'label': 'hg19', 'value': 'hg19'},
                                {'label': 'hg38', 'value': 'hg38'},
                                {'label': 'mm9', 'value': 'mm9'},
                                {'label': 'mm10', 'value': 'mm10'},
                            ],
                            value = 'hg19',
                        ),

                        ], className = 'four columns'),

                    html.Div([

                        html.Label('Combining Replicates', style = {'font-weight':'bold'}),
                        dcc.RadioItems(
                            id = 'avg',
                            options=[
                                {'label': 'Mean', 'value': 'mean'},
                                {'label': 'Median', 'value': 'median'},
                            ],
                            value = 'median',
                            labelStyle={'display': 'inline-block'}
                            ),

                        ], className = 'four columns'),

                    ], className = 'row'),

                html.Label(id = 'perturbation-range-show', children = 'Characteristic Perturbation Length: None', style = {'font-weight':'bold'}),
                dcc.Slider(
                    id = 'range',
                    min=1,
                    max=300,
                    step=1,
                    value = 7,
                    ),

                html.Label('Perturbation Profile', style = {'font-weight':'bold'}),
                dcc.Graph(id='perturbation-profile', animate=False),

                ], style = {'border-style':'solid','border-width':'2px', 'border-color':'#DCDCDC','border-radius':'10px','border-spacing':'15px','padding':'10px'})

            ], className = 'six columns'),

        html.Div([

            html.H3('Advanced Parameters'),

            html.Button(id = 'advanced-parameters-button', children = 'Show (+)', n_clicks = 0),

            html.Br(),

            html.Div(id = 'advanced-parameters-container', children = [

                html.H5('Deconvolution Parameters'),

                html.Div([

                    html.Div([

                        html.Label(id = 'scale-show', children = 'Scale: 1 bp', style = {'font-weight':'bold'}),
                        dcc.Slider(
                            id = 'scale',
                            min=1,
                            max=50,
                            step=1,
                            value = 1
                            ),

                        html.Label('Gamma List (i.e. 1,2,3,4,5)', style = {'font-weight':'bold'}),
                        dcc.Input(
                            id = 'gamma_list',
                            placeholder='Use Default ...',
                            type='text',
                            value = ''
                            ),

                        ], className = 'six columns'),

                    html.Div([

                        html.Label(id = 'limit-show', children = 'Limit: 25 bps', style = {'font-weight':'bold'}),
                        dcc.Slider(
                            id = 'limit',
                            min=1,
                            max=500,
                            step=1,
                            value = 25
                            ),

                        html.Label('Use Rapid Mode', style = {'font-weight':'bold'}),
                        dcc.RadioItems(
                            id = 'rapid',
                            options=[
                                {'label': 'Yes', 'value': True},
                                {'label': 'No', 'value': False},
                            ],
                            value = False,
                            labelStyle={'display': 'inline-block'}
                            ),

                        ], className = 'six columns'),

                    ], className = 'row'),

                html.Hr(),

                html.H5('Simulation Parameters'),

                html.Div([

                    html.Div([

                        html.Label('Null Distribution Simulation', style = {'font-weight':'bold'}),
                        dcc.Dropdown(
                            id = 'sim_type',
                            options=[
                                {'label': 'Gaussian', 'value': 'gaussian'},
                                {'label': 'Laplace', 'value': 'laplace'},
                                {'label': 'Negative Control sgRNAs', 'value': 'negative_control'}
                            ],
                            value = 'gaussian',
                        ),

                        html.Label(id = 'simulation-n-show', children = 'Number of Simulations: 100', style = {'font-weight':'bold'}),
                        dcc.Slider(
                            id = 'sim_n',
                            min=10,
                            max=1000,
                            step=10,
                            value = 100,
                            ),

                        ], className = 'six columns'),

                    html.Div([

                        html.Label('Effect Size', style = {'font-weight':'bold'}),
                        dcc.Input(
                            id = 'effect_size',
                            type='number',
                            value = 1,
                            ),

                        ], className = 'six columns')

                    ], className = 'row'),

                ]),

            html.Hr(),
            
            html.H3('Parameters Help'),

            html.Button(id = 'parameters-help-button', children = 'Show (+)', n_clicks = 0),

            html.Br(),

            html.Div(id = 'parameters-help-container', children = [

                html.Label('Combining Replicates', style = {'font-weight':'bold'}),
                html.Div(['The method for combining deconvolution scores across biological replicates.']),

                html.Label('Characteristic Perturbation Length', style = {'font-weight':'bold'}),
                html.Div(['The average perturbation length of the CRISPR screening modality used.']),

                html.Label('Scale', style = {'font-weight':'bold'}),
                html.Div(['Determines the resolution of output. Scale of 1 outputs scores for single bps, while scale of 20 outputs scores for groups of 20 bps.']),

                html.Label('Limit', style = {'font-weight':'bold'}),
                html.Div(['Maximum perturbation length (from center) of the CRISPR screening modality used.']),

                html.Label('Gamma List', style = {'font-weight':'bold'}),
                html.Div(['List of gammas to use for deconvolution computation.']),

                html.Label('Use Rapid Mode', style = {'font-weight':'bold'}),
                html.Div(['Rapid mode aggregates null distributions across all bps. Speeds up computation at the cost of eliminating density-aware statistical tests.']),

                html.Label('Null Distribution Simulation', style = {'font-weight':'bold'}),
                html.Div(['Method for determining null L2FC distribution to undergo deconvolution simulations. Gaussian and Laplace parameterize observed and negative control sgRNAs. Negative Control sgRNAs samples from only negative control sgRNA L2FC scores.']),

                html.Label('Number of Simulations', style = {'font-weight':'bold'}),
                html.Div(['Number of deconvolution simulations to construct null distribution.']),

                html.Label('Effect Size', style = {'font-weight':'bold'}),
                html.Div(['The effect size of interest for estimating statistical power across CRISPR tiling screen.']),

                ]),

            ], className = 'six columns'),

        ], className = 'row', style = {'display':'none'}),

    html.Div(id = 'deconvolution-notification-container-total', children = [html.H2('Please Complete Steps 1 and 2')], style = {'display':'none'}),

    html.Div(id = 'deconvolution-container-total', children = [

        html.Div(id = 'custom-loading-states-1',
            children = [

            html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

            ], style = {'display':'none'}),

        # dcc.Interval(id='common-interval-1', interval=1000000),

        html.H3('Run CRISPR-SURF'),

        html.Button("Let's SURF!", id='deconvolution-button', n_clicks = 0),

        html.Hr(),

        html.Div([

            html.Div(id = 'deconvolution-container', children = [

                html.Div([

                    html.Div([

                        html.Div([

                            html.Div([

                                html.Label('Chr', style = {'font-weight':'bold'}),
                                dcc.Dropdown(
                                    id = 'chr',
                                    options=[],
                                    value = None
                                ),

                                ], className = 'three columns'),

                            html.Div([

                                html.Label('Start', style = {'font-weight':'bold'}),
                                dcc.Input(
                                    id = 'start',
                                    type = 'text',
                                    placeholder='Enter Start Coordinate ...',
                                    ),

                                ], className = 'three columns', style = {'text-align':'center'}),

                            html.Div([

                                html.Label('Stop', style = {'font-weight':'bold'}),
                                dcc.Input(
                                    id = 'stop',
                                    type = 'text',
                                    placeholder='Enter Stop Coordinate ...',
                                    ),

                                ], className = 'three columns', style = {'text-align':'center'}),

                            html.Div([

                                html.Br(),
                                html.Br(),

                                html.Button(id = 'update-graph-button1', children = 'Update Graph'),

                                ], className = 'three columns', style = {'text-align':'left'}),

                            ], className = 'row'),

                        ], className = 'eight columns'),

                    ], className = 'row'),
                
                html.Div([

                    html.Div([

                        html.Label(id = 'replicate-correlation-show', children = 'Replicate Correlation: 0.8', style = {'font-weight':'bold'}),
                        dcc.Slider(
                            id = 'replicate-correlation',
                            min=0.5,
                            max=1.0,
                            step=0.01,
                            value = 0.8,
                            ),

                        ], className = 'six columns')

                    ], className = 'row'),

                ], style = {'display': 'none'}),

            dcc.Graph(id='deconvolution-plot', animate=True),

            ])

        ], style = {'display':'none'}),

    html.Div(id = 'significance-notification-container-total', children = [], style = {'display':'none'}),

    html.Div(id = 'significance-container-total', children = [

        html.Div(id = 'custom-loading-states-2',
            children = [

            html.Div(id = 'custom-loading-state1', className = '_dash-loading-callback_custom', children = ['Loading...', html.Center(children=[html.Div(id = 'custom-loading-state2', className = 'loader', style = {'display':'block'})])],  style = {'display':'block'})

            ], style = {'display':'none'}),

        # dcc.Interval(id='common-interval-2', interval=1000000),

        html.H3('Find Regions'),

        html.Button(id = 'significance-button', children = 'Find Regions', n_clicks = 0),

        html.Label(id = 'time-estimate', children = 'Time Estimate: NA', style = {'font-weight':'bold'}),

        html.Hr(),

        html.Div([

            html.Div(id = 'significance-container', children = [

                html.Div([

                    html.Div([

                        html.Div([

                            html.Div([

                                html.Label('Chr', style = {'font-weight':'bold'}),
                                dcc.Dropdown(
                                    id = 'chr2',
                                    options=[],
                                    value = None
                                ),

                                ], className = 'three columns'),

                            html.Div([

                                html.Label('Start', style = {'font-weight':'bold'}),
                                dcc.Input(
                                    id = 'start2',
                                    type = 'text',
                                    placeholder='Enter Start Coordinate ...',
                                    ),

                                ], className = 'three columns', style = {'text-align':'center'}),

                            html.Div([

                                html.Label('Stop', style = {'font-weight':'bold'}),
                                dcc.Input(
                                    id = 'stop2',
                                    type = 'text',
                                    placeholder='Enter Stop Coordinate ...',
                                    ),

                                ], className = 'three columns', style = {'text-align':'center'}),

                            html.Div([

                                html.Br(),
                                html.Br(),

                                html.Button(id = 'update-graph-button2', children = 'Update Graph'),

                                ], className = 'three columns', style = {'text-align':'left'}),

                            ], className = 'row'),

                        ], className = 'eight columns'),

                    ], className = 'row'),
                
                html.Div([

                    html.Div([

                        html.Label(id = 'fdr-show', children = 'FDR < 0.05', style = {'font-weight':'bold'}),
                        dcc.Slider(
                            id = 'fdr',
                            min=0.01,
                            max=0.25,
                            step=0.01,
                            value = 0.05,
                            ),

                        ], className = 'six columns')

                    ], className = 'row'),

                ], style = {'display': 'none'}),

            dcc.Graph(id='significance-plot', animate=True),

            ])

        ], style = {'display':'none'}),

    html.Div(id = 'download-notification-container-total', children = [], style = {'display':'none'}),

    html.Div(
        id = 'download-container-total',
        children = [

        html.Div([

            html.Label('Experiment Title', style = {'font-weight':'bold', 'padding-right':'10px'}),
            dcc.Input(id = 'title-input', value = '', type = 'text', placeholder = 'Optional'),

            html.Label('Experiment Description', style = {'font-weight':'bold', 'padding-right':'10px'}),
            dcc.Input(id = 'description-input', value = '', type = 'text', placeholder = 'Optional'),

            ]),

        html.Br(),

        html.Div([

            html.A(
                html.Button('Download Files'),
                id='download-total',
                download = "surf-outputs.zip",
                href="",
                target="_blank",
                n_clicks = 0,
                style = {'font-weight':'bold', 'font-size':'100%', 'text-align':'center'}
                ),

            ]),

        ], style = {'display':'none'}),


    html.Br(),
    html.Br(),

    html.Div(id = 'tmp1', style = {'display':'none'}),
    html.Div(id = 'tmp2', style = {'display':'none'}),
    html.Div(id = 'tmp3', style = {'display':'none'}),
    html.Div(id = 'tmp4', style = {'display':'none'}),
    html.Div(id = 'tmp5', style = {'display':'none'}),
    html.Div(id = 'tmp6', style = {'display':'none'}),
    html.Div(id = 'tmp7', style = {'display':'none'}),
    html.Div(id = 'tmp8', style = {'display':'none'}),
    html.Div(id = 'tmp9', style = {'display':'none'}),
    html.Div(id = 'tmp10', style = {'display':'none'}),

    ])

# Pseudo tabs code
@app.callback(
    Output('upload-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    if value == 'upload':
        return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(
    Output('parameters-notification-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    if value == 'parameters':
        if os.path.isfile(fname):
            return {'display':'none'}
        else:
            return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(
    Output('parameters-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    if value == 'parameters':
        if os.path.isfile(fname):
            return {'display':'block'}
        else:
            return {'display':'none'}
    else:
        return {'display':'none'}

@app.callback(
    Output('deconvolution-notification-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    if value == 'deconvolution':
        if os.path.isfile(fname):
            return {'display':'none'}
        else:
            return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(
    Output('deconvolution-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    if value == 'deconvolution':
        if os.path.isfile(fname):
            return {'display':'block'}
        else:
            return {'display':'none'}
    else:
        return {'display':'none'}

@app.callback(
    Output('significance-notification-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'significance':
        if os.path.isfile(fname):
            if data_dict['computed']:
                return {'display':'none'}
            else:
                return {'display':'block'}
        else:
            return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(
    Output('significance-notification-container-total', 'children'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'significance':
        if os.path.isfile(fname):
            if data_dict['computed']:
                return html.Div()
            else:
                return html.Div([html.H2("Please complete Step 3")])
        else:
            return html.Div([html.H2("Please complete Steps 1, 2, and 3")])
    else:
        return html.Div()

@app.callback(
    Output('significance-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'significance':
        if os.path.isfile(fname):
            if data_dict['computed']:
                return {'display':'block'}
            else:
                return {'display':'none'}
        else:
            return {'display':'none'}
    else:
        return {'display':'none'}

@app.callback(
    Output('download-notification-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'download':
        if os.path.isfile(fname):
            if data_dict['computed']:
                if data_dict2['computed']:
                    return {'display':'none'}
                else:
                    return {'display':'block'}
            else:
                return {'display':'block'}
        else:
            return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(
    Output('download-notification-container-total', 'children'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'download':
        if os.path.isfile(fname):
            if data_dict['computed']:
                if data_dict2['computed']:
                    return html.Div()
                else:
                    return html.Div([html.H2("Please complete Step 4")])
            else:
                return html.Div([html.H2("Please complete Steps 3 and 4")])
        else:
            return html.Div([html.H2("Please complete Steps 1, 2, 3, and 4")])
    else:
        return html.Div()

@app.callback(
    Output('download-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    fname = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    if value == 'download':
        if os.path.isfile(fname):
            if data_dict['computed']:
                if data_dict2['computed']:
                    return {'display':'block'}
                else:
                    return {'display':'none'}
            else:
                return {'display':'none'}
        else:
            return {'display':'none'}
    else:
        return {'display':'none'}

# CALLBACKS FOR UPLOAD
@app.callback(Output('upload-log', 'children'),
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename'),
               Input('example-data-button', 'n_clicks')],
               state = [State('url', 'pathname')])

def download_file(contents, filename, n_clicks, pathname):

    if pathname:

        UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
        RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

        try:
            sb.call('rm %s/sgRNAs_summary_table.csv' % UPLOADS_FOLDER, shell = True)
        except:
            pass

        # hack
        json_good = False
        while not json_good:
            with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
                json_string = f.readline().strip()
                try:
                    param_dict = json.loads(json_string)
                    json_good = True
                except:
                    pass

        if n_clicks > param_dict['example_data_clicks']:

            sb.call('cp /SURF/exampleDataset/sgRNAs_summary_table.csv %s/sgRNAs_summary_table.csv' % (UPLOADS_FOLDER), shell = True)

            sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
            df = pd.read_csv(sgRNA_summary)

            json_good = False
            while not json_good:
                with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
                    json_string = f.readline().strip()
                    try:
                        data_dict = json.loads(json_string)
                        json_good = True
                    except:
                        pass

            json_good = False
            while not json_good:
                with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
                    json_string = f.readline().strip()
                    try:
                        data_dict2 = json.loads(json_string)
                        json_good = True
                    except:
                        pass

            param_dict['example_data_clicks'] = n_clicks
            data_dict['chr'] = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()][0]
            data_dict2['chr'] = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()][0]

            with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
                new_json_string = json.dumps(param_dict)
                f.write(new_json_string + '\n')

            with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
                new_json_string = json.dumps(data_dict)
                f.write(new_json_string + '\n')

            with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
                new_json_string = json.dumps(data_dict)
                f.write(new_json_string + '\n')

            return 'Upload Status: Complete'

        if contents is not None:
            df = parse_contents(contents, filename)
            if df is not None:

                if all(x in df.columns for x in ['Chr', 'Start', 'Stop', 'Perturbation_Index', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type', 'Log2FC_Replicate1']):
                    df.to_csv(UPLOADS_FOLDER + '/sgRNAs_summary_table.csv', index = False)

                    with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
                        json_string = f.readline().strip()
                        data_dict = json.loads(json_string)

                    with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
                        json_string = f.readline().strip()
                        data_dict2 = json.loads(json_string)

                    data_dict['chr'] = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()][0]
                    data_dict2['chr'] = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()][0]

                    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
                        new_json_string = json.dumps(data_dict)
                        f.write(new_json_string + '\n')

                    with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
                        new_json_string = json.dumps(data_dict)
                        f.write(new_json_string + '\n')

                    return 'Upload Status: Complete'

                else:
                    return 'Upload Status: Incomplete (Required Columns - Chr, Start, Stop, Perturbation_Index, sgRNA_Sequence, Strand, sgRNA_Type, Log2FC_Replicate1, etc.)'

            else:
                return 'Upload Status: Incomplete (Please upload .csv format)'

        return 'Upload Status: Incomplete'


    return 'Upload Status: Incomplete'

@app.callback(Output('sgRNA-table-view', 'rows'),
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename'),
               Input('example-data-button', 'n_clicks')],
               state = [State('url', 'pathname')])

def update_output(contents, filename, n_clicks, pathname):

    # if pathname:

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if n_clicks > param_dict['example_data_clicks']:

        try:
            time.sleep(0.1)
            sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
            df = pd.read_csv(sgRNA_summary)
            return df.to_dict('records')

        except:
            return [{}]

    if contents is not None:
        df = parse_contents(contents, filename)
        if df is not None:
            if all(x in df.columns for x in ['Chr', 'Start', 'Stop', 'Perturbation_Index', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type', 'Log2FC_Replicate1']):
                return df.to_dict('records')

            else:
                try:
                    time.sleep(0.1)
                    sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
                    df = pd.read_csv(sgRNA_summary)
                    return df.to_dict('records')

                except:
                    return [{}]

        else:
            try:
                time.sleep(0.1)
                sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
                df = pd.read_csv(sgRNA_summary)
                return df.to_dict('records')

            except:
                return [{}]
    else:

        try:
            time.sleep(0.1)
            sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
            df = pd.read_csv(sgRNA_summary)
            return df.to_dict('records')

        except:
            return [{}]

@app.callback(Output('sgRNA-table-view', 'columns'),
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename'),
               Input('example-data-button', 'n_clicks')],
               state = [State('url', 'pathname')])
def update_output(contents, filename, n_clicks, pathname):

    # if pathname:

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if n_clicks > param_dict['example_data_clicks']:

        try:
            time.sleep(0.1)
            sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
            df = pd.read_csv(sgRNA_summary)
            return df.columns
            
        except:
            return []

    if contents is not None:
        df = parse_contents(contents, filename)
        if df is not None:
            if all(x in df.columns for x in ['Chr', 'Start', 'Stop', 'Perturbation_Index', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type', 'Log2FC_Replicate1']):
                return df.columns

            else:
                try:
                    time.sleep(0.1)
                    sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
                    df = pd.read_csv(sgRNA_summary)
                    return df.columns
                    
                except:
                    return []

        else:
            try:
                time.sleep(0.1)
                sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
                df = pd.read_csv(sgRNA_summary)
                return df.columns
                
            except:
                return []
    else:

        try:
            time.sleep(0.1)
            sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
            df = pd.read_csv(sgRNA_summary)
            return df.columns
            
        except:
            return []


# CALLBACKS FOR PARAMETERS
@app.callback(
    Output('advanced-parameters-button', 'children'),
    [Input('advanced-parameters-button', 'n_clicks')])

def update_advanced_parameter_button_children(n_clicks):
    if n_clicks%2 == 0:
        return 'Show (+)'
    else:
        return 'Hide (-)'

@app.callback(
    Output('advanced-parameters-container', 'style'),
    [Input('advanced-parameters-button', 'n_clicks')])

def update_advanced_parameter_container(n_clicks):

    if n_clicks%2 == 0:
        return {'border-style':'solid','border-width':'2px', 'border-color':'#DCDCDC','border-radius':'10px','border-spacing':'15px','padding':'10px', 'display':'none'}
    else:
        return {'border-style':'solid','border-width':'2px', 'border-color':'#DCDCDC','border-radius':'10px','border-spacing':'15px','padding':'10px', 'display':'block'}

@app.callback(
    Output('parameters-help-button', 'children'),
    [Input('parameters-help-button', 'n_clicks')])

def update_help_parameter_button_children(n_clicks):
    if n_clicks%2 == 0:
        return 'Show (+)'
    else:
        return 'Hide (-)'

@app.callback(
    Output('parameters-help-container', 'style'),
    [Input('parameters-help-button', 'n_clicks')])

def update_help_parameter_container(n_clicks):
    if n_clicks%2 == 0:
        return {'border-style':'solid','border-width':'2px', 'border-color':'#DCDCDC','border-radius':'10px','border-spacing':'15px','padding':'10px', 'display':'none'}
    else:
        return {'border-style':'solid','border-width':'2px', 'border-color':'#DCDCDC','border-radius':'10px','border-spacing':'15px','padding':'10px', 'display':'block'}

@app.callback(Output('range', 'value'),
              [Input('pert', 'value'),
              Input('tabs', 'value')],
              state = [State('range', 'value'),
              State('url', 'pathname')])

def update_perturbation_range_label(perturbation_type, tabs, original_range, pathname):

    if tabs == 'parameters':

        if perturbation_type == 'cas9' or perturbation_type == 'cpf1':
            return 7

        elif perturbation_type == 'crispri' or perturbation_type == 'crispra':
            return 200

        elif perturbation_type == 'be':
            return 5

    else:
        return original_range


@app.callback(Output('perturbation-range-show', 'children'),
              [Input('range', 'value')],
              state = [State('url', 'pathname')])

def update_perturbation_range_label(perturbation_range, pathname):

    return 'Characteristic Perturbation Length: %s bps' % int(perturbation_range)

@app.callback(Output('perturbation-profile', 'figure'),
              [Input('range', 'value')])

def update_perturbation_profile_figure(perturbation_range):

    try:
        scale = 1
        limit =  perturbation_range*3
        perturbation_profile = gaussian_pattern(characteristic_perturbation_range = perturbation_range, scale = scale, limit = limit)

        trace = [go.Scatter(
            x= range(-limit, limit + 1),
            y = perturbation_profile,
            fill='tozeroy'
            )]

    except:
        trace = []

    return {
        'data': trace,
        'layout': go.Layout(
            height = 250,
            autosize = True,
            margin=dict(l=50,r=50,t=0),
            hovermode='closest',
            xaxis = dict(showgrid = False, zeroline=False, showline=True, title = 'Base-Pairs (bps)', range = [-800,800]),
            yaxis = dict(showgrid = False, zeroline=False, title = ''),
        )
    }

@app.callback(Output('scale', 'value'),
              [Input('range', 'value')],
              state = [State('tabs', 'value'),
              State('scale', 'value'),
              State('url', 'pathname')])

def update_perturbation_range_label(perturbation_range, tab, original_scale, pathname):

    if tab == 'parameters':
        if perturbation_range >= 10:
            return int(float(perturbation_range)/10.0)
        else:
            return 1
    else:
        return original_scale

@app.callback(Output('scale-show', 'children'),
              [Input('scale', 'value')],
              state = [State('url', 'pathname')])

def update_perturbation_range_label(scale, pathname):

    return 'Scale: %s' % int(scale)

@app.callback(Output('limit', 'value'),
              [Input('pert', 'value')],
              state = [State('url', 'pathname')])

def update_perturbation_range_label(perturbation_type, pathname):

    if perturbation_type == 'cas9' or perturbation_type == 'cpf1':
        return 25

    elif perturbation_type == 'crispri' or perturbation_type == 'crispra':
        return 300

    elif perturbation_type == 'be':
        return 10

@app.callback(Output('limit-show', 'children'),
              [Input('limit', 'value')],
              state = [State('url', 'pathname')])

def update_perturbation_range_label(limit, pathname):

    return 'Limit: %s' % int(limit)

@app.callback(Output('simulation-n-show', 'children'),
              [Input('sim_n', 'value')],
              state = [State('url', 'pathname')])

def update_simulation_n(simulation_n, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    data_dict['sim_n'] = simulation_n

    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
        new_json_string = json.dumps(data_dict)
        f.write(new_json_string + '\n')

    return 'Number of Simulations: %s' % int(simulation_n)

@app.callback(Output('tmp4', 'style'),
              [Input('sim_type', 'value')],
              state = [State('url', 'pathname')])

def update_simulation_n(simulation_type, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    data_dict['sim_type'] = simulation_type

    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
        new_json_string = json.dumps(data_dict)
        f.write(new_json_string + '\n')

    return {'display':'none'}

@app.callback(Output('tmp5', 'style'),
              [Input('effect_size', 'value')],
              state = [State('url', 'pathname')])

def update_effect_size(effect_size, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    data_dict['effect_size'] = effect_size

    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
        new_json_string = json.dumps(data_dict)
        f.write(new_json_string + '\n')

    return {'display':'none'}

### CALLBACKS FOR DECONVOLUTION ###
@app.callback(
    Output('tmp2', 'style'),
    [Input('deconvolution-button', 'n_clicks')],
    state=[State('range', 'value'),
    State('scale', 'value'),
    State('limit', 'value'),
    State('gamma_list', 'value'),
    State('avg', 'value'),
    State('sim_type', 'value'),
    State('sim_n', 'value'),
    State('effect_size', 'value'),
    State('rapid', 'value'),
    State('genome', 'value'),
    State('url', 'pathname')])

def perform_deconvolution(n_clicks, range_val, scale_val, limit_val, gamma_list, avg, sim_type, sim_n, effect_size, rapid, genome, pathname):

    if n_clicks > 0:

        UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
        RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

        # hack
        json_good = False
        while not json_good:
            with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
                json_string = f.readline().strip()
                try:
                    data_dict = json.loads(json_string)
                    json_good = True
                except:
                    pass

        perturbation_profile = gaussian_pattern(characteristic_perturbation_range = range_val, scale = scale_val, limit = limit_val)

        if gamma_list == '':

            if range_val <= 50:
                gamma_list = [0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
            else:
                gamma_list = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]

        else:
            try:
                gamma_list = [float(x) for x in gamma_list.split(',')]

            except:
                if range_val <= 50:
                    gamma_list = [0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
                else:
                    gamma_list = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]


        # Grab dataframe to set up deconvolution
        sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
        df = pd.read_csv(sgRNA_summary)
        replicates = len([x for x in df.columns.tolist() if 'Log2FC_Replicate' in x])

        # Set up inputs
        sgRNA_indices = [int(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Perturbation_Index']).flatten().tolist()]
        chromosomes = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()]

        observations = []
        for i in range(1, replicates + 1):
            observations.append([float(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Log2FC_Replicate' + str(i)]).flatten().tolist()])

        data_dict['perturbation_profile'] = perturbation_profile
        data_dict['sgRNA_indices'] = sgRNA_indices
        data_dict['chromosomes'] = chromosomes
        data_dict['observations'] = observations
        data_dict['gamma_list'] = gamma_list
        data_dict['scale'] = scale_val
        data_dict['limit'] = limit_val
        data_dict['avg'] = avg
        data_dict['sim_type'] = sim_type
        data_dict['sim_n'] = sim_n
        data_dict['effect_size'] = effect_size
        data_dict['rapid'] = rapid
        data_dict['genome'] = genome

        print 'SCALE VALUE %s' % scale_val

        with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
            new_json_string = json.dumps(data_dict)
            f.write(new_json_string + '\n')

        sb.Popen('python /SURF/SURF_deconvolution_webapp.py -uploads_dir %s -results_dir %s' % (UPLOADS_FOLDER, RESULTS_FOLDER), shell = True)

    return {'display': 'none'}

@app.callback(
    Output('chr', 'options'),
    [Input('deconvolution-plot', 'figure')],
    state = [State('url', 'pathname')])

def update_chr(figure, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    unique_chr = set(data_dict['chromosomes'])

    return [{'label':entry, 'value':entry} for entry in unique_chr]

@app.callback(
    Output('custom-loading-states-1', 'style'),
    [Input('deconvolution-button', 'n_clicks'),
    Input('deconvolution-container', 'style')],
    state = [State('url', 'pathname')])

def update_container(n_clicks, deconvolution_container, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'CUSTOM LOADING STATE'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # First time analysis
    if n_clicks == param_dict['checkbutton1']:
        print 'CUSTOM LOADING STATE PASS'
        return {'display': 'block'}

    else:
        print 'CUSTOM LOADING FAIL'
        return {'display': 'none'}

@app.callback(
    Output('deconvolution-container', 'style'),
    [Input('deconvolution-plot', 'figure')],
    state = [State('url', 'pathname')])

def smoothing_container(fig_update, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'DECONVOLUTION CONTAINER'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if param_dict['deconvolution-clicks'] > 0:
        print 'DECONVOLUTION CONTAINER PASS'
        return {'display': 'block'}

    else:
        print 'DECONVOLUTION CONTAINER FAIL'
        return {'display': 'none'}

# @app.callback(
#     Output('deconvolution-button', 'style'),
#     [Input('deconvolution-plot', 'figure')],
#     state = [State('deconvolution-button', 'n_clicks'),
#     State('url', 'pathname')])

# def button_container(fig_update, n_clicks, pathname):

#     UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
#     RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

#     if n_clicks > 0:
#         return {'display': 'none'}

#     else:
#         return {'display': 'block'}

@app.callback(
    Output('common-interval-1', 'interval'),
    [Input('deconvolution-button', 'n_clicks'),
    Input('deconvolution-container', 'style')],
    state = [State('url', 'pathname')])

def update_container(n_clicks, deconvolution_container, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'COMMON INTERVAL 1'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # First time analysis
    if n_clicks == param_dict['checkbutton1']:
        print 'COMMON INTERVAL 1 PASS'
        return 10000

    else:
        print 'COMMON INTERVAL 1 FAIL'
        return 2000000000

@app.callback(Output('replicate-correlation-show', 'children'),
              [Input('replicate-correlation', 'value'),
              Input('deconvolution-plot', 'figure')],
              state = [State('url', 'pathname')])

def update_perturbation_range_label(replicate_correlation, figure, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    gamma_range = data_dict['gamma_range']

    gamma_list = [float(x[0]) for x in gamma_range]
    correlation_list = [float(x[1]) for x in gamma_range]

    gamma_index_use = min(range(len(correlation_list)), key=lambda i: abs(correlation_list[i] - float(replicate_correlation)))
    gamma_use = gamma_list[gamma_index_use]

    data_dict['gamma_use'] = gamma_use

    with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
        new_json_string = json.dumps(data_dict)
        f.write(new_json_string + '\n')

    return 'Replicate Correlation: %s (Gamma: %s)' % (replicate_correlation, gamma_use)

@app.callback(Output('chr', 'value'),
              [Input('sgRNA-table-view', 'rows')],
              state = [State('url', 'pathname')])

def initialize_chr(table, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    return data_dict['chr']

@app.callback(
    Output('deconvolution-plot', 'figure'),
    [Input('update-graph-button1', 'n_clicks'),
    Input('start','type')],
    state=[State('avg', 'value'),
    State('scale', 'value'),
    State('chr', 'value'),
    State('start', 'value'),
    State('stop', 'value'),
    State('replicate-correlation', 'value'),
    State('deconvolution-button', 'n_clicks'),
    State('custom-loading-states-1', 'style'),
    State('url', 'pathname')],
    events=[Event('common-interval-1', 'interval')])

def update_deconvolution_plot(update_graph_clicks, chrom_opt, avg, scale_val, chrom, start, stop, correlation_use, deconvolution_n_clicks, loading_style, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'DECONVOLUTION PLOT'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if os.path.exists(RESULTS_FOLDER + '/flag1.txt'):

        print 'DECONVOLUTION PLOT PASS'

        print 'PASS 1 -.-.--.--.---.-.-.-.-.-.-.---.-.-.-.-.-.-.-.-.-.'

        sb.call('rm %s/flag1.txt' % RESULTS_FOLDER, shell = True)

        fig = tools.make_subplots(rows=2, cols=1, specs=[[{}], [{}]],
                                  shared_xaxes=True, shared_yaxes=True,
                                  vertical_spacing=0.1)

        sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
        df = pd.read_csv(sgRNA_summary)
        replicates = sorted([x for x in df.columns.tolist() if 'Log2FC_Replicate' in x])

        chrom = data_dict['chr']
        df = df.loc[df['Chr'] == chrom]

        df = df.loc[df['sgRNA_Type'] != 'negative_control']
        df = df.dropna(axis=0, how='any')

        for replicate in replicates:

            fig.append_trace(go.Scatter(
                x=df['Perturbation_Index'],
                y=df[replicate],
                text=replicate,
                mode='markers',
                opacity=0.5,
                marker={
                    'size': 6,
                    'line': {'width': 0.3, 'color': 'white'}
                },
                name= replicate,
                # yaxis = 'y2'
            ), 1, 1)

        gammas2betas = data_dict['gammas2betas']
        gamma_range = data_dict['gamma_range']

        gamma_list = [float(x[0]) for x in gamma_range]
        correlation_list = [float(x[1]) for x in gamma_range]

        gamma_index_use = min(range(len(correlation_list)), key=lambda i: abs(correlation_list[i] - float(correlation_use)))
        gamma_use = gamma_list[gamma_index_use]

        deconvolved_signal = {}
        for i in [x for x in gammas2betas.keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr') and (x != 'p'))]:

            for j in range(len(gammas2betas[i][str(gamma_use)])):

                if j not in deconvolved_signal:
                    deconvolved_signal[j] = []

                deconvolved_signal[j].append(gammas2betas[i][str(gamma_use)][j])

        # Create mean or median profile
        if avg == 'mean':
            y_p = [np.mean(deconvolved_signal[x]) for x in deconvolved_signal]

        elif avg == 'median':
            y_p = [np.median(deconvolved_signal[x]) for x in deconvolved_signal]

        # Filter for chrom
        indices_filt = []
        y_p_filt = []
        for chrom_i, coord_i, value_i in zip(gammas2betas['chr'], gammas2betas['indices'], y_p):
            if chrom_i == chrom:
                indices_filt.append(coord_i)
                y_p_filt.append(value_i)

        # Boundaries of inference
        diff_vec = [0] + list(np.diff(indices_filt))
        mult_vec = [0] + [float(y_p_filt[i])*float(y_p_filt[i+1]) for i in range(len(y_p_filt[:-1]))]

        boundaries = [0]
        for index, value in enumerate(diff_vec):
            try:
                if int(value) > int(scale_val):
                    boundaries.append(index)
            except:
                pass

        boundaries = sorted(boundaries)

        boundaries.append(len(indices_filt))

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1], y=[0] + y_p_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                yaxis = 'y2',
                showlegend=False,
                line=dict(color='rgb(169,169,169)')), 2, 1)

        fig['layout'].update(
            height=600,
            autosize = True,
            xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[min(indices_filt)-1000, max(indices_filt)+1000]},
            yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Log2FC Scores', 'autorange':True},
            yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution'},
            hovermode='closest',
            title='Raw Scores and Deconvolution')

        param_dict['deconvolution-clicks'] += 1
        param_dict['checkbutton1'] = deconvolution_n_clicks + 1
        data_dict['computed'] = True

        with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
            new_json_string = json.dumps(param_dict)
            f.write(new_json_string + '\n')

        with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
            new_json_string = json.dumps(data_dict)
            f.write(new_json_string + '\n')

        return fig

    elif (update_graph_clicks > 0) and (loading_style == {'display': 'none'}):

        print 'PASS 2 -.-.--.--.---.-.-.-.-.-.-.---.-.-.-.-.-.-.-.-.-.'

        fig = tools.make_subplots(rows=2, cols=1, specs=[[{}], [{}]],
                                  shared_xaxes=True, shared_yaxes=True,
                                  vertical_spacing=0.1)

        sgRNA_summary = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'
        df = pd.read_csv(sgRNA_summary)
        replicates = sorted([x for x in df.columns.tolist() if 'Log2FC_Replicate' in x])

        # chrom = data_dict['chr']
        df = df.loc[df['Chr'] == chrom]

        df = df.loc[df['sgRNA_Type'] != 'negative_control']
        df = df.dropna(axis=0, how='any')

        for replicate in replicates:

            fig.append_trace(go.Scatter(
                x=df['Perturbation_Index'],
                y=df[replicate],
                text=replicate,
                mode='markers',
                opacity=0.5,
                marker={
                    'size': 6,
                    'line': {'width': 0.3, 'color': 'white'}
                },
                name= replicate,
                # yaxis = 'y1'
            ), 1, 1)

        gammas2betas = data_dict['gammas2betas']
        gamma_range = data_dict['gamma_range']

        gamma_list = [float(x[0]) for x in gamma_range]
        correlation_list = [float(x[1]) for x in gamma_range]

        gamma_index_use = min(range(len(correlation_list)), key=lambda i: abs(correlation_list[i] - float(correlation_use)))
        gamma_use = gamma_list[gamma_index_use]

        deconvolved_signal = {}
        for i in [x for x in gammas2betas.keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr') and (x != 'p'))]:

            for j in range(len(gammas2betas[i][str(gamma_use)])):

                if j not in deconvolved_signal:
                    deconvolved_signal[j] = []

                deconvolved_signal[j].append(gammas2betas[i][str(gamma_use)][j])

        # Create mean or median profile
        if avg == 'mean':
            y_p = [np.mean(deconvolved_signal[x]) for x in deconvolved_signal]

        elif avg == 'median':
            y_p = [np.median(deconvolved_signal[x]) for x in deconvolved_signal]

        # Filter for chrom
        indices_filt = []
        y_p_filt = []
        for chrom_i, coord_i, value_i in zip(gammas2betas['chr'], gammas2betas['indices'], y_p):
            if chrom_i == chrom:
                indices_filt.append(coord_i)
                y_p_filt.append(value_i)

        # Boundaries of inference
        diff_vec = [0] + list(np.diff(indices_filt))
        mult_vec = [0] + [float(y_p_filt[i])*float(y_p_filt[i+1]) for i in range(len(y_p_filt[:-1]))]

        boundaries = [0]
        for index, value in enumerate(diff_vec):
            try:
                if int(value) > int(scale_val):
                    boundaries.append(index)
            except:
                pass

        boundaries = sorted(boundaries)

        boundaries.append(len(indices_filt))

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1], y=[0] + y_p_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                yaxis = 'y2',
                showlegend=False,
                line=dict(color='rgb(169,169,169)')), 2, 1)

        try:
            start = int(start)
            stop = int(stop)
            fig['layout'].update(
                height=600,
                autosize = True,
                xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[start, stop]},
                yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Log2FC Scores', 'autorange':True},
                yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
                hovermode='closest',
                title='Raw Scores and Deconvolution')

        except:
            fig['layout'].update(
                height=600,
                autosize = True,
                xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False},
                yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Log2FC Scores', 'autorange':True},
                yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution'},
                hovermode='closest',
                title='Raw Scores and Deconvolution')

        param_dict['deconvolution-clicks'] += 1
        param_dict['checkbutton1'] = deconvolution_n_clicks + 1
        data_dict['computed'] = True

        with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
            new_json_string = json.dumps(param_dict)
            f.write(new_json_string + '\n')

        with open(UPLOADS_FOLDER + '/data.json', 'w') as f:
            new_json_string = json.dumps(data_dict)
            f.write(new_json_string + '\n')

        return fig

    print 'DECONVOLUTION PLOT FAIL'

### CALLBACKS FOR SIGNIFICANCE ###
@app.callback(
    Output('time-estimate', 'children'),
    [Input('deconvolution-plot', 'figure'),
    Input('sim_n', 'value')],
    state = [State('url', 'pathname')])

def update_time_estimate(figure, sim_n, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    time_per_gamma = data_dict['time_per_gamma']

    estimate_simulation_time = 2.0*float(time_per_gamma)*float(sim_n)

    if estimate_simulation_time >= 1:
        return 'Estimated Time: %s Minutes' % "{0:.2f}".format(estimate_simulation_time)

    else:
        return 'Estimated Time: %s Seconds' % "{0:.2f}".format(estimate_simulation_time*60.0)


@app.callback(
    Output('chr2', 'options'),
    [Input('significance-plot', 'figure')],
    state = [State('url', 'pathname')])

def update_chr(figure, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    unique_chr = set(data_dict['chromosomes'])

    return [{'label':entry, 'value':entry} for entry in unique_chr]

@app.callback(
    Output('tmp3', 'style'),
    [Input('significance-button', 'n_clicks')],
    state=[State('url', 'pathname')])

def find_regions(n_clicks, pathname): #n_clicks, pert_range, scale, limit, gamma_list, pathname):

    if n_clicks > 0:

        UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
        RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

        sb.Popen('python /SURF/SURF_significance_webapp.py -uploads_dir %s -results_dir %s' % (UPLOADS_FOLDER, RESULTS_FOLDER), shell = True)

    return {'display': 'none'}

@app.callback(
    Output('custom-loading-states-2', 'style'),
    [Input('significance-button', 'n_clicks'),
    Input('significance-container', 'style')],
    state = [State('url', 'pathname')])

def update_container(n_clicks, significance_container, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'CUSTOM LOADING STATE 2'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # First time analysis
    if n_clicks == param_dict['checkbutton2']:
        print 'CUSTOM LOADING STATE 2 PASS'
        return {'display': 'block'}

    else:
        print 'CUSTOM LOADING STATE 2 FAIL'
        return {'display': 'none'}

@app.callback(
    Output('significance-container', 'style'),
    [Input('significance-plot', 'figure')],
    state = [State('url', 'pathname')])

def significance_container(fig_update, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'SIGNIFICANCE CONTAINER'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    if param_dict['significance-clicks'] > 0:
        print 'SIGNIFICANCE CONTAINER PASS'
        return {'display': 'block'}

    else:
        print 'SIGNIFICANCE CONTAINER FAIL'
        return {'display': 'none'}

@app.callback(
    Output('common-interval-2', 'interval'),
    [Input('significance-button', 'n_clicks'),
    Input('significance-container', 'style')],
    state = [State('url', 'pathname')])

def update_container(n_clicks, significance_container, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'COMMON INTERVAL 2'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # First time analysis
    if n_clicks == param_dict['checkbutton2']:
        print 'COMMON INTERVAL 2 PASS'
        return 10000

    else:
        print 'COMMON INTERVAL 2 FAIL'
        return 2000000000

@app.callback(Output('fdr-show', 'children'),
              [Input('fdr', 'value'),
              Input('significance-plot', 'figure')],
              state = [State('url', 'pathname')])

def update_fdr_label(fdr, figure, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    data_dict2['fdr'] = fdr

    with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
            new_json_string = json.dumps(data_dict2)
            f.write(new_json_string + '\n')

    return 'FDR < %s' % fdr

@app.callback(Output('chr2', 'value'),
              [Input('sgRNA-table-view', 'rows')],
              state = [State('url', 'pathname')])

def initialize_chr(table, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    return data_dict2['chr']

@app.callback(
    Output('significance-plot', 'figure'),
    [Input('update-graph-button2', 'n_clicks'),
    Input('start2','type')],
    state=[State('scale', 'value'),
    State('chr2', 'value'),
    State('start2', 'value'),
    State('stop2', 'value'),
    State('fdr', 'value'),
    State('significance-button', 'n_clicks'),
    State('custom-loading-states-2', 'style'),
    State('url', 'pathname')],
    events=[Event('common-interval-2', 'interval')])

def update_significance_plot(update_graph_clicks, chrom_opt, scale_val, chrom, start, stop, fdr, significance_n_clicks, loading_style, pathname):

    UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
    RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

    print 'SIGNIFICANCE PLOT'
    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/params.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                param_dict = json.loads(json_string)
                json_good = True
            except:
                pass

    # hack
    json_good = False
    while not json_good:
        with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
            json_string = f.readline().strip()
            try:
                data_dict2 = json.loads(json_string)
                json_good = True
            except:
                pass

    if os.path.exists(RESULTS_FOLDER + '/flag2.txt'):

        print 'SIGNIFICANCE PLOT PASS'

        print 'PASS 1 -.-.--.--.---.-.-.-.-.-.-.---.-.-.-.-.-.-.-.-.-.'

        sb.call('rm %s/flag2.txt' % RESULTS_FOLDER, shell = True)

        fig = tools.make_subplots(rows=2, cols=1, specs=[[{}], [{}]],
                                  shared_xaxes=True, shared_yaxes=True,
                                  vertical_spacing=0.1)

        gammas2betas = data_dict2['gammas2betas']

        # Filter for chrom
        indices_filt = []
        y_p_filt = []
        power_filt = []

        for chrom_i, coord_i, value_i, power_i in zip(gammas2betas['chr'], gammas2betas['indices'], gammas2betas['combined'], gammas2betas['power']):
            if chrom_i == data_dict2['chr']:
                indices_filt.append(int(coord_i))
                y_p_filt.append(float(value_i))
                power_filt.append(float(power_i))

        # Boundaries of inference
        diff_vec = [0] + list(np.diff(indices_filt))

        boundaries = [0]
        for index, value in enumerate(diff_vec):
            try:
                if int(value) > int(scale_val):
                    boundaries.append(index)
            except:
                pass

        boundaries = sorted(boundaries)

        boundaries.append(len(indices_filt))

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1],
                y=[0] + power_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                showlegend=False,
                # yaxis = 'y2',
                line=dict(color='rgb(255,165,0)')), 1, 1)

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1],
                y=[0] + y_p_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                showlegend=False,
                yaxis = 'y2',
                line=dict(color='rgb(169,169,169)')), 2, 1)

        # Find indices with significant padj-values
        significant_boundary_indices = []
        padj_list = [float(x) for x in gammas2betas['padj']]
        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]
            significant_indices = [1 if x < fdr else 0 for x in padj_list[start_index:stop_index]]
            significant_boundary_indices_tmp = [j + start_index for j, k in enumerate(np.diff([0] + significant_indices)) if k != 0]

            if len(significant_boundary_indices_tmp)%2 == 1:
                significant_boundary_indices_tmp.append(stop_index - 1)

            significant_boundary_indices += significant_boundary_indices_tmp

        for i, j in zip(significant_boundary_indices[0::2], significant_boundary_indices[1::2]):

            boundary_start = int(i)
            boundary_stop = int(j)

            genomic_boundary_start = int(gammas2betas['indices'][boundary_start])
            genomic_boundary_stop = int(gammas2betas['indices'][boundary_stop])

            signal_mean = np.mean(gammas2betas['combined'][boundary_start:boundary_stop])

            if signal_mean > np.median(gammas2betas['combined']):
                fig.append_trace(go.Scatter(
                    x=[genomic_boundary_start, genomic_boundary_stop],
                    y=[max(y_p_filt) + 0.01, max(y_p_filt) + 0.01],
                    mode = 'lines',
                    # fill='tozeroy',
                    showlegend=False,
                    yaxis = 'y2',
                    line=dict(color='rgb(255,0,0)', width = 5)), 2, 1)
            else:
                fig.append_trace(go.Scatter(
                    x=[genomic_boundary_start, genomic_boundary_stop],
                    y=[min(y_p_filt) - 0.01, min(y_p_filt) - 0.01],
                    mode = 'lines',
                    # fill='tozeroy',
                    showlegend=False,
                    yaxis = 'y2',
                    line=dict(color='rgb(30,144,255)', width = 5)), 2, 1)

        fig['layout'].update(
            height=600,
            autosize = True,
            xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[min(indices_filt)-1000, max(indices_filt)+1000]},
            yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Statistical Power', 'autorange':True},
            yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
            hovermode='closest',
            title='Finding Regions')

        param_dict['significance-clicks'] += 1
        param_dict['checkbutton2'] = significance_n_clicks + 1
        data_dict2['computed'] = True

        with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
            new_json_string = json.dumps(param_dict)
            f.write(new_json_string + '\n')

        with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
            new_json_string = json.dumps(data_dict2)
            f.write(new_json_string + '\n')

        return fig

    elif (update_graph_clicks > 0) and (loading_style == {'display': 'none'}):

        print 'PASS 2 -.-.--.--.---.-.-.-.-.-.-.---.-.-.-.-.-.-.-.-.-.'

        fig = tools.make_subplots(rows=2, cols=1, specs=[[{}], [{}]],
                                  shared_xaxes=True, shared_yaxes=True,
                                  vertical_spacing=0.1)

        gammas2betas = data_dict2['gammas2betas']

        # Filter for chrom
        indices_filt = []
        y_p_filt = []
        power_filt = []

        for chrom_i, coord_i, value_i, power_i in zip(gammas2betas['chr'], gammas2betas['indices'], gammas2betas['combined'], gammas2betas['power']):
            if chrom_i == chrom:
                indices_filt.append(int(coord_i))
                y_p_filt.append(float(value_i))
                power_filt.append(float(power_i))

        # Boundaries of inference
        diff_vec = [0] + list(np.diff(indices_filt))

        boundaries = [0]
        for index, value in enumerate(diff_vec):
            try:
                if int(value) > int(scale_val):
                    boundaries.append(index)
            except:
                pass

        boundaries = sorted(boundaries)

        boundaries.append(len(indices_filt))

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1],
                y=[0] + power_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                showlegend=False,
                # yaxis = 'y2',
                line=dict(color='rgb(255,165,0)')), 1, 1)

        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]

            fig.append_trace(go.Scatter(
                x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1],
                y=[0] + y_p_filt[start_index:stop_index] + [0],
                mode = 'lines',
                fill='tozeroy',
                showlegend=False,
                yaxis = 'y2',
                line=dict(color='rgb(169,169,169)')), 2, 1)

        # Find indices with significant padj-values
        significant_boundary_indices = []
        padj_list = [float(x) for x in gammas2betas['padj']]
        for i in range(len(boundaries) - 1):
            start_index, stop_index = boundaries[i], boundaries[i + 1]
            significant_indices = [1 if x < fdr else 0 for x in padj_list[start_index:stop_index]]
            significant_boundary_indices_tmp = [j + start_index for j, k in enumerate(np.diff([0] + significant_indices)) if k != 0]

            if len(significant_boundary_indices_tmp)%2 == 1:
                significant_boundary_indices_tmp.append(stop_index - 1)

            significant_boundary_indices += significant_boundary_indices_tmp

        for i, j in zip(significant_boundary_indices[0::2], significant_boundary_indices[1::2]):

            boundary_start = int(i)
            boundary_stop = int(j)

            genomic_boundary_start = int(gammas2betas['indices'][boundary_start])
            genomic_boundary_stop = int(gammas2betas['indices'][boundary_stop])

            signal_mean = np.mean(gammas2betas['combined'][boundary_start:boundary_stop])

            if signal_mean > np.median(gammas2betas['combined']):
                fig.append_trace(go.Scatter(
                    x=[genomic_boundary_start, genomic_boundary_stop],
                    y=[max(y_p_filt) + 0.01, max(y_p_filt) + 0.01],
                    mode = 'lines',
                    # fill='tozeroy',
                    showlegend=False,
                    yaxis = 'y2',
                    line=dict(color='rgb(255,0,0)', width = 5)), 2, 1)
            else:
                fig.append_trace(go.Scatter(
                    x=[genomic_boundary_start, genomic_boundary_stop],
                    y=[min(y_p_filt) - 0.01, min(y_p_filt) - 0.01],
                    mode = 'lines',
                    # fill='tozeroy',
                    showlegend=False,
                    yaxis = 'y2',
                    line=dict(color='rgb(30,144,255)', width = 5)), 2, 1)

        try:
            start = int(start)
            stop = int(stop)
            fig['layout'].update(
                height=600,
                autosize = True,
                xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[start, stop]},
                yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Statistical Power', 'autorange':True},
                yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
                hovermode='closest',
                title='Finding Regions')

        except:
            fig['layout'].update(
                height=600,
                autosize = True,
                xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[min(indices_filt)-1000, max(indices_filt)+1000]},
                yaxis={'domain':[0, 0.5], 'showgrid': False, 'title': 'Statistical Power', 'autorange':True},
                yaxis2={'domain':[0.5, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
                hovermode='closest',
                title='Finding Regions')

        param_dict['significance-clicks'] += 1
        param_dict['checkbutton2'] = significance_n_clicks + 1
        data_dict2['computed'] = True

        with open(UPLOADS_FOLDER + '/params.json', 'w') as f:
            new_json_string = json.dumps(param_dict)
            f.write(new_json_string + '\n')

        with open(UPLOADS_FOLDER + '/data2.json', 'w') as f:
            new_json_string = json.dumps(data_dict2)
            f.write(new_json_string + '\n')

        return fig

    print 'SIGNIFICANCE PLOT FAIL'

### CALLBACKS FOR DOWNLOAD ###
@app.callback(
    Output('download-total', 'href'),
    [Input('download-total', 'n_clicks'),
    Input('url', 'pathname')],
    state = [State('title-input', 'value'),
    State('description-input', 'value')])

def zip_dir(n_clicks, pathname, title_input, description_input):

    if pathname:

        UPLOADS_FOLDER = app.server.config['UPLOADS_FOLDER'] + '/' + str(pathname).split('/')[-1]
        RESULTS_FOLDER = app.server.config['RESULTS_FOLDER'] + '/' + str(pathname).split('/')[-1]

        if n_clicks > 0:

            json_file = {'title': title_input,
                        'description': description_input}

            with open('%s/surf.json' % (RESULTS_FOLDER), 'w') as f:
                json_string = json.dumps(json_file)
                f.write(json_string + '\n')

            full_path = RESULTS_FOLDER + '/' + 'surf-outputs.zip'

            return '/dash/urldownload%s' % full_path

        else:
            full_path = RESULTS_FOLDER + '/' + 'surf-outputs.zip'
            return '/dash/urldownload%s' % full_path

@app.server.route('/dash/urldownload/tmp/RESULTS_FOLDER/<directory>/surf-outputs.zip')
def generate_report_url(directory):

    overview_folder = 'surf-outputs'
    results_folder = 'SURF_result'

    RESULTS_FOLDER = '/tmp/RESULTS_FOLDER/%s' % directory
    UPLOADS_FOLDER = '/tmp/UPLOADS_FOLDER/%s' % directory

    proc = sb.Popen('python /SURF/SURF_download_webapp.py -uploads_dir %s -results_dir %s' % (UPLOADS_FOLDER, RESULTS_FOLDER), shell = True)
    proc.wait()

    if os.path.exists('%s/%s' % (RESULTS_FOLDER, overview_folder)):
        sb.call('rm -r %s/%s' % (RESULTS_FOLDER, overview_folder), shell = True)

    sb.call('mkdir %s/%s' % (RESULTS_FOLDER, overview_folder), shell = True)
    sb.call('mkdir %s/%s/%s' % (RESULTS_FOLDER, overview_folder, results_folder), shell = True)

    sgRNA_summary_table = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv'

    sb.call('cp %s %s/%s/' % (sgRNA_summary_table, RESULTS_FOLDER, overview_folder), shell = True)

    sb.call('cp %s/surf.json %s/%s/' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder), shell = True)
    sb.call('cp -r %s/*bed %s/%s/%s' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder, results_folder), shell = True)
    sb.call('cp -r %s/*bedgraph %s/%s/%s' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder, results_folder), shell = True)
    sb.call('cp -r %s/*csv %s/%s/%s' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder, results_folder), shell = True)
    sb.call('cp -r %s/*xml %s/%s/%s' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder, results_folder), shell = True)
    # sb.call('cp -r %s/*log %s/%s/%s' % (RESULTS_FOLDER, RESULTS_FOLDER, overview_folder, results_folder), shell = True)

    proc = sb.Popen('pushd %s && zip -r %s.zip %s && popd' % (RESULTS_FOLDER, overview_folder, overview_folder), shell=True, executable='/bin/bash')
    proc.wait()

    return send_file('/tmp/RESULTS_FOLDER/%s/surf-outputs.zip' % (directory), attachment_filename = 'surf-outputs.zip', as_attachment = True)

### Precomputed Layout ###
app2.layout = html.Div([

    dcc.Location(id='url', refresh=False),

    html.Img(src='data:image/png;base64,{}'.format(crisprsurf_logo_image), width = '50%'),
    html.H2('CRISPR Screening Uncharacterized Region Function'),

    dcc.Tabs(
        tabs=[
            {'label': 'Visualize Datasets', 'value': 'dataset'},
            {'label': 'Download IGV Session', 'value': 'significance'}
        ],
        value='dataset',
        id='tabs',
        style = {'font-weight':'bold'}
    ),

    html.Div(id = 'dataset-container-total', children = [

        html.H3('Choose Dataset'),
        dcc.Dropdown(
            id = 'dataset',
            options=[
                {'label': 'Canver et al. 2015 BCL11A CRISPR-Cas9', 'value': 'Canver_2015'},
                {'label': 'Fulco et al. 2016 MYC CRISPRi', 'value': 'Fulco_2016'},
                {'label': 'Simeonov and Gowen et al. 2017 IL2RA CRISPRa', 'value': 'Simeonov_Gowen_2017'},
                {'label': 'Hsu et al. 2018 BCL11A CRISPRi', 'value': 'Hsu_CRISPRi_2018'},
                {'label': 'Hsu et al. 2018 BCL11A CRISPR-Cas9', 'value': 'Hsu_Cas9_2018'},
            ],
            value = 'Canver_2015',
        ),

        html.Hr(),

        html.Div(id = 'deconvolution-container', children = [

            html.Div([

                html.Div([

                    html.Div([

                        html.Div([

                            html.Label('Chr', style = {'font-weight':'bold'}),
                            dcc.Dropdown(
                                id = 'chr',
                                options=[],
                                value = None
                            ),

                            ], className = 'three columns'),

                        html.Div([

                            html.Label('Start', style = {'font-weight':'bold'}),
                            dcc.Input(
                                id = 'start',
                                type = 'text',
                                placeholder='Enter Start Coordinate ...',
                                ),

                            ], className = 'three columns', style = {'text-align':'center'}),

                        html.Div([

                            html.Label('Stop', style = {'font-weight':'bold'}),
                            dcc.Input(
                                id = 'stop',
                                type = 'text',
                                placeholder='Enter Stop Coordinate ...',
                                ),

                            ], className = 'three columns', style = {'text-align':'center'}),

                        html.Div([

                            html.Br(),
                            # html.Br(),

                            html.Button(id = 'update-graph-button1', children = 'Update Graph'),

                            ], className = 'three columns', style = {'text-align':'left'}),

                        ], className = 'row'),

                    ], className = 'eight columns'),

                ], className = 'row'),

            ]),

        dcc.Graph(id='deconvolution-plot', animate=False),

        ]),

    html.Div(id = 'significance-container-total', children = [

        html.H3('View Significant Regions'),

        html.Div(id = 'significance-container', children = [

            html.A(
                html.Button('Download IGV Session'),
                id='download-total',
                download = "igv-session.zip",
                href="",
                target="_blank",
                n_clicks = 0,
                style = {'font-weight':'bold', 'font-size':'100%', 'text-align':'center'}
                ),

            dt.DataTable(id='regions-table-view', filterable=True, sortable=True, column_widths = [250]*100000, rows=[{}], columns = ['FDR','Chr','Start','Stop','Direction','Deconvolved_Signal_Area','Deconvolved_Signal_Mean','Padj_Mean','Supporting_sgRNAs']),

            ]),

        ], style = {'display':'none'}),

    html.Br(),
    html.Br(),

    html.Div(id = 'tmp11', style = {'display':'none'}),
    html.Div(id = 'tmp12', style = {'display':'none'}),
    html.Div(id = 'tmp13', style = {'display':'none'}),
    html.Div(id = 'tmp14', style = {'display':'none'}),
    html.Div(id = 'tmp15', style = {'display':'none'}),
    html.Div(id = 'tmp16', style = {'display':'none'}),
    html.Div(id = 'tmp17', style = {'display':'none'}),
    html.Div(id = 'tmp18', style = {'display':'none'}),
    html.Div(id = 'tmp19', style = {'display':'none'}),
    html.Div(id = 'tmp20', style = {'display':'none'}),

    ])

# Pseudo tabs code
@app2.callback(
    Output('dataset-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    if value == 'dataset':
        return {'display':'block'}
    else:
        return {'display':'none'}

@app2.callback(
    Output('significance-container-total', 'style'),
    [Input('tabs', 'value')],
    state = [State('url', 'pathname')])

def display_content(value, pathname):

    if value == 'significance':
        return {'display':'block'}
    else:
        return {'display':'none'}

### CALLBACKS FOR DECONVOLUTION ###
@app2.callback(
    Output('chr', 'options'),
    [Input('dataset', 'value')])

def update_chr(dataset):

    sgRNA_summary = '/SURF/precomputed/%s/SURF_result/sgRNAs_summary_table_updated.csv' % (dataset)
    df = pd.read_csv(sgRNA_summary)

    df = df.loc[df['sgRNA_Type'] != 'negative_control']
    df = df.dropna(axis=0, how='any')

    unique_chr = [x for x in df['Chr'].unique()]

    return [{'label':entry, 'value':entry} for entry in unique_chr]

@app2.callback(Output('chr', 'value'),
              [Input('dataset', 'value')])

def initialize_chr(dataset):

    sgRNA_summary = '/SURF/precomputed/%s/SURF_result/sgRNAs_summary_table_updated.csv' % (dataset)
    df = pd.read_csv(sgRNA_summary)

    df = df.loc[df['sgRNA_Type'] != 'negative_control']
    df = df.dropna(axis=0, how='any')

    unique_chr = [x for x in df['Chr'].unique()]

    return unique_chr[0]

@app2.callback(
    Output('deconvolution-plot', 'figure'),
    [Input('dataset', 'value'),
    Input('update-graph-button1', 'n_clicks')],
    state = [
    State('tabs', 'value'),
    State('chr', 'value'),
    State('start', 'value'),
    State('stop', 'value')])

def update_deconvolution_plot(dataset, update_graph_clicks, tab, chrom, start, stop):

    fig = tools.make_subplots(rows=3, cols=1, specs=[[{}], [{}], [{}]],
                              shared_xaxes=True, shared_yaxes=True,
                              vertical_spacing=0.1)

    sgRNA_summary = '/SURF/precomputed/%s/SURF_result/sgRNAs_summary_table_updated.csv' % (dataset)
    df = pd.read_csv(sgRNA_summary)
    replicates = sorted([x for x in df.columns.tolist() if 'Log2FC_Replicate' in x])

    df = df.loc[df['sgRNA_Type'] != 'negative_control']
    df = df.dropna(axis=0, how='any')

    if chrom not in df['Chr'].unique():
        unique_chr = [x for x in df['Chr'].unique()]
        chrom = unique_chr[0]

    df = df.loc[df['Chr'] == chrom]

    # fig.append_trace(go.Scattergl(
    #     x = np.random.randn(1000),
    #     y = np.random.randn(1000),
    #     mode = 'markers',
    #     ), 1, 1)


    for replicate in replicates:

        fig.append_trace(go.Scattergl(
            x=df['Perturbation_Index'],
            y=df[replicate],
            text=replicate,
            mode='markers',
            opacity=0.5,
            marker={
                'size': 6,
                'line': {'width': 0.3, 'color': 'white'}
            },
            name= replicate,
            # yaxis = 'y1'
        ), 1, 1)

    beta_profile = '/SURF/precomputed/%s/SURF_result/beta_profile.csv' % (dataset)
    df2 = pd.read_csv(beta_profile)

    # Filter for chrom
    indices_filt = []
    y_p_filt = []
    power_filt = []

    for chrom_i, coord_i, value_i, power_i in zip(df2['Chr'], df2['Index'], df2['Beta'], df2['Statistical_Power']):
        if chrom_i == chrom:
            indices_filt.append(int(coord_i))
            y_p_filt.append(float(value_i))
            power_filt.append(float(power_i))

    # Boundaries of inference
    parameters_json = '/SURF/precomputed/%s/%s.json' % (dataset, dataset)

    with open(parameters_json, 'r') as f:
        json_string = f.readline().strip()
        param_dict = json.loads(json_string)

    scale_val = int(param_dict['scale'])

    diff_vec = [0] + list(np.diff(indices_filt))

    boundaries = [0]
    for index, value in enumerate(diff_vec):
        try:
            if int(value) > int(scale_val):
                boundaries.append(index)
        except:
            pass

    boundaries = sorted(boundaries)

    boundaries.append(len(indices_filt))

    for i in range(len(boundaries) - 1):
        start_index, stop_index = boundaries[i], boundaries[i + 1]

        fig.append_trace(go.Scatter(
            x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1],
            y=[0] + power_filt[start_index:stop_index] + [0],
            mode = 'lines',
            fill='tozeroy',
            showlegend=False,
            # yaxis = 'y2',
            line=dict(color='rgb(255,165,0)')), 2, 1)

    for i in range(len(boundaries) - 1):
        start_index, stop_index = boundaries[i], boundaries[i + 1]

        fig.append_trace(go.Scatter(
            x=[indices_filt[start_index] - 1] + indices_filt[start_index:stop_index] + [indices_filt[stop_index - 1] + 1], y=[0] + y_p_filt[start_index:stop_index] + [0],
            mode = 'lines',
            fill='tozeroy',
            yaxis = 'y4',
            showlegend=False,
            line=dict(color='rgb(169,169,169)')), 3, 1)

    # Find indices with significant padj-values
    fdr = float(param_dict['fdr'])
    significant_boundary_indices = []
    padj_list = []
    for padj in df2['Pval_adj.']:
        try:
            padj_list.append(float(padj))
        except:
            padj_list.append(0)

    for i in range(len(boundaries) - 1):
        start_index, stop_index = boundaries[i], boundaries[i + 1]
        significant_indices = [1 if x < fdr else 0 for x in padj_list[start_index:stop_index]]
        significant_boundary_indices_tmp = [j + start_index for j, k in enumerate(np.diff([0] + significant_indices)) if k != 0]

        if len(significant_boundary_indices_tmp)%2 == 1:
            significant_boundary_indices_tmp.append(stop_index - 1)

        significant_boundary_indices += significant_boundary_indices_tmp

    for i, j in zip(significant_boundary_indices[0::2], significant_boundary_indices[1::2]):

        boundary_start = int(i)
        boundary_stop = int(j)

        genomic_boundary_start = int(df2['Index'][boundary_start])
        genomic_boundary_stop = int(df2['Index'][boundary_stop])

        signal_mean = np.mean(df2['Beta'][boundary_start:boundary_stop])

        if signal_mean > np.median(df2['Beta']):
            fig.append_trace(go.Scatter(
                x=[genomic_boundary_start, genomic_boundary_stop],
                y=[max(y_p_filt) + 0.01, max(y_p_filt) + 0.01],
                mode = 'lines',
                # fill='tozeroy',
                showlegend=False,
                yaxis = 'y2',
                line=dict(color='rgb(255,0,0)', width = 5)), 3, 1)
        else:
            fig.append_trace(go.Scatter(
                x=[genomic_boundary_start, genomic_boundary_stop],
                y=[min(y_p_filt) - 0.01, min(y_p_filt) - 0.01],
                mode = 'lines',
                # fill='tozeroy',
                showlegend=False,
                yaxis = 'y2',
                line=dict(color='rgb(30,144,255)', width = 5)), 3, 1)

    try:
        start = int(start)
        stop = int(stop)
        fig['layout'].update(
            height=800,
            autosize = True,
            xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False, 'range':[start, stop]},
            yaxis={'domain':[0, 0.4], 'showgrid': False, 'title': 'Log2FC Scores', 'autorange':True},
            yaxis2={'domain':[0.45, 0.55], 'showgrid': False, 'title': 'Statistical Power', 'autorange':True},
            yaxis3={'domain':[0.6, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
            hovermode='closest',
            title='Raw Scores and Deconvolution')

    except:
        fig['layout'].update(
            height=800,
            autosize = True,
            xaxis={'type': 'linear', 'title': 'Genomic Coordinate', 'showgrid': False},
            yaxis={'domain':[0, 0.4], 'showgrid': False, 'title': 'Log2FC Scores', 'autorange':True},
            yaxis2={'domain':[0.45, 0.55], 'showgrid': False, 'title': 'Statistical Power', 'autorange':True},
            yaxis3={'domain':[0.6, 1], 'showgrid': False, 'title': 'Deconvolution', 'autorange':True},
            hovermode='closest',
            title='Raw Scores and Deconvolution')

    return fig

### CALLBACKS FOR SIGNIFICANCE VIEW ###
@app2.callback(Output('regions-table-view', 'rows'),
              [Input('dataset', 'value')])

def update_output(dataset):

    sig_regions = '/SURF/precomputed/%s/SURF_result/significant_regions.csv' % (dataset)

    largest_row = 0
    with open(sig_regions, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            if len(line) > largest_row:
                largest_row = len(line)

    df = pd.read_csv(sig_regions, names = range(largest_row))
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(0))

    df = df[['FDR','Chr','Start','Stop','Direction','Deconvolved_Signal_Area','Deconvolved_Signal_Mean','Padj_Mean','Supporting_sgRNAs']]
    return df.to_dict('records')

@app2.callback(Output('regions-table-view', 'columns'),
              [Input('dataset', 'value')])
def update_output(dataset):

    sig_regions = '/SURF/precomputed/%s/SURF_result/significant_regions.csv' % (dataset)

    largest_row = 0
    with open(sig_regions, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            if len(line) > largest_row:
                largest_row = len(line)

    df = pd.read_csv(sig_regions, names = range(largest_row))
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(0))

    df = df[['FDR','Chr','Start','Stop','Direction','Deconvolved_Signal_Area','Deconvolved_Signal_Mean','Padj_Mean','Supporting_sgRNAs']]
    return df.columns

@app2.callback(
    Output('download-total', 'href'),
    [Input('dataset', 'value')])

def update_link(dataset):

    igv_sesh = '/SURF/precomputed/%s/SURF_result/igv_session.zip' % (dataset)
    return igv_sesh

@app2.server.route('/SURF/precomputed/<directory>/SURF_result/igv_session.zip')
def generate_report_url2(directory):

    return send_file('/SURF/precomputed/%s/SURF_result/igv_session.zip' % (directory), attachment_filename = 'igv_session.zip', as_attachment = True)

def main():
    app.run_server(debug = True, processes = 5, port = 9992, host = '0.0.0.0')

if __name__ == '__main__':
	main()
