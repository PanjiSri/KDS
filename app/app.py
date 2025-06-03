import dash
from dash import html, dcc, Input, Output, State, dash_table
import base64
import io
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import traceback
import urllib.parse

from app.data_processing import *

app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server

app.layout = html.Div(
    style={
        'maxWidth': '800px',
        'margin': '50px auto',
        'padding': '30px',
        'fontFamily': '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
        'boxShadow': '0px 4px 20px rgba(0,0,0,0.05)',
        'borderRadius': '12px',
        'backgroundColor': '#fff'
    },
    children=[
        html.H1(
            'Analisis Struktur Populasi RAD-Seq',
            style={
                'fontSize': '28px',
                'fontWeight': '600',
                'color': '#2c3e50',
                'marginBottom': '15px',
                'textAlign': 'center',
                'letterSpacing': '-0.5px'
            }
        ),
        
        dcc.Tabs(id='analysis-tabs', value='pca-tab', children=[
            dcc.Tab(label='Analisis PCA', value='pca-tab', children=[
                html.Div([
                    dcc.Store(id='store-vcf-raw-data', storage_type='memory'),
                    
                    dcc.Upload(
                        id='upload-vcf',
                        children=html.Div([
                            html.Div(
                                '📁',
                                style={'fontSize': '38px', 'marginBottom': '12px', 'color': '#3498db'}
                            ),
                            html.Div(
                                'Unggah VCF',
                                style={'fontSize': '15px', 'color': '#555', 'fontWeight': '500'}
                            ),
                            html.Div(
                                '(.vcf atau .vcf.gz)',
                                style={'fontSize': '12px', 'color': '#7f8c8d', 'marginTop': '5px'}
                            )
                        ]),
                        style={
                            'width': '100%',
                            'height': '150px',
                            'borderWidth': '2px',
                            'borderStyle': 'dashed',
                            'borderRadius': '10px',
                            'borderColor': '#bdc3c7',
                            'backgroundColor': '#f8f9fa',
                            'cursor': 'pointer',
                            'display': 'flex',
                            'flexDirection': 'column',
                            'alignItems': 'center',
                            'justifyContent': 'center',
                            'transition': 'all 0.3s ease',
                            'marginBottom': '20px',
                            'marginTop': '20px'
                        },
                        multiple=False
                    ),
                    
                    html.Div(id='output-vcf-upload-status', style={'textAlign': 'center', 'minHeight': '20px', 'marginBottom': '20px'}),

                    html.Div([
                        html.H4("Pengaturan Analisis", style={'marginTop': '20px', 'marginBottom': '15px', 'fontSize': '18px', 'color': '#34495e'}),
                        html.Div([
                            html.Label("Ambang Batas MAF (0.00 - 0.50):", style={'fontWeight': '500', 'fontSize': '14px', 'color': '#555'}),
                            dcc.Input(id='maf-threshold-input', type='number', value=0.05, min=0.00, max=0.50, step=0.01, style={'width': '100%', 'padding': '8px', 'borderRadius': '4px', 'border': '1px solid #ccc', 'boxSizing': 'border-box', 'marginTop': '5px'}),
                        ], style={'marginBottom': '15px'}),
                        
                        html.Div([
                            html.Label("Ambang Batas Data Hilang per SNP (0.00 - 1.00):", style={'fontWeight': '500', 'fontSize': '14px', 'color': '#555'}),
                            dcc.Input(id='snp-missing-threshold-input', type='number', value=0.2, min=0.00, max=1.00, step=0.05, style={'width': '100%', 'padding': '8px', 'borderRadius': '4px', 'border': '1px solid #ccc', 'boxSizing': 'border-box', 'marginTop': '5px'}),
                        ], style={'marginBottom': '15px'}),

                        html.Div([
                            html.Label("Ambang Batas Data Hilang per Individu (0.00 - 1.00):", style={'fontWeight': '500', 'fontSize': '14px', 'color': '#555'}),
                            dcc.Input(id='ind-missing-threshold-input', type='number', value=0.2, min=0.00, max=1.00, step=0.05, style={'width': '100%', 'padding': '8px', 'borderRadius': '4px', 'border': '1px solid #ccc', 'boxSizing': 'border-box', 'marginTop': '5px'}),
                        ], style={'marginBottom': '15px'}),

                        html.Div([
                            html.Label("Jumlah Komponen PCA (1 - 10):", style={'fontWeight': '500', 'fontSize': '14px', 'color': '#555'}),
                            dcc.Input(id='n-pca-components-input', type='number', value=3, min=1, max=10, step=1, style={'width': '100%', 'padding': '8px', 'borderRadius': '4px', 'border': '1px solid #ccc', 'boxSizing': 'border-box', 'marginTop': '5px'}),
                        ], style={'marginBottom': '15px'}),
                    ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'borderRadius': '8px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)', 'marginBottom': '20px'}),

                    html.Button(
                        'Mulai Analisis',
                        id='start-analysis-button',
                        n_clicks=0,
                        style={
                            'width': '100%',
                            'padding': '15px',
                            'fontSize': '17px',
                            'fontWeight': '600',
                            'cursor': 'pointer',
                            'backgroundColor': '#2c3e50',
                            'color': 'white',
                            'border': 'none',
                            'borderRadius': '8px',
                            'transition': 'all 0.2s ease',
                            'boxShadow': '0px 2px 5px rgba(0,0,0,0.1)'
                        }
                    ),
                    html.Div(id='analysis-error-output', style={'color': 'red', 'marginTop': '15px', 'textAlign': 'center', 'fontWeight': 'bold'}),
                    
                    dcc.Store(id='store-pca-results-json'),

                    dcc.Loading(
                        id="loading-analysis-output",
                        type="circle",
                        children=[
                            html.Div(id='analysis-results-container', style={'marginTop': '30px'})
                        ],
                        overlay_style={'visibility':'hidden'}
                    )
                ])
            ]),
            
            dcc.Tab(label='Analisis FST', value='fst-tab', children=[
                html.Div([
                    dcc.Store(id='store-pooled-data', storage_type='memory'),
                    dcc.Store(id='store-fst-results', storage_type='memory'),
                    
                    dcc.Upload(
                        id='upload-pooled-data',
                        children=html.Div([
                            html.Div(
                                '📊',
                                style={'fontSize': '38px', 'marginBottom': '12px', 'color': '#e74c3c'}
                            ),
                            html.Div(
                                'Unggah Data Pooled',
                                style={'fontSize': '15px', 'color': '#555', 'fontWeight': '500'}
                            ),
                            html.Div(
                                '(data.auto)',
                                style={'fontSize': '12px', 'color': '#7f8c8d', 'marginTop': '5px'}
                            )
                        ]),
                        style={
                            'width': '100%',
                            'height': '150px',
                            'borderWidth': '2px',
                            'borderStyle': 'dashed',
                            'borderRadius': '10px',
                            'borderColor': '#bdc3c7',
                            'backgroundColor': '#f8f9fa',
                            'cursor': 'pointer',
                            'display': 'flex',
                            'flexDirection': 'column',
                            'alignItems': 'center',
                            'justifyContent': 'center',
                            'transition': 'all 0.3s ease',
                            'marginBottom': '20px',
                            'marginTop': '20px'
                        },
                        multiple=False
                    ),
                    
                    html.Div(id='output-pooled-upload-status', style={'textAlign': 'center', 'minHeight': '20px', 'marginBottom': '20px'}),
                    
                    html.Div([
                        html.H4("Pengaturan Analisis FST", style={'marginTop': '20px', 'marginBottom': '15px', 'fontSize': '18px', 'color': '#34495e'}),
                        html.Div([
                            html.Label("Kedalaman Minimal Per Pool (1-50):", style={'fontWeight': '500', 'fontSize': '14px', 'color': '#555'}),
                            dcc.Input(id='min-depth-input', type='number', value=10, min=1, max=50, step=1, style={'width': '100%', 'padding': '8px', 'borderRadius': '4px', 'border': '1px solid #ccc', 'boxSizing': 'border-box', 'marginTop': '5px'}),
                        ], style={'marginBottom': '15px'}),
                    ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'borderRadius': '8px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)', 'marginBottom': '20px'}),
                    
                    html.Button(
                        'Hitung FST',
                        id='calculate-fst-button',
                        n_clicks=0,
                        style={
                            'width': '100%',
                            'padding': '15px',
                            'fontSize': '17px',
                            'fontWeight': '600',
                            'cursor': 'pointer',
                            'backgroundColor': '#e74c3c',
                            'color': 'white',
                            'border': 'none',
                            'borderRadius': '8px',
                            'transition': 'all 0.2s ease',
                            'boxShadow': '0px 2px 5px rgba(0,0,0,0.1)'
                        }
                    ),
                    
                    html.Div(id='fst-error-output', style={'color': 'red', 'marginTop': '15px', 'textAlign': 'center', 'fontWeight': 'bold'}),
                    
                    dcc.Loading(
                        id="loading-fst-output",
                        type="circle",
                        children=[
                            html.Div(id='fst-results-container', style={'marginTop': '30px'})
                        ]
                    )
                ])
            ])
        ])
    ]
)

app.clientside_callback(
    """
    function(n) {
        var upload = document.getElementById('upload-vcf');
        if (upload) {
            upload.onmouseover = function() {
                this.style.borderColor = '#3498db';
                this.style.backgroundColor = '#eaf4fb';
            };
            upload.onmouseout = function() {
                this.style.borderColor = '#bdc3c7';
                this.style.backgroundColor = '#f8f9fa';
            };
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('upload-vcf', 'style', allow_duplicate=True),
    Input('upload-vcf', 'id'),
    prevent_initial_call=True
)

app.clientside_callback(
    """
    function(n) {
        var upload = document.getElementById('upload-pooled-data');
        if (upload) {
            upload.onmouseover = function() {
                this.style.borderColor = '#e74c3c';
                this.style.backgroundColor = '#fbeaea';
            };
            upload.onmouseout = function() {
                this.style.borderColor = '#bdc3c7';
                this.style.backgroundColor = '#f8f9fa';
            };
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('upload-pooled-data', 'style', allow_duplicate=True),
    Input('upload-pooled-data', 'id'),
    prevent_initial_call=True
)

app.clientside_callback(
    """
    function(n) {
        var button = document.getElementById('start-analysis-button');
        if (button) {
            button.onmouseover = function() {
                this.style.backgroundColor = '#34495e'; 
                this.style.transform = 'translateY(-1px)';
            };
            button.onmouseout = function() {
                this.style.backgroundColor = '#2c3e50';
                this.style.transform = 'translateY(0)';
            };
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('start-analysis-button', 'style', allow_duplicate=True),
    Input('start-analysis-button', 'id'),
    prevent_initial_call=True
)

app.clientside_callback(
    """
    function(n) {
        var button = document.getElementById('calculate-fst-button');
        if (button) {
            button.onmouseover = function() {
                this.style.backgroundColor = '#c0392b'; 
                this.style.transform = 'translateY(-1px)';
            };
            button.onmouseout = function() {
                this.style.backgroundColor = '#e74c3c';
                this.style.transform = 'translateY(0)';
            };
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('calculate-fst-button', 'style', allow_duplicate=True),
    Input('calculate-fst-button', 'id'),
    prevent_initial_call=True
)

@app.callback(
    [Output('store-vcf-raw-data', 'data'),
     Output('output-vcf-upload-status', 'children'),
     Output('analysis-results-container', 'children'),
     Output('analysis-error-output', 'children')],
    Input('upload-vcf', 'contents'),
    State('upload-vcf', 'filename'),
    prevent_initial_call=True
)
def update_vcf_store_and_clear_results(contents, filename):
    if contents is not None:
        data_summary, message = parse_vcf_to_json_summary(contents, filename)
        if data_summary:
            return (
                {'filename': data_summary['filename'], 'vcf_contents_base64': data_summary['vcf_contents_base64']},
                html.Div(
                    f"✓ Berkas '{data_summary['filename']}' siap dianalisis. ({message})",
                    style={'fontSize': '14px', 'color': '#27ae60'}
                ),
                None,
                None
            )
        else:
            return (
                None,
                html.Div(
                    f"✗ Gagal memuat berkas '{filename}': {message}",
                    style={'fontSize': '14px', 'color': '#c0392b'}
                ),
                None,
                None
            )
    return dash.no_update, dash.no_update, None, None

@app.callback(
    [Output('store-pooled-data', 'data'),
     Output('output-pooled-upload-status', 'children'),
     Output('fst-results-container', 'children', allow_duplicate=True),
     Output('fst-error-output', 'children', allow_duplicate=True)],
    Input('upload-pooled-data', 'contents'),
    State('upload-pooled-data', 'filename'),
    prevent_initial_call=True
)
def update_pooled_store_and_clear_results(contents, filename):
    if contents is not None:
        pooled_data, message = parse_pooled_data(contents, filename)
        if pooled_data:
            return (
                {'pooled_df_json': pooled_data},
                html.Div(
                    f"✓ {message}",
                    style={'fontSize': '14px', 'color': '#27ae60'}
                ),
                None,
                None
            )
        else:
            return (
                None,
                html.Div(
                    f"✗ Gagal memuat berkas '{filename}': {message}",
                    style={'fontSize': '14px', 'color': '#c0392b'}
                ),
                None,
                None
            )
    return dash.no_update, dash.no_update, None, None

@app.callback(
    [Output('store-pca-results-json', 'data'),
     Output('analysis-error-output', 'children', allow_duplicate=True)],
    Input('start-analysis-button', 'n_clicks'),
    [State('store-vcf-raw-data', 'data'),
     State('maf-threshold-input', 'value'),
     State('snp-missing-threshold-input', 'value'),
     State('ind-missing-threshold-input', 'value'),
     State('n-pca-components-input', 'value')],
    prevent_initial_call=True
)
def run_analysis_pipeline_callback(n_clicks, vcf_raw_data,
                                   maf_thresh, snp_miss_thresh,
                                   ind_miss_thresh, n_pca_components):
    if not vcf_raw_data or 'vcf_contents_base64' not in vcf_raw_data or 'filename' not in vcf_raw_data:
        return None, "Silakan unggah berkas VCF terlebih dahulu."

    vcf_contents_base64 = vcf_raw_data['vcf_contents_base64']
    filename = vcf_raw_data['filename']
    
    error_messages = []
    if maf_thresh is None or not (0.0 <= maf_thresh <= 0.5):
        error_messages.append("Ambang batas MAF harus antara 0.00 dan 0.50.")
    if snp_miss_thresh is None or not (0.0 <= snp_miss_thresh <= 1.0):
        error_messages.append("Ambang batas data hilang SNP harus antara 0.00 dan 1.00.")
    if ind_miss_thresh is None or not (0.0 <= ind_miss_thresh <= 1.0):
        error_messages.append("Ambang batas data hilang individu harus antara 0.00 dan 1.00.")
    if n_pca_components is None or not (1 <= n_pca_components <= 10):
        error_messages.append("Jumlah komponen PCA harus antara 1 dan 10.")
    
    if error_messages:
        return None, " ".join(error_messages)

    try:
        results = trigger_analysis_pipeline(vcf_contents_base64, filename,
                                            maf_thresh=float(maf_thresh),
                                            snp_miss_thresh=float(snp_miss_thresh),
                                            ind_miss_thresh=float(ind_miss_thresh),
                                            n_pca_components=int(n_pca_components))
        return results, None
    except Exception as e:
        error_message = f"Terjadi kesalahan: {str(e)}"
        if "Pipeline analisis gagal" in str(e) or "Kesalahan saat" in str(e) or "Tidak ditemukan" in str(e) or "Tidak ada" in str(e) or "Data tidak cukup" in str(e):
            error_message = str(e)
        return None, error_message

@app.callback(
    [Output('store-fst-results', 'data'),
     Output('fst-error-output', 'children', allow_duplicate=True)],
    Input('calculate-fst-button', 'n_clicks'),
    [State('store-pooled-data', 'data'),
     State('min-depth-input', 'value')],
    prevent_initial_call=True
)
def calculate_fst_callback(n_clicks, pooled_data, min_depth):
    if not pooled_data or 'pooled_df_json' not in pooled_data:
        return None, "Silakan unggah berkas pooled data terlebih dahulu."
    
    if min_depth is None or not (1 <= min_depth <= 50):
        return None, "Minimum depth harus antara 1 dan 50."
    
    print(f"DEBUG: calculate_fst_callback triggered. min_depth from input: {min_depth} (type: {type(min_depth)})")
    try:
        results = analyze_fst_from_pooled_data(pooled_data['pooled_df_json'], min_depth=int(min_depth))
        print("DEBUG: calculate_fst_callback: analyze_fst_from_pooled_data successful.")
        return results, None
    except Exception as e:
        print("ERROR in calculate_fst_callback during FST analysis:")
        print(f"Exception type: {type(e).__name__}")
        print(f"Exception message: {str(e)}")
        print("Full Traceback:")
        traceback.print_exc()
        
        error_message = str(e)
        if "Kesalahan selama analisis FST:" in error_message or "Kesalahan saat memparsing" in error_message :
            user_facing_error = error_message
        else:
            user_facing_error = f"Kesalahan selama analisis FST: {error_message}"
        return None, user_facing_error

@app.callback(
    Output('analysis-results-container', 'children', allow_duplicate=True),
    Input('store-pca-results-json', 'data'),
    prevent_initial_call=True
)
def display_analysis_results(pca_results_json):
    if not pca_results_json:
        return None

    pca_coords_df_json = pca_results_json.get('pca_coords_df_json')
    variance_explained = pca_results_json.get('variance_explained')
    analysis_summary = pca_results_json.get('analysis_summary')

    pca_plot_div = html.Div("Gagal memuat plot PCA.")
    pca_df = None
    fig = None 

    if pca_coords_df_json and variance_explained:
        try:
            pca_df = pd.read_json(io.StringIO(pca_coords_df_json), orient='split')
            
            has_pc1 = 'PC1' in pca_df.columns and len(variance_explained) >= 1
            has_pc2 = 'PC2' in pca_df.columns and len(variance_explained) >= 2
            has_pc3 = 'PC3' in pca_df.columns and len(variance_explained) >= 3

            plot_title_prefix = "<b>Plot Analisis Komponen Utama (PCA)</b>"
            plot_font = dict(family='-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif')
            plot_paper_bgcolor = 'rgba(0,0,0,0)'
            plot_bgcolor = '#f8f9fa'
            plot_height = 550 

            if has_pc3:
                pc1_var = variance_explained[0]*100
                pc2_var = variance_explained[1]*100
                pc3_var = variance_explained[2]*100
                
                fig = px.scatter_3d(pca_df, x='PC1', y='PC2', z='PC3', hover_name='Sample',
                                   labels={'PC1': f'PC1 ({pc1_var:.2f}%)',
                                           'PC2': f'PC2 ({pc2_var:.2f}%)',
                                           'PC3': f'PC3 ({pc3_var:.2f}%)'})
                fig.update_layout(
                    title_text=f"{plot_title_prefix} - 3D", title_x=0.5, font=plot_font,
                    scene=dict(
                        xaxis_title=f'PC1 ({pc1_var:.2f}%)',
                        yaxis_title=f'PC2 ({pc2_var:.2f}%)',
                        zaxis_title=f'PC3 ({pc3_var:.2f}%)',
                        xaxis=dict(backgroundcolor=plot_bgcolor, gridcolor='#e0e0e0', zeroline=False),
                        yaxis=dict(backgroundcolor=plot_bgcolor, gridcolor='#e0e0e0', zeroline=False),
                        zaxis=dict(backgroundcolor=plot_bgcolor, gridcolor='#e0e0e0', zeroline=False)
                    ),
                    margin=dict(l=0, r=0, b=0, t=40), 
                    height=plot_height,
                    paper_bgcolor=plot_paper_bgcolor
                )
            elif has_pc2: 
                pc1_var = variance_explained[0]*100
                pc2_var = variance_explained[1]*100
                
                fig = px.scatter(pca_df, x='PC1', y='PC2', hover_name='Sample',
                                 labels={'PC1': f'PC1 ({pc1_var:.2f}%)',
                                         'PC2': f'PC2 ({pc2_var:.2f}%)'})
                fig.update_layout(
                    title_text=f"{plot_title_prefix} - 2D", title_x=0.5, font=plot_font,
                    paper_bgcolor=plot_paper_bgcolor, plot_bgcolor=plot_bgcolor,
                    xaxis=dict(gridcolor='#e0e0e0', zeroline=False),
                    yaxis=dict(gridcolor='#e0e0e0', zeroline=False),
                    height=500 
                )
            elif has_pc1: 
                 pc1_var = variance_explained[0]*100
                 pca_df['dummy_y'] = 0 
                 fig = px.scatter(pca_df, x='PC1', y='dummy_y', hover_name='Sample',
                                 labels={'PC1': f'PC1 ({pc1_var:.2f}%)', 'dummy_y': ''})
                 fig.update_yaxes(visible=False, showticklabels=False, zeroline=False, showgrid=False)
                 fig.update_xaxes(zeroline=False)
                 fig.update_layout(
                    title_text=f"{plot_title_prefix} - 1D (PC1)", title_x=0.5, font=plot_font,
                    paper_bgcolor=plot_paper_bgcolor, plot_bgcolor=plot_bgcolor,
                    height=300 
                )
            else:
                pca_plot_div = html.Div("Data PCA tidak memiliki komponen yang cukup untuk diplot (minimal PC1 diperlukan).")
            
            if fig:
                 pca_plot_div = dcc.Graph(figure=fig)

        except Exception as e:
            pca_plot_div = html.Div(f"Kesalahan saat membuat plot PCA: {str(e)}")
            pca_df = None 

    download_pca_link_div = html.Div() 
    if pca_df is not None and not pca_df.empty:
        try:
            pc_cols = [col for col in pca_df.columns if col.startswith('PC')]
            if 'Sample' in pca_df.columns and pc_cols:
                df_for_pca_file = pca_df[['Sample'] + pc_cols]
                pca_data_string = df_for_pca_file.to_csv(sep=' ', index=False, header=True)
                
                download_pca_link_div = html.Div(
                    html.A(
                        'Download Hasil PCA (.pca)',
                        id='download-pca-link',
                        download="pca_results.pca",
                        href=f"data:text/plain;charset=utf-8,{urllib.parse.quote(pca_data_string)}",
                        target="_blank",
                        style={
                            'display': 'inline-block',
                            'padding': '10px 20px',
                            'backgroundColor': '#3498db',
                            'color': 'white',
                            'textDecoration': 'none',
                            'borderRadius': '5px',
                            'fontSize': '14px',
                            'fontWeight': '500'
                        }
                    ), style={'textAlign': 'center', 'marginTop': '20px'}
                )
        except Exception as e:
            download_pca_link_div = html.Div(f"Gagal membuat link download PCA: {str(e)}")


    variance_table_div = html.Div("Gagal memuat tabel varians.")
    if variance_explained:
        try:
            rows = []
            max_pcs_to_show = len(variance_explained)
            for i in range(max_pcs_to_show):
                rows.append(html.Tr([
                    html.Td(f"PC{i+1}", style={'padding': '8px', 'border': '1px solid #ddd'}),
                    html.Td(f"{variance_explained[i]*100:.2f}%", style={'padding': '8px', 'border': '1px solid #ddd', 'textAlign': 'right'})
                ]))
            
            variance_table_div = html.Div([
                html.H4("Varians yang Dijelaskan", style={'textAlign': 'left', 'marginTop': '25px', 'marginBottom': '10px', 'fontSize': '18px', 'color': '#34495e'}),
                html.Table([
                    html.Thead(html.Tr([
                        html.Th("Komponen", style={'padding': '10px', 'border': '1px solid #ddd', 'backgroundColor': '#ecf0f1', 'textAlign': 'left'}),
                        html.Th("Persentase Varians", style={'padding': '10px', 'border': '1px solid #ddd', 'backgroundColor': '#ecf0f1', 'textAlign': 'right'})
                    ])),
                    html.Tbody(rows)
                ], style={'margin': '0', 'borderCollapse': 'collapse', 'width': 'auto', 'fontSize': '14px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)'}),
            ])
        except Exception as e:
            variance_table_div = html.Div(f"Kesalahan saat membuat tabel varians: {str(e)}")

    summary_stats_div = html.Div("Gagal memuat statistik ringkasan.")
    if analysis_summary:
        try:
            summary_stats_div = html.Div([
                html.H4("Ringkasan Kontrol Kualitas", style={'textAlign': 'left', 'marginTop': '25px', 'marginBottom': '10px', 'fontSize': '18px', 'color': '#34495e'}),
                html.P(f"Jumlah sampel awal: {analysis_summary.get('samples_original', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'5px'}),
                html.P(f"Jumlah sampel setelah QC: {analysis_summary.get('samples_after_qc', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'5px'}),
                html.P(f"Jumlah SNP awal: {analysis_summary.get('snps_original', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'5px'}),
                html.P(f"Jumlah SNP setelah QC: {analysis_summary.get('snps_after_qc', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'15px'}),
            ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'borderRadius': '8px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)'})
        except Exception as e:
            summary_stats_div = html.Div(f"Kesalahan saat menampilkan statistik ringkasan: {str(e)}")
            
    return html.Div([
        html.Hr(style={'borderColor': '#ecf0f1', 'margin': '30px 0'}),
        pca_plot_div,
        download_pca_link_div,
        html.Div([
            variance_table_div,
            summary_stats_div
        ], style={'display': 'grid', 'gridTemplateColumns': 'minmax(250px, auto) 1fr', 'gap': '30px', 'alignItems': 'start', 'marginTop': '20px'})
    ], style={'padding': '0 20px'})

@app.callback(
    Output('fst-results-container', 'children'),
    Input('store-fst-results', 'data'),
    prevent_initial_call=True
)
def display_fst_results(fst_results):
    if not fst_results:
        return None
    
    try:
        fst_matrix_df = pd.read_json(io.StringIO(fst_results['fst_matrix']), orient='split')
        
        heatmap_fig = go.Figure(data=go.Heatmap(
            z=fst_matrix_df.values,
            x=fst_matrix_df.columns,
            y=fst_matrix_df.index,
            colorscale='RdBu_r',
            colorbar=dict(title="FST"),
            text=np.round(fst_matrix_df.values, 3),
            texttemplate="%{text}",
            textfont={"size": 10},
            hovertemplate="Pool 1: %{y}<br>Pool 2: %{x}<br>FST: %{z:.4f}<extra></extra>"
        ))
        
        heatmap_fig.update_layout(
            title=dict(text="<b>Heatmap FST Antar Pool</b>", x=0.5, font=dict(size=20)),
            xaxis=dict(tickangle=-45),
            height=600,
            margin=dict(l=100, r=100, t=100, b=100),
            font=dict(family='-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif')
        )
        
        download_link = html.A(
            'Download Matrix FST (CSV)',
            id='download-fst-csv',
            download="fst_matrix.csv",
            href="data:text/csv;charset=utf-8," + urllib.parse.quote(fst_matrix_df.to_csv()), # Menggunakan urllib.parse.quote untuk FST juga
            target="_blank",
            style={
                'display': 'inline-block',
                'padding': '10px 20px',
                'backgroundColor': '#3498db',
                'color': 'white',
                'textDecoration': 'none',
                'borderRadius': '5px',
                'fontSize': '14px',
                'fontWeight': '500'
            }
        )

        fst_summary_data = fst_results.get('fst_summary_stats')
        fst_value_distribution_data = fst_results.get('fst_value_distribution_table_data')

        summary_div = html.Div("Gagal memuat ringkasan analisis FST.")
        if fst_summary_data:
            try:
                summary_div = html.Div([
                    html.H4("Ringkasan Analisis FST", style={'textAlign': 'left', 'marginTop': '25px', 'marginBottom': '10px', 'fontSize': '18px', 'color': '#34495e'}),
                    html.P(f"Jumlah pool dianalisis: {fst_summary_data.get('num_pools', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'5px'}),
                    html.P(f"Jumlah SNP dalam data input: {fst_summary_data.get('num_snps_input', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'5px'}),
                    html.P(f"Ambang batas kedalaman minimal per pool: {fst_summary_data.get('min_depth_filter', 'N/A')}", style={'fontSize': '14px', 'marginBottom':'15px'}),
                ], style={'padding': '20px', 'backgroundColor': '#f8f9fa', 'borderRadius': '8px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)'})
            except Exception as e:
                summary_div = html.Div(f"Kesalahan saat menampilkan ringkasan FST: {str(e)}")

        distribution_table_div = html.Div("Gagal memuat tabel distribusi nilai FST.")
        if fst_value_distribution_data:
            try:
                rows = []
                for item in fst_value_distribution_data:
                    rows.append(html.Tr([
                        html.Td(item['Statistic'], style={'padding': '8px', 'border': '1px solid #ddd'}),
                        html.Td(str(item['Value']), style={'padding': '8px', 'border': '1px solid #ddd', 'textAlign': 'right'})
                    ]))
                
                distribution_table_div = html.Div([
                    html.H4("Distribusi Nilai FST Pairwise", style={'textAlign': 'left', 'marginTop': '25px', 'marginBottom': '10px', 'fontSize': '18px', 'color': '#34495e'}),
                    html.Table([
                        html.Thead(html.Tr([
                            html.Th("Statistik", style={'padding': '10px', 'border': '1px solid #ddd', 'backgroundColor': '#ecf0f1', 'textAlign': 'left'}),
                            html.Th("Nilai", style={'padding': '10px', 'border': '1px solid #ddd', 'backgroundColor': '#ecf0f1', 'textAlign': 'right'})
                        ])),
                        html.Tbody(rows)
                    ], style={'margin': '0', 'borderCollapse': 'collapse', 'width': 'auto', 'fontSize': '14px', 'boxShadow': '0px 1px 3px rgba(0,0,0,0.05)'}),
                ])
            except Exception as e:
                distribution_table_div = html.Div(f"Kesalahan saat membuat tabel distribusi FST: {str(e)}")
        
        return html.Div([
            html.Hr(style={'borderColor': '#ecf0f1', 'margin': '30px 0'}),
            dcc.Graph(figure=heatmap_fig),
            html.Div([ 
                distribution_table_div,
                summary_div
            ], style={'display': 'grid', 'gridTemplateColumns': 'minmax(300px, auto) 1fr', 'gap': '30px', 'alignItems': 'start', 'marginTop': '20px'}),
            html.Div(download_link, style={'textAlign': 'center', 'marginTop': '30px', 'marginBottom': '20px'})
        ], style={'padding': '0 20px'})
        
    except Exception as e:
        print(f"Error displaying FST results: {e}")
        traceback.print_exc() 
        return html.Div(f"Kesalahan saat menampilkan hasil FST: {str(e)}", 
                       style={'color': 'red', 'textAlign': 'center', 'marginTop': '20px'})


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)