import dash
from dash import html, dcc, Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
import base64
import io
import pandas as pd
import numpy as np

from data_processing import (
    parse_vcf_to_json,     
    parse_pca_to_json,
    parse_admixture_to_json,
    parse_metadata_to_json
)

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(children=[
    html.H1(children='Zosterops borbonicus Population Structure Dashboard'),
    html.Div(children='A web dashboard for analyzing population structure based on RADSeq data.'),

    dcc.Store(id='store-vcf-data', storage_type='memory'),
    dcc.Store(id='store-pca-data', storage_type='memory'),
    dcc.Store(id='store-admixture-data', storage_type='memory'),
    dcc.Store(id='store-metadata-data', storage_type='memory'),

    html.Div(id='upload-section', children=[
        html.H3("Upload Data Files"),
        dcc.Upload(
            id='upload-vcf',
            children=html.Div(['Drag and Drop or ', html.A('Select VCF File')]),
            style={'width': '90%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px auto'},
            multiple=False
        ),
        html.Div(id='output-vcf-upload-status'),

        # PCA
        dcc.Upload(id='upload-pca', children=html.Div(['Drag and Drop or ', html.A('Select PCA File (.evec, .pca, or .csv)')]), style={'width': '90%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px auto'}, multiple=False),
        html.Div(id='output-pca-upload-status'),

        # ADMIXTURE
        dcc.Upload(id='upload-admixture', children=html.Div(['Drag and Drop or ', html.A('Select ADMIXTURE File (.Q)')]), style={'width': '90%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px auto'}, multiple=False),
        html.Div(id='output-admixture-upload-status'),

        # METADATA
        dcc.Upload(id='upload-metadata', children=html.Div(['Drag and Drop or ', html.A('Select Metadata File (.csv)')]), style={'width': '90%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px auto'}, multiple=False),
        html.Div(id='output-metadata-upload-status'),
    ]),

    html.Div(id='visualization-section', children=[
        html.H3("Visualizations"),
        html.Div(id='vcf-summary-output'),
        html.Div(id='pca-plot-output'),
        html.Div(id='admixture-plot-output'),
        html.Div(id='metadata-summary-output')
    ]),
])

@app.callback(
    [Output('store-vcf-data', 'data'), Output('output-vcf-upload-status', 'children')],
    Input('upload-vcf', 'contents'),
    State('upload-vcf', 'filename')
)
def update_vcf_store(contents, filename):
    if contents is not None:
        data, message = parse_vcf_to_json(contents, filename)
        return data, message
    return None, "Upload a VCF file."

@app.callback(
    [Output('store-pca-data', 'data'), Output('output-pca-upload-status', 'children')],
    Input('upload-pca', 'contents'),
    State('upload-pca', 'filename')
)
def update_pca_store(contents, filename):
    if contents is not None:
        data, message = parse_pca_to_json(contents, filename)
        return data, message
    return None, "Upload a PCA file."

@app.callback(
    [Output('store-admixture-data', 'data'), Output('output-admixture-upload-status', 'children')],
    Input('upload-admixture', 'contents'),
    State('upload-admixture', 'filename')
)
def update_admixture_store(contents, filename):
    if contents is not None:
        data, message = parse_admixture_to_json(contents, filename)
        return data, message
    return None, "Upload an ADMIXTURE file."

@app.callback(
    [Output('store-metadata-data', 'data'), Output('output-metadata-upload-status', 'children')],
    Input('upload-metadata', 'contents'),
    State('upload-metadata', 'filename')
)
def update_metadata_store(contents, filename):
    if contents is not None:
        data, message = parse_metadata_to_json(contents, filename)
        return data, message
    return None, "Upload a Metadata file."

@app.callback(
    Output('vcf-summary-output', 'children'),
    Input('store-vcf-data', 'data')
)
def display_vcf_summary(vcf_data):
    if vcf_data is not None:
        try:
            return [
                html.H4("VCF Data Summary:"),
                html.P(f"Number of samples: {len(vcf_data.get('samples', []))}"),
                html.P(f"Total variants: {vcf_data.get('total_variants', 'Unknown')}"),
                html.P(f"Samples: {', '.join(vcf_data.get('samples', [])[:10])}{'...' if len(vcf_data.get('samples', [])) > 10 else ''}"),
                html.Hr()
            ]
        except Exception as e:
            return f"Error displaying VCF summary: {e}"
    return ""

@app.callback(
    Output('metadata-summary-output', 'children'),
    Input('store-metadata-data', 'data')
)
def display_metadata_summary(jsonified_metadata):
    if jsonified_metadata is not None:
        try:
            df = pd.read_json(io.StringIO(jsonified_metadata), orient='split')
            return [
                html.H4("Metadata Summary:"),
                html.P(f"Number of samples: {len(df)}"),
                html.P(f"Columns: {', '.join(df.columns)}"),
                html.Hr()
            ]
        except Exception as e:
            return f"Error displaying metadata: {e}"
    return ""

@app.callback(
    Output('pca-plot-output', 'children'),
    [Input('store-pca-data', 'data'), Input('store-metadata-data', 'data')]
)
def create_pca_plot(pca_json, metadata_json):
    if pca_json is None:
        return html.Div("Upload PCA data to see visualization")
    
    try:
        pca_df = pd.read_json(io.StringIO(pca_json), orient='split')
        
        if pca_df.shape[1] < 2:
            return html.Div("Error: PCA data should have at least 2 principal components")
        
        if pca_df.iloc[:, 0].dtype == 'object':
            pca_df.index = pca_df.iloc[:, 0]
            pca_df = pca_df.iloc[:, 1:]
        
        pca_df.columns = [f'PC{i+1}' for i in range(pca_df.shape[1])]
        
        fig = go.Figure()
        
        if metadata_json is not None:
            try:
                metadata_df = pd.read_json(io.StringIO(metadata_json), orient='split')
                if len(set(pca_df.index) & set(metadata_df.index)) > 0:
                    merged_df = pca_df.merge(metadata_df, left_index=True, right_index=True, how='left')
                    color_col = None
                    for col in metadata_df.columns:
                        if metadata_df[col].dtype == 'object':
                            color_col = col
                            break
                    
                    if color_col:
                        for group in merged_df[color_col].unique():
                            if pd.notna(group):
                                group_data = merged_df[merged_df[color_col] == group]
                                fig.add_trace(go.Scatter(
                                    x=group_data['PC1'],
                                    y=group_data['PC2'],
                                    mode='markers',
                                    name=str(group),
                                    text=group_data.index,
                                    marker=dict(size=8)
                                ))
                    else:
                        fig.add_trace(go.Scatter(
                            x=pca_df['PC1'],
                            y=pca_df['PC2'],
                            mode='markers',
                            text=pca_df.index,
                            marker=dict(size=8)
                        ))
                else:
                    fig.add_trace(go.Scatter(
                        x=pca_df['PC1'],
                        y=pca_df['PC2'],
                        mode='markers',
                        text=pca_df.index,
                        marker=dict(size=8)
                    ))
            except:
                fig.add_trace(go.Scatter(
                    x=pca_df['PC1'],
                    y=pca_df['PC2'],
                    mode='markers',
                    text=pca_df.index,
                    marker=dict(size=8)
                ))
        else:
            fig.add_trace(go.Scatter(
                x=pca_df['PC1'],
                y=pca_df['PC2'],
                mode='markers',
                text=pca_df.index,
                marker=dict(size=8)
            ))
        
        fig.update_layout(
            title="PCA Plot",
            xaxis_title="PC1",
            yaxis_title="PC2",
            hovermode='closest',
            height=500
        )
        
        return [
            html.H4("Principal Component Analysis"),
            dcc.Graph(figure=fig),
            html.Hr()
        ]
        
    except Exception as e:
        return html.Div(f"Error creating PCA plot: {str(e)}")

@app.callback(
    Output('admixture-plot-output', 'children'),
    [Input('store-admixture-data', 'data'), Input('store-vcf-data', 'data'), Input('store-metadata-data', 'data')]
)
def create_admixture_plot(admixture_json, vcf_data, metadata_json):
    if admixture_json is None:
        return html.Div("Upload ADMIXTURE data to see visualization")
    
    try:
        admix_df = pd.read_json(io.StringIO(admixture_json), orient='split')
        
        if vcf_data and 'samples' in vcf_data:
            sample_names = vcf_data['samples']
            if len(sample_names) == len(admix_df):
                admix_df.index = sample_names
        
        admix_df.columns = [f'K{i+1}' for i in range(admix_df.shape[1])]
        
        fig = go.Figure()
        
        for col in admix_df.columns:
            fig.add_trace(go.Bar(
                name=col,
                x=list(range(len(admix_df))),
                y=admix_df[col],
                text=[f'{v:.3f}' for v in admix_df[col]],
                textposition='auto',
            ))
        
        fig.update_layout(
            title=f"ADMIXTURE Results (K={admix_df.shape[1]})",
            xaxis_title="Samples",
            yaxis_title="Ancestry Proportion",
            barmode='stack',
            height=400,
            xaxis=dict(
                tickmode='array',
                tickvals=list(range(len(admix_df))),
                ticktext=[str(idx)[:10] for idx in admix_df.index],
                tickangle=-45
            )
        )
        
        return [
            html.H4("ADMIXTURE Analysis"),
            dcc.Graph(figure=fig),
            html.Hr()
        ]
        
    except Exception as e:
        return html.Div(f"Error creating ADMIXTURE plot: {str(e)}")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)