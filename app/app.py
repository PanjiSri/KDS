import dash
from dash import html, dcc
from dash.dependencies import Input, Output

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(children=[
    html.H1(children='Zosterops borbonicus Population Structure Dashboard'),

    html.Div(children='''
        A web dashboard for analyzing population structure based on RADSeq data.
    '''),

    html.Div(id='upload-section', children=[
        html.H3("Upload Data Files"),
    ]),

    html.Div(id='visualization-section', children=[
        html.H3("Visualizations"),
        html.Div(id='map-plot-placeholder'),
        html.Div(id='pca-plot-placeholder'),
        html.Div(id='admixture-plot-placeholder'),
    ]),

])

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)