import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import dash
from dash import dcc, html, dash_table, ctx
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from sqlalchemy import create_engine
import dash_bio as dashbio
import logging

# --- Logging Configuration ---
logging.basicConfig(level=logging.ERROR)

# --- Database Configuration ---
DEFAULT_DB_PATH = "quantms.db"

# --- Initialize Dash App ---
app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server

# --- Custom CSS Styles ---
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Custom styles
CUSTOM_STYLES = {
    'header': {
        'backgroundColor': '#2c3e50',
        'color': 'white',
        'padding': '1.5rem',
        'marginBottom': '2rem',
        'borderRadius': '4px'
    },
    'input': {
        'padding': '8px',
        'borderRadius': '4px',
        'border': '1px solid #ddd',
        'width': '100%'
    },
    'label': {
        'fontWeight': '600',
        'marginBottom': '4px',
        'color': '#2c3e50'
    },
    'card': {
        'backgroundColor': 'white',
        'borderRadius': '4px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
        'padding': '1.5rem',
        'marginBottom': '2rem'
    },
    'graph-container': {
        'backgroundColor': 'white',
        'borderRadius': '4px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
        'padding': '1rem',
        'marginBottom': '2rem'
    },
    'table-container': {
        'backgroundColor': 'white',
        'borderRadius': '4px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
        'padding': '1rem',
        'marginBottom': '2rem'
    }
}

# --- Helper Function ---
def create_empty_figure(title: str = "No data available"):
    return {
        'data': [],
        'layout': {
            'xaxis': {'visible': False},
            'yaxis': {'visible': False},
            'annotations': [{
                'text': title,
                'xref': 'paper',
                'yref': 'paper',
                'showarrow': False,
                'font': {'size': 16}
            }],
            'template': 'plotly_white',
            'plot_bgcolor': 'rgba(0,0,0,0)',
            'paper_bgcolor': 'rgba(0,0,0,0)',
        }
    }

def create_small_figure(peptide, aa_before, aa_after, start, end):
    # Wrap the peptide sequence every 30 characters for better display
    def wrap_sequence(seq, width=30):
        return '<br>'.join([seq[i:i+width] for i in range(0, len(seq), width)])

    wrapped_peptide = wrap_sequence(peptide)
    fig = go.Figure()
    fig.add_shape(
        type="rect",
        x0=0.5, y0=0.7, x1=3.5, y1=1.3,
        fillcolor="lightblue",
        opacity=0.2,
        line_width=0
    )
    fig.add_trace(go.Scatter(
        x=[1, 3],
        y=[1, 1],
        mode='markers+text',
        text=[f"Pos: {start}", f"Pos: {end}"],
        textposition="bottom center",
        marker=dict(
            size=25,
            color=['#2ca02c', '#d62728'],
            line=dict(width=2, color='DarkSlateGrey')
        ),
        textfont=dict(size=12, color='black'),
        hoverinfo='text',
        hovertext=f"Change: {aa_before}→{aa_after}",
        showlegend=False
    ))
    fig.add_annotation(
        x=3, y=1, ax=1, ay=1,
        xref="x", yref="y",
        axref="x", ayref="y",
        text="",
        showarrow=True,
        arrowhead=3,
        arrowsize=1.5,
        arrowwidth=2,
        arrowcolor='#636363'
    )
    fig.add_annotation(
        x=2, y=1.2,
        text=f"<b>{wrapped_peptide}</b>",
        showarrow=False,
        font=dict(size=9, color='black'),
        align='center'
    )
    fig.add_annotation(
        x=2, y=0.85,
        text=f"<b>{aa_before} → {aa_after}</b>",
        showarrow=False,
        font=dict(size=14, color='#d62728'),
        align='center'
    )
    fig.add_annotation(
        x=2, y=1.5,
        text="<b>Peptide Modification</b>",
        showarrow=False,
        font=dict(size=16, color='black'),
        align='center'
    )
    fig.update_layout(
        xaxis=dict(visible=False, range=[0, 4]),
        yaxis=dict(visible=False, range=[0.5, 1.6]),
        margin=dict(l=20, r=20, t=40, b=20),
        height=200,
        width=400,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        hovermode='closest'
    )
    return fig

# --- App Layout ---
app.layout = html.Div([
    # Header
    html.Div([
        html.H1("Proteomics Data Visualization", style={'marginBottom': '0.5rem'}),
        html.P("Interactive visualization of mass spectrometry data",
              style={'color': '#bdc3c7', 'marginBottom': '0'}),
    ], style=CUSTOM_STYLES['header']),

    # Main content container
    html.Div([
        # Sample selection
        html.Div([
            html.Div([
                html.Label("Sample ID", style=CUSTOM_STYLES['label']),
                dcc.Dropdown(
                    id='sample-id',
                    options=[],  # Will be populated dynamically
                    value=None,  # Will be set by callback
                    style=CUSTOM_STYLES['input']
                ),
            ], className="twelve columns", style={'padding': '0 15px'}),
        ], className="row", style={'marginBottom': '2rem'}),

        # First row of graphs
        html.Div([
            html.Div([
                dcc.Graph(
                    id='mz-rt',
                    style={'height': '400px'},
                    config={'displayModeBar': True}
                ),
            ], className="six columns", style=CUSTOM_STYLES['graph-container']),

            html.Div([
                dcc.Graph(
                    id='rt-intensity',
                    style={'height': '400px'},
                    config={'displayModeBar': True}
                ),
            ], className="six columns", style=CUSTOM_STYLES['graph-container']),
        ], className="row"),

        # Second row of graphs
        html.Div([
            html.Div([
                dcc.Graph(
                    id='peptide-length',
                    style={'height': '400px'},
                    config={'displayModeBar': True}
                ),
            ], className="six columns", style=CUSTOM_STYLES['graph-container']),

            html.Div([
                dcc.Graph(
                    id='charge-distribution',
                    style={'height': '400px'},
                    config={'displayModeBar': True}
                ),
            ], className="six columns", style=CUSTOM_STYLES['graph-container']),
        ], className="row"),

        # Filters above the table
        html.Div([
            html.Div([
                html.Label("MSGF+ Score Cutoff", style=CUSTOM_STYLES['label']),
                dcc.Input(
                    id='msgf-cutoff',
                    type='number',
                    value=-10,
                    style=CUSTOM_STYLES['input']
                ),
            ], className="three columns", style={'padding': '0 15px'}),

            html.Div([
                html.Label("Percolator Score Cutoff", style=CUSTOM_STYLES['label']),
                dcc.Input(
                    id='percolator-cutoff',
                    type='number',
                    value=-10,
                    style=CUSTOM_STYLES['input']
                ),
            ], className="three columns", style={'padding': '0 15px'}),

            html.Div([
                html.Label("Q-value Cutoff", style=CUSTOM_STYLES['label']),
                dcc.Input(
                    id='qvalue-cutoff',
                    type='number',
                    value=0.05,
                    step=0.01,
                    style=CUSTOM_STYLES['input']
                ),
            ], className="three columns", style={'padding': '0 15px'}),

            html.Div([
                html.Label("Update Data", style={'visibility': 'hidden'}),
                html.Button(
                    'Apply Filters',
                    id='submit-button',
                    n_clicks=0,
                    style={
                        'width': '100%',
                        'padding': '10px',
                        'backgroundColor': '#3498db',
                        'color': 'white',
                        'border': 'none',
                        'borderRadius': '4px',
                        'cursor': 'pointer'
                    }
                )
            ], className="three columns", style={'padding': '0 15px'}),
        ], className="row", style={'marginBottom': '2rem'}),

        # Data table
        html.Div([
            dash_table.DataTable(
                id='data-table',
                columns=[],
                data=[],
                page_size=10,
                filter_action='native',
                sort_action='native',
                export_format='csv',
                row_selectable='multi',
                style_table={'overflowX': 'auto'},
                style_header={
                    'backgroundColor': '#f8f9fa',
                    'fontWeight': 'bold',
                    'border': '1px solid #dee2e6'
                },
                style_cell={
                    'padding': '8px',
                    'textAlign': 'left',
                    'border': '1px solid #dee2e6',
                    'fontFamily': 'Arial, sans-serif',
                    'minWidth': '80px', 'width': '120px', 'maxWidth': '180px',
                    'whiteSpace': 'normal',
                    'height': 'auto'
                },
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    }
                ],
                tooltip_data=[],
                tooltip_duration=None
            ),
        ], style=CUSTOM_STYLES['table-container']),

        # Sequence viewer section
        html.Div([
            html.H3("Selected Peptide Sequence Viewer", style={'marginBottom': '1rem', 'color': '#2c3e50'}),
            html.Div(
                id='sequence-viewer',
                style={
                    'padding': '20px',
                    'backgroundColor': 'white',
                    'borderRadius': '4px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                }
            )
        ]),
    ], style={'padding': '0 2rem', 'maxWidth': '1400px', 'margin': '0 auto'})
])

# --- Callback to Update Graphs and Table ---
@app.callback(
    Output('mz-rt', 'figure'),
    Output('rt-intensity', 'figure'),
    Output('peptide-length', 'figure'),
    Output('charge-distribution', 'figure'),
    Output('data-table', 'columns'),
    Output('data-table', 'data'),
    Output('data-table', 'tooltip_data'),
    Input('submit-button', 'n_clicks'),
    Input('sample-id', 'value'),
    Input('msgf-cutoff', 'value'),
    Input('percolator-cutoff', 'value'),
    Input('qvalue-cutoff', 'value')
)
def update_graphs(n_clicks, sample_id, msgf_cutoff, percolator_cutoff, qvalue_cutoff):
    if n_clicks is None:
        raise PreventUpdate

    try:
        engine = create_engine(f'sqlite:///{DEFAULT_DB_PATH}')
        df_combined = pd.read_sql("SELECT * FROM combined_score", engine)
        df = df_combined[df_combined['sample_name'] == sample_id].copy()

        if df.empty:
            empty = create_empty_figure("No data for sample")
            return empty, empty, empty, empty, [], [], []

        df_filtered = df.copy()

        try:
            df_filtered = df_filtered[df_filtered['msgfplus_score'] >= float(msgf_cutoff)]
        except:
            pass
        try:
            df_filtered = df_filtered[df_filtered['percolator_score'] >= float(percolator_cutoff)]
        except:
            pass
        try:
            df_filtered = df_filtered[df_filtered['qvalue_score'] <= float(qvalue_cutoff)]
        except:
            pass

        # Create mz vs rt plot with professional styling
        fig_mz_rt = px.scatter(
            df_filtered,
            x='rt',
            y='mz',
            color='charge',
            title='m/z vs Retention Time',
            labels={'rt': 'Retention Time (min)', 'mz': 'm/z'},
            hover_data=['sequence', 'charge', 'msgfplus_score']
        )
        fig_mz_rt.update_layout(
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            font=dict(family="Arial, sans-serif", size=12, color="#2c3e50"),
            title_font=dict(size=18, color="#2c3e50"),
            legend_title_text='Charge State',
            margin=dict(l=50, r=50, t=60, b=50)
        )
        fig_mz_rt.update_traces(
            marker=dict(size=8, opacity=0.7, line=dict(width=0.5, color='DarkSlateGrey')),
            selector=dict(mode='markers')
        )

        # Create peptide length distribution plot
        df_filtered['sequence'] = df_filtered['sequence'].astype(str)
        df_filtered['peptide_length'] = df_filtered['sequence'].str.len()
        fig_pep_len = px.violin(
            df_filtered,
            y='peptide_length',
            title='Peptide Length Distribution',
            labels={'peptide_length': 'Peptide Length (AA)'}
        )
        fig_pep_len.update_layout(
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            font=dict(family="Arial, sans-serif", size=12, color="#2c3e50"),
            title_font=dict(size=18, color="#2c3e50"),
            margin=dict(l=50, r=50, t=60, b=50)
        )
        fig_pep_len.update_traces(
            fillcolor='#3498db',
            line_color='#2c3e50',
            opacity=0.7
        )

        # Create charge state distribution plot
        df_filtered['charge'] = pd.to_numeric(df_filtered['charge'], errors='coerce')
        fig_charge = go.Figure()
        fig_charge.add_trace(go.Histogram(
            x=df_filtered['charge'],
            name='Charge States',
            marker_color='#3498db',
            opacity=0.7
        ))
        fig_charge.update_layout(
            title='Charge State Distribution',
            xaxis_title='Charge',
            yaxis_title='Frequency',
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            font=dict(family="Arial, sans-serif", size=12, color="#2c3e50"),
            title_font=dict(size=18, color="#2c3e50"),
            margin=dict(l=50, r=50, t=60, b=50),
            bargap=0.1
        )

        # Create RT vs Intensity plot
        df_ms_info = pd.read_sql(f"SELECT * FROM ms_info WHERE sample_name = '{sample_id}'", engine)
        fig_rt_intensity = px.scatter(
            df_ms_info,
            x='Retention_Time',
            y='Base_Peak_Intensity',
            color='Charge',
            title='Base Peak Intensity vs Retention Time',
            labels={'Retention_Time': 'Retention Time (min)', 'Base_Peak_Intensity': 'Base Peak Intensity'},
            hover_data=['MSLevel', 'Charge']
        )
        fig_rt_intensity.update_layout(
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            font=dict(family="Arial, sans-serif", size=12, color="#2c3e50"),
            title_font=dict(size=18, color="#2c3e50"),
            legend_title_text='Charge State',
            margin=dict(l=50, r=50, t=60, b=50)
        )
        fig_rt_intensity.update_traces(
            marker=dict(size=8, opacity=0.7, line=dict(width=0.5, color='DarkSlateGrey')),
            selector=dict(mode='markers')
        )

        # Prepare data table
        columns = [{"name": i, "id": i} for i in df_filtered.columns]
        data = df_filtered.to_dict('records')

        # Create tooltips for the table
        tooltip_data = [{
            column: {'value': str(value), 'type': 'markdown'}
            for column, value in row.items()
        } for row in data]

        return fig_mz_rt, fig_rt_intensity, fig_pep_len, fig_charge, columns, data, tooltip_data

    except Exception as e:
        logging.error(f"Error in callback: {e}")
        empty = create_empty_figure("Error loading data")
        return empty, empty, empty, empty, [], [], []

# --- Callback to Display Sequence Viewer ---
@app.callback(
    Output('sequence-viewer', 'children'),
    Input('data-table', 'derived_virtual_data'),
    Input('data-table', 'derived_virtual_selected_rows')
)
def display_selected_sequence(rows, selected_rows):
    if not rows or not selected_rows:
        return html.Div([
            html.P("Select one or more rows in the table to view peptide sequences.",
                  style={'color': '#7f8c8d', 'textAlign': 'center'})
        ], style={'padding': '2rem'})

    try:
        figures = []

        for idx, i in enumerate(selected_rows):
            row = rows[i]
            seq = str(row.get('sequence', ''))
            aa_before = str(row.get('aa_before', 'X'))
            aa_after = str(row.get('aa_after', 'Y'))
            start = str(row.get('start', 1))
            end = str(row.get('end', len(seq)))

            # Safely format the qvalue score
            qvalue_score = row.get('qvalue_score', 'N/A')
            try:
                qvalue_display = f"{float(qvalue_score):.4f}" if qvalue_score != 'N/A' else 'N/A'
            except:
                qvalue_display = str(qvalue_score)

            # Extract and format protein name from 'accessions'
            protein = row.get('accessions', 'N/A')
            if not protein or not isinstance(protein, str) or protein.strip() == '':
                protein = 'N/A'
            else:
                protein = protein.strip()
                if ';' in protein:
                    protein = protein.split(';')[0].strip()
                if '|' in protein:
                    parts = protein.split('|')
                    if len(parts) > 2 and parts[2]:
                        protein = parts[2].strip()
                    else:
                        protein = parts[-1].strip() if len(parts) > 1 else parts[0].strip()
                if not protein:
                    protein = 'N/A'

            figures.append(html.Div([
                html.Div([
                    html.H5(f"Sequence {idx + 1}", style={"marginBottom": "10px", "color": "#2c3e50"}),
                    html.Div([
                        html.Span("Protein: ", style={"fontWeight": "bold"}),
                        html.Span(protein),
                        html.Span(" | ", style={"margin": "0 5px"}),
                        html.Span("Score: ", style={"fontWeight": "bold"}),
                        html.Span(str(row.get('msgfplus_score', 'N/A'))),
                        html.Span(" | ", style={"margin": "0 5px"}),
                        html.Span("Q-value: ", style={"fontWeight": "bold"}),
                        html.Span(qvalue_display)
                    ], style={"marginBottom": "15px", "color": "#34495e"})
                ], style={"marginBottom": "15px"}),

                html.Div([
                    dashbio.SequenceViewer(
                        id=f'sequence-viewer-{idx}',
                        sequence=seq,
                        showLineNumbers=True,
                        wrapAminoAcids=True,
                        charsPerLine=60,
                        coverage=[{
                            'start': int(start) if start.isdigit() else 1,
                            'end': int(end) if end.isdigit() else len(seq),
                            'color': '#e74c3c',
                            'bgcolor': '#fadbd8'
                        }]
                    ),
                ], style={"marginBottom": "15px"}),

                dcc.Graph(
                    figure=create_small_figure(seq, aa_before, aa_after, start, end),
                    config={'displayModeBar': False}
                )
            ], style={
                'marginBottom': '30px',
                'padding': '20px',
                'backgroundColor': 'white',
                'borderRadius': '4px',
                'boxShadow': '0 1px 3px rgba(0,0,0,0.1)'
            }))

        return figures

    except Exception as e:
        logging.error(f"Error displaying sequence: {e}")
        return html.Div([
            html.P("An error occurred while displaying sequences.",
                  style={'color': '#e74c3c', 'textAlign': 'center'})
        ], style={'padding': '2rem'})

# --- Callback to populate sample filter dynamically ---
@app.callback(
    Output('sample-id', 'options'),
    Output('sample-id', 'value'),
    Input('sample-id', 'id')  # Dummy input to trigger on page load
)
def populate_sample_dropdown(_):
    try:
        engine = create_engine(f'sqlite:///{DEFAULT_DB_PATH}')
        df = pd.read_sql("SELECT DISTINCT sample_name FROM combined_score", engine)
        options = [
            {'label': str(s), 'value': str(s)}
            for s in df['sample_name'].dropna().unique()
        ]
        value = options[0]['value'] if options else None
        return options, value
    except Exception as e:
        logging.error(f"Error loading sample names: {e}")
        return [], None

if __name__ == '__main__':
    app.run(debug=True, port=8051)