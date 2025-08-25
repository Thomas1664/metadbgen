import pandas as pd
from sklearn.manifold import TSNE
import plotly.express as px
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import sys
from io import StringIO


dist_mat, meta, plot_html, plot_svg = sys.argv[1:]


if __name__ == '__main__':
    with open(dist_mat, 'r') as f:
        buffer = f.read()
        buffer = '\n'.join(buffer.split('\n', 1)[1:])
        buffer = buffer.replace('  ', ' ')
        buffer = buffer.replace('  ', ' ')
        buffer = buffer.replace(' \n', '\n')
        with open('dist-mat-test.csv', 'w') as o:
            o.write(buffer)
    
    mat = pd.read_csv(StringIO(buffer), sep=' ', index_col=0, header=None)
    mat.columns = mat.index

    df = pd.read_csv(meta, index_col='accession', sep='\t')
    merged = df[df.index.isin(mat.index)]

    tsne = TSNE(metric='precomputed', perplexity=30,  init='random', random_state=42)
    coords = tsne.fit_transform(mat.to_numpy())

    # DataFrame for plotting
    df_plot = pd.DataFrame({
        'Accession': merged.index,
        'Family': merged['family'],
        'Species': merged['organism_name_sanitized'],
        't-SNE1': coords[:, 0],
        't-SNE2': coords[:, 1]
    })

    # Custom colors
    unique_families = df_plot['Family'].unique()
    palette = px.colors.qualitative.Light24[:len(unique_families)]
    hex_colors = [mcolors.to_hex(color) for color in palette]
    color_discrete_map = {family: hex_colors[i] for i, family in enumerate(unique_families)}
    color_discrete_map['Lagenidiales incertae sedis'] = '#E0E0E0'
    color_discrete_map['Anisolpidiaceae'] = '#611800'
    color_discrete_map['Olpidiopsidaceae'] = '#BDE500'


    # Plot
    fig = px.scatter(
        df_plot,
        x='t-SNE1',
        y='t-SNE2',
        color='Family',
        hover_name='Accession',
        hover_data={'Species': True, 't-SNE1': False, 't-SNE2': False},
        title=None,
        color_discrete_map=color_discrete_map,
        width=800,
        height=600,
    )

    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='black')))
    fig.update_layout(
        legend=dict(
            itemsizing='trace'
        ),
        margin=dict(b=120)  # Make room at the bottom for the legend
        )
    fig.write_image(plot_svg, format='svg')
    fig.write_html(plot_html)
