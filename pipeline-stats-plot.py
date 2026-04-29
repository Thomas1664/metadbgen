import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors

BASE = 'data/cox2_th'

stage_files = {
    'Download': f'{BASE}/seq-meta.csv',
    'Trim': f'{BASE}/primer-trim/metadata.tsv',
    'Clean DB': f'{BASE}/seq-meta-filtered.tsv',
    'Dereplicate': f'{BASE}/dereplicate/metadata.tsv',
    'Quality Filter': f'{BASE}/quality-filter/metadata.tsv',
    'Dereplicate<br>(Alignment)': f'{BASE}/alignment/metadata.tsv',
    'Mislabel<br>Detection': f'{BASE}/output/cox2_th-meta-clean.tsv',
}

dfs = {}
for stage, path in stage_files.items():
    sep = '\t' if path.endswith('.tsv') else ','
    dfs[stage] = pd.read_csv(path, sep=sep)

stage_names = list(stage_files.keys())
n_stages = len(stage_names)

# Collect all families, sort by total sequence count descending so large families anchor the bottom
family_totals = (
    pd.concat([df[['family']] for df in dfs.values()])
    .groupby('family')
    .size()
    .sort_values(ascending=True)  # ascending so largest is at top of stack → legend order
)
all_families = family_totals.index.tolist()
n_families = len(all_families)

seq_mat = np.zeros((n_families, n_stages), dtype=int)
spe_mat = np.zeros((n_families, n_stages), dtype=int)

for j, stage in enumerate(stage_names):
    df = dfs[stage]
    seq_grp = df.groupby('family').size()
    spe_grp = df.groupby('family')['organism_name_sanitized'].nunique()
    for i, fam in enumerate(all_families):
        seq_mat[i, j] = seq_grp.get(fam, 0)
        spe_mat[i, j] = spe_grp.get(fam, 0)

# 21 families — use Light24 palette (24 colors)
palette = plotly.colors.qualitative.Light24
colors = [palette[i % len(palette)] for i in range(n_families)]

fig = go.Figure()

for i, (fam, color) in enumerate(zip(all_families, colors)):
    # Sequences bar (left bar in each stage group)
    fig.add_trace(go.Bar(
        name=fam,
        x=[stage_names, ['Sequences'] * n_stages],
        y=seq_mat[i],
        marker_color=color,
        legendgroup=fam,
        showlegend=True,
        hovertemplate=f'<b>{fam}</b><br>Stage: %{{x}}<br>Sequences: %{{y}}<extra></extra>',
    ))
    # Species bar (right bar in each stage group)
    fig.add_trace(go.Bar(
        name=fam,
        x=[stage_names, ['Species'] * n_stages],
        y=spe_mat[i],
        marker_color=color,
        legendgroup=fam,
        showlegend=False,
        hovertemplate=f'<b>{fam}</b><br>Stage: %{{x}}<br>Species: %{{y}}<extra></extra>',
    ))

# Annotate totals above each bar
seq_totals = seq_mat.sum(axis=0)
spe_totals = spe_mat.sum(axis=0)

annotations = []
for j, stage in enumerate(stage_names):
    annotations.append(dict(
        x=[stage, 'Sequences'],
        y=seq_totals[j],
        text=str(seq_totals[j]),
        showarrow=False,
        yanchor='bottom',
        font=dict(size=11),
    ))
    annotations.append(dict(
        x=[stage, 'Species'],
        y=spe_totals[j],
        text=str(spe_totals[j]),
        showarrow=False,
        yanchor='bottom',
        font=dict(size=11),
    ))

fig.update_layout(
    barmode='stack',
    title=dict(text='Sequences and Species per Family across Pipeline Stages', font=dict(size=16)),
    xaxis=dict(title='Pipeline Stage'),
    yaxis=dict(title='Count'),
    legend=dict(title='Family', traceorder='normal'),
    height=650,
    width=1300,
    annotations=annotations,
    plot_bgcolor='white',
    paper_bgcolor='white',
)
fig.update_xaxes(showgrid=False)
fig.update_yaxes(gridcolor='lightgrey')

fig.write_html(f'{BASE}/pipeline-stats.html')
print('Saved pipeline-stats.html')
