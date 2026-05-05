import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.colors
from src.database import sanitize_org_name

parser = argparse.ArgumentParser()
parser.add_argument('--base', default='data/cox2_th', help='Base data directory')
parser.add_argument('--rank', default='family', help='Taxonomic rank to group by')
parser.add_argument('--out', default=None, help='Output file path (without extension); defaults to <base>/pipeline-stats')
args = parser.parse_args()

BASE = args.base
RANK = args.rank
OUT = args.out if args.out is not None else f'{BASE}/pipeline-stats'

stage_files = {
    'Download': f'{BASE}/seq-meta.csv',
    'Trim': f'{BASE}/primer-trim/metadata.tsv',
    'Clean DB': f'{BASE}/seq-meta-filtered.tsv',
    'Dereplicate': f'{BASE}/dereplicate/metadata.tsv',
    'Quality Filter': f'{BASE}/quality-filter/metadata.tsv',
    'Trim<br>(Alignment)': f'{BASE}/alignment/metadata-trim.tsv',
    'Dereplicate<br>(Alignment)': f'{BASE}/alignment/metadata.tsv',
    'Gappiness': f'{BASE}/gappiness/metadata.tsv',
    'Mislabel<br>Detection': f'{BASE}/output/cox2_th-meta-clean.tsv',
}

dfs = {}
for stage, path in stage_files.items():
    sep = '\t' if path.endswith('.tsv') else ','
    df = pd.read_csv(path, sep=sep)
    for _, row in df.iterrows():
        if pd.isna(row["family"]):
            row["family"] = row["genus"]
        if pd.isna(row["order"]):
            row["order"] = row["family"]
        if pd.isna(row["class"]):
            row["class"] = row["order"]
        if pd.isna(row["phylum"]):
            row["phylum"] = row["class"]
        if pd.isna(row["kingdom"]):
            row["kingdom"] = row["phylum"]
    dfs[stage] = df

stage_names = list(stage_files.keys())
n_stages = len(stage_names)

# Collect all families, sort by total sequence count descending so large families anchor the bottom
family_totals = (
    pd.concat([df[[RANK]] for df in dfs.values()])
    .groupby(RANK)
    .size()
    .sort_values(ascending=True)  # ascending so largest is at top of stack → legend order
)
all_families = family_totals.index.tolist()
n_families = len(all_families)

seq_mat = np.zeros((n_families, n_stages), dtype=int)
spe_mat = np.zeros((n_families, n_stages), dtype=int)
non_on_mat = np.zeros((n_families, n_stages), dtype=int)

for j, stage in enumerate(stage_names):
    df = dfs[stage]
    if 'organism_name_sanitized' not in df.columns:
        df = df.copy()
        df['organism_name_sanitized'] = df['organism_name'].apply(sanitize_org_name)
    seq_grp = df.groupby(RANK).size()
    spe_grp = df.groupby(RANK)['tax_id'].nunique()
    non_on_grp = df[df['organism_name'] == df['organism_name_sanitized']].groupby(RANK)['tax_id'].nunique()
    for i, fam in enumerate(all_families):
        seq_mat[i, j] = seq_grp.get(fam, 0)
        spe_mat[i, j] = spe_grp.get(fam, 0)
        non_on_mat[i, j] = non_on_grp.get(fam, 0)

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
    # Species bar (middle bar in each stage group)
    fig.add_trace(go.Bar(
        name=fam,
        x=[stage_names, ['Species'] * n_stages],
        y=spe_mat[i],
        marker_color=color,
        legendgroup=fam,
        showlegend=False,
        hovertemplate=f'<b>{fam}</b><br>Stage: %{{x}}<br>Species: %{{y}}<extra></extra>',
    ))
    # Non-ON species bar (right bar in each stage group)
    fig.add_trace(go.Bar(
        name=fam,
        x=[stage_names, ['Non-ON Species'] * n_stages],
        y=non_on_mat[i],
        marker_color=color,
        legendgroup=fam,
        showlegend=False,
        hovertemplate=f'<b>{fam}</b><br>Stage: %{{x}}<br>Non-ON Species: %{{y}}<extra></extra>',
    ))

# Annotate totals above each bar
seq_totals = seq_mat.sum(axis=0)
spe_totals = spe_mat.sum(axis=0)
non_on_totals = non_on_mat.sum(axis=0)

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
    annotations.append(dict(
        x=[stage, 'Non-ON Species'],
        y=non_on_totals[j],
        text=str(non_on_totals[j]),
        showarrow=False,
        yanchor='bottom',
        font=dict(size=11),
    ))

fig.update_layout(
    barmode='stack',
    title=dict(text=f'Sequences, Species and Non-ON Species per {RANK.title()} across Pipeline Stages', font=dict(size=16)),
    xaxis=dict(title='Pipeline Stage'),
    yaxis=dict(title='Count'),
    legend=dict(title=RANK.title(), traceorder='normal'),
    height=650,
    width=1300,
    annotations=annotations,
    plot_bgcolor='white',
    paper_bgcolor='white',
)
fig.update_xaxes(showgrid=False)
fig.update_yaxes(gridcolor='lightgrey')

fig.write_html(f'{OUT}.html')
fig.write_image(f'{OUT}.svg', format='svg')
print(f'Saved {OUT}.html and {OUT}.svg')
