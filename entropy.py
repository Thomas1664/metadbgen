import sys
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

tsv_path, html_out, svg_out = sys.argv[1:]

df = pd.read_csv(tsv_path, sep="\t")

window = 10
df["RollingMean"] = df["Entropy"].rolling(window, center=True, min_periods=1).mean()

fig = make_subplots(rows=2, cols=1, row_heights=[0.75, 0.25], shared_xaxes=True, vertical_spacing=0.08)

fig.add_trace(
    go.Bar(
        x=df["Site"],
        y=df["Entropy"],
        marker_color=df["Entropy"],
        marker_colorscale="Viridis",
        marker_showscale=True,
        marker_colorbar=dict(title="Entropy", len=0.7, y=0.65),
        name="Entropy",
        hovertemplate="Site %{x}<br>Entropy: %{y:.3f}<extra></extra>",
    ),
    row=1, col=1,
)

fig.add_trace(
    go.Scatter(
        x=df["Site"],
        y=df["RollingMean"],
        mode="lines",
        line=dict(color="red", width=2),
        name=f"Rolling mean (w={window})",
        hovertemplate="Site %{x}<br>Mean: %{y:.3f}<extra></extra>",
    ),
    row=1, col=1,
)

fig.add_trace(
    go.Heatmap(
        z=[df["Entropy"].values],
        x=df["Site"],
        colorscale="Viridis",
        showscale=False,
        hovertemplate="Site %{x}<br>Entropy: %{z:.3f}<extra></extra>",
        name="",
    ),
    row=2, col=1,
)

fig.update_layout(
    title="Shannon entropy along the alignment",
    xaxis2_title="Alignment position",
    yaxis_title="Shannon entropy (bits)",
    yaxis2=dict(showticklabels=False, showgrid=False),
    legend=dict(orientation="h", y=1.05),
    template="plotly_white",
    height=550,
)

fig.write_html(html_out)
fig.write_image(svg_out)
