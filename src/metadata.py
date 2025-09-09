import pandas as pd
import plotly.express as px
import pycountry


def parse_country(df: pd.DataFrame):
    df['country'] = df['geo_loc_name'].apply(lambda x: x.split(':', 1)[0] if type(x) == str else x)
    df['country'] = df['country'].str.replace('USA', 'United States of America')
    return df

def plot_country_heatmap(df: pd.DataFrame, html_out_file: str, svg_out_file: str):
    # Count occurrences
    country_counts = df["country"].value_counts().reset_index()
    country_counts.columns = ["country", "count"]

    # Plotly expects ISO Alpha-3 codes, so let's convert
    # Plotly can handle country names, but ISO3 is safer
    def country_to_iso3(name):
        try:
            return pycountry.countries.lookup(name).alpha_3
        except LookupError:
            return None

    country_counts["iso_alpha"] = country_counts["country"].apply(country_to_iso3)

    # Drop unknowns
    country_counts = country_counts.dropna(subset=["iso_alpha"])

    # Make interactive choropleth
    fig = px.choropleth(
        country_counts,
        locations="iso_alpha",
        color="count",
        hover_name="country",
        color_continuous_scale="Reds",
        projection="natural earth"
    )
    fig.write_html(html_out_file)
    fig.write_image(svg_out_file, 'svg')
