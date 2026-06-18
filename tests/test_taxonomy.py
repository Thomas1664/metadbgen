import pandas as pd
from src.taxonomy import update_species_name


def make_df(genera, species):
    return pd.DataFrame({'genus': genera, 'species': species})


def test_mixed_rows_only_bracketed_updated():
    df = make_df(
        ['Newgenus', 'Stablegenus'],
        ['[Oldgenus] epithet', 'Stablegenus epithet'],
    )
    result = update_species_name(df)
    assert result['species'].iloc[0] == 'Newgenus epithet'
    assert result['species'].iloc[1] == 'Stablegenus epithet'


def test_does_not_mutate_input():
    df = make_df(['Newgenus'], ['[Oldgenus] species'])
    original_species = df['species'].iloc[0]
    update_species_name(df)
    assert df['species'].iloc[0] == original_species
