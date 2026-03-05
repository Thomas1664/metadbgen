import io
import pytest
import pandas as pd
from src import database

def test_prefer_refseq():
    df = pd.read_csv('tests/prefer-refseq.csv')
    res_df = database.prefer_refseq(df)
    assert len(res_df) == 5
    assert res_df['accession'].isin(['NC_001.1', 'AT001.1', 'AT002.1', 'SY_001.1', 'AT003.1']).all()

def test_sanitize_org_name():
    assert database.sanitize_org_name('Phytophthora taxon juncus') == 'Phytophthora juncus'
    
