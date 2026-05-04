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
    assert database.sanitize_org_name('Phytophthora x juncus') == 'Phytophthora juncus'
    assert database.sanitize_org_name('Phytophthora cf. juncus') == 'Phytophthora juncus'
    assert database.sanitize_org_name('Phytophthora aff. juncus') == 'Phytophthora juncus'
    assert database.sanitize_org_name('Phytophthora v juncus') == 'Phytophthora juncus'
    assert database.sanitize_org_name('Phytophthora irregulare juncus') == 'Phytophthora irregulare'
    assert database.sanitize_org_name('Phytophthora juncus') == 'Phytophthora juncus'
    assert database.sanitize_org_name('Phytophthora') == 'Phytophthora'


def _make_clust_results(rows):
    """rows: list of (type, cluster_id, query_length, similarity, strand, align_res, query_label)"""
    return pd.DataFrame(rows, columns=['type', 'cluster_id', 'query_length', 'similarity', 'strand', 'align_res', 'query_label'])

def _make_metadata(rows, species_col=None):
    """rows: list of (accession, tax_id, trimed_seq_len)"""
    df = pd.DataFrame(rows, columns=['accession', 'tax_id', 'trimed_seq_len'])
    df['organism_name'] = species_col if species_col is not None else 'Genus species'
    df['organism_name_sanitized'] = df['organism_name'].apply(database.sanitize_org_name)
    return df

def test_dereplicate_keeps_longest_sequence():
    # Two sequences same cluster + tax_id: ACC1 (len=500) and ACC2 (len=300) → ACC1 kept
    metadata = _make_metadata([
        ('ACC1', 1, 500),
        ('ACC2', 1, 300),
    ])
    clust = _make_clust_results([
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC1'),
        ('H', 'C1', 300, 1.0, '+', '8I', 'ACC2'),
    ])
    result = database.dereplicate(clust, metadata)
    assert list(result['accession']) == ['ACC1']

def test_dereplicate_no_duplicate_cluster_tax_id():
    metadata = _make_metadata([
        ('ACC1', 1, 500),
        ('ACC2', 1, 400),
        ('ACC3', 2, 300),
        ('ACC4', 2, 200),
        ('ACC5', 1, 500),
    ])
    clust = _make_clust_results([
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC1'),
        ('H', 'C1', 400, 1.0, '+', '8I', 'ACC2'),
        ('H', 'C1', 300, 1.0, '+', '8I', 'ACC3'),
        ('H', 'C1', 200, 1.0, '+', '8I', 'ACC4'),
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC5'),
    ])
    result = database.dereplicate(clust, metadata)
    assert result.duplicated(subset=['cluster_id', 'tax_id']).sum() == 0

def test_dereplicate_multiple_clusters_and_tax_ids():
    # C1/tax1: ACC1(500) > ACC2(300) → ACC1
    # C1/tax2: ACC3(400)             → ACC3
    # C2/tax1: ACC4(600)             → ACC4
    metadata = _make_metadata([
        ('ACC1', 1, 500),
        ('ACC2', 1, 300),
        ('ACC3', 2, 400),
        ('ACC4', 1, 600),
    ])
    clust = _make_clust_results([
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC1'),
        ('H', 'C1', 300, 1.0, '+', '8I', 'ACC2'),
        ('H', 'C1', 400, 1.0, '+', '8I', 'ACC3'),
        ('H', 'C2', 600, 1.0, '+', '8I', 'ACC4'),
    ])
    result = database.dereplicate(clust, metadata)
    assert set(result['accession']) == {'ACC1', 'ACC3', 'ACC4'}

def test_dereplicate_single_row_cluster():
    metadata = _make_metadata([('ACC1', 1, 500)])
    clust = _make_clust_results([('S', 'C1', 500, 1.0, '+', '8I', 'ACC1')])
    result = database.dereplicate(clust, metadata)
    assert len(result) == 1
    assert result['accession'].iloc[0] == 'ACC1'

def test_dereplicate_drops_columns():
    metadata = _make_metadata([('ACC1', 1, 500)])
    clust = _make_clust_results([('S', 'C1', 500, 1.0, '+', '8I', 'ACC1')])
    result = database.dereplicate(clust, metadata)
    dropped = {'type', 'query_length', 'similarity', 'align_res', 'query_label', 'strand'}
    assert dropped.isdisjoint(result.columns)

def test_dereplicate_removes_on_when_non_on_in_same_cluster():
    # ACC1 is non-ON, ACC2 is ON; same cluster → ACC2 removed even though it's longer
    metadata = _make_metadata(
        [('ACC1', 1, 500), ('ACC2', 2, 600)],
        species_col=['Genus species', 'Genus aff. species'],
    )
    clust = _make_clust_results([
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC1'),
        ('H', 'C1', 600, 1.0, '+', '8I', 'ACC2'),
    ])
    result = database.dereplicate(clust, metadata)
    assert list(result['accession']) == ['ACC1']

def test_dereplicate_keeps_on_when_no_non_on_in_cluster():
    # Both sequences are ON; no non-ON in their cluster → keep them, dedup by tax_id keeps longest
    metadata = _make_metadata(
        [('ACC1', 1, 500), ('ACC2', 1, 300)],
        species_col=['Genus aff. species', 'Genus cf. species'],
    )
    clust = _make_clust_results([
        ('H', 'C1', 500, 1.0, '+', '8I', 'ACC1'),
        ('H', 'C1', 300, 1.0, '+', '8I', 'ACC2'),
    ])
    result = database.dereplicate(clust, metadata)
    assert list(result['accession']) == ['ACC1']
