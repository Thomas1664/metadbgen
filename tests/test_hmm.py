import pandas as pd
from src.hmm import (
    TBLOUT_COLUMNS,
    compute_score_threshold,
    filter_metadata_by_accession,
    parse_hmm_tblout,
    passing_hits,
    read_metadata,
    read_threshold_file,
    split_by_strand,
    write_threshold_file,
)


def write_tblout(path, rows):
    """Write a synthetic nhmmer --tblout file.

    Each row is a dict overriding a set of realistic defaults for the 16
    whitespace-separated tblout columns; only the fields a test cares about
    need to be specified.
    """
    defaults = {
        "target_name": "seq1", "target_accession": "-", "query_name": "profile",
        "query_accession": "-", "hmm_from": 1, "hmm_to": 100, "ali_from": 1,
        "ali_to": 100, "env_from": 1, "env_to": 100, "sq_len": 100,
        "strand": "+", "e_value": "1e-50", "score": 100.0, "bias": 0.0,
        "description": "-",
    }
    with open(path, 'w') as f:
        f.write("# synthetic tblout fixture\n")
        for row in rows:
            merged = {**defaults, **row}
            f.write(' '.join(str(merged[c]) for c in TBLOUT_COLUMNS) + '\n')


def test_parse_hmm_tblout_keeps_best_score_per_target(tmp_path):
    """A target hit on both strands should be reduced to its single best-scoring row."""
    path = tmp_path / 'scores.tsv'
    write_tblout(path, [
        {"target_name": "seqA", "strand": "-", "score": 30.0},
        {"target_name": "seqA", "strand": "+", "score": 55.0},
        {"target_name": "seqB", "strand": "+", "score": 10.0},
    ])

    df = parse_hmm_tblout(str(path))

    assert sorted(df['target_name']) == ['seqA', 'seqB']
    row = df[df['target_name'] == 'seqA'].iloc[0]
    assert row['score'] == 55.0
    assert row['strand'] == '+'


def test_compute_score_threshold(tmp_path):
    path = tmp_path / 'refseq-scores.tsv'
    write_tblout(path, [
        {"target_name": "ref1", "score": 70.0},
        {"target_name": "ref2", "score": 50.0},
        {"target_name": "ref3", "score": 60.0},
    ])

    min_score, threshold = compute_score_threshold(str(path), margin=10.0)

    assert min_score == 50.0
    assert threshold == 40.0


def test_threshold_file_roundtrip(tmp_path):
    path = tmp_path / 'score-threshold.txt'
    write_threshold_file(str(path), 42.5)

    assert read_threshold_file(str(path)) == 42.5


def test_passing_hits_filters_by_threshold(tmp_path):
    path = tmp_path / 'query-scores.tsv'
    write_tblout(path, [
        {"target_name": "good1", "score": 60.0},
        {"target_name": "good2", "score": 40.0},
        {"target_name": "bad1", "score": 39.9},
    ])

    good = passing_hits(str(path), threshold=40.0)

    assert sorted(good['target_name']) == ['good1', 'good2']


def test_split_by_strand():
    hits = pd.DataFrame({
        'target_name': ['a', 'b', 'c'],
        'strand': ['+', '-', '+'],
    })

    plus, minus = split_by_strand(hits)

    assert sorted(plus) == ['a', 'c']
    assert sorted(minus) == ['b']


def test_read_metadata_detects_csv_separator(tmp_path):
    path = tmp_path / 'meta.csv'
    path.write_text('accession,species\nACC1,Foo bar\nACC2,Baz qux\n')

    df = read_metadata(str(path))

    assert list(df['accession']) == ['ACC1', 'ACC2']


def test_read_metadata_detects_tsv_separator(tmp_path):
    path = tmp_path / 'meta.tsv'
    path.write_text('accession\tspecies\nACC1\tFoo bar\nACC2\tBaz qux\n')

    df = read_metadata(str(path))

    assert list(df['accession']) == ['ACC1', 'ACC2']


def test_filter_metadata_by_accession():
    meta_df = pd.DataFrame({'accession': ['ACC1', 'ACC2', 'ACC3'], 'value': [1, 2, 3]})

    filtered = filter_metadata_by_accession(meta_df, ['ACC1', 'ACC3'])

    assert sorted(filtered['accession']) == ['ACC1', 'ACC3']
