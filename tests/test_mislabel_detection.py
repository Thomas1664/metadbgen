import pytest
import pandas as pd
from src.mislabel_detection import validate_at_rank


def write_dist_matrix(path, seq_ids, matrix):
    """Write a distance matrix in gotree tab-separated format (count line + rows)."""
    n = len(seq_ids)
    with open(path, 'w') as f:
        f.write(f'{n}\n')
        for i, sid in enumerate(seq_ids):
            row = '\t'.join(f'{v:.6f}' for v in matrix[i])
            f.write(f'{sid}\t{row}\n')


def make_tax_df(seq_ids, genus_labels):
    return pd.DataFrame({'seq_id': seq_ids, 'genus': genus_labels})


# ---------------------------------------------------------------------------
# Basic classification
# ---------------------------------------------------------------------------

def test_correct_labels_not_flagged(tmp_path):
    """Sequences close to same-genus neighbours are not flagged."""
    seq_ids = ['A', 'B', 'C', 'D']
    matrix = [
        [0.00, 0.02, 0.15, 0.15],
        [0.02, 0.00, 0.15, 0.15],
        [0.15, 0.15, 0.00, 0.02],
        [0.15, 0.15, 0.02, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    tax_df = make_tax_df(['A', 'B', 'C', 'D'], ['GenuX', 'GenuX', 'GenuY', 'GenuY'])

    result = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.10,
                               min_neighbors=1, confidence_threshold=0.6)

    assert len(result) == 4
    assert not result['is_mislabeled'].any()


def test_mislabeled_sequence_flagged(tmp_path):
    """Sequence surrounded by different-genus neighbours is flagged as mislabeled."""
    seq_ids = ['A', 'B', 'C', 'D']
    # A is labeled GenuX but sits close to B, C, D which are all GenuY
    matrix = [
        [0.00, 0.02, 0.03, 0.03],
        [0.02, 0.00, 0.04, 0.04],
        [0.03, 0.04, 0.00, 0.02],
        [0.03, 0.04, 0.02, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    tax_df = make_tax_df(['A', 'B', 'C', 'D'], ['GenuX', 'GenuY', 'GenuY', 'GenuY'])

    result = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.10,
                               min_neighbors=2, confidence_threshold=0.6)

    a = result[result['seq_id'] == 'A'].iloc[0]
    assert a['is_mislabeled']
    assert a['predicted_label'] == 'GenuY'
    assert a['score'] == pytest.approx(1.0)

    # The correctly-labeled sequences must not be flagged
    others = result[result['seq_id'] != 'A']
    assert not others['is_mislabeled'].any()


# ---------------------------------------------------------------------------
# Not enough neighbours
# ---------------------------------------------------------------------------

def test_insufficient_neighbours_yields_none(tmp_path):
    """Sequence with fewer in-radius neighbours than min_neighbors gets None predictions."""
    seq_ids = ['A', 'B', 'C']
    # A is far from everything; B and C are close to each other
    matrix = [
        [0.00, 0.50, 0.50],
        [0.50, 0.00, 0.02],
        [0.50, 0.02, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    tax_df = make_tax_df(['A', 'B', 'C'], ['GenuX', 'GenuX', 'GenuX'])

    result = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.10,
                               min_neighbors=2, confidence_threshold=0.6)

    a = result[result['seq_id'] == 'A'].iloc[0]
    assert pd.isna(a['predicted_label'])
    assert pd.isna(a['score'])
    assert not a['is_mislabeled']


# ---------------------------------------------------------------------------
# Self-exclusion
# ---------------------------------------------------------------------------

def test_self_not_counted_as_neighbour(tmp_path):
    """The diagonal (distance 0) must not be used as a neighbour."""
    seq_ids = ['A', 'B']
    # Neither A nor B is within max_radius of the other
    matrix = [
        [0.00, 0.50],
        [0.50, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    tax_df = make_tax_df(['A', 'B'], ['GenuX', 'GenuX'])

    # If self (distance 0) were counted, min_neighbors=1 would pass — it must not
    result = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.10,
                               min_neighbors=1, confidence_threshold=0.6)

    for _, row in result.iterrows():
        assert pd.isna(row['predicted_label'])
        assert pd.isna(row['score'])
        assert not row['is_mislabeled']


# ---------------------------------------------------------------------------
# Sequences missing from lookup
# ---------------------------------------------------------------------------

def test_missing_from_lookup_not_in_output(tmp_path):
    """Sequences absent from tax_df must not appear as rows in the output."""
    seq_ids = ['A', 'B', 'C']
    matrix = [
        [0.00, 0.02, 0.02],
        [0.02, 0.00, 0.02],
        [0.02, 0.02, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    # C is absent from tax_df
    tax_df = make_tax_df(['A', 'B'], ['GenuX', 'GenuX'])

    result = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.10,
                               min_neighbors=1, confidence_threshold=0.6)

    assert set(result['seq_id'].tolist()) == {'A', 'B'}


def test_missing_from_lookup_not_counted_as_neighbour(tmp_path):
    """Sequences absent from tax_df must not occupy neighbour slots."""
    # E is in the matrix but not in tax_df.
    # A is close to B, C, D, and E (all within radius).
    # With correct exclusion A has 3 labelled neighbours; without it, 4.
    seq_ids = ['A', 'B', 'C', 'D', 'E']
    matrix = [
        [0.00, 0.01, 0.01, 0.01, 0.005],
        [0.01, 0.00, 0.01, 0.01, 0.005],
        [0.01, 0.01, 0.00, 0.01, 0.005],
        [0.01, 0.01, 0.01, 0.00, 0.005],
        [0.005, 0.005, 0.005, 0.005, 0.00],
    ]
    dist_file = tmp_path / 'dist.tsv'
    write_dist_matrix(dist_file, seq_ids, matrix)
    tax_df = make_tax_df(['A', 'B', 'C', 'D'], ['GenuX', 'GenuX', 'GenuX', 'GenuX'])

    # min_neighbors=4 — only reachable if E (unlabelled) were counted
    result_strict = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.05,
                                      min_neighbors=4, confidence_threshold=0.6)
    a_strict = result_strict[result_strict['seq_id'] == 'A'].iloc[0]
    assert pd.isna(a_strict['predicted_label']), \
        "E must not count as a neighbour — A should have only 3 labelled neighbours"

    # min_neighbors=3 — exactly 3 labelled neighbours available
    result_lenient = validate_at_rank(tax_df, str(dist_file), k=10, max_radius=0.05,
                                       min_neighbors=3, confidence_threshold=0.6)
    a_lenient = result_lenient[result_lenient['seq_id'] == 'A'].iloc[0]
    assert a_lenient['predicted_label'] == 'GenuX'
    assert not a_lenient['is_mislabeled']