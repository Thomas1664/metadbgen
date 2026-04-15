import pandas as pd
import numpy as np
from collections import Counter


def validate_at_rank(
    tax_df: pd.DataFrame,
    dist_mat_path: str,
    k: int = 10,
    max_radius: float = 0.05,
    min_neighbors: int = 5,
    confidence_threshold: float = 0.8,
    rank_level: str = 'genus',
) -> pd.DataFrame:
    label_lookup: dict[str, str] = dict(
        zip(tax_df['seq_id'].to_list(), tax_df[rank_level].to_list())
    )

    seq_ids: list[str] = []
    predictions: list[str | None] = []
    scores: list[float | None] = []
    is_mislabeled: list[bool] = []

    # Format: line 1 is the sequence count; each data row is tab-separated
    # with the seq_id as the first field followed by distances. There is no
    # column-header row — column order equals row order.
    #
    # First pass: collect all seq_ids so we can build the valid_mask before
    # processing any distances. Sequences absent from label_lookup are excluded
    # both as subjects AND as potential neighbours.
    with open(dist_mat_path) as f:
        f.readline()  # skip count line
        col_ids = np.array([line[:line.index('\t')] for line in f])

    valid_mask = np.array([sid in label_lookup for sid in col_ids])
    valid_col_ids = col_ids[valid_mask]
    # Cumulative count of valid entries; lets us find a row's position within
    # valid_col_ids in O(1) per row.
    valid_cumsum = np.cumsum(valid_mask)

    # Second pass: stream distances one row at a time
    with open(dist_mat_path) as f:
        f.readline()  # skip count line
        for row_idx, line in enumerate(f):
            if not valid_mask[row_idx]:
                continue
            seq_id = col_ids[row_idx]
            actual_label = label_lookup[seq_id]

            # Skip the leading seq_id field; parse distances directly
            distances = np.fromstring(line[line.index('\t') + 1:], dtype=np.float32, sep='\t')

            # Restrict to columns with labels, then exclude self
            valid_distances = distances[valid_mask]
            self_pos = int(valid_cumsum[row_idx]) - 1
            neighbor_mask = np.ones(len(valid_col_ids), dtype=bool)
            neighbor_mask[self_pos] = False

            candidate_distances = valid_distances[neighbor_mask]
            candidate_ids = valid_col_ids[neighbor_mask]

            # Filter by radius
            radius_mask = candidate_distances <= max_radius
            neighbor_distances = candidate_distances[radius_mask]
            neighbor_ids = candidate_ids[radius_mask]

            seq_ids.append(seq_id)

            if len(neighbor_ids) < min_neighbors:
                predictions.append(None)
                scores.append(None)
                is_mislabeled.append(False)
                continue

            top_k = np.argsort(neighbor_distances)[:k]
            neighbor_labels = [label_lookup[sid] for sid in neighbor_ids[top_k]]

            label_counts = Counter(neighbor_labels)
            predicted_label, vote_count = label_counts.most_common(1)[0]
            confidence = vote_count / len(neighbor_labels)

            predictions.append(predicted_label)
            scores.append(float(confidence))
            is_mislabeled.append(
                actual_label != predicted_label and confidence >= confidence_threshold
            )

    return pd.DataFrame({
        'seq_id': seq_ids,
        'predicted_label': predictions,
        'score': scores,
        'is_mislabeled': is_mislabeled,
    })
