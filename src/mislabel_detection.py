import pandas as pd
import numpy as np
import dendropy
from collections import Counter


def _tree_distances(nwk_data: str) -> tuple[np.ndarray, list[str]]:
    tree: dendropy.Tree = dendropy.Tree.get(data=nwk_data, schema='newick')
    pdm = tree.phylogenetic_distance_matrix()
    taxa = tree.taxon_namespace
    n = len(taxa)
    dist_matrix = np.zeros((n, n))
    taxon_list = list(taxa)
    
    for i, taxon1 in enumerate(taxon_list):
        for j, taxon2 in enumerate(taxon_list):
            dist_matrix[i, j] = pdm(taxon1, taxon2)
    # dendropy replaces '_' with ' '
    return dist_matrix, [t.label.replace(' ', '_') for t in taxon_list]

def tree_distances(tree_path: str):
    with open(tree_path, 'r') as nwk:
        nwk_data = nwk.read()
    dist_mat, seq_ids = _tree_distances(nwk_data)
    return pd.DataFrame(dist_mat, index=seq_ids, columns=seq_ids)


def validate_at_rank(tax_df: pd.DataFrame, dist_mat: np.ndarray, k=10, max_radius=0.05, min_neighbors=5, confidence_threshold=0.8, rank_level='genus'):
    labels = tax_df[rank_level].values
    n_seqs = len(labels)
    seq_ids = tax_df.index
    
    all_predictions: list[float] = []
    all_scores: list[float] = []
    sus: list[str] = []
    
    for i in range(n_seqs):
        actual_label: str = labels[i]
        seq_id = seq_ids[i]
        # Get distances to all other sequences
        distances = dist_mat[i, :]
        # Get neighbor indices
        neighbor_indices = np.where((np.arange(n_seqs) != i) & (distances <= max_radius))[0]
        neighbor_distances = distances[neighbor_indices]
        # Sort by distance and take k nearest
        sorted_idx = np.argsort(neighbor_distances)[:k]
        nearest_neighbors = neighbor_indices[sorted_idx]
        neighbor_labels = labels[nearest_neighbors]
        # Check if we have enough neighbors
        if len(nearest_neighbors) < min_neighbors:
            all_predictions.append(np.nan)
            all_scores.append(np.nan)
            continue
        # Majority vote
        label_counts = Counter(neighbor_labels)
        predicted_label, vote_count = label_counts.most_common(1)[0]
        confidence = vote_count / len(neighbor_labels)
        all_predictions.append(predicted_label)
        all_scores.append(confidence)
        # Categorize
        if actual_label != predicted_label and confidence >= confidence_threshold:
            sus.append(seq_id)
    return pd.DataFrame({'seq_id': seq_ids, 'predicted_label': all_predictions, 'score': all_scores, 'is_mislabeled': tax_df.index.isin(sus)})
