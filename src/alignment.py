from Bio import AlignIO
import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices


def _compute_gappiness(alignment):
    alignment_array = np.array([list(record.seq) for record in alignment])
    num_sequences = alignment_array.shape[0]
    alignment_length = alignment_array.shape[1]
    gappiness_scores = []

    for _, sequence in enumerate(alignment_array):
        gpp = 0
        for column_index in range(alignment_length):
            column = alignment_array[:, column_index]
            num_residues_in_column = np.sum(column != '-')
            is_residue_in_this_sequence = sequence[column_index] != '-'
            majority_have_gaps = num_residues_in_column < (num_sequences / 2)

            if is_residue_in_this_sequence and majority_have_gaps:
                penalty = num_sequences - num_residues_in_column
                gpp += penalty

        normalized_gpp = gpp / (alignment_length * num_sequences)
        gappiness_scores.append(normalized_gpp)
    return gappiness_scores

def _compute_gap_penalty(matrix: np.ndarray):
    scores = matrix.flatten()  # get all values
    return np.min(scores) - np.std(scores)


def _compute_sequence_weights(alignment, matrix_name="NUC.4.4"):
    substitution_matrix = substitution_matrices.load(matrix_name)
    gap_penalty = _compute_gap_penalty(substitution_matrix)

    # Convert alignment to NumPy array of characters
    alignment_array = np.array([list(record.seq.upper()) for record in alignment])
    n_seqs, alignment_length = alignment_array.shape

    # --- Step 1: Map characters to integer indices ---
    chars = sorted(set("".join("".join(seq) for seq in alignment_array)))
    if '-' not in chars:
        chars.append('-')
    char_to_idx = {c: i for i, c in enumerate(chars)}

    # Build substitution matrix as a NumPy array
    subst_array = np.full((len(chars), len(chars)), gap_penalty, dtype=float)
    for a, a_idx in char_to_idx.items():
        for b, b_idx in char_to_idx.items():
            try:
                subst_array[a_idx, b_idx] = substitution_matrix[a][b]
            except (KeyError, IndexError):
                try:
                    subst_array[a_idx, b_idx] = substitution_matrix[b][a]
                except (KeyError, IndexError):
                    subst_array[a_idx, b_idx] = gap_penalty

    # --- Step 2: Convert alignment to integer matrix ---
    int_alignment = np.vectorize(char_to_idx.get)(alignment_array)

    # --- Step 3: Compute pairwise scores with broadcasting ---
    weights = np.zeros(n_seqs, dtype=float)

    for i in range(n_seqs):
        # Compare seq i to all others at once (broadcast over positions)
        seq_i = int_alignment[i]  # shape (L,)
        others = int_alignment    # shape (N, L)
        # Look up scores for all pairwise char comparisons
        scores = subst_array[seq_i[:, None], others.T]  # shape (L, N)
        # Remove self-comparison (i vs i)
        scores[:, i] = 0
        # Sum across positions and other sequences
        weights[i] = scores.sum()
    return weights.tolist()

def calculate_gappiness(msa_path, substitution_matrix='NUC.4.4'):
    alignment = AlignIO.read(msa_path, 'fasta')
    gappiness = _compute_gappiness(alignment)
    weights = _compute_sequence_weights(alignment, substitution_matrix)

    # Normalize
    norm_weights = (weights - np.min(weights)) / (np.max(weights) - np.min(weights))
    norm_gpp = (gappiness - np.min(gappiness)) / (np.max(gappiness) - np.min(gappiness))

    df = pd.DataFrame({
        'accession': [record.id for record in alignment],
        'weight': weights,
        'weight_normalized': norm_weights,
        'gappiness': gappiness,
        'gappiness_normalized': norm_gpp
    })
    return df.sort_values('gappiness_normalized', ascending=False)

def filter_by_gappiness(df: pd.DataFrame, gappiness_treshold=0.1):
    return df[df['gappiness'].lt(gappiness_treshold)]
