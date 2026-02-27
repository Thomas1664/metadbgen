from Bio import AlignIO
from Bio import Align
import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices


def _compute_gappiness(alignment: Align.MultipleSeqAlignment) -> np.ndarray:
    # construct numpy array
    seq_strings = [str(record.seq) for record in alignment]
    ali_mat = np.array([list(seq) for seq in seq_strings])
    num_seqs, aln_len = ali_mat.shape
    
    # mask all gaps
    residue_mask = (ali_mat != '-')
    residues_per_column: np.ndarray = np.sum(residue_mask, axis=0)
    
    # only consider columns with more gaps than residues
    more_gaps_than_residues_mask = residues_per_column < (num_seqs / 2)
    residues_per_column = np.where(more_gaps_than_residues_mask, residues_per_column, 0)
    
    # expand residue counts per column for each sequence
    # shape: num_seqs, aln_length
    residues_per_column_mat = np.broadcast_to(residues_per_column, (num_seqs, aln_len))
    # replace gaps with 0 using mask
    residues_per_column_mat_filtered = np.where(residue_mask, residues_per_column_mat, 0)
    
    # source: https://stackoverflow.com/questions/34755244/apply-bincount-to-each-row-of-a-2d-numpy-array
    N: int = num_seqs + 1 # +1 because we use base-1 ranges
    # encode residues per row because bincount only works with 1D arrays
    id = residues_per_column_mat_filtered + (N * np.arange(residues_per_column_mat_filtered.shape[0]))[:, None]
    # count how often each bin shows up per row
    # remove first column because there are no columns with 0 residues
    residue_hist = np.bincount(id.ravel(), minlength=N * residues_per_column_mat_filtered.shape[0]).reshape(-1, N)[:,1:]
    
    n_vals = np.arange(1, num_seqs + 1)
    weights = (num_seqs - n_vals)
    gappiness: np.ndarray = (residue_hist * weights).sum(axis=1) / (aln_len * num_seqs)
    return gappiness

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
    # FIXME: Weights disabled because they're unused & take a really long time.
    # also delete wight computation above.
    #weights = _compute_sequence_weights(alignment, substitution_matrix)

    # Normalize
    #norm_weights = (weights - np.min(weights)) / (np.max(weights) - np.min(weights))
    norm_gpp = (gappiness - np.min(gappiness)) / (np.max(gappiness) - np.min(gappiness))

    df = pd.DataFrame({
        'accession': [record.id for record in alignment],
        #'weight': weights,
        #'weight_normalized': norm_weights,
        'gappiness': gappiness,
        'gappiness_normalized': norm_gpp
    })
    return df.sort_values('gappiness_normalized', ascending=False)

def filter_by_gappiness(df: pd.DataFrame, gappiness_threshold=0.1):
    return df[df['gappiness'].lt(gappiness_threshold)]
