import pandas as pd
import numpy as np
import dendropy
from collections import Counter
from src import taxonomy

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
