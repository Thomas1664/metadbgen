from collections.abc import Iterable
from io import TextIOWrapper
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def build_search_query(keywords: list[str], organism: str):
    keywords = [f'"{keyword}"' for keyword in keywords]
    joined_keywords = ' OR '.join(keywords)
    return f'({joined_keywords}) AND "{organism}"[organism]'

def remove_if_contains_substr(df: pd.DataFrame, col: str, substr: str):
    return df[~df[col].str.contains(substr, case=False, regex=False)]

def sanitize_org_name(org_name: str) -> str:
    org_name_arr = org_name.split(' ')
    to_remove = ['x', 'var.', 'var', 'aff.', 'aff', 'cf.', 'cf']
    for s in to_remove:
        try:
            org_name_arr.remove(s)
        except ValueError:
            pass
    return ' '.join(org_name_arr[:2])

def fix_tax_id(row: pd.Series, trimmed: pd.DataFrame):
    if row['organism_name_sanitized'] == row['organism_name']:
        return row['tax_id']
    if row['organism_name_sanitized'] in trimmed['organism_name'].unique():
        new_tax_id = trimmed.loc[trimmed['organism_name'] == row['organism_name_sanitized']]['tax_id'].iloc[0]
        return new_tax_id
    return row['tax_id']

def prefer_refseq(trimmed: pd.DataFrame):
    species_groups = trimmed.groupby('organism_name_sanitized')
    for _, species in species_groups:
        # accessions of refseq sequences contain an underscore
        # see https://www.ncbi.nlm.nih.gov/refseq/about/
        refseq = species[species['accession'].str.contains('_')]
        if not refseq.empty:
            non_refseq = species[~species['accession'].str.contains('_')]
            if not non_refseq.empty:
                trimmed.drop(trimmed[trimmed['accession'].isin(non_refseq['accession'])].index, inplace=True)
    return trimmed

def filter_data(df: pd.DataFrame):
    # remove unverified sequences
    filtered = df.drop(df[df['title'].str.lower().str.startswith('unverified')].index).copy()
    # remove sequences with unknown species (i.e. sp.)
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'sp.')
    # remove sequences from uncultured samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'uncultured')
    # remove environmental samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'environmental')
    # remove nom. inval. species
    filtered = remove_if_contains_substr(filtered, 'title', 'nom. inval.')
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'nom. inval.')

    removed = df.merge(filtered, how='outer', indicator=True).query('_merge == "left_only"').drop(columns='_merge')

    # Remove sequences that are too short
    trimmed = filtered.copy()
    # sanitize organism name
    trimmed['organism_name_sanitized'] = trimmed['organism_name']
    trimmed['organism_name_sanitized'] = trimmed['organism_name_sanitized'].apply(sanitize_org_name)
    # sanitize tax_id
    trimmed['tax_id_sanitized'] = trimmed['tax_id']
    trimmed['tax_id_sanitized'] = trimmed.apply(lambda row: fix_tax_id(row, trimmed), axis=1)

    trimmed = prefer_refseq(trimmed)
    return trimmed, removed

def write_seq(file: TextIOWrapper, accession: str, seq: SeqRecord | str):
    file.write(f">{accession}\n{str(seq)}\n")

def write_fasta(records: Iterable[SeqRecord] | str, out_file: TextIOWrapper | str, accessions: Iterable[str] | None=None):
    if isinstance(records, str):
        records = SeqIO.parse(records, 'fasta')
    if isinstance(accessions, pd.Series):
        accessions = accessions.values
    if accessions is not None:
        filter_fn = filter(lambda x: x.id in accessions, records)
    else:
        filter_fn = filter(lambda _: True, records)
    if isinstance(out_file, str):
        with open(out_file, 'w') as f:
            for record in filter_fn:
                write_seq(f, record.id, record.seq)
    else:
        for record in filter_fn:
            write_seq(out_file, record.id, record.seq)

def dereplicate(clust_results: pd.DataFrame, metadata: pd.DataFrame):
    merged = metadata.merge(clust_results, how='left', left_on='accession', right_on='query_label')
    df_sorted = merged.sort_values(by=['cluster_id', 'tax_id', 'trimed_seq_len'], ascending=[True, True, False])

    def filter_func(group: pd.DataFrame):
        # only keep identical sequences if they are from different species
        # Removes identical (but potentially shorter) sequences from the same species
        # ==> Remove redundant sequences that don't give any new information
        # PRE: each group is sorted (desc) by trimed_seq_len
        if group.shape[0] == 1:
            return group
        subgroup = group.groupby('tax_id', group_keys=False)
        return subgroup.head(1)

    df_filtered = df_sorted.copy().groupby('cluster_id', group_keys=False).apply(filter_func)
    return df_filtered.drop(columns=['type', 'query_length', 'similarity', 'align_res', 'query_label', 'strand'])


def import_clust_results(path: str):
    clust_res = pd.read_csv(path, sep='\t', header=None)
    clust_res.columns = ['type', 'cluster_id', 'query_length', 'similarity', 'strand', 0, 1, 'align_res', 'query_label', 'centroid_label']
    clust_res = clust_res[['type', 'cluster_id', 'query_length', 'similarity', 'strand', 'align_res', 'query_label']]
    return clust_res[(clust_res['type'] == 'S') | (clust_res['type'] == 'H')]

def make_qiime_db(meta_df: pd.DataFrame, tax_out_path: str):
    with open(tax_out_path, "w") as tax_out:
        for _, row in meta_df.iterrows():
            if pd.isna(row["family"]):
                row["family"] = row["genus"]
            if pd.isna(row["order"]):
                row["order"] = row["family"]
            if pd.isna(row["class"]):
                row["class"] = row["order"]
            if pd.isna(row["phylum"]):
                row["phylum"] = row["class"]
            if pd.isna(row["kingdom"]):
                row["kingdom"] = row["phylum"]

            taxonomy = f'd__{row["domain"]};k__{row["kingdom"]};p__{row["phylum"]};c__{row["class"]};o__{row["order"]};f__{row["family"]};g__{row["genus"]};s__{row["organism_name_sanitized"]}'
            taxonomy = taxonomy.replace(" ", "+")
            print(f'{row["accession"]}\t{taxonomy}', file=tax_out)

def export_sintax(msa_path: str, tax_path: str, out_db: str):
    tax = pd.read_csv(tax_path, sep='\t', header=None)
    tax[1] = tax[1].str.replace(';', ',')
    tax[1] = tax[1].str.replace('__', ':') + ';'

    records = SeqIO.parse(msa_path, 'fasta')
    with open(out_db, 'w') as fasta:
        for record in records:
            print(f'>{record.id};tax={tax.loc[tax[0] == record.id][1].values[0]}', file=fasta)
            print(str(record.seq).upper(), file=fasta)


def seq_base_composition(seq: str):
    seq = seq.upper()
    seq_arr = np.frombuffer(seq.encode('ascii'), dtype='S1')
    bases = np.array([b'A', b'T', b'G', b'C'])
    # Vectorized comparison: create a 2D boolean matrix, then sum along axis 1
    counts_arr = (seq_arr[:, None] == bases).sum(axis=0, dtype=np.int32)
    #perc_arr = (counts_arr / seq_arr.size) * 100
    counts = dict(zip(bases.astype(str), counts_arr))
    counts['other'] = int(seq_arr.size - counts_arr.sum())
    counts['seq_len'] = seq_arr.size
    return counts
