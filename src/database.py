from collections.abc import Iterable
from io import TextIOWrapper
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import re


def build_search_query(keywords: list[str], organism: str):
    keywords = [f'"{keyword}"' for keyword in keywords]
    joined_keywords = ' OR '.join(keywords)
    # HOTFIX: For reproducability, limit release date
    return f'({joined_keywords}) AND "{organism}"[organism] AND ("0001/01/01"[PDAT] : "2026/04/15"[PDAT])'

def remove_if_contains_substr(df: pd.DataFrame, col: str, substr: str):
    return df[~df[col].str.contains(substr, case=False, regex=False)]

def sanitize_org_name(org_name: str) -> str:
    words = org_name.split(' ')
    explicit_modifiers = {'taxon'}
    # TODO: first check if name matches binomen: "^[A-Z][a-z\-]+ [a-z\-]+$"gmu
    # also check 3rd word when removing second word
    if len(words) > 2:
        second = words[1]
        epitet_pat = re.compile(r"^[a-z\-]+$")
        if second in explicit_modifiers or len(second) == 1 or epitet_pat.match(second) is None:
            words = [words[0]] + words[2:]

    return ' '.join(words[:2])

def fix_tax_id(row: pd.Series, trimmed: pd.DataFrame):
    if row['organism_name_sanitized'] == row['organism_name']:
        return row['tax_id']
    if row['organism_name_sanitized'] in trimmed['organism_name'].unique():
        new_tax_id = trimmed.loc[trimmed['organism_name'] == row['organism_name_sanitized']]['tax_id'].iloc[0]
        return new_tax_id
    return row['tax_id']

def prefer_refseq(trimmed: pd.DataFrame):
    # does this row’s group contain ANY refseq?
    # accessions of refseq sequences contain an underscore
    # see https://www.ncbi.nlm.nih.gov/refseq/about/
    has_refseq = trimmed.groupby('organism_name_sanitized')['accession'] \
                        .transform(lambda s: s.str.contains('_').any())
    # is THIS row refseq?
    is_refseq = trimmed['accession'].str.contains('_', regex=False)
    # if group has refseq, keep only refseq rows
    # if group has no refseq, keep all rows
    mask = (~has_refseq) | (is_refseq)
    return trimmed[mask].copy()


def filter_data(df: pd.DataFrame, remove_on_species: bool = True):
    # remove unverified sequences
    filtered = df.drop(df[df['title'].str.lower().str.startswith('unverified')].index).copy()
    filtered = filtered.drop_duplicates(subset=['accession'])
    # remove sequences with unknown species (i.e. sp.)
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'sp.')
    # remove sequences from uncultured samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'uncultured')
    # remove environmental samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'environmental')
    # remove nom. inval. species
    filtered = remove_if_contains_substr(filtered, 'title', 'nom. inval.')
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'nom. inval.')
    # filter "imposters"
    filtered = remove_if_contains_substr(filtered, 'title', '-like ')
    filtered = remove_if_contains_substr(filtered, 'title', 'predicted')

    removed = df.merge(filtered, how='outer', indicator=True).query('_merge == "left_only"').drop(columns='_merge')

    trimmed = filtered.copy()
    # sanitize organism name
    trimmed['organism_name_sanitized'] = trimmed['organism_name']
    trimmed['organism_name_sanitized'] = trimmed['organism_name_sanitized'].apply(sanitize_org_name)
    trimmed.index = trimmed['accession']
    if remove_on_species:
        is_on_mask = trimmed['organism_name'] != trimmed['organism_name_sanitized']
        non_on_seqs = trimmed[~is_on_mask]
        non_on_species = non_on_seqs['organism_name'].unique()
        # ON species that don't have a non-ON species name
        # e.g. Globisporangium aff. irregulare while there is no Globisporangium irregulare
        on_native_mask = is_on_mask & ~trimmed['organism_name_sanitized'].isin(non_on_species)
        on_native = trimmed[on_native_mask]
        trimmed = trimmed.filter([*non_on_seqs.index, *on_native.index], axis='index')
    trimmed = trimmed.dropna(subset='genus')
    trimmed = prefer_refseq(trimmed)
    return trimmed, removed

def write_seq(file: TextIOWrapper, accession: str, seq: SeqRecord | str):
    file.write(f">{accession}\n{str(seq)}\n")

def write_fasta(in_fasta: str, out_fasta: str, accessions: Iterable[str] | pd.Series):
    if not isinstance(accessions, pd.Series):
        accessions = pd.Series(accessions)
    acc_list_file = f'{out_fasta}.txt'
    accessions.to_csv(acc_list_file, header=False, index=False)
    subprocess.run(['seqkit', 'grep', '-f', acc_list_file, '-o', out_fasta, in_fasta], check=True)

def dereplicate(clust_results: pd.DataFrame, metadata: pd.DataFrame, remove_on_species: bool = True):
    merged = metadata.merge(clust_results, how='left', left_on='accession', right_on='query_label', suffixes=('_x', ''))
    if remove_on_species:
        # Remove ON sequences from clusters that also contain non-ON sequences.
        is_on = merged['organism_name'] != merged['organism_name_sanitized']
        clusters_with_non_on = merged.loc[~is_on, 'cluster_id'].unique()
        # remove rows if row is ON and cluster has non-ON sequences
        merged = merged[~(is_on & merged['cluster_id'].isin(clusters_with_non_on))]
    # Sort so that the longest sequence per (cluster_id, tax_id) comes first,
    # then drop_duplicates keeps exactly one row per (cluster_id, tax_id) — the longest.
    df_sorted = merged.sort_values(by=['cluster_id', 'tax_id', 'trimed_seq_len', 'accession'], ascending=[True, True, False, True])
    df_filtered = df_sorted.drop_duplicates(subset=['cluster_id', 'tax_id'], keep='first')
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

            taxonomy = f'd__{row["domain"]};k__{row["kingdom"]};p__{row["phylum"]};c__{row["class"]};o__{row["order"]};f__{row["family"]};g__{row["genus"]};s__{row["species"]}'
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
