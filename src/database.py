import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def remove_if_contains_substr(df: pd.DataFrame, col: str, substr: str):
    return df[~df[col].str.contains(substr, case=False, regex=False)]

def filter_short_cox2(df: pd.DataFrame):
    def filter_fn(row: pd.Series):
        start = int(str(row['cox2_start']).lstrip('<'))
        end = int(str(row['cox2_end']).lstrip('>'))
        return end - start < 200
    return df[~df.apply(filter_fn, axis=1)]

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
    filtered = df.drop(df[df['title'].str.lower().str.startswith('unverified')].index)
    # remove sequences with unknown cox2 region
    # TODO: Possible annotation errors?
    filtered = filtered.dropna(subset=['cox2_start', 'cox2_end'])
    # remove sequences with unknown species (i.e. sp.)
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'sp.')
    # remove sequences from uncultured samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'uncultured')
    # remove environmental samples
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'environmental')
    # remove nom. inval. species
    filtered = remove_if_contains_substr(filtered, 'title', 'nom. inval.')
    filtered = remove_if_contains_substr(filtered, 'organism_name', 'nom. inval.')
    # Remove sequences that are too short
    trimmed = filter_short_cox2(filtered.copy())
    # sanitize organism name
    trimmed['organism_name_sanitized'] = trimmed['organism_name']
    trimmed['organism_name_sanitized'] = trimmed['organism_name_sanitized'].apply(sanitize_org_name)
    # sanitize tax_id
    trimmed['tax_id_sanitized'] = trimmed['tax_id']
    trimmed['tax_id_sanitized'] = trimmed.apply(lambda row: fix_tax_id(row, trimmed), axis=1)

    trimmed = prefer_refseq(trimmed)
    return trimmed

def _trim_seq(cox2_start: str, cox2_end: str, seq: str) -> str:
    start = int(str(cox2_start).lstrip('<'))
    end = int(str(cox2_end).lstrip('>'))
    assert(start < end)
    return seq[start:end]

def trim(df: pd.DataFrame):
    df['sequence'] = df.apply(lambda row: _trim_seq(row['cox2_start'], row['cox2_end'], row['sequence']), axis=1)
    return df

def reverse_complement(df: pd.DataFrame):
    df['sequence'] = df.apply(lambda row: str(Seq(row['sequence']).reverse_complement()) if row['cox2_strand'] == -1 else row['sequence'], axis=1)
    return df

def export_fasta_no_tax(df: pd.DataFrame, file: str):
    with open(file, 'w') as out:
        for _, row in df.iterrows():
            print(f">{row['accession']}", file=out)
            print(row['sequence'], file=out)

def import_fasta_no_tax(df: pd.DataFrame, file: str):
    rows = [{'accession': record.id, 'sequence': record.seq} for record in SeqIO.parse(file, 'fasta')]
    seq = pd.DataFrame(rows)
    df = df.drop(columns=['sequence'])
    return df.merge(seq, how='inner', on='accession')
