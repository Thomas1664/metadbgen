# Helper functions for validating the integrety of the database.
import pandas as pd
from Bio import SeqIO


def metadata(df: pd.DataFrame):
    assert not df.empty, "Metadata is empty."
    assert df['accession'].is_unique, "Metadata contains duplicate accessions."

def fasta(file_name: str):
    records = {record.id: str(record.seq) for record in SeqIO.parse(file_name, 'fasta')}
    accessions = pd.Series(records.keys())
    assert not accessions.empty, f"Fasta file {file_name} is empty."
    assert accessions.is_unique, f"Fasta file {file_name} contains duplicate accessions."

def fasta_and_metadata(fasta: str, meta_file: str):
    records = {record.id: str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}
    if meta_file.endswith('.csv'):
        df = pd.read_csv(meta_file)
    else:
        df = pd.read_csv(meta_file, sep='\t')
    # same numer of sequences
    assert len(records) == df.shape[0], f"{len(records)} sequences in fasta ({fasta}), but {df.shape[0]} sequences in metadata ({meta_file})."
    # same sequences
    assert set(records.keys()) == set(df['accession']), f"Metadata ({meta_file}) and fasta ({fasta}) contain different sequences."

def ingroup_vs_outgroup(meta_ingroup: pd.DataFrame, meta_outgroup: pd.DataFrame):
    # no overlap between in-group and outgroup sequences
    assert set(meta_ingroup['accession']).isdisjoint(set(meta_outgroup['accession']))
