import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import subprocess
import io
import sys
import matplotlib.pyplot as plt
import logomaker as lm
import numpy as np


def parse_cutadapt_info(info_file: str, out_path: str):
    df = pd.read_csv(info_file, sep='\t', header=None)
    df = df.drop(columns=[8, 9, 10]) # Quality columns are empty
    df.columns = ['accession', 'errors', 'start', 'end', 'left_seq', 'primer_seq', 'right_seq', 'primer_name', 'reverse_complement']
    df = df[df['errors'] != -1] # Discard untrimmed sequences
    df['accession'] = df['accession'].str.replace(' rc', '')
    df['reverse_complement'] = df['reverse_complement'].astype(bool)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['primer_size'] = df['end'] - df['start']

    primer_map = {"1;1": "forward", "1;2": "reverse"}
    df["primer_direction"] = df["primer_name"].map(primer_map)

    # Pivot so forward/reverse get their own columns
    wide = df.pivot_table(
        index="accession",
        columns="primer_direction",
        values=["errors", "start", "end", "left_seq", "primer_seq", "right_seq", "primer_name", 'primer_size'],
        aggfunc="first"
    )
    # Flatten the MultiIndex column names
    wide.columns = [f"{col[0]}_{col[1]}" for col in wide.columns]
    # Reset index so accession is a normal column
    wide = wide.reset_index()
    # Keep accession first, then all forward, then all reverse
    ordered_cols = ["accession", "start_forward", "end_forward", "errors_forward", 'primer_size_forward', "primer_name_forward", "primer_seq_forward", "left_seq_forward", "right_seq_forward", "start_reverse", "end_reverse", 'primer_size_reverse', "errors_reverse", "primer_name_reverse", "primer_seq_reverse", "left_seq_reverse", "right_seq_reverse"]
    wide = wide[ordered_cols]
    wide = wide.merge(df[['accession', 'reverse_complement']], on='accession', how='inner')
    wide = wide.drop_duplicates(subset="accession")
    wide.to_csv(out_path, sep='\t', index=False)


def trim_msa(cutadapt_path: str, meta_path: str, msa_path: str, fasta_out: str):
    cutadapt = pd.read_csv(cutadapt_path, sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    cut_filtered = cutadapt[cutadapt['accession'].isin(meta['accession'])]
    
    references = cut_filtered.dropna(subset=['primer_seq_forward', 'primer_seq_reverse'])
    accessions = references['accession'].tolist()

    alignment = AlignIO.read(msa_path, 'fasta')
    firsts = []
    lasts = []
    if len(accessions) == 0:
        # No references available; skip trimming
        with open(fasta_out, 'w') as fasta:
            for record in alignment:
                #TODO: check again for min_size?
                ungapped_seq = str(record.seq).replace('-', '').upper()
                fasta.write(f">{record.id}\n{ungapped_seq}\n")
        return
    
    for sequence in alignment:
        if sequence.id not in accessions:
            continue
        seq = str(sequence.seq)
        # find first non-gap
        first = next(i for i, base in enumerate(seq) if base != "-")
        # find last non-gap
        last = len(seq) - next(i for i, base in enumerate(reversed(seq)) if base != "-") - 1
        firsts.append(first)
        lasts.append(last)
    firsts = pd.Series(firsts)
    lasts = pd.Series(lasts)
    # trim at most frequent position among the reference sequences
    trim_start = int(firsts.value_counts(sort=True, ascending=False).index[0])
    trim_end = int(lasts.value_counts(sort=True, ascending=False).index[0])

    # slice alignment: [rows, columns]
    msa_trimmed = alignment[:, trim_start:trim_end + 1]
    # write fasta file
    with open(fasta_out, 'w') as fasta:
        for record in msa_trimmed:
            #TODO: check again for min_size?
            ungapped_seq = str(record.seq).replace('-', '').upper()
            fasta.write(f">{record.id}\n{ungapped_seq}\n")


def _write_primer_fasta(info_path: str, meta_path: str, fasta_path: str, primer_size: int, direction: str):
    info = pd.read_csv(info_path, sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    cut_filtered = info[info['accession'].isin(meta['accession'])]

    with open(fasta_path, 'w') as out:
        for _, row in cut_filtered.iterrows():
            if pd.isna(row[f'primer_size_{direction}']):
                continue
            if row[f'primer_size_{direction}'] > primer_size:
                continue
            seq = str(row[f'primer_seq_{direction}']).upper()
            # make reverse primer reverse complement because cutadapt returns sequence in direction of amplicon
            if direction == 'reverse':
                seq = str(Seq(seq).reverse_complement())
            out.write(f">{row['accession']}\n{seq}\n")
    return meta.shape[0]

def _align_primers(fasta: str):
    result = subprocess.run(['mafft', '--auto', fasta], stdout=subprocess.PIPE, stderr=None, check=True, text=True)
    return [str(record.seq).upper() for record in SeqIO.parse(io.StringIO(result.stdout), 'fasta')]

def _draw_logo(seqs: list[str], logo_path: str, total_seqs: int):
    lengths = [len(seq) for seq in seqs]
    counts_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='-N')
    counts_df = counts_df.map(lambda x: x / total_seqs)
    lm.Logo(counts_df, color_scheme='classic')
    plt.xticks(range(1, np.max(lengths), 2))
    plt.xlabel("Base positions")
    plt.ylabel("Proportions")
    plt.ylim((0, 1))
    plt.savefig(logo_path, format='svg', bbox_inches='tight')

def logo_plot(info_path: str, meta_path: str, fasta_path: str, primer_size: int, logo_path: str, direction: str):
    total_seqs = _write_primer_fasta(info_path, meta_path, fasta_path, primer_size, direction)
    seqs = _align_primers(fasta_path)
    _draw_logo(seqs, logo_path, total_seqs)
