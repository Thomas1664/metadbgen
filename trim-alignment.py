import sys
from src import primer, alignment, database
import pandas as pd
from Bio import SeqIO


msa, cutadapt_info, meta, fasta_tmp, fasta_out, out_path, gap_treshold, min_seq_size = sys.argv[1:]

primer.trim_msa(cutadapt_info, meta, msa, fasta_tmp, int(min_seq_size))

gap_df = alignment.calculate_gappiness(msa).sort_values('gappiness_normalized', ascending=False)
gap_df.to_csv(f'{out_path}/msa-gap-scores.csv', index=False)
# Re-order sequences in alignment based on gappiness scores, desc
msa_records = {record.id: str(record.seq) for record in SeqIO.parse(msa, 'fasta')}
with open(msa, 'w') as f:
    for _, row in gap_df.iterrows():
        f.write(f'>{row["accession"]}\n{msa_records[row["accession"]]}\n')

gap_df = alignment.filter_by_gappiness(gap_df, float(gap_treshold))
meta_df = pd.read_csv(meta, sep='\t')
meta_df = meta_df[meta_df['accession'].isin(gap_df['accession'])]
meta_df.to_csv(f'{out_path}/metadata.tsv', sep='\t', index=False)
database.write_fasta(fasta_tmp, fasta_out, meta_df['accession'])
