import sys
from src import alignment, database
import pandas as pd
from Bio import SeqIO


msa, meta, fasta, meta_out, fasta_out, gap_out, gap_treshold = sys.argv[1:]

gap_df = alignment.calculate_gappiness(msa).sort_values('gappiness_normalized', ascending=False)
gap_df.to_csv(gap_out, index=False)
# Re-order sequences in alignment based on gappiness scores, desc
msa_records = {record.id: str(record.seq) for record in SeqIO.parse(msa, 'fasta')}
with open(msa, 'w') as f:
    for _, row in gap_df.iterrows():
        f.write(f'>{row["accession"]}\n{msa_records[row["accession"]]}\n')

gap_df = alignment.filter_by_gappiness(gap_df, float(gap_treshold))
meta_df = pd.read_csv(meta, sep='\t')
meta_df = meta_df[meta_df['accession'].isin(gap_df['accession'])]
meta_df.to_csv(meta_out, sep='\t', index=False)
database.write_fasta(fasta, fasta_out, meta_df['accession'])
