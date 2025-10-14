import sys
from src import primer, alignment, database
import pandas as pd
from Bio import SeqIO


msa, cutadapt_info, meta, fasta_out, meta_out, gap_treshold = sys.argv[1:]

primer.trim_msa(cutadapt_info, meta, msa, fasta_out)

gap_df = alignment.calculate_gappiness(msa)
gap_df = alignment.filter_by_gappiness(gap_df, float(gap_treshold))
meta_df = pd.read_csv(meta, sep='\t')
meta_df = meta_df[meta_df['accession'].isin(gap_df['accession'])]
meta_df.to_csv(meta_out, sep='\t', index=False)
records = [record for record in SeqIO.parse(fasta_out, 'fasta') if record.id in meta_df['accession'].values]
database.write_fasta(records, fasta_out)
