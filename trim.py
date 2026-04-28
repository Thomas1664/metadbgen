import pandas as pd
from Bio import SeqIO
from src import primer
import sys


in_meta, fasta, out_meta, cutadapt_info, cutadapt_out = sys.argv[1:]
meta = pd.read_csv(in_meta)
# calculate length of trimmed sequences
rows = pd.DataFrame([{'accession': record.id, 'trimed_seq_len': len(str(record.seq))} for record in SeqIO.parse(fasta, 'fasta')])
# filter metadata by trimmed sequences & add sequence length to metadata
meta = meta.merge(rows, how='inner', on='accession')
meta.to_csv(out_meta, sep='\t', index=False)
primer.parse_cutadapt_info(cutadapt_info, cutadapt_out)
