import pandas as pd
from src import database
import sys


in_meta, in_fasta, clust, out_meta, out_fasta = sys.argv[1:]

meta = pd.read_csv(in_meta, sep='\t')
clust_df = database.import_clust_results(clust)
df = database.dereplicate(clust_df, meta)
acc_list_file = f'{out_fasta}.txt'
database.write_fasta(in_fasta, out_fasta, df['accession'])
df.to_csv(out_meta, index=False, sep='\t')
