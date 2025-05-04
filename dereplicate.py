import pandas as pd
from src import database
import sys


in_meta, clust, out_meta, out_fasta = sys.argv[1:]

meta = pd.read_csv(in_meta, sep='\t')
clust = database.import_clust_results(clust)
df = database.dereplicate(clust, meta)
database.export_fasta_no_tax(df, out_fasta)
df.to_csv(out_meta, index=False, sep='\t')
