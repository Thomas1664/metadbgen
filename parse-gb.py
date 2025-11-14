import csv
import sys
from src import genbank_parser
import pandas as pd
import json

gb_dir, out_file, fasta_out = sys.argv[1:]

with open(out_file, 'w', newline='') as csv_file, open(fasta_out, 'w') as fasta_file:
    csv_writer = csv.writer(csv_file)
    genbank_parser.dataclass_to_csv_header(csv_writer, genbank_parser.ParsedRecord)
    genbank_parser.parse_genbank_file(gb_dir, csv_writer, fasta_file)

df = pd.read_csv(out_file, index_col='accession')
df['metadata'] = df['metadata'].apply(json.loads)
new_df = pd.DataFrame(df['metadata'].tolist(), index=df.index)

# Remove single-element lists
def unwrap_single_element(x):
    if isinstance(x, list) and len(x) == 1:
        return x[0]
    return x

# Apply it across the whole DataFrame
new_df = new_df.map(unwrap_single_element)

df = df.drop(columns='metadata').join(new_df)
df.to_csv(out_file)
