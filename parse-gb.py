import csv
import sys
from src import genbank_parser, taxonomy
import pandas as pd
import json

gb_dir, missing_seqs, tax_db, tax_data, out_file, fasta_out = sys.argv[1:]

sequence_lookup = genbank_parser.load_sequence_lookup(missing_seqs)

with open(f"{out_file}.tmp", 'w', newline='') as csv_file, open(fasta_out, 'w') as fasta_file:
    csv_writer = csv.writer(csv_file)
    genbank_parser.dataclass_to_csv_header(csv_writer, genbank_parser.ParsedRecord)
    genbank_parser.parse_genbank_file(gb_dir, csv_writer, fasta_file, sequence_lookup=sequence_lookup)

print("done parsing gb")
df = pd.read_csv(f"{out_file}.tmp", index_col='accession')
df = df[~df.index.duplicated(keep='first')]
df['metadata'] = df['metadata'].apply(json.loads)
new_df = pd.DataFrame(df['metadata'].tolist(), index=df.index)

# Remove single-element lists
def unwrap_single_element(x):
    if isinstance(x, list) and len(x) == 1:
        return x[0]
    return x
new_df = new_df.map(unwrap_single_element)

df = df.drop(columns='metadata').join(new_df)
df = taxonomy.update_tax_ids(tax_db, df, index_col='accession')
tax_df = pd.read_csv(tax_data)
df = taxonomy.merge_taxonomy(df, tax_df)
df.to_csv(out_file, index=False)
