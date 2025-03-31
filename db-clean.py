import argparse
import pandas as pd
from src import database


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean the COX2 database.")
    parser.add_argument("--in-metadata", type=str, required=True, help="Input metadata file")
    parser.add_argument("--out-metadata", type=str, required=True, help="Output metadata file")
    parser.add_argument("--out-fasta", type=str, required=True, help="Output FASTA file")
    args = parser.parse_args()

    df = pd.read_csv(args.in_metadata)
    filtered = database.filter_data(df)
    trimmed = database.trim(filtered)
    rev_comp = database.reverse_complement(trimmed)
    rev_comp.to_csv(args.out_metadata, sep='\t', index=False)
    database.export_fasta_no_tax(rev_comp, args.out_fasta)
