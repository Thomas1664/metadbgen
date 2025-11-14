import argparse
import pandas as pd
from src import database


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean the COX2 database.")
    parser.add_argument("--in-metadata", type=str, required=True, help="Input metadata file")
    parser.add_argument("--in-fasta", type=str, required=True, help="Input fasta file")
    parser.add_argument("--out-metadata", type=str, required=True, help="Output metadata file")
    parser.add_argument("--out-fasta", type=str, required=True, help="Output FASTA file")
    parser.add_argument("--out-metadata-removed", type=str, required=True, help="Output metadata file containing removed sequences")
    parser.add_argument("--out-fasta-removed", type=str, required=True, help="Output FASTA file containing removed sequences")
    args = parser.parse_args()

    df = pd.read_csv(args.in_metadata)
    filtered, removed = database.filter_data(df)
    filtered.to_csv(args.out_metadata, sep='\t', index=False)
    removed.to_csv(args.out_metadata_removed, sep='\t', index=False)
    database.write_fasta(args.in_fasta, args.out_fasta, filtered['accession'])
    database.write_fasta(args.in_fasta, args.out_fasta_removed, removed['accession'])
