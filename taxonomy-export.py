from src import database
import sys

if __name__ == '__main__':
    meta_in, tax_in, fasta_out, meta_out, tax_out = sys.argv[1:]
    meta_and_tax = database.merge_taxonomy(meta_in, tax_in)
    database.make_qiime_db(meta_and_tax, fasta_out, tax_out)
    meta_and_tax.to_csv(meta_out, sep='\t', index=False)
