from src import database, taxonomy
import sys

if __name__ == '__main__':
    meta_in, tax_in, backup_tax, meta_out, taxonomy_meta, tax_out = sys.argv[1:]
    taxonomy.filter_taxonomy(meta_in, tax_in, backup_tax, taxonomy_meta)
    meta_and_tax = taxonomy.merge_taxonomy(meta_in, taxonomy_meta)
    database.make_qiime_db(meta_and_tax, tax_out)
    meta_and_tax.to_csv(meta_out, sep='\t', index=False)
