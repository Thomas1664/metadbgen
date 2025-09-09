from src import database, metadata, taxonomy
import sys
import asyncio

if __name__ == '__main__':
    meta_in, tax_in, backup_tax, meta_out, taxonomy_meta, tax_out, host_tax = sys.argv[1:]
    taxonomy.filter_taxonomy(meta_in, tax_in, backup_tax, taxonomy_meta)
    meta_and_tax = taxonomy.merge_taxonomy(meta_in, taxonomy_meta)
    meta_and_tax = metadata.parse_country(meta_and_tax)
    database.make_qiime_db(meta_and_tax, tax_out)
    meta_and_tax.to_csv(meta_out, sep='\t', index=False)
    asyncio.run(taxonomy.from_organsim_names(meta_and_tax['host'], f'{host_tax}.{{}}.xml', host_tax))
