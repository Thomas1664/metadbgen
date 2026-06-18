from src import database, metadata, taxonomy
import sys
import asyncio
import pandas as pd

if __name__ == '__main__':
    meta_in, tax_db, meta_out, tax_out, host_tax = sys.argv[1:]
    meta_and_tax = pd.read_csv(meta_in, sep='\t')
    meta_and_tax = metadata.parse_country(meta_and_tax)
    database.make_qiime_db(meta_and_tax, tax_out)
    meta_and_tax.to_csv(meta_out, sep='\t', index=False)
    if 'host' in meta_and_tax:
        asyncio.run(taxonomy.from_organsim_names(tax_db, meta_and_tax['host'], host_tax))
