from collections.abc import Sequence
from itertools import batched
import pandas as pd
import asyncio
import os
import glob
import shutil
from src import request


missing_df = pd.read_csv(snakemake.input.missing)
tmp_dir = snakemake.params.tmp_dir
            
def build_query(organisms: Sequence[str], exclude_keywords: list[str]):
    exclude = " NOT ".join([f'"{keyword}"' for keyword in exclude_keywords])
    q_organisms = " OR ".join([f'"{organism}"[organism]' for organism in organisms])
    return f'({q_organisms}) NOT {exclude} AND ("0001/01/01"[PDAT] : "2026/04/15"[PDAT])'

async def download_all():
    client = request.HttpClient()
    os.makedirs(tmp_dir, exist_ok=True)
    for batch, species in enumerate(batched(missing_df["species"].tolist(), 10)):
        if os.path.exists(f'{tmp_dir}/{batch}.gb'):
            continue
        query = build_query(species, snakemake.params.keywords)
        print(query)
        res = await client.esearch(query)
        if res.result_count == 0:
            continue
        await client.efetch_all(res, f'{tmp_dir}/{batch}.gb')
        for split_file in glob.glob(f'{tmp_dir}/{batch}_*.gb'):
            os.remove(split_file)
    gb_files = sorted(glob.glob(f'{tmp_dir}/*.gb'))
    with open(snakemake.output.gb_file, 'wb') as out_f:
        for gb_file in gb_files:
            with open(gb_file, 'rb') as f:
                shutil.copyfileobj(f, out_f)

asyncio.run(download_all())
