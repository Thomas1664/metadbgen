from asyncio import Semaphore
import asyncio
from collections import defaultdict
from dataclasses import dataclass
import random
from urllib.parse import urlencode
from httpx_retries import Retry, RetryTransport
from pydantic import BaseModel, Field
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from bs4 import BeautifulSoup
from src import request, taxonomy, database
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker as lm
import numpy as np
import pandas as pd
import csv
import httpx
import io
import re
import subprocess
import time
import threading


def parse_cutadapt_info(info_file: str, out_path: str):
    df = pd.read_csv(info_file, sep='\t', header=None)
    df = df.drop(columns=[8, 9, 10]) # Quality columns are empty
    df.columns = ['accession', 'errors', 'start', 'end', 'left_seq', 'primer_seq', 'right_seq', 'primer_name', 'reverse_complement']
    df = df[df['errors'] != -1] # Discard untrimmed sequences
    df['accession'] = df['accession'].str.replace(' rc', '')
    df['reverse_complement'] = df['reverse_complement'].astype(bool)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['primer_size'] = df['end'] - df['start']

    primer_map = {"1;1": "forward", "1;2": "reverse"}
    df["primer_direction"] = df["primer_name"].map(primer_map)

    # Pivot so forward/reverse get their own columns
    wide = df.pivot_table(
        index="accession",
        columns="primer_direction",
        values=["errors", "start", "end", "left_seq", "primer_seq", "right_seq", "primer_name", 'primer_size'],
        aggfunc="first"
    )
    # Flatten the MultiIndex column names
    wide.columns = [f"{col[0]}_{col[1]}" for col in wide.columns]
    # Reset index so accession is a normal column
    wide = wide.reset_index()
    # Keep accession first, then all forward, then all reverse
    ordered_cols = ["accession", "start_forward", "end_forward", "errors_forward", 'primer_size_forward', "primer_name_forward", "primer_seq_forward", "left_seq_forward", "right_seq_forward", "start_reverse", "end_reverse", 'primer_size_reverse', "errors_reverse", "primer_name_reverse", "primer_seq_reverse", "left_seq_reverse", "right_seq_reverse"]
    wide = wide[ordered_cols]
    wide = wide.merge(df[['accession', 'reverse_complement']], on='accession', how='inner')
    wide = wide.drop_duplicates(subset="accession")
    wide.to_csv(out_path, sep='\t', index=False)


def trim_msa(cutadapt_path: str, meta_path: str, msa_path: str, fasta_out: str):
    cutadapt = pd.read_csv(cutadapt_path, sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    cut_filtered = cutadapt[cutadapt['accession'].isin(meta['accession'])]
    
    references = cut_filtered.dropna(subset=['primer_seq_forward', 'primer_seq_reverse'])
    accessions = references['accession'].tolist()

    alignment = AlignIO.read(msa_path, 'fasta')
    firsts = []
    lasts = []
    if len(accessions) == 0:
        # No references available; skip trimming
        with open(fasta_out, 'w') as fasta:
            for record in alignment:
                #TODO: check again for min_size?
                ungapped_seq = str(record.seq).replace('-', '').upper()
                database.write_seq(fasta, record.id, ungapped_seq)
        return
    
    for sequence in alignment:
        if sequence.id not in accessions:
            continue
        seq = str(sequence.seq)
        # find first non-gap
        first = next(i for i, base in enumerate(seq) if base != "-")
        # find last non-gap
        last = len(seq) - next(i for i, base in enumerate(reversed(seq)) if base != "-") - 1
        firsts.append(first)
        lasts.append(last)
    firsts = pd.Series(firsts)
    lasts = pd.Series(lasts)
    # trim at most frequent position among the reference sequences
    trim_start = int(firsts.value_counts(sort=True, ascending=False).index[0])
    trim_end = int(lasts.value_counts(sort=True, ascending=False).index[0])

    # slice alignment: [rows, columns]
    msa_trimmed = alignment[:, trim_start:trim_end + 1]
    # write fasta file
    with open(fasta_out, 'w') as fasta:
        for record in msa_trimmed:
            #TODO: check again for min_size?
            ungapped_seq = str(record.seq).replace('-', '').upper()
            database.write_seq(fasta, record.id, ungapped_seq)


def _write_primer_fasta(info_path: str, meta_path: str, fasta_path: str, primer_size: int, direction: str):
    info = pd.read_csv(info_path, sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    cut_filtered = info[info['accession'].isin(meta['accession'])]

    with open(fasta_path, 'w') as out:
        for _, row in cut_filtered.iterrows():
            if pd.isna(row[f'primer_size_{direction}']):
                continue
            if row[f'primer_size_{direction}'] > primer_size:
                continue
            seq = str(row[f'primer_seq_{direction}']).upper()
            # make reverse primer reverse complement because cutadapt returns sequence in direction of amplicon
            if direction == 'reverse':
                seq = str(Seq(seq).reverse_complement())
            database.write_seq(out, row['accession'], seq)
    return meta.shape[0]

def _align_primers(fasta: str):
    result = subprocess.run(['mafft', '--auto', fasta], stdout=subprocess.PIPE, stderr=None, check=True, text=True)
    return [str(record.seq).upper() for record in SeqIO.parse(io.StringIO(result.stdout), 'fasta')]

def _draw_logo(seqs: list[str], logo_path: str, total_seqs: int):
    lengths = [len(seq) for seq in seqs]
    counts_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='-N')
    counts_df = counts_df.map(lambda x: x / total_seqs)
    ax = plt.subplot()
    lm.Logo(counts_df, color_scheme='classic', ax=ax)
    plt.xticks(range(1, np.max(lengths), 2))
    plt.xlabel("Base positions")
    plt.ylabel("Proportions")
    plt.ylim((0, 1))
    plt.savefig(logo_path, format='svg', bbox_inches='tight')

def logo_plot(info_path: str, meta_path: str, fasta_path: str, primer_size: int, logo_path: str, direction: str):
    total_seqs = _write_primer_fasta(info_path, meta_path, fasta_path, primer_size, direction)
    seqs = _align_primers(fasta_path)
    _draw_logo(seqs, logo_path, total_seqs)


@dataclass
class BlastJobInfo:
    id: str
    time: str

class PrimerBlastResult(BaseModel):
    accession: str
    uid: int
    product_size: int = Field(validation_alias='product_length')
    forward_template_start: int = Field(validation_alias='template_1_start')
    forward_template_end: int = Field(validation_alias='template_1_end')
    forward_template_sequence: str = Field(validation_alias='template_1_seq')
    reverse_template_start: int = Field(validation_alias='template_2_start')
    reverse_template_end: int = Field(validation_alias='template_2_end')
    reverse_template_sequence: str = Field(validation_alias='template_2_seq')

headers = {
  'accept': '*/*',
  'accept-language': 'de,de-DE;q=0.9,en;q=0.8,en-GB;q=0.7,en-US;q=0.6',
  'cache-control': 'no-cache',
  'origin': 'https://www.ncbi.nlm.nih.gov',
  'pragma': 'no-cache',
  'priority': 'u=1, i',
  'referer': 'https://www.ncbi.nlm.nih.gov/',
  'sec-ch-ua': '"Chromium";v="140", "Not=A?Brand";v="24", "Microsoft Edge";v="140"',
  'sec-ch-ua-mobile': '?0',
  'sec-ch-ua-platform': '"macOS"',
  'sec-fetch-dest': 'empty',
  'sec-fetch-mode': 'no-cors',
  'sec-fetch-site': 'cross-site',
  'sec-fetch-storage-access': 'active',
  'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/140.0.0.0 Safari/537.36 Edg/140.0.0.0'
}


def primer_blast_client():
    return httpx.AsyncClient(
        timeout=httpx.Timeout(30, read=30),
        transport=RetryTransport(retry=Retry(total=5, backoff_factor=0.5)),
        headers=headers
    )

class PrimerBlastJob():
    def __init__(self, client: httpx.AsyncClient, primer_f: str, primer_r: str, tax_ids: list[int], db='core_nt') -> None:
        self._data = [
            ('PRIMER_LEFT_INPUT', primer_f),
            ('PRIMER_RIGHT_INPUT', primer_r),
            ('PRIMER_SPECIFICITY_DATABASE', db),
            ('SEARCH_SPECIFIC_PRIMER', 'on'),
            ('NUM_TARGETS_WITH_PRIMERS', '10000'),
            ('HITSIZE', '100000')
        ]
        for tax_id in tax_ids:
            self._data.append(('ORGANISM', str(tax_id)))
        self._client = client
    
    async def submit(self):
        response = await self._client.post('https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi', content=urlencode(self._data))
        response.raise_for_status()
        self._job_info = self._parse_form_submission(response.text)
        if self._job_info is None:
            raise AttributeError()
        print(f'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time={self._job_info.time}&job_key={self._job_info.id}')
    
    @staticmethod
    def _parse_form_submission(html: str):
        soup = BeautifulSoup(html, "html.parser")
        time_input = soup.find("input", {"name": "ctg_time"})
        job_key_input = soup.find("input", {"name": "job_key"})
        if time_input is None or job_key_input is None:
            # Error handling
            error_html = soup.find('p', {'class': 'error'})
            if error_html is not None:
                error_msg: str = error_html.text
                print(error_msg)
                if error_msg.find('not present in selected database') != -1:
                    return None
            raise ValueError('Failed to parse response: Cannot extract job info!')
        return BlastJobInfo(job_key_input["value"], time_input['value'])

    async def _poll_request(self):
        response = await self._client.get('https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi', params={
            'ctg_time': self._job_info.time,
            'job_key': self._job_info.id
        })
        response.raise_for_status()
        return response
    
    @staticmethod
    def _is_job_done(response: str):
        return response.find('Detailed primer reports') != -1
    
    async def _poll_for_results(self, delay: int):
        response = await self._poll_request()
        while response.status_code == 200 and not self._is_job_done(response.text):
            response = await self._poll_request()
            time.sleep(delay)
        return response.text

    async def get_results(self, delay: int):
        html = await self._poll_for_results(delay)
        return self._parse_results(html)

    @staticmethod
    def _parse_results(html: str):
        def parse_line(line):
            parts = line.split()
            if len(parts) < 4:
                return None
            label = " ".join(parts[:-3])  # handle labels with spaces
            return {
                "label": label,
                "start": int(parts[-3]),
                "sequence": parts[-2],
                "end": int(parts[-1]),
            }
        
        bs = BeautifulSoup(html, "html.parser")
        res_accessions = bs.find_all("a", {"target": "new_entrez"})

        rows: list[PrimerBlastResult] = []

        for accession in res_accessions:
            hit_pre = accession.find_next('pre').text
            product_length = re.search(r"product length = (\d+)", hit_pre).group(1)
            lines = [parse_line(l) for l in hit_pre.splitlines() if l.strip() and "product length" not in l]
            parsed = [l for l in lines if l]

            # Handle duplicates
            label_counter = defaultdict(int)
            row = {"product_length": product_length}

            for entry in parsed:
                base_label = entry["label"].replace(" ", "_").lower()
                label_counter[base_label] += 1
                label = f"{base_label}_{label_counter[base_label]}"  # make unique
                    
                row[f"{label}_start"] = entry["start"]
                row[f"{label}_seq"] = entry["sequence"]
                row[f"{label}_end"] = entry["end"]
            
            # accession
            row['accession'] = str(accession.text)
            acc_link = accession.attrs['href']
            match = re.search(r"id=(\d+)", acc_link)
            row['uid'] = match.group(1)

            parsed_row = PrimerBlastResult(**row)
            rows.append(parsed_row)
        return rows

def export_primer_blast_results(rows: list[dict], file_name: str):
    if len(rows) > 0:
        with open(file_name, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

async def run_primer_blast(client: httpx.AsyncClient, primer_f: str, primer_r: str, organisms: list[int], db='core_nt'):
    try:
        job = PrimerBlastJob(client, primer_f, primer_r, organisms, db)
        await job.submit()
    except AttributeError:
        return
    return await job.get_results(15)

async def find_off_targets(primer_f: str, primer_r: str, target_organism: str, out_file: str):
    client = request.HttpClient()
    tax_id = await taxonomy.get_taxid_by_name(target_organism, client)
    tax_xml = await client.efetch_by_id('taxonomy', str(tax_id))
    tax_info = taxonomy.parse_xml(tax_xml)[0]
    # lineage is ordered from top to bottom, but we want to start from the bottom
    lineage = list(reversed(tax_info.lineage))
    taxon = taxonomy.Taxon(TaxId=tax_info.tax_id, ScientificName=tax_info.name, Rank=tax_info.rank)
    
    field_names = list(PrimerBlastResult.model_json_schema(by_alias=False)['properties'].keys())
    with open(out_file, 'w') as fp:
        csv_writer = csv.DictWriter(fp, fieldnames=field_names)
        csv_writer.writeheader()

        parent_taxids = [str(tax.tax_id) for tax in lineage]
        parents = await taxonomy.by_taxids(parent_taxids)
        sisters: list[int] = []

        # get sister tax ids
        for parent, exclude_child in zip(lineage, [taxon, *lineage[:-1]]):
            parent_data = next((item for item in parents if item['query'][0] == str(parent.tax_id)), None)
            if parent_data is None:
                raise ValueError(f'Missing results for query {parent.tax_id}')
            sisters.extend(taxonomy.get_sisters(parent_data, exclude_child.tax_id))

        # Allow max 10 concurrent jobs
        rate_limit = Semaphore(10)
        tasks = []
        client = primer_blast_client()
        results = []
        failed_tax_ids = []
        failed_tax_ids_lock = threading.Semaphore(1)
        try:
            async def rate_limited(item):
                 async with rate_limit:
                    try:
                        return await run_primer_blast(client, primer_f, primer_r, [item])
                    except (httpx.ConnectTimeout, httpx.RemoteProtocolError):
                        with failed_tax_ids_lock:
                            failed_tax_ids.append(item)

            for i, sister in enumerate(sisters):
                # Random jitter added to the base delay
                jitter = random.uniform(-0.3, 0.3)   # +/- 300ms wiggle
                delay = max(0, i * 1.0 + jitter)     # base 1s per job, jittered
                # Staggered launch
                await asyncio.sleep(delay - (i - 1) if i > 0 else delay)
                task = asyncio.create_task(rate_limited(sister))
                tasks.append(task)
            results = await asyncio.gather(*tasks)
            for item in failed_tax_ids:
                results.append(await run_primer_blast(client, primer_f, primer_r, [item]))
        except Exception as e:
            for sister in results:
                if sister is not None:
                    csv_writer.writerows(map(lambda x: x.model_dump(), sister))
            raise e
        for sister in results:
            if sister is not None:
                csv_writer.writerows(map(lambda x: x.model_dump(), sister))

async def get_outgroup_gb(input_path: str, output_path: str):
    df = pd.read_csv(input_path)
    uids = df['uid'].astype(int).astype(str)
    client = request.HttpClient()
    with open(output_path, 'w') as f:
        await client.efetch_by_ids('nuccore', uids.to_list(), lambda x: f.write(x), ret_type='gb')
