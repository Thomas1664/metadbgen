from asyncio import Semaphore
import asyncio
from collections import defaultdict
from dataclasses import dataclass
from itertools import filterfalse
import random
from httpx_retries import Retry, RetryTransport
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential
from pydantic import AliasChoices, BaseModel, Field
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from bs4 import BeautifulSoup
from src import request, taxonomy, database, cache
from src.taxonomy import SisterTaxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker as lm
import numpy as np
import pandas as pd
import csv
import httpx
import re
import subprocess
import json


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


def trim_msa(cutadapt_path: str, meta_path: str, msa_path: str, fasta_out: str, meta_out: str, min_seq_size: int, align_format='fasta'):
    cutadapt = pd.read_csv(cutadapt_path, sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    cut_filtered = cutadapt[cutadapt['accession'].isin(meta['accession'])]

    references = cut_filtered.dropna(subset=['primer_seq_forward', 'primer_seq_reverse'])
    accessions = set(references['accession'])

    alignment = AlignIO.read(msa_path, align_format)
    accessions_written = []
    if len(accessions) == 0:
        # No references available; skip trimming
        with open(fasta_out, 'w') as fasta:
            for record in alignment:
                ungapped_seq = str(record.seq).replace('-', '').upper()
                if len(ungapped_seq) >= min_seq_size:
                    database.write_seq(fasta, record.id, ungapped_seq)
                    accessions_written.append(record.id)
        meta = meta[meta['accession'].isin(accessions_written)]
        meta.to_csv(meta_out, sep='\t', index=False)
        return

    ref_seqs = [str(seq.seq) for seq in alignment if seq.id in accessions]
    aln_len = len(ref_seqs[0])
    matrix = np.frombuffer(''.join(ref_seqs).encode('ascii'), dtype='S1').reshape(len(ref_seqs), aln_len)
    firsts_arr = np.argmax(matrix != b'-', axis=1)
    lasts_arr = aln_len - 1 - np.argmax(matrix[:, ::-1] != b'-', axis=1)
    # trim at most frequent position among the reference sequences
    trim_start = int(np.bincount(firsts_arr).argmax())
    trim_end = int(np.bincount(lasts_arr).argmax())
    # slice alignment: [rows, columns]
    msa_trimmed = alignment[:, trim_start:trim_end + 1]
    # write fasta file
    with open(fasta_out, 'w') as fasta:
        for record in msa_trimmed:
            ungapped_seq = str(record.seq).replace('-', '').upper()
            if len(ungapped_seq) >= min_seq_size:
                database.write_seq(fasta, record.id, ungapped_seq)
                accessions_written.append(record.id)
    meta = meta[meta['accession'].isin(accessions_written)]
    meta.to_csv(meta_out, sep='\t', index=False)


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

def _align_primers(fasta: str, tmp_out_file: str):
    subprocess.run(['muscle', '-super5', fasta, '-output', tmp_out_file, '-nt'], check=True, text=True)
    return [str(record.seq).upper() for record in SeqIO.parse(tmp_out_file, 'fasta')]

def _draw_logo(seqs: list[str], logo_path: str, total_seqs: int, primer_seq: str = ''):
    lengths = [len(seq) for seq in seqs]
    counts_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='-N')
    counts_df = counts_df.map(lambda x: x / total_seqs)
    width_in = 4.8  # 340px at 96 CSS px/in (SVG web standard)
    height_in = width_in * (4.8 / 6.4)  # preserve default matplotlib aspect ratio
    fig, ax = plt.subplots(figsize=(width_in, height_in))
    lm.Logo(counts_df, color_scheme='classic', ax=ax)
    if primer_seq:
        ax.set_xticks(range(len(primer_seq)))
        ax.set_xticklabels(list(primer_seq.upper()), fontfamily='monospace', fontsize=12)
        ax.set_xlabel("Reference primer", fontsize=14)
    else:
        ax.set_xticks(range(1, np.max(lengths), 2))
        ax.set_xlabel("Base positions", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylabel("Base frequency", fontsize=14)
    ax.set_ylim((0, 1))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0%}'))
    with plt.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(logo_path, format='svg', bbox_inches='tight')
    plt.close(fig)

def logo_plot(info_path: str, meta_path: str, fasta_path: str, primer_size: int, logo_path: str, direction: str, primer_seq: str = ''):
    total_seqs = _write_primer_fasta(info_path, meta_path, fasta_path, primer_size, direction)
    seqs = _align_primers(fasta_path, f'{fasta_path}.tmp')
    _draw_logo(seqs, logo_path, total_seqs, primer_seq)


@dataclass
class BlastJobInfo:
    id: str
    time: str

class PrimerBlastResult(BaseModel):
    accession: str
    uid: int
    product_size: int = Field(validation_alias=AliasChoices('product_length', 'product_size'))
    forward_template_start: int = Field(validation_alias=AliasChoices('template_1_start', 'forward_template_start'))
    forward_template_end: int = Field(validation_alias=AliasChoices('template_1_end', 'forward_template_end'))
    forward_template_sequence: str = Field(validation_alias=AliasChoices('template_1_seq', 'forward_template_sequence'))
    reverse_template_start: int = Field(validation_alias=AliasChoices('template_2_start', 'reverse_template_start'))
    reverse_template_end: int = Field(validation_alias=AliasChoices('template_2_end', 'reverse_template_end'))
    reverse_template_sequence: str = Field(validation_alias=AliasChoices('template_2_seq', 'reverse_template_sequence'))

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
    def __init__(self, client: httpx.AsyncClient, primer_f: str, primer_r: str, tax_id: int, db='core_nt') -> None:
        self._data = {
            'PRIMER_LEFT_INPUT': primer_f,
            'PRIMER_RIGHT_INPUT': primer_r,
            'PRIMER_SPECIFICITY_DATABASE': db,
            'SEARCH_SPECIFIC_PRIMER': 'on',
            'NUM_TARGETS_WITH_PRIMERS': '10000',
            'HITSIZE': '100000',
            # These parameters allow for mor flexibility. Usually, Primer-BLAST is very strict about
            # mismatches in the 3' region but our experience has shown that such sequences will be amplified anyways.
            # This allows finding more unrelated outgroups with more mismatches in the primer region.
            'TOTAL_PRIMER_SPECIFICITY_MISMATCH': '4', # 4th option in selection, actually 5!
            'PRIMER_3END_SPECIFICITY_MISMATCH': '4', # 4th option in selection, actually 5!
            'TOTAL_MISMATCH_IGNORE': '5',
            'WORD_SIZE': '6',
            'CON_ANEAL_OLIGO': '50.0',
            'CON_DNTPS': '0.6',
            'DIVA_CATIONS': '1.5',
            'EVALUE': '30000',
            'GC_CLAMP': '0',
            'LOW_COMPLEXITY_FILTER': 'on',
            'MAX_CANDIDATE_PRIMER': '500',
            'MAX_INTRON_SIZE': '1000000',
            'MAX_TARGET_PER_TEMPLATE': '100',
            'MAX_TARGET_SIZE': '4000',
            'MIN_INTRON_SIZE': '1000',
            'MISMATCH_REGION_LENGTH': '5',
            'MONO_CATIONS': '50.0',
            'NUM_DIFFS': '7',
            'NUM_OPTS_DIFFS': '2',
            'NUM_TARGETS': '20',
            'OVERLAP_3END': '4',
            'OVERLAP_5END': '7',
            'POLYX': '5',
            'PRIMER_INTERNAL_OLIGO_MAX_GC': '80.0',
            'PRIMER_INTERNAL_OLIGO_MAX_SIZE': '27',
            'PRIMER_INTERNAL_OLIGO_MAX_TM': '63.0',
            'PRIMER_INTERNAL_OLIGO_MIN_GC': '20.0',
            'PRIMER_INTERNAL_OLIGO_MIN_SIZE': '18',
            'PRIMER_INTERNAL_OLIGO_MIN_TM': '57.0',
            'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT': '50',
            'PRIMER_INTERNAL_OLIGO_OPT_SIZE': '20',
            'PRIMER_INTERNAL_OLIGO_OPT_TM': '60.0',
            'PRIMER_MAX_DIFF_TM': '3',
            'PRIMER_MAX_END_GC': '5',
            'PRIMER_MAX_END_STABILITY': '9',
            'PRIMER_MAX_GC': '80.0',
            'PRIMER_MAX_HAIRPIN_TH': '24.0',
            'PRIMER_MAX_SELF_ANY_TH': '45.0',
            'PRIMER_MAX_SELF_END_TH': '35.0',
            'PRIMER_MAX_SIZE': '25',
            'PRIMER_MAX_TEMPLATE_MISPRIMING': '12.00',
            'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': '40.00',
            'PRIMER_MAX_TM': '63.0',
            'PRIMER_MIN_GC': '20.0',
            'PRIMER_MIN_SIZE': '15',
            'PRIMER_MIN_TM': '57.0',
            'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
            'PRIMER_NUM_RETURN': '10',
            'PRIMER_ON_SPLICE_SITE': '0',
            'PRIMER_OPT_SIZE': '20',
            'PRIMER_OPT_TM': '60.0',
            'PRIMER_PAIR_MAX_COMPL_ANY': '8.00',
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': '45.0',
            'PRIMER_PAIR_MAX_COMPL_END': '3.00',
            'PRIMER_PAIR_MAX_COMPL_END_TH': '35.0',
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': '24.00',
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': '70.00',
            'PRIMER_PRODUCT_MAX': '1000',
            'PRIMER_PRODUCT_MIN': '70',
            'SALT_FORMULAR': '1',
            'SEARCHMODE': '0',
            'SELF_ANY': '8.00',
            'SELF_END': '3.00',
            'SHOW_SVIEWER': 'on',
            'SPLICE_SITE_OVERLAP_3END': '4',
            'SPLICE_SITE_OVERLAP_3END_MAX': '8',
            'SPLICE_SITE_OVERLAP_5END': '7',
            'TM_METHOD': '1',
            'UNGAPPED_BLAST': 'on',
            'USER_TYPE': '2',
            'CMD': 'request',
        }
        self._data['ORGANISM'] = str(tax_id)
        self._client = client

    async def submit(self):
        response = await self._client.post('https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi', data=self._data)
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
            await asyncio.sleep(delay)
        return response.text

    async def get_results(self, delay: int):
        html = await self._poll_for_results(delay)
        return self._parse_results(html)

    @staticmethod
    def _parse_results(html: str):
        def parse_line(line: str):
            if line.startswith('Features associated with this product'):
                return None
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

        for accession in filterfalse(lambda accession: accession.parent.name == 'pre', res_accessions):
            try:
                hit_pre = accession.find_next('pre').text
            except AttributeError as ex:
                print(accession.text)
                raise ex
            product_length = re.search(r"product length = (\d+)", hit_pre).group(1)
            lines = [parse_line(l) for l in hit_pre.splitlines() if l.strip() and "product length" not in l and "Features flanking this product" not in l and "bp at" not in l]
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
            if not match:
                match = re.search(r"nucleotide/(\d+)\?", acc_link)
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

def blast_cache_key(_, primer_f, primer_r, organism, db="core_nt"):
    return json.dumps((primer_f, primer_r, organism, db))

@retry(
    retry=retry_if_exception_type((httpx.ConnectTimeout, httpx.RemoteProtocolError, httpx.ReadError)),
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=5, max=40),
    reraise=True
)
@cache.decorator(key_func=blast_cache_key, model_cls=PrimerBlastResult)
async def run_primer_blast(client: httpx.AsyncClient, primer_f: str, primer_r: str, organism: int, db='core_nt'):
    try:
        job = PrimerBlastJob(client, primer_f, primer_r, organism, db)
        await job.submit()
    except AttributeError:
        return
    return await job.get_results(15)


async def _off_targets_in_sisters(sisters: set[SisterTaxon], primer_f: str, primer_r: str):
    # Allow max 10 concurrent jobs
    rate_limit = Semaphore(10)
    tasks = []
    client = primer_blast_client()
    results: list[tuple[SisterTaxon, list[PrimerBlastResult] | None]] = []
    async def rate_limited(item: SisterTaxon):
        async with rate_limit:
            await asyncio.sleep(random.uniform(5, 20))
            return (item, await run_primer_blast(client, primer_f, primer_r, item.taxon.tax_id))

    tasks = [asyncio.create_task(rate_limited(sister)) for sister in sisters]
    try:
        results = await asyncio.gather(*tasks)
    finally:
        for t in tasks:
            t.cancel()
    return results


async def find_off_targets(primer_f: str, primer_r: str, target_organism: str, out_file: str, use_cache: bool, cache_path: str, tax_db_path: str):
    if use_cache:
        print(f'Using cache at {cache_path}')
        cache.cache_var.set(cache.AsyncSQLiteCache(cache_path))
    sisters = taxonomy.get_sister_taxa(tax_db_path, target_organism)
    print(sisters)

    with open(out_file, 'w') as fp:
        # Write CSV header
        field_names = list(PrimerBlastResult.model_json_schema(by_alias=False)['properties'].keys())
        csv_writer = csv.DictWriter(fp, fieldnames=field_names)
        csv_writer.writeheader()

        # Run Primer-BLAST on all sister taxa
        results = await _off_targets_in_sisters(sisters, primer_f, primer_r)
        try_below: list[SisterTaxon] = []
        for item, sister in results:
            if sister is not None and len(sister) != 0:
                csv_writer.writerows(map(lambda x: x.model_dump(), sister))
            else:
                # Try one level below if we didn't get any results.
                # This is a limitation of Primer-BLAST if it needs to go through many target sequences.
                if item.parent.rank == 'domain':
                    print(f'Nothing found in {item.taxon.name}(taxid {item.taxon.tax_id}); trying one level below...')
                    try_below.append(item)

        if len(try_below) == 0:
            return # Nothing to do

        # Run Primer-BLAST on children of sister taxa
        children = taxonomy.get_sister_childs(tax_db_path, try_below)
        results = await _off_targets_in_sisters(children, primer_f, primer_r)
        for item, sister in results:
            if sister is not None and len(sister) != 0:
                csv_writer.writerows(map(lambda x: x.model_dump(), sister))


async def get_outgroup_gb(input_path: str, output_path: str):
    df = pd.read_csv(input_path)
    uids = df['uid'].astype(int).astype(str)
    client = request.HttpClient()
    with open(output_path, 'w') as f:
        await client.efetch_by_ids('nuccore', uids.to_list(), lambda x: f.write(x), ret_type='gb')
