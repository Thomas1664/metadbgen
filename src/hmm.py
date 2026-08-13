from collections.abc import Sequence
from itertools import filterfalse
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
from src import request
import glob
import os

TBLOUT_COLUMNS = [
    "target_name", "target_accession", "query_name", "query_accession",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "sq_len", "strand", "e_value", "score", "bias", "description",
]


# Relevant GenBank feature annotations containing possible keywords
# Source: https://www.insdc.org/submitting-standards/feature-table/#7.2
GENBANK_FEATURES = {
    'CDS': ['gene', 'operon', 'product'],
    'gene': ['gene', 'operon', 'product'],
    'misc_feature': ['gene', 'product'],
    'mRNA': ['gene', 'operon', 'product'],
    'operon': ['operon'],
    'prim_transcript': ['gene', 'operon'],
    'propeptide': ['gene', 'product'],
    'rRNA': ['gene', 'operon', 'product'],
}


def parse_hmm_tblout(path: str) -> pd.DataFrame:
    """Parse an nhmmer --tblout file, keeping only the best-scoring hit per target."""
    try:
        df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, engine='python', names=TBLOUT_COLUMNS)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=TBLOUT_COLUMNS)
    return df.sort_values('score', ascending=False).drop_duplicates(subset=['target_name'], keep='first')


def compute_score_threshold(tblout_path: str, margin: float) -> tuple[float, float]:
    """Return (min_score, threshold) for a reference-set tblout; threshold = min_score - margin."""
    df = parse_hmm_tblout(tblout_path)
    min_score = df['score'].min()
    return min_score, min_score - margin


def write_threshold_file(path: str, threshold: float) -> None:
    with open(path, 'w') as f:
        f.write(f"{threshold}\n")


def read_threshold_file(path: str) -> float:
    with open(path) as f:
        return float(f.read().strip())


def passing_hits(tblout_path: str, threshold: float) -> pd.DataFrame:
    """Rows from a tblout whose score meets the given threshold."""
    df = parse_hmm_tblout(tblout_path)
    return df[df['score'].ge(threshold)]


def split_by_strand(hits: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Split hit target_name accessions into (plus_strand, minus_strand) based on the 'strand' column."""
    plus = hits.loc[hits['strand'] == '+', 'target_name']
    minus = hits.loc[hits['strand'] == '-', 'target_name']
    return plus, minus


def _is_hit_in_feature(hit: pd.Series, feature: SeqFeature) -> bool:
    feat_start = int(feature.location.start)
    feat_end = int(feature.location.end)
    return feat_start <= hit['ali_from'] <= feat_end or feat_start <= hit['ali_to'] <= feat_end


def find_additional_keywords(hmm_hits: pd.DataFrame, gb_records: Sequence[SeqRecord], existing_keywords: Sequence[str]) -> set[str]:
    existing_keywords = list(map(lambda kw: kw.lower(), existing_keywords))
    keywords: set[str] = set()
    for record in filter(lambda record: record.id in hmm_hits['target_name'].values, gb_records):
        # Find HMM hit location
        hit: pd.Series = hmm_hits.loc[hmm_hits['target_name'] == record.id].iloc[0]
        features_in_hit = filter(
            lambda feature: feature.type in GENBANK_FEATURES.keys() and _is_hit_in_feature(hit, feature),
            record.features
        )
        for feature in features_in_hit:
            # Add feature name to keyword list
            for field in GENBANK_FEATURES[feature.type]:
                keywords.update(feature.qualifiers.get(field, []))
    return set(filterfalse(lambda x: x in existing_keywords, set(x.lower() for x in keywords)))


def build_query_additional_keywords(new_kw: Sequence[str], exclude_kw: Sequence[str], organism: str):
    exclude = " NOT ".join([f'"{keyword}"' for keyword in exclude_kw])
    q_new_kw = " OR ".join([f'"{keyword}"' for keyword in new_kw])
    return f'({q_new_kw}) NOT {exclude} AND "{organism}"[Organism] AND ("0001/01/01"[PDAT] : "2026/04/15"[PDAT])'


async def download_additional_keywords(query: str, out_dir: str):
    client = request.HttpClient()
    eserch_res = await client.esearch(query)
    if eserch_res.result_count == 0:
        print(f'No additional records found using query: {query}')
        with open(f'{out_dir}/input-sequences.gb', 'w') as f:
            f.write('')
        return
    await client.efetch_all(eserch_res, f'{out_dir}/input-sequences.gb')
    for split_file in glob.glob(f'{out_dir}/input-sequences_*.gb'):
        os.remove(split_file)


def read_metadata(path: str) -> pd.DataFrame:
    """Read a metadata table, guessing comma vs tab separator from the file extension."""
    sep = ',' if str(path).endswith('.csv') else '\t'
    try:
        return pd.read_csv(path, sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['accession'])


def filter_metadata_by_accession(meta_df: pd.DataFrame, accessions) -> pd.DataFrame:
    return meta_df[meta_df['accession'].isin(accessions)]
