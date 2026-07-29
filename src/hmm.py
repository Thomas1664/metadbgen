import pandas as pd

TBLOUT_COLUMNS = [
    "target_name", "target_accession", "query_name", "query_accession",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "sq_len", "strand", "e_value", "score", "bias", "description",
]


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


def read_metadata(path: str) -> pd.DataFrame:
    """Read a metadata table, guessing comma vs tab separator from the file extension."""
    sep = ',' if str(path).endswith('.csv') else '\t'
    try:
        return pd.read_csv(path, sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['accession'])


def filter_metadata_by_accession(meta_df: pd.DataFrame, accessions) -> pd.DataFrame:
    return meta_df[meta_df['accession'].isin(accessions)]
