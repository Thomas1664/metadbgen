from typing import Any
from src import validate
import pandas as pd
import pytest
import yaml

@pytest.fixture
def config() -> dict[str, Any]:
    with open('config.yml') as config_f:
        return yaml.safe_load(config_f)

@pytest.fixture
def base_paths(config: dict[str, Any]):
    base_paths = [f"{config['base_dir']}/{config['db_name']}"]
    if config['find_outgroups']:
        base_paths.append(f"{config['base_dir']}/{config['db_name']}/outgroups")
    return base_paths

FILES_TO_TEST = [
    ('sequences.fasta', 'seq-meta.csv'),
    ('seq-filtered.fasta', 'seq-meta-filtered.tsv'),
    ('primer-trim/sequences.fasta', 'primer-trim/metadata.tsv'),
    ('dereplicate/sequences.fasta', 'dereplicate/metadata.tsv'),
    ('quality-filter/sequences.fasta', 'quality-filter/metadata.tsv'),
    ('output/{db}-qiime.fa', 'output/{db}-meta.tsv')
]

@pytest.mark.parametrize("fasta_t, meta_t", FILES_TO_TEST)
def test_workflow(base_paths, config, fasta_t: str, meta_t: str):
    for base_path in base_paths:
        fasta = f"{base_path}/" + fasta_t.format(db=config['db_name'])
        meta = f"{base_path}/" + meta_t.format(db=config['db_name'])
        
        validate.fasta(fasta)
        if meta.endswith('.csv'):
            df = pd.read_csv(meta)
        else:
            df = pd.read_csv(meta, sep='\t')
        validate.metadata(df)
        validate.fasta_and_metadata(fasta, meta)
