from src import validate
import pandas as pd
import snakemake


def test_workflow():
    config = snakemake.load_configfile('config.yml')
    files_to_test = [
        ('sequences.fasta', 'seq-meta.csv'),
        ('seq-filtered.fasta', 'seq-meta-filtered.tsv'),
        ('primer-trim/sequences.fasta', 'primer-trim/metadata.tsv'),
        ('dereplicate/sequences.fasta', 'dereplicate/metadata.tsv'),
        ('quality-filter/sequences.fasta', 'quality-filter/metadata.tsv'),
        (f'output/{config["db_name"]}-qiime.fa', f'output/{config["db_name"]}-meta.tsv')
    ]
    base_paths = [f"{config['base_dir']}/{config['db_name']}"]
    if config['find_outgroups']:
        base_paths.append(f"{config['base_dir']}/{config['db_name']}/outgroups")

    for base_path in base_paths:
        for fasta, meta in files_to_test:
            fasta = f"{base_path}/{fasta}"
            validate.fasta(fasta)
            
            meta = f"{base_path}/{meta}"
            if meta.endswith('.csv'):
                df = pd.read_csv(meta)
            else:
                df = pd.read_csv(meta, sep='\t')
            validate.metadata(df)
            validate.fasta_and_metadata(fasta, meta)
