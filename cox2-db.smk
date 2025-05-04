import glob
import os
import asyncio
from src import request


gb_prefix = config['gb_prefix']
db_prefix = f"{config['base_dir']}/{config['db_name']}"

rule all:
    input:
        gb_prefix,
        f"{db_prefix}/taxonomy-filtered.tsv",
        f"{db_prefix}/seq-meta.csv",
        f"{db_prefix}/seq-meta-filtered.csv",
        f"{db_prefix}/seq-filtered.fasta",
        f"{db_prefix}/seq-trimmed.fasta",
        f"{db_prefix}/seq-meta-trimmed.csv",
        f"{db_prefix}/clust-msa.fa",
        f"{db_prefix}/clust-results.tsv",
        f"{db_prefix}/db-clust.csv",
        f"{db_prefix}/db-clust.fasta",
        f"{db_prefix}/{config['db_name']}-meta.csv",
        f"{db_prefix}/{config['db_name']}-qiime.fa",
        f"{db_prefix}/{config['db_name']}-qiime.tsv",


rule download:
    output:
        gb_files = directory(gb_prefix)
    params:
        query = config['query_str']
    run:
        os.makedirs(gb_prefix, exist_ok=True)
        asyncio.run(request.fetch_sequences(params.query, f'{gb_prefix}/cox2.gb'))

rule taxonomy:
    input:
        script = "taxonomy.py"
    params:
        tax_file = config['tax_file'],
    output:
        out_tax = f"{db_prefix}/taxonomy.tsv",
        filtered_tax = f"{db_prefix}/taxonomy-filtered.tsv"
    shell:
        """
        python {input.script} {params.tax_file} {output.out_tax} {output.filtered_tax}
        """

rule parse_metadata:
    input:
        gb_files = rules.download.output.gb_files,
        script = "parse-gb.py"
    output:
        meta_csv = f"{db_prefix}/seq-meta.csv"
    shell:
        """
        python {input.script} {input.gb_files} {output.meta_csv}
        """
        
rule clean_db:
    input:
        meta = rules.parse_metadata.output.meta_csv,
        script = "db-clean.py"
    output:
        meta = f"{db_prefix}/seq-meta-filtered.csv",
        fasta = f"{db_prefix}/seq-filtered.fasta",
    shell:
        """
        python {input.script} --in-metadata {input.meta} --out-metadata {output.meta} --out-fasta {output.fasta}
        """

rule trim:
    input:
        fasta = rules.clean_db.output.fasta,
        meta = rules.clean_db.output.meta,
        script = 'trim.py'
    output:
        fasta_f = f"{db_prefix}/seq-trimmed_f.fasta",
        fasta_tmp = f"{db_prefix}/seq-trimmed_tmp.fasta",
        fasta = f"{db_prefix}/seq-trimmed.fasta",
        meta = f"{db_prefix}/seq-meta-trimmed.csv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = config['primer_r']['seq'],
        primer_r_len = config['primer_r']['min_len'],
        min_seq_len = config['min_seq_len']
    shell:
        """
        cutadapt -g {params.primer_f_seq} --revcomp --overlap {params.primer_f_len} -e 0.2 -o {output.fasta_f} {input.fasta}
        cutadapt -a {params.primer_r_seq} --revcomp --overlap {params.primer_r_len} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} -o {output.fasta_tmp} {output.fasta_f}
        cutadapt -g TTCAAGATCC --revcomp --overlap 8 -e 0.0 -o {output.fasta} {output.fasta_tmp}
        python {input.script} {input.meta} {output.fasta} {output.meta}
        """

rule dereplicate:
    input:
        fasta = rules.trim.output.fasta,
    output:
        msa = f"{db_prefix}/clust-msa.fa",
        clust = f"{db_prefix}/clust-results.tsv"
    shell:
        """
        vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}
        """

rule dereplicate_filter:
    input:
        meta = rules.trim.output.meta,
        clust = rules.dereplicate.output.clust,
        script = 'dereplicate.py'
    output:
        meta = f"{db_prefix}/db-clust.csv",
        fasta = f"{db_prefix}/db-clust.fasta",
    shell:
        """
        python {input.script} {input.meta} {input.clust} {output.meta} {output.fasta}
        """

rule export_db:
    input:
        meta = rules.dereplicate_filter.output.meta,
        tax = rules.taxonomy.output.filtered_tax,
        script = 'taxonomy-export.py'
    output:
        meta = f"{db_prefix}/{config['db_name']}-meta.csv",
        fasta = f"{db_prefix}/{config['db_name']}-qiime.fa",
        tax = f"{db_prefix}/{config['db_name']}-qiime.tsv",
    shell:
        """
        python {input.script} {input.meta} {input.tax} {output.fasta} {output.meta} {output.tax}
        """
