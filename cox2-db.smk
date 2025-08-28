import glob
import os
import asyncio
from src import request, taxonomy, primer


gb_prefix = config['gb_prefix']
db_prefix = f"{config['base_dir']}/{config['db_name']}"

rule all:
    input:
        f"{gb_prefix}/{config['db_name']}.gb",
        f"{db_prefix}/taxonomy_result.xml",
        f"{db_prefix}/taxonomy.tsv",
        f"{db_prefix}/taxonomy-filtered.tsv",
        f"{db_prefix}/taxonomy-backup.csv",
        f"{db_prefix}/seq-meta.csv",
        f"{db_prefix}/seq-meta-filtered.csv",
        f"{db_prefix}/seq-meta-filtered-removed.csv",
        f"{db_prefix}/seq-filtered.fasta",
        f"{db_prefix}/seq-filtered-removed.fasta",
        f"{db_prefix}/seq-trimmed.fasta",
        f"{db_prefix}/seq-trimmed-removed.fasta",
        f"{db_prefix}/seq-meta-trimmed.csv",
        f"{db_prefix}/seq-meta-trimmed-removed.csv",
        f"{db_prefix}/clust-msa.fa",
        f"{db_prefix}/clust-results.tsv",
        f"{db_prefix}/db-clust.csv",
        f"{db_prefix}/db-clust.fasta",
        f"{db_prefix}/db-clust-msa.fasta",
        f"{db_prefix}/{config['db_name']}-meta.csv",
        f"{db_prefix}/{config['db_name']}-qiime.fa",
        f"{db_prefix}/{config['db_name']}-qiime.tsv",
        f"{db_prefix}/{config['db_name']}-sintax.fasta",
        f"{db_prefix}/primer-f.svg",
        f"{db_prefix}/primer-r.svg",


rule download:
    output:
        gb_file = f"{gb_prefix}/{config['db_name']}.gb",
    params:
        query = config['query_str']
    run:
        os.makedirs(gb_prefix, exist_ok=True)
        asyncio.run(request.fetch_sequences(params.query, output.gb_file))

rule taxonomy:
    params:
        tax_query = config['tax_query'],
        backup = config['tax_backup'],
    output:
        tax_file = f"{db_prefix}/taxonomy_result.xml",
        out_tax = f"{db_prefix}/taxonomy.tsv",
        backup = f"{db_prefix}/taxonomy-backup.csv",
    run:
        asyncio.run(taxonomy.download_taxonomy(params.tax_query, output.tax_file))
        taxonomy.parse_taxonomy(output.tax_file, output.out_tax)
        taxonomy.parse_col_tree(params.backup, output.backup)

rule parse_metadata:
    input:
        gb_file = rules.download.output.gb_file,
        script = "parse-gb.py"
    output:
        meta_csv = f"{db_prefix}/seq-meta.csv"
    shell:
        """
        python {input.script} {input.gb_file} {output.meta_csv}
        """
        
rule clean_db:
    input:
        meta = rules.parse_metadata.output.meta_csv,
        script = "db-clean.py"
    output:
        meta = f"{db_prefix}/seq-meta-filtered.csv",
        meta_removed = f"{db_prefix}/seq-meta-filtered-removed.csv",
        fasta = f"{db_prefix}/seq-filtered.fasta",
        fasta_removed = f"{db_prefix}/seq-filtered-removed.fasta",
    shell:
        """
        python {input.script} --in-metadata {input.meta} --out-metadata {output.meta} --out-fasta {output.fasta} --out-metadata-removed {output.meta_removed} --out-fasta-removed {output.fasta_removed}
        """


if "outer_primer_f" in config or "outer_primer_r" in config:
    rule pcr_outer:
        input:
            rules.clean_db.output.fasta
        output:
            f"{db_prefix}/nested.fasta"
        params:
            f=config.get("outer_primer_f", None),
            r=config.get("outer_primer_r", None),
            flags_f = lambda wc: (
                (f"-g '{config['outer_primer_f']['seq']};min_overlap={config['outer_primer_f']['min_len']}'" if config.get("outer_primer_f") else "")
            ),
            flags_r = lambda wc: (
                (f"-a '{config['outer_primer_r']['seq']};min_overlap={config['outer_primer_r']['min_len']}'" if config.get("outer_primer_r") else "")
            )
        shell:
            "cutadapt {params.flags_f} {params.flags_r} -o {output} {input}"


rule trim:
    input:
        fasta = lambda wildcards: (
            f"{db_prefix}/nested.fasta"
            if "outer_primer_f" in config or "outer_primer_r" in config
            else rules.clean_db.output.fasta
        ),
        meta = rules.clean_db.output.meta,
        script = 'trim.py'
    output:
        fasta = f"{db_prefix}/seq-trimmed.fasta",
        meta = f"{db_prefix}/seq-meta-trimmed.csv",
        info = f"{db_prefix}/cutadapt-info.tsv",
        info_out = f"{db_prefix}/cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = config['primer_r']['seq'],
        primer_r_len = config['primer_r']['min_len'],
        min_seq_len = config['min_seq_len'],
        max_seq_len = config['max_seq_len'],
    threads: workflow.cores
    shell:
        """
        cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --overlap {params.primer_f_len} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info} -o {output.fasta} --cores {threads} {input.fasta}
        awk -F'\t' '$2 != -1' {output.info} > tmp && mv tmp {output.info}
        python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}
        """


rule trim_removed:
    input:
        fasta = rules.clean_db.output.fasta_removed,
        meta = rules.clean_db.output.meta_removed,
        script = 'trim.py',
    output:
        fasta = f"{db_prefix}/seq-trimmed-removed.fasta",
        meta = f"{db_prefix}/seq-meta-trimmed-removed.csv",
        info = f"{db_prefix}/removed-cutadapt-info.tsv",
        info_out = f"{db_prefix}/removed-cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = config['primer_r']['seq'],
        primer_r_len = config['primer_r']['min_len'],
        min_seq_len = config['min_seq_len'],
        max_seq_len = config['max_seq_len'],
    threads: workflow.cores
    shell:
        """
        cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --overlap {params.primer_f_len} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info} -o {output.fasta} --cores {threads} {input.fasta}
        awk -F'\t' '$2 != -1' {output.info} > tmp && mv tmp {output.info}
        python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}
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

rule trim_alignment:
    input:
        fasta = rules.dereplicate_filter.output.fasta,
        meta = rules.dereplicate_filter.output.meta,
        cutadapt_info = rules.trim.output.info_out,
        script = "trim-alignment.py"
    output:
        msa = f"{db_prefix}/db-clust-msa.fasta",
        fasta = f"{db_prefix}/{config['db_name']}-qiime.fa",
    params:
        gap_open = config['gap_open']
    threads: workflow.cores
    shell:
        """
        mafft --op {params.gap_open} --thread {threads} {input.fasta} > {output.msa}
        python {input.script} {output.msa} {input.cutadapt_info} {input.meta} {output.fasta}
        """

rule export_db:
    input:
        meta = rules.dereplicate_filter.output.meta,
        tax = rules.taxonomy.output.out_tax,
        backup = rules.taxonomy.output.backup,
        script = 'taxonomy-export.py'
    output:
        meta = f"{db_prefix}/{config['db_name']}-meta.csv",
        tax = f"{db_prefix}/{config['db_name']}-qiime.tsv",
        tax_meta = f"{db_prefix}/taxonomy-filtered.tsv"
    shell:
        """
        python {input.script} {input.meta} {input.tax} {input.backup} {output.meta} {output.tax_meta} {output.tax}
        """

rule export_sintax:
    input:
        fasta = rules.trim_alignment.output.fasta,
        tax = rules.export_db.output.tax,
        script = "export-sintax.py"
    output:
        msa = f"{db_prefix}/{config['db_name']}-msa.fasta",
        db = f"{db_prefix}/{config['db_name']}-sintax.fasta",
    params:
        gap_open = config['gap_open']
    shell:
        """
        mafft --op {params.gap_open} --thread {threads} {input.fasta} > {output.msa}
        python {input.script} {output.msa} {input.tax} {output.db}
        """

rule plot:
    input:
        meta = rules.export_db.output.meta,
        info = rules.trim.output.info_out,
    output:
        fasta_f = f"{db_prefix}/primer-f.fasta",
        fasta_r = f"{db_prefix}/primer-r.fasta",
        plot_f = f"{db_prefix}/primer-f.svg",
        plot_r = f"{db_prefix}/primer-r.svg",
    params:
        primer_f_len = len(config['primer_f']['seq']),
        primer_r_len = len(config['primer_r']['seq']),
    run:
        primer.logo_plot(input.info, input.meta, output.fasta_f, params.primer_f_len, output.plot_f, 'forward')
        primer.logo_plot(input.info, input.meta, output.fasta_r, params.primer_r_len, output.plot_r, 'reverse')
