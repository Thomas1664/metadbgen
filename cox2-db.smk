import glob
import os
import asyncio
from src import metadata, request, taxonomy, primer
from Bio.Seq import Seq
import pandas as pd
import csv


gb_prefix = config['gb_prefix']
db_prefix = f"{config['base_dir']}/{config['db_name']}"
sources = [db_prefix, f'{db_prefix}/outgroups']
print(expand("{src}/seq-meta.csv", src=sources))

rule all:
    input:
        f"{db_prefix}/taxonomy-filtered.tsv",
        f"{db_prefix}/taxonomy-backup.csv",
        expand("{src}/seq-meta.csv", src=sources),
        expand("{src}/seq-meta-filtered-removed.csv", src=sources),
        expand("{src}/seq-filtered-removed.fasta", src=sources),
        expand("{src}/seq-trimmed-removed.fasta", src=sources),
        expand("{src}/seq-meta-trimmed-removed.csv", src=sources),
        expand("{src}/clust-msa.fa", src=sources),
        expand("{src}/clust-results.tsv", src=sources),
        expand("{src}/db-clust-msa.fasta", src=sources),
        expand(f"{{src}}/{config['db_name']}-meta.tsv",src=sources),
        expand(f"{{src}}/{config['db_name']}-qiime.fa",src=sources),
        expand(f"{{src}}/{config['db_name']}-qiime.tsv",src=sources),
        expand(f"{{src}}/{config['db_name']}-sintax.fasta", src=sources),
        expand("{src}/primer-f.svg", src=sources),
        expand("{src}/primer-r.svg", src=sources),
        expand("{src}/country-heat.html", src=sources),


rule download:
    output:
        gb_file = f"{db_prefix}/{config['db_name']}.gb",
    params:
        query = config['query_str']
    run:
        os.makedirs(gb_prefix, exist_ok=True)
        asyncio.run(request.fetch_sequences(params.query, output.gb_file))

rule primer_blast_outgroups:
    output:
        out = f"{db_prefix}/outgroups/primer-blast.csv",
    params:
        primer_f = config['primer_f']['seq'],
        primer_r = config['primer_r']['seq'],
        organism = config['organism']
    run:
        asyncio.run(primer.find_off_targets(params.primer_f, params.primer_r, params.organism, output.out))

rule download_outgroups:
    input:
        in_file = rules.primer_blast_outgroups.output.out
    output:
        gb_file = f"{db_prefix}/outgroups/{config['db_name']}.gb",
    run:
        asyncio.run(primer.get_outgroup_gb(input.in_file, output.gb_file))

rule taxonomy:
    params:
        organism = config['organism'],
        backup = config['tax_backup'],
    output:
        tax_file = f"{db_prefix}/taxonomy_result.xml",
        out_tax = f"{db_prefix}/taxonomy.csv",
        backup = f"{db_prefix}/taxonomy-backup.csv",
    run:
        asyncio.run(taxonomy.download_taxonomy(params.organism, output.tax_file))
        with open(output.tax_file, 'r') as f:
            taxa = taxonomy.parse_xml(f.read())
        with open(output.out_tax, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            taxonomy.to_csv(taxa, csv_writer, header=True)
        taxonomy.parse_col_tree(params.backup, output.backup)

rule parse_metadata:
    input:
        gb_file = f"{{src}}/{config['db_name']}.gb",
        script = "parse-gb.py"
    output:
        meta_csv = "{src}/seq-meta.csv"
    shell:
        """
        python {input.script} {input.gb_file} {output.meta_csv}
        """


rule clean_db:
    input:
        meta = rules.parse_metadata.output.meta_csv,
        script = "db-clean.py"
    output:
        meta = "{src}/seq-meta-filtered.csv",
        meta_removed = "{src}/seq-meta-filtered-removed.csv",
        fasta = "{src}/seq-filtered.fasta",
        fasta_removed = "{src}/seq-filtered-removed.fasta",
    shell:
        """
        python {input.script} --in-metadata {input.meta} --out-metadata {output.meta} --out-fasta {output.fasta} --out-metadata-removed {output.meta_removed} --out-fasta-removed {output.fasta_removed}
        """


if "outer_primer_f" in config or "outer_primer_r" in config:
    rule pcr_outer:
        input:
            rules.clean_db.output.fasta
        output:
            "{src}/nested.fasta"
        params:
            flags_f = lambda wc: (
                (f"-g '{config['outer_primer_f']['seq']};min_overlap={config['outer_primer_f']['min_len']}'" if config.get("outer_primer_f") else "")
            ),
            flags_r = lambda wc: (
                (f"-a '{str(Seq(config['outer_primer_r']['seq']).reverse_complement())};min_overlap={config['outer_primer_r']['min_len']}'" if config.get("outer_primer_r") else "")
            )
        shell:
            "cutadapt {params.flags_f} {params.flags_r} -o {output} {input}"


rule trim:
    input:
        fasta = lambda wildcards: (
            "{src}/nested.fasta"
            if "outer_primer_f" in config or "outer_primer_r" in config
            else rules.clean_db.output.fasta
        ),
        meta = rules.clean_db.output.meta,
        script = 'trim.py'
    output:
        fasta = "{src}/seq-trimmed.fasta",
        meta = "{src}/seq-meta-trimmed.csv",
        info = "{src}/cutadapt-info.tsv",
        info_out = "{src}/cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = str(Seq(config['primer_r']['seq']).reverse_complement()),
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
        fasta = "{src}/seq-trimmed-removed.fasta",
        meta = "{src}/seq-meta-trimmed-removed.csv",
        info = "{src}/removed-cutadapt-info.tsv",
        info_out = "{src}/removed-cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = str(Seq(config['primer_r']['seq']).reverse_complement()),
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
        msa = "{src}/clust-msa.fa",
        clust = "{src}/clust-results.tsv"
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
        meta = "{src}/db-clust.csv",
        fasta = "{src}/db-clust.fasta",
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
        msa = "{src}/db-clust-msa.fasta",
        fasta = f"{{src}}/{config['db_name']}-qiime.fa",
        meta = "{src}/db-clust-msa.csv",
    params:
        gap_open = config['gap_open'],
        gap_treshold = config['gappiness_treshold'],
    threads: workflow.cores
    shell:
        """
        mafft --op {params.gap_open} --thread {threads} {input.fasta} > {output.msa}
        python {input.script} {output.msa} {input.cutadapt_info} {input.meta} {output.fasta} {output.meta} {params.gap_treshold}
        """

rule taxonomy_outgroup:
    input:
        meta = f"{db_prefix}/outgroups/db-clust-msa.csv",
    output:
        tax = f"{db_prefix}/outgroups/taxonomy.csv",
    run:
        meta_df = pd.read_csv(input.meta, sep='\t')
        tax_ids = meta_df['tax_id_sanitized'].astype(int).astype(str).tolist()
        client = request.HttpClient()
        asyncio.run(taxonomy.fetch_by_tax_ids(tax_ids, client, output.tax))


rule export_db:
    input:
        meta = rules.trim_alignment.output.meta,
        tax = "{src}/taxonomy.csv",
        backup = rules.taxonomy.output.backup,
        script = 'taxonomy-export.py'
    output:
        meta = f"{{src}}/{config['db_name']}-meta.tsv",
        tax = f"{{src}}/{config['db_name']}-qiime.tsv",
        tax_meta = "{src}/taxonomy-filtered.tsv",
    params:
        # It's not guaranteed this file will be always present
        # because it might be possible that none of the organisms in the database has a host.
        # Therefore, this file is not an output file.
        host_tax = "{src}/host-tax.csv",
    shell:
        """
        python {input.script} {input.meta} {input.tax} {input.backup} {output.meta} {output.tax_meta} {output.tax} {params.host_tax}
        """

rule export_sintax:
    input:
        fasta = rules.trim_alignment.output.fasta,
        tax = rules.export_db.output.tax,
        script = "export-sintax.py"
    output:
        msa = f"{{src}}/{config['db_name']}-msa.fasta",
        db = f"{{src}}/{config['db_name']}-sintax.fasta",
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
        fasta_f = "{src}/primer-f.fasta",
        fasta_r = "{src}/primer-r.fasta",
        plot_f = "{src}/primer-f.svg",
        plot_r = "{src}/primer-r.svg",
        geo_html = "{src}/country-heat.html",
        geo_svg = "{src}/country-heat.svg",
    params:
        primer_f_len = len(config['primer_f']['seq']),
        primer_r_len = len(config['primer_r']['seq']),
    run:
        primer.logo_plot(input.info, input.meta, output.fasta_f, params.primer_f_len, output.plot_f, 'forward')
        primer.logo_plot(input.info, input.meta, output.fasta_r, params.primer_r_len, output.plot_r, 'reverse')
        meta = pd.read_csv(input.meta, sep='\t')
        metadata.plot_country_heatmap(meta, output.geo_html, output.geo_svg)
