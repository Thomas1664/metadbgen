import glob
import os
import asyncio
from src import metadata, request, taxonomy, primer, database, mislabel_detection
from src.ncbi_taxonomy import download_and_build_db
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import csv
import io
from concurrent.futures import ThreadPoolExecutor


db_prefix = f"{config['base_dir']}/{config['db_name']}"

sources = [db_prefix]
if config['find_outgroups']:
    sources.append(f'{db_prefix}/outgroups')

rule all:
    input:
        # Taxonomy
        f"{db_prefix}/taxonomy/parsed.csv",

        # Metadata plots
        expand("{src}/primer-f.fasta", src=sources),
        expand("{src}/primer-r.fasta", src=sources),
        expand("{src}/primer-f.svg", src=sources),
        expand("{src}/primer-r.svg", src=sources),
        expand("{src}/metadata/country-heat.html", src=sources),
        expand("{src}/metadata/country-heat.svg", src=sources),

        # Dereplicate
        expand("{src}/dereplicate/clust-msa.fa", src=sources),

        # Tree
        expand(f"{{src}}/tree/{config['db_name']}-tree-dist.csv", src=sources),

        # Misslabel detection
        expand(f"{{src}}/output/{config['db_name']}-meta-clean.tsv", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-qiime-clean.tsv", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-qiime-clean.fa", src=sources),
        expand(f"{{src}}/mislabel-detection.csv", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-accessions.txt", src=sources),

        # Outputs
        expand(f"{{src}}/output/{config['db_name']}-meta.tsv", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-qiime.fa", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-qiime.tsv", src=sources),
        expand(f"{{src}}/output/{config['db_name']}-sintax.fasta", src=sources),

        # Primer trimming
        expand("{src}/primer-trim/removed-sequences.fasta", src=sources),
        expand("{src}/primer-trim/removed-metadata.tsv", src=sources),
        expand("{src}/primer-trim/removed-cutadapt-info.tsv", src=sources),
        expand("{src}/primer-trim/removed-cutadapt-info-parsed.tsv", src=sources),

        expand("{src}/seq-meta-filtered-removed.tsv", src=sources),
        expand("{src}/seq-filtered-removed.fasta", src=sources),
        expand("{src}/seq-filtered-removed.txt", src=sources),
        expand("{src}/alignment/msa.fasta", src=sources),


rule download:
    retries: 5
    output:
        gb_file = f"{db_prefix}/input-sequences.gb",
    params:
        keywords = config['keywords'],
        organism = config['organism']
    run:
        query = database.build_search_query(params.keywords, params.organism)
        with ThreadPoolExecutor() as pool:
            future = pool.submit(lambda: asyncio.run(request.fetch_sequences(query, output.gb_file)))
            future.result()

rule primer_blast_outgroups:
    retries: 5
    output:
        out = f"{db_prefix}/outgroups/primer-blast.csv",
    params:
        primer_f = config['primer_f']['seq'],
        primer_r = config['primer_r']['seq'],
        organism = config['organism'],
        use_cache = config['primer_blast']['enable_cache'],
        cache_path = config['primer_blast']['cache_file'],
    run:
        with ThreadPoolExecutor() as pool:
            future = pool.submit(lambda: asyncio.run(primer.find_off_targets(params.primer_f, params.primer_r, params.organism, output.out, params.use_cache, params.cache_path)))
            future.result()


rule tax_db_download:
    output:
        tax_db = f"{db_prefix}/taxonomy/tax_db.sqlite",
    run:
        download_and_build_db(db_path=output.tax_db, taxdump_dir=f"{db_prefix}/ncbi_tax_dump")


rule download_outgroups:
    retries: 5
    input:
        in_file = rules.primer_blast_outgroups.output.out,
        tax_db = rules.tax_db_download.output.tax_db,
    output:
        gb_file = f"{db_prefix}/outgroups/input-sequences.gb",
    run:
        with ThreadPoolExecutor() as pool:
            future = pool.submit(lambda: asyncio.run(primer.get_outgroup_gb(input.in_file, output.gb_file)))
            future.result()


rule taxonomy:
    input:
        tax_db = rules.tax_db_download.output.tax_db,
    params:
        organism = config['organism'],
    output:
        tax = f"{db_prefix}/taxonomy/parsed.csv",
    run:
        with ThreadPoolExecutor() as pool:
            future = pool.submit(lambda: asyncio.run(taxonomy.download_taxonomy(params.organism, input.tax_db)))
            tax_data = future.result()
        with open(output.tax, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            taxonomy.to_csv(tax_data, csv_writer, header=True)

rule parse_metadata:
    input:
        gb_file = f"{{src}}/input-sequences.gb",
        tax_db = rules.tax_db_download.output.tax_db,
        script = "parse-gb.py"
    output:
        meta_csv = "{src}/seq-meta.csv",
        fasta_tmp = "{src}/sequences.fasta.tmp",
        fasta = "{src}/sequences.fasta",
    shell:
        """
        python {input.script} {input.gb_file} {input.tax_db} {output.meta_csv} {output.fasta_tmp}
        seqkit rmdup -n -o {output.fasta} {output.fasta_tmp}
        """


if config['find_outgroups']:
    rule outgroup_tax:
        input:
            meta = f"{db_prefix}/outgroups/seq-meta.csv",
            tax_db = rules.tax_db_download.output.tax_db,
        output:
            tax = f"{db_prefix}/outgroups/taxonomy/parsed.csv",
        run:
            with taxonomy.NCBITaxonomyDB(input.tax_db) as db:
                tax_ids = pd.read_csv(input.meta)['tax_id'].unique().tolist()
                tax_data = db.get_taxon_data_bulk(tax_ids)
            with open(output.tax, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                taxonomy.to_csv(tax_data, csv_writer, header=True)


rule clean_db:
    input:
        meta = rules.parse_metadata.output.meta_csv,
        fasta = rules.parse_metadata.output.fasta,
        tax_db = rules.tax_db_download.output.tax_db,
        tax = "{src}/taxonomy/parsed.csv",
        script = "db-clean.py"
    output:
        meta = "{src}/seq-meta-filtered.tsv",
        meta_removed = "{src}/seq-meta-filtered-removed.tsv",
        fasta = "{src}/seq-filtered.fasta",
        fasta_removed = "{src}/seq-filtered-removed.fasta",
        accessions = "{src}/seq-filtered.txt",
        accessions_removed = "{src}/seq-filtered-removed.txt",
    params:
        organism = lambda wildcards: f"--organism '{config['organism']}'" if wildcards.src == db_prefix else ""
    shell:
        """
        python {input.script} \\
            --in-metadata {input.meta} \\
            --in-fasta {input.fasta} \\
            --in-tax {input.tax} \\
            --out-metadata {output.meta} \\
            --out-fasta {output.accessions} \\
            --out-metadata-removed {output.meta_removed} \\
            --out-fasta-removed {output.accessions_removed} \\
            --tax-db {input.tax_db} {params.organism}
        seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}
        seqkit grep -f {output.accessions_removed} -o {output.fasta_removed} {input.fasta}
        """


if "outer_primer_f" in config or "outer_primer_r" in config:
    rule pcr_outer:
        input:
            rules.clean_db.output.fasta
        output:
            "{src}/primer-trim/nested.fasta"
        params:
            flags_f = lambda wc: (
                (f"-g '{config['outer_primer_f']['seq']};min_overlap={config['outer_primer_f']['min_len']}'" if config.get("outer_primer_f") else "")
            ),
            flags_r = lambda wc: (
                (f"-a '{str(Seq(config['outer_primer_r']['seq']).reverse_complement())};min_overlap={config['outer_primer_r']['min_len']}'" if config.get("outer_primer_r") else "")
            )
        run:
            os.makedirs(f'{wildcards.src}/primer-trim', exist_ok=True)
            shell("cutadapt {params.flags_f} {params.flags_r} -o {output} {input}")


rule trim:
    input:
        fasta = lambda wildcards: (
            "{src}/primer-trim/nested.fasta"
            if "outer_primer_f" in config or "outer_primer_r" in config
            else rules.clean_db.output.fasta
        ),
        meta = rules.clean_db.output.meta,
        script = "trim.py",
    output:
        fasta = "{src}/primer-trim/sequences.fasta",
        meta = "{src}/primer-trim/metadata.tsv",
        info_tmp = "{src}/primer-trim/cutadapt-info-tmp.tsv",
        info = "{src}/primer-trim/cutadapt-info.tsv",
        info_out = "{src}/primer-trim/cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = str(Seq(config['primer_r']['seq']).reverse_complement()),
        primer_r_len = config['primer_r']['min_len'],
        min_seq_len = config['min_seq_len'],
        max_seq_len = config['max_seq_len'],
    threads: workflow.cores / 2
    run:
        os.makedirs(f'{wildcards.src}/primer-trim', exist_ok=True)
        max_seq_len = pd.read_csv(io.StringIO(shell("seqkit stats -T {input.fasta}", read=True)), sep='\t')['max_len'][0]
        buffer_size = max(max_seq_len * 2, 1048576)
        shell("cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --overlap {params.primer_f_len} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info_tmp} -o {output.fasta} --cores {threads} --buffer-size " + str(buffer_size) + " {input.fasta}")
        shell("python filter_cutadapt_info.py {output.info_tmp} {output.info}")
        shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


rule trim_removed:
    input:
        fasta = rules.clean_db.output.fasta_removed,
        meta = rules.clean_db.output.meta_removed,
        script = "trim.py",
    output:
        fasta = "{src}/primer-trim/removed-sequences.fasta",
        meta = "{src}/primer-trim/removed-metadata.tsv",
        info_tmp = "{src}/primer-trim/removed-cutadapt-info-tmp.tsv",
        info = "{src}/primer-trim/removed-cutadapt-info.tsv",
        info_out = "{src}/primer-trim/removed-cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f']['seq'],
        primer_f_len = config['primer_f']['min_len'],
        primer_r_seq = str(Seq(config['primer_r']['seq']).reverse_complement()),
        primer_r_len = config['primer_r']['min_len'],
        min_seq_len = config['min_seq_len'],
        max_seq_len = config['max_seq_len'],
    threads: workflow.cores / 2
    run:
        shell("cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --overlap {params.primer_f_len} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info_tmp} -o {output.fasta} --cores {threads} --buffer-size 104857600 {input.fasta}")
        shell("python filter_cutadapt_info.py {output.info_tmp} {output.info}")
        shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


rule dereplicate:
    input:
        meta = rules.trim.output.meta,
        fasta = rules.trim.output.fasta,
        script = 'dereplicate.py',
    output:
        msa = "{src}/dereplicate/clust-msa.fa",
        clust = "{src}/dereplicate/clust-results.tsv",
        meta = "{src}/dereplicate/metadata.tsv",
        fasta = "{src}/dereplicate/sequences.fasta",
    shell:
        """
        vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}
        python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta}
        """


rule sequence_quality_filter:
    input:
        fasta = rules.dereplicate.output.fasta,
        meta = rules.dereplicate.output.meta,
    params:
        max_ambiguous_bases = config['max_ambiguous_bases']
    output:
        meta = "{src}/quality-filter/metadata.tsv",
        fasta = "{src}/quality-filter/sequences.fasta",
    run:
        composition = {}
        records = list(SeqIO.parse(input.fasta, 'fasta'))
        for record in records:
            composition[record.id] = database.seq_base_composition(str(record.seq))
        df = pd.DataFrame.from_dict(composition, orient='index')
        df = df[df['other'] <= params.max_ambiguous_bases]
        meta = pd.read_csv(input.meta, sep='\t')
        meta = meta[meta['accession'].isin(df.index)]
        meta.to_csv(output.meta, sep='\t', index=False)
        database.write_fasta(input.fasta, output.fasta, meta['accession'])

rule trim_alignment:
    input:
        fasta = rules.sequence_quality_filter.output.fasta,
        meta = rules.sequence_quality_filter.output.meta,
        cutadapt_info = rules.trim.output.info_out,
        script = "trim-alignment.py"
    output:
        msa = "{src}/alignment/msa.fasta",
        gap_scores = "{src}/alignment/msa-gap-scores.csv",
        fasta_tmp = "{src}/alignment/seq-trimmed.fasta",
        fasta = "{src}/alignment/sequences.fasta",
        meta = "{src}/alignment/metadata-trim.tsv",
    params:
        gap_threshold = config['gappiness_threshold'],
        min_length = config['min_seq_len'],
    threads: workflow.cores / 2
    run:
        os.makedirs(f'{wildcards.src}/alignment', exist_ok=True)
        os.makedirs(f'{wildcards.src}/output', exist_ok=True)
        shell("muscle -super5 {input.fasta} -nt -threads {threads} -output {output.msa}")
        shell("python {input.script} {output.msa} {input.cutadapt_info} {input.meta} {output.fasta_tmp} {output.fasta} {output.meta} {wildcards.src}/alignment {params.gap_threshold} {params.min_length}")


rule dereplicate_alignment:
    input:
        meta = rules.trim_alignment.output.meta,
        fasta = rules.trim_alignment.output.fasta,
        script = 'dereplicate.py',
    output:
        msa = "{src}/alignment/clust-msa.fa",
        clust = "{src}/alignment/clust-results.tsv",
        meta = "{src}/alignment/metadata.tsv",
        fasta = f"{{src}}/output/{config['db_name']}-qiime.fa",
    shell:
        """
        vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}
        python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta}
        """


rule export_db:
    input:
        meta = rules.dereplicate_alignment.output.meta,
        tax_db = rules.tax_db_download.output.tax_db,
        script = 'taxonomy-export.py'
    output:
        meta = f"{{src}}/output/{config['db_name']}-meta.tsv",
        tax = f"{{src}}/output/{config['db_name']}-qiime.tsv",
    params:
        # It's not guaranteed this file will be always present
        # because it might be possible that none of the organisms in the database has a host.
        # Therefore, this file is not an output file.
        host_tax = "{src}/metadata/host-tax.csv",
    run:
        os.makedirs(f'{wildcards.src}/metadata', exist_ok=True)
        shell("python {input.script} {input.meta} {input.tax_db} {output.meta} {output.tax} {params.host_tax}")


rule tree:
    input: rules.dereplicate_alignment.output.fasta,
    params:
        max_iterations = config['tree']['max_iterations']
    output: 
        msa = f"{{src}}/tree/{config['db_name']}-msa.fasta",
        tree = f"{{src}}/tree/{config['db_name']}-msa.fasta.treefile",
        tree_dist = f"{{src}}/tree/{config['db_name']}-tree-dist.csv",
    threads: workflow.cores / 2
    run: 
        shell("muscle -super5 {input} -nt -threads {threads} -output {output.msa}")
        shell("iqtree -s {output.msa} -m 'GTR+G' --fast -n {params.max_iterations} -nt {threads} --redo")
        shell("gotree matrix -i {output.tree} -o {output.dist_mat} -t {threads}")


rule mislabel_detection:
    input:
        tree_dist = rules.tree.output.tree_dist,
        tax = rules.export_db.output.tax,
        meta = rules.export_db.output.meta,
        fasta = rules.dereplicate_alignment.output.fasta
    params:
        mislabel = config['mislabel_detection'],
    output:
        meta = f"{{src}}/output/{config['db_name']}-meta-clean.tsv",
        tax = f"{{src}}/output/{config['db_name']}-qiime-clean.tsv",
        fasta = f"{{src}}/output/{config['db_name']}-qiime-clean.fa",
        mislabel = f"{{src}}/mislabel-detection.csv",
        accessions = f"{{src}}/output/{config['db_name']}-accessions.txt",
        removed = f"{{src}}/mislabeled-sequences.fasta",
    run:
        tax_df = taxonomy.load_tax_qiime(input.tax).rename_axis('seq_id').reset_index()
        mislabel_df = mislabel_detection.validate_at_rank(tax_df, input.tree_dist, **params.mislabel)
        seq_ids_to_keep = mislabel_df[~mislabel_df['is_mislabeled'].values]['seq_id']
        meta_df = pd.read_csv(input.meta, sep='\t')
        meta_df = meta_df[meta_df['accession'].isin(seq_ids_to_keep)]
        meta_df.to_csv(output.meta, sep='\t', index=False)
        qiime_df = pd.read_csv(input.tax, sep='\t', header=None)
        qiime_df = qiime_df[qiime_df[0].isin(seq_ids_to_keep)]
        qiime_df.to_csv(output.tax, sep='\t', index=False, header=False)
        seq_ids_to_keep.to_csv(output.accessions, header=False, index=False)
        mislabel_df.to_csv(output.mislabel, index=False)
        shell("seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}")
        shell("seqkit grep -v -f {output.accessions} -o {output.removed} {input.fasta}")


rule export_sintax:
    input:
        fasta = rules.mislabel_detection.output.fasta,
        tax = rules.mislabel_detection.output.tax,
        script = "export-sintax.py"
    output:
        msa = f"{{src}}/{config['db_name']}-msa.fasta",
        db = f"{{src}}/output/{config['db_name']}-sintax.fasta",
    threads: workflow.cores / 2
    run:
        shell("muscle -super5 {input.fasta} -nt -threads {threads} -output {output.msa}")
        shell("python {input.script} {output.msa} {input.tax} {output.db}")

rule plot:
    input:
        meta = rules.mislabel_detection.output.meta,
        info = rules.trim.output.info_out,
    output:
        fasta_f = "{src}/primer-f.fasta",
        fasta_r = "{src}/primer-r.fasta",
        plot_f = "{src}/primer-f.svg",
        plot_r = "{src}/primer-r.svg",
        geo_html = "{src}/metadata/country-heat.html",
        geo_svg = "{src}/metadata/country-heat.svg",
    params:
        primer_f_len = len(config['primer_f']['seq']),
        primer_r_len = len(config['primer_r']['seq']),
    run:
        primer.logo_plot(input.info, input.meta, output.fasta_f, params.primer_f_len, output.plot_f, 'forward')
        primer.logo_plot(input.info, input.meta, output.fasta_r, params.primer_r_len, output.plot_r, 'reverse')
        meta = pd.read_csv(input.meta, sep='\t')
        metadata.plot_country_heatmap(meta, output.geo_html, output.geo_svg)
