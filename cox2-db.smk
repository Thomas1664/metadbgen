import glob
import os
import asyncio
from src import metadata, request, taxonomy, primer, database, mislabel_detection, genbank_parser
from src.ncbi_taxonomy import download_and_build_db
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import csv
import io
import shutil
from concurrent.futures import ThreadPoolExecutor


db_prefix = f"{config['base_dir']}/{config['db_name']}"

sources = [db_prefix]
if config['find_outgroups']:
    sources.append(f'{db_prefix}/outgroups')
if config.get('find_missing_species', False):
    sources.append(f'{db_prefix}/round2')

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
        primer_f = config['primer_f'],
        primer_r = config['primer_r'],
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
    wildcard_constraints:
        src = f"(?!{db_prefix}/round2$).*",
    input:
        gb_file = f"{{src}}/input-sequences.gb",
        tax_db = rules.tax_db_download.output.tax_db,
        tax = "{src}/taxonomy/parsed.csv",
        script = "parse-gb.py"
    output:
        meta_csv = "{src}/seq-meta.csv",
        fasta_tmp = "{src}/sequences.fasta.tmp",
        fasta = "{src}/sequences.fasta",
    shell:
        """
        python {input.script} {input.gb_file} {input.tax_db} {input.tax} {output.meta_csv} {output.fasta_tmp}
        seqkit rmdup -n -o {output.fasta} {output.fasta_tmp}
        """


if config['find_outgroups']:
    rule outgroup_tax:
        input:
            gb_file = f"{db_prefix}/outgroups/input-sequences.gb",
            tax_db = rules.tax_db_download.output.tax_db,
        output:
            tax = f"{db_prefix}/outgroups/taxonomy/parsed.csv",
        run:
            raw_tax_ids = list({
                int(genbank_parser.parse_taxid(r.features))
                for r in SeqIO.parse(input.gb_file, 'genbank')
                if genbank_parser.parse_taxid(r.features) is not None
            })
            with taxonomy.NCBITaxonomyDB(input.tax_db) as db:
                resolved = db.resolve_tax_ids(raw_tax_ids)
                tax_ids = list({t.tax_id for t in resolved.values()})
                tax_data = db.get_taxon_data_bulk(tax_ids)
            with open(output.tax, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                taxonomy.to_csv(tax_data, csv_writer, header=True)


if "outer_primer_f" in config or "outer_primer_r" in config:
    rule pcr_outer:
        input:
            rules.parse_metadata.output.fasta
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
            else rules.parse_metadata.output.fasta
        ),
        meta = rules.parse_metadata.output.meta_csv,
        script = "trim.py",
    output:
        fasta = "{src}/primer-trim/sequences.fasta",
        meta = "{src}/primer-trim/metadata.tsv",
        info_tmp = "{src}/primer-trim/cutadapt-info-tmp.tsv",
        info = "{src}/primer-trim/cutadapt-info.tsv",
        info_out = "{src}/primer-trim/cutadapt-info-parsed.tsv",
    params:
        primer_f_seq = config['primer_f'],
        primer_r_seq = str(Seq(config['primer_r']).reverse_complement()),
        min_primer_overlap = config['min_primer_overlap'],
        min_seq_len = config['min_seq_len'],
        max_seq_len = config['max_seq_len'],
    threads: workflow.cores / 2
    run:
        os.makedirs(f'{wildcards.src}/primer-trim', exist_ok=True)
        max_seq_len = pd.read_csv(io.StringIO(shell("seqkit stats -T {input.fasta}", read=True)), sep='\t')['max_len'][0]
        buffer_size = max(max_seq_len * 2, 1048576)
        shell("cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --no-indel --overlap {params.min_primer_overlap} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info_tmp} -o {output.fasta} --cores {threads} --buffer-size " + str(buffer_size) + " {input.fasta}")
        shell("python filter_cutadapt_info.py {output.info_tmp} {output.info}")
        shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


rule clean_db:
    input:
        meta = rules.trim.output.meta,
        fasta = rules.trim.output.fasta,
        tax_db = rules.tax_db_download.output.tax_db,
        script = "db-clean.py"
    output:
        meta = "{src}/seq-meta-filtered.tsv",
        meta_removed = "{src}/seq-meta-filtered-removed.tsv",
        fasta = "{src}/seq-filtered.fasta",
        fasta_removed = "{src}/seq-filtered-removed.fasta",
        accessions = "{src}/seq-filtered.txt",
        accessions_removed = "{src}/seq-filtered-removed.txt",
    params:
        organism = lambda wildcards: "" if wildcards.src == f"{db_prefix}/outgroups" else f"--organism '{config['organism']}'",
        no_remove_on = "--no-remove-on-species" if not config.get('remove_on_species', True) else "",
        no_prefer_refseq = "--no-prefer-refseq" if not config.get('prefer_refseq', True) else "",
    shell:
        """
        python {input.script} \\
            --in-metadata {input.meta} \\
            --in-fasta {input.fasta} \\
            --out-metadata {output.meta} \\
            --out-fasta {output.accessions} \\
            --out-metadata-removed {output.meta_removed} \\
            --out-fasta-removed {output.accessions_removed} \\
            --tax-db {input.tax_db} {params.organism} {params.no_remove_on} {params.no_prefer_refseq}
        seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}
        seqkit grep -f {output.accessions_removed} -o {output.fasta_removed} {input.fasta}
        """


rule dereplicate:
    input:
        meta = rules.clean_db.output.meta,
        fasta = rules.clean_db.output.fasta,
        script = 'dereplicate.py',
    output:
        msa = "{src}/dereplicate/clust-msa.fa",
        clust = "{src}/dereplicate/clust-results.tsv",
        meta = "{src}/dereplicate/metadata.tsv",
        fasta = "{src}/dereplicate/sequences.fasta",
    params:
        remove_on = str(config.get('remove_on_species', True)),
    shell:
        """
        vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}
        python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta} {params.remove_on}
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
    output:
        msa = "{src}/alignment/msa.fasta",
        fasta = "{src}/alignment/sequences-trim.fasta",
        meta = "{src}/alignment/metadata-trim.tsv",
    params:
        min_length = config['min_seq_len'],
    threads: workflow.cores / 2
    run:
        os.makedirs(f'{wildcards.src}/alignment', exist_ok=True)
        shell("muscle -super5 {input.fasta} -nt -threads {threads} -output {output.msa}")
        primer.trim_msa(input.cutadapt_info, input.meta, output.msa, output.fasta, output.meta, int(params.min_length))


rule dereplicate_alignment:
    input:
        meta = rules.trim_alignment.output.meta,
        fasta = rules.trim_alignment.output.fasta,
        script = 'dereplicate.py',
    output:
        msa = "{src}/alignment/clust-msa.fa",
        clust = "{src}/alignment/clust-results.tsv",
        meta = "{src}/alignment/metadata.tsv",
        fasta = "{src}/alignment/sequences.fasta",
    params:
        remove_on = str(config.get('remove_on_species', True)),
    shell:
        """
        vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}
        python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta} {params.remove_on}
        """

rule gappiness:
    input:
        fasta = rules.dereplicate_alignment.output.fasta,
        meta = rules.dereplicate_alignment.output.meta,
        script = "gappiness.py",
    output:
        msa = "{src}/gappiness/msa.fasta",
        gap_scores = "{src}/gappiness/msa-gap-scores.csv",
        fasta = "{src}/gappiness/sequences.fasta",
        meta = "{src}/gappiness/metadata.tsv",
    params:
        gap_threshold = config['gappiness_threshold'],
    threads: workflow.cores / 2
    run:
        os.makedirs(f'{wildcards.src}/gappiness', exist_ok=True)
        shell("muscle -super5 {input.fasta} -nt -threads {threads} -output {output.msa}")
        shell("python {input.script} {output.msa} {input.meta} {input.fasta} {output.meta} {output.fasta} {output.gap_scores} {params.gap_threshold}")


rule export_db:
    input:
        meta = rules.gappiness.output.meta,
        fasta = rules.gappiness.output.fasta,
        tax_db = rules.tax_db_download.output.tax_db,
        script = 'taxonomy-export.py'
    output:
        meta = f"{{src}}/output/{config['db_name']}-meta.tsv",
        tax = f"{{src}}/output/{config['db_name']}-qiime.tsv",
        fasta = f"{{src}}/output/{config['db_name']}-qiime.fa",
    params:
        # It's not guaranteed this file will be always present
        # because it might be possible that none of the organisms in the database has a host.
        # Therefore, this file is not an output file.
        host_tax = "{src}/metadata/host-tax.csv",
    run:
        os.makedirs(f'{wildcards.src}/metadata', exist_ok=True)
        shell("python {input.script} {input.meta} {input.tax_db} {output.meta} {output.tax} {params.host_tax}")
        shutil.copy2(input.fasta, output.fasta)


rule tree:
    input: rules.gappiness.output.fasta,
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
        shell("gotree matrix -i {output.tree} -o {output.tree_dist} -t {threads}")


rule mislabel_detection:
    input:
        tree_dist = rules.tree.output.tree_dist,
        tax = rules.export_db.output.tax,
        meta = rules.export_db.output.meta,
        fasta = rules.export_db.output.fasta,
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
        primer_f_len = len(config['primer_f']),
        primer_r_len = len(config['primer_r']),
    run:
        primer.logo_plot(input.info, input.meta, output.fasta_f, params.primer_f_len, output.plot_f, 'forward')
        primer.logo_plot(input.info, input.meta, output.fasta_r, params.primer_r_len, output.plot_r, 'reverse')
        meta = pd.read_csv(input.meta, sep='\t')
        metadata.plot_country_heatmap(meta, output.geo_html, output.geo_svg)


if config.get('find_missing_species', False):

    rule round2_taxonomy:
        input:
            f"{db_prefix}/taxonomy/parsed.csv",
        output:
            f"{db_prefix}/round2/taxonomy/parsed.csv",
        shell:
            "cp {input} {output}"

    rule find_missing_species:
        input:
            tax = f"{db_prefix}/taxonomy/parsed.csv",
            meta = f"{db_prefix}/output/{config['db_name']}-meta-clean.tsv",
        output:
            missing = f"{db_prefix}/round2/missing-species.csv",
        run:
            tax_df = pd.read_csv(input.tax)
            tax_df['species_sanitized'] = tax_df['species'].apply(database.sanitize_org_name)
            binom_df = tax_df[
                (tax_df['species'] == tax_df['species_sanitized']) &
                (~tax_df['species'].str.endswith('sp.')) &
                (~tax_df['species'].str.startswith('uncultured'))
            ]
            db_tax_df = pd.read_csv(input.meta, sep='\t')
            missing_df = binom_df[~binom_df['species'].isin(db_tax_df['species'])]
            missing_df[['species', 'tax_id']].drop_duplicates().to_csv(output.missing, index=False)


    rule download_missing:
        retries: 5
        input:
            missing = f"{db_prefix}/round2/missing-species.csv",
        output:
            gb_file = f"{db_prefix}/round2/input-sequences.gb",
        params:
            keywords = config['keywords'],
        run:
            missing_df = pd.read_csv(input.missing)
            tmp_dir = f"{db_prefix}/round2/tmp"

            async def download_all():
                client = request.HttpClient()
                os.makedirs(tmp_dir, exist_ok=True)
                for _, row in missing_df.iterrows():
                    query = database.build_search_query(params.keywords, row['species'])
                    res = await client.esearch(query)
                    if res.result_count == 0:
                        continue
                    await client.efetch_all(res, f'{tmp_dir}/{row["tax_id"]}.gb')
                    for split_file in glob.glob(f'{tmp_dir}/{row["tax_id"]}_*.gb'):
                        os.remove(split_file)
                gb_files = sorted(glob.glob(f'{tmp_dir}/*.gb'))
                with open(output.gb_file, 'w') as out_f:
                    for gb_file in gb_files:
                        with open(gb_file) as f:
                            out_f.write(f.read())

            with ThreadPoolExecutor() as pool:
                future = pool.submit(lambda: asyncio.run(download_all()))
                future.result()


    rule parse_missing:
        input:
            gb_file = f"{db_prefix}/round2/input-sequences.gb",
            tax_db = rules.tax_db_download.output.tax_db,
            tax = f"{db_prefix}/taxonomy/parsed.csv",
            script = "parse-gb.py",
        output:
            meta_csv = f"{db_prefix}/round2/seq-meta.csv",
            fasta_tmp = f"{db_prefix}/round2/sequences.fasta.tmp",
            fasta = f"{db_prefix}/round2/sequences.fasta",
        shell:
            """
            python {input.script} {input.gb_file} {input.tax_db} {input.tax} {output.meta_csv} {output.fasta_tmp}
            seqkit rmdup -n -o {output.fasta} {output.fasta_tmp}
            """


    rule trim_missing:
        input:
            fasta = f"{db_prefix}/round2/sequences.fasta",
            meta = f"{db_prefix}/round2/seq-meta.csv",
            script = "trim.py",
        output:
            fasta = f"{db_prefix}/round2/primer-trim-new/sequences.fasta",
            meta = f"{db_prefix}/round2/primer-trim-new/metadata.tsv",
            info_tmp = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-tmp.tsv",
            info = f"{db_prefix}/round2/primer-trim-new/cutadapt-info.tsv",
            info_out = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-parsed.tsv",
        params:
            primer_f_seq = config['primer_f'],
            primer_r_seq = str(Seq(config['primer_r']).reverse_complement()),
            min_primer_overlap = config['min_primer_overlap'],
            min_seq_len = config['min_seq_len'],
            max_seq_len = config['max_seq_len'],
        threads: workflow.cores / 2
        run:
            os.makedirs(f"{db_prefix}/round2/primer-trim-new", exist_ok=True)
            stats = pd.read_csv(io.StringIO(shell("seqkit stats -T {input.fasta}", read=True)), sep='\t')
            try:
                seq_max_len = int(stats['max_len'][0])
            except (ValueError, TypeError):
                seq_max_len = 0
            buffer_size = max(seq_max_len * 2, 1048576)
            shell("cutadapt -a {params.primer_f_seq}...{params.primer_r_seq} --revcomp --no-indel --overlap {params.min_primer_overlap} -e 0.2 --discard-untrimmed --minimum-length {params.min_seq_len} --maximum-length {params.max_seq_len} --info-file {output.info_tmp} -o {output.fasta} --cores {threads} --buffer-size " + str(buffer_size) + " {input.fasta}")
            shell("python filter_cutadapt_info.py {output.info_tmp} {output.info}")
            shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


    rule combine_sequences:
        """
        Merges round-1 primer-trim outputs with newly trimmed missing-species sequences.
        Downstream wildcard rules use {src}/primer-trim/... paths; since this rule's
        outputs are concrete (no wildcards), Snakemake prefers them over the trim rule
        when building the round2 source tree.
        """
        input:
            fasta_orig = f"{db_prefix}/primer-trim/sequences.fasta",
            meta_orig = f"{db_prefix}/primer-trim/metadata.tsv",
            info_orig = f"{db_prefix}/primer-trim/cutadapt-info-parsed.tsv",
            fasta_new = f"{db_prefix}/round2/primer-trim-new/sequences.fasta",
            meta_new = f"{db_prefix}/round2/primer-trim-new/metadata.tsv",
            info_new = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-parsed.tsv",
        output:
            fasta = f"{db_prefix}/round2/primer-trim/sequences.fasta",
            meta = f"{db_prefix}/round2/primer-trim/metadata.tsv",
            info_out = f"{db_prefix}/round2/primer-trim/cutadapt-info-parsed.tsv",
            info_tmp = f"{db_prefix}/round2/primer-trim/cutadapt-info-tmp.tsv",
            info = f"{db_prefix}/round2/primer-trim/cutadapt-info.tsv",
        run:
            os.makedirs(f"{db_prefix}/round2/primer-trim", exist_ok=True)

            orig_meta = pd.read_csv(input.meta_orig, sep='\t')
            existing_accessions = set(orig_meta['accession'])

            new_meta = pd.read_csv(input.meta_new, sep='\t') if os.path.getsize(input.meta_new) > 0 else pd.DataFrame(columns=orig_meta.columns)
            new_meta_filtered = new_meta[~new_meta['accession'].isin(existing_accessions)]
            new_accessions = set(new_meta_filtered['accession'])

            pd.concat([orig_meta, new_meta_filtered], ignore_index=True).to_csv(output.meta, sep='\t', index=False)

            shutil.copy2(input.fasta_orig, output.fasta)
            if new_accessions:
                acc_file = f"{db_prefix}/round2/primer-trim/.new_acc.txt"
                with open(acc_file, 'w') as f:
                    f.write('\n'.join(new_accessions))
                shell("seqkit grep -f " + acc_file + " {input.fasta_new} >> {output.fasta}")
                os.remove(acc_file)

            orig_info = pd.read_csv(input.info_orig, sep='\t')
            new_info = pd.read_csv(input.info_new, sep='\t') if os.path.getsize(input.info_new) > 0 else pd.DataFrame(columns=orig_info.columns)
            new_info_filtered = new_info[new_info['accession'].isin(new_accessions)]
            pd.concat([orig_info, new_info_filtered], ignore_index=True).to_csv(output.info_out, sep='\t', index=False)

            shutil.copy2(output.info_out, output.info_tmp)
            shutil.copy2(output.info_out, output.info)
