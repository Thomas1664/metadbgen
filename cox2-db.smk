import glob
import os
import re
import asyncio
from src import metadata, request, taxonomy, primer, database, mislabel_detection, genbank_parser, hmm
from src.ncbi_taxonomy import download_and_build_db
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import csv
import io
import shutil
from concurrent.futures import ThreadPoolExecutor
from itertools import batched


db_prefix = f"{config['base_dir']}/{config['db_name']}"
find_missing_species = config.get('find_missing_species', False)
find_no_primers = config.get('find_no_primers', False)
autofetch_additional_keywords = config.get('autofetch_additional_keywords', False)

sources = [db_prefix]
if config['find_outgroups']:
    sources.append(f'{db_prefix}/outgroups')
if find_missing_species or find_no_primers:
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

        # Entropy plot
        expand("{src}/output/msa-entropy.tsv", src=sources),
        expand("{src}/output/msa-entropy.html", src=sources),
        expand("{src}/output/msa-entropy.svg", src=sources),


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


rule tax_db_download:
    output:
        tax_db = f"{db_prefix}/taxonomy/tax_db.sqlite",
    run:
        download_and_build_db(db_path=output.tax_db, taxdump_dir=f"{db_prefix}/ncbi_tax_dump")


rule primer_blast_outgroups:
    retries: 5
    input:
        tax_db = rules.tax_db_download.output.tax_db
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
            future = pool.submit(lambda: asyncio.run(primer.find_off_targets(params.primer_f, params.primer_r, params.organism, output.out, params.use_cache, params.cache_path, input.tax_db)))
            future.result()


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
        tax_df = pd.read_csv(output.tax)
        taxonomy.update_species_name(tax_df).to_csv(output.tax, index=False)

rule fetch_missing_sequences:
    retries: 5
    input:
        gb_file = "{src}/input-sequences.gb",
    output:
        fasta = "{src}/missing-sequences.fasta",
    run:
        with ThreadPoolExecutor() as pool:
            future = pool.submit(lambda: asyncio.run(genbank_parser.fetch_missing_sequences(input.gb_file, output.fasta)))
            future.result()


rule parse_metadata:
    wildcard_constraints:
        src = f"(?!{db_prefix}/round2(/new-keywords)?$).*",
    input:
        gb_file = "{src}/input-sequences.gb",
        missing_seqs = "{src}/missing-sequences.fasta",
        tax_db = rules.tax_db_download.output.tax_db,
        tax = "{src}/taxonomy/parsed.csv",
        script = "parse-gb.py"
    output:
        meta_csv = "{src}/seq-meta.csv",
        fasta_tmp = "{src}/sequences.fasta.tmp",
        fasta = "{src}/sequences.fasta",
    run:
        shell("python {input.script} {input.gb_file} {input.missing_seqs} {input.tax_db} {input.tax} {output.meta_csv} {output.fasta_tmp}")
        shell("seqkit rmdup -n -o {output.fasta} {output.fasta_tmp}")


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
                resolved_ids = list({t.tax_id for t in resolved.values()})
                species_tax = db.find_parent_at_rank(tax_ids=resolved_ids, rank='species')
                tax_ids = list({t.tax_id for t in species_tax.values() if t is not None})
                tax_data = db.get_taxon_data_bulk(tax_ids)
            with open(output.tax, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                taxonomy.to_csv(tax_data, csv_writer, header=True)
            tax_df = pd.read_csv(output.tax)
            taxonomy.update_species_name(tax_df).to_csv(output.tax, index=False)


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
    wildcard_constraints:
        src = f"(?!{re.escape(db_prefix)}/round2$).+"
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
        info_untrimmed = "{src}/primer-trim/cutadapt-info-untrimmed.tsv",
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
        shell("python filter_cutadapt_info.py {output.info_tmp} {output.info} {output.info_untrimmed}")
        shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


rule clean_db:
    input:
        meta = "{src}/primer-trim/metadata.tsv",
        fasta = "{src}/primer-trim/sequences.fasta",
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
        organism = lambda wildcards: "" if wildcards.src == f"{db_prefix}/outgroups" else f'--organism "{config['organism']}"',
        no_remove_on = "--no-remove-on-species" if not config.get('remove_on_species', True) else "",
        no_prefer_refseq = "--no-prefer-refseq" if not config.get('prefer_refseq', True) else "",
    run:
        shell("python {input.script} --in-metadata {input.meta} --in-fasta {input.fasta} --out-metadata {output.meta} --out-fasta {output.accessions} --out-metadata-removed {output.meta_removed} --out-fasta-removed {output.accessions_removed}  --tax-db {input.tax_db} {params.organism} {params.no_remove_on} {params.no_prefer_refseq}")
        shell("seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}")
        shell("seqkit grep -f {output.accessions_removed} -o {output.fasta_removed} {input.fasta}")


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
    run:
        shell("vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}")
        shell("python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta} {params.remove_on}")


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
        cutadapt_info = "{src}/primer-trim/cutadapt-info-parsed.tsv",
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
    run:
        shell("vsearch --cluster_fast {input.fasta} -id 1.0 --strand both --msaout {output.msa} --clusterout_id --uc {output.clust}")
        shell("python {input.script} {input.meta} {input.fasta} {output.clust} {output.meta} {output.fasta} {params.remove_on}")

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
        shell('iqtree3 -s {output.msa} -m "GTR+G" --fast -n {params.max_iterations} -nt {threads} --redo')
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
        mislabel = "{src}/mislabel-detection.csv",
        accessions = f"{{src}}/output/{config['db_name']}-accessions.txt",
        removed = "{src}/mislabeled-sequences.fasta",
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
        tax = rules.mislabel_detection.output.meta,
    output:
        msa = f"{{src}}/{config['db_name']}-msa.fasta",
        db = f"{{src}}/output/{config['db_name']}-sintax.fasta",
    threads: workflow.cores / 2
    run:
        shell("muscle -super5 {input.fasta} -nt -threads {threads} -output {output.msa}")
        database.export_sintax(output.msa, input.tax, output.db)

rule msa_entropy:
    input:
        msa = rules.export_sintax.output.db,
        script = "entropy.py"
    output:
        tsv  = "{src}/output/msa-entropy.tsv",
        html = "{src}/output/msa-entropy.html",
        svg  = "{src}/output/msa-entropy.svg",
    threads: workflow.cores / 2
    run:
        shell("goalign compute entropy -g -i {input.msa} -t {threads} > {output.tsv}")
        shell("python {input.script} {output.tsv} {output.html} {output.svg}")


rule plot:
    input:
        meta = rules.mislabel_detection.output.meta,
        info = "{src}/primer-trim/cutadapt-info-parsed.tsv",
    output:
        fasta_f = "{src}/primer-f.fasta",
        fasta_r = "{src}/primer-r.fasta",
        plot_f = "{src}/primer-f.svg",
        plot_r = "{src}/primer-r.svg",
        geo_html = "{src}/metadata/country-heat.html",
        geo_svg = "{src}/metadata/country-heat.svg",
    params:
        primer_f = config['primer_f'],
        primer_r = config['primer_r'],
    run:
        primer.logo_plot(input.info, input.meta, output.fasta_f, len(params.primer_f), output.plot_f, 'forward', params.primer_f)
        primer.logo_plot(input.info, input.meta, output.fasta_r, len(params.primer_r), output.plot_r, 'reverse', params.primer_r)
        meta = pd.read_csv(input.meta, sep='\t')
        metadata.plot_country_heatmap(meta, output.geo_html, output.geo_svg)


if find_missing_species or find_no_primers:
    # alignment.smk builds an HMM once from round-1's already-finalized
    # reference set (round2's own reference doesn't exist yet — it's built
    # downstream of trim_missing). hmm-query-filter.smk then scores the newly
    # parsed round2 candidates against that HMM and keeps only the sequences
    # (and metadata rows) that plausibly belong, so trim_missing only runs
    # cutadapt on those. Both modules are scoped to {db_prefix}/round2 via
    # prefix:, reproducing the same {db_prefix}/round2/hmm/... layout as
    # before.
    module alignment:
        snakefile: "alignment.smk"
        config: {**config, 'ref_fasta': f"{db_prefix}/output/{config['db_name']}-qiime-clean.fa"}
        prefix: f"{db_prefix}/round2"

    use rule * from alignment as alignment_*


if find_missing_species:
    # combine_sequences produces round2/primer-trim/ directly, bypassing cutadapt.
    # Without this declaration Snakemake may prefer the wildcard 'trim' rule, which
    # would run cutadapt on new sequences only — silently dropping round-1 sequences
    # (already trimmed = no primers = discarded as untrimmed by cutadapt).
    ruleorder: combine_sequences > trim
    # parse_missing handles round2 and round2/new-keywords (via its own {src}
    # wildcard_constraints); the more general parse_metadata rule must not run
    # for those paths, hence the explicit ordering as a backstop.
    ruleorder: parse_missing > parse_metadata

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
            keywords = config["keywords"]
        run:
            missing_df = pd.read_csv(input.missing)
            tmp_dir = f"{db_prefix}/round2/tmp"
            
            def build_query(organisms: list[str], exclude_keywords: list[str]):
                exclude = " NOT ".join([f'"{keyword}"' for keyword in exclude_keywords])
                q_organisms = " OR ".join([f'"{organism}"[organism]' for organism in organisms])
                return f'({q_organisms}) NOT {exclude} AND ("0001/01/01"[PDAT] : "2026/04/15"[PDAT])'

            async def download_all():
                client = request.HttpClient()
                os.makedirs(tmp_dir, exist_ok=True)
                for batch, species in enumerate(batched(missing_df["species"].tolist(), 10)):
                    query = build_query(species, params.keywords)
                    res = await client.esearch(query)
                    if res.result_count == 0:
                        continue
                    await client.efetch_all(res, f'{tmp_dir}/{batch}.gb')
                    for split_file in glob.glob(f'{tmp_dir}/{batch}_*.gb'):
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
        wildcard_constraints:
            src = f"{db_prefix}/round2(/new-keywords)?",
        input:
            r1_meta = f"{db_prefix}/seq-meta.csv",
            gb_file = "{src}/input-sequences.gb",
            missing_seqs = "{src}/missing-sequences.fasta",
            tax_db = rules.tax_db_download.output.tax_db,
            tax = f"{db_prefix}/taxonomy/parsed.csv",
            script = "parse-gb.py",
        output:
            meta_csv = "{src}/seq-meta.csv",
            r1_accs = "{src}/round1_accs.txt",
            fasta_tmp = "{src}/sequences.fasta.tmp",
            fasta = "{src}/sequences.fasta",
        run:
            shell("python {input.script} {input.gb_file} {input.missing_seqs} {input.tax_db} {input.tax} {output.meta_csv} {output.fasta_tmp}")
            pd.read_csv(input.r1_meta)['accession'].to_csv(output.r1_accs, header=False, index=False)
            shell("seqkit rmdup -n {output.fasta_tmp} | seqkit grep -v -f {output.r1_accs} -o {output.fasta}")


    if autofetch_additional_keywords:
        module hmm_filter_new_keywords:
            snakefile: "hmm-query-filter.smk"
            config:
                {**config,
                'hmm': rules.alignment_hmm_profile.output[0],
                'threshold_file': rules.alignment_score_db_seq.output.threshold_file,
                'query_fasta': f"{db_prefix}/round2/sequences.fasta",
                'query_meta': f"{db_prefix}/round2/seq-meta.csv"}
            prefix: f"{db_prefix}/round2/new-keywords/r2-hits"
        use rule * from hmm_filter_new_keywords as hmm_filter_new_keywords*
        
        rule download_new_keywords:
            retries: 5
            input:
                gb_file = f"{db_prefix}/round2/input-sequences.gb",
                hmm_res = f"{db_prefix}/round2/new-keywords/r2-hits/hmm/sequence-scores.tsv",
                hmm_meta = f"{db_prefix}/round2/new-keywords/r2-hits/hmm/metadata-hmm-filter.tsv",
            output:
                gb_file = f"{db_prefix}/round2/new-keywords/input-sequences.gb",
            params:
                orig_keywords = config["keywords"],
                organism = config["organism"]
            run:
                tmp_dir = f"{db_prefix}/round2/new-keywords"
                os.makedirs(tmp_dir, exist_ok=True)
                meta_df = hmm.read_metadata(input.hmm_meta)
                hmm_hits = hmm.parse_hmm_tblout(input.hmm_res)
                hmm_hits = hmm_hits[hmm_hits['target_name'].isin(meta_df['accession'])]
                
                gb_records = SeqIO.parse(input.gb_file, format='genbank')
                new_keywords = hmm.find_additional_keywords(hmm_hits, gb_records, params.orig_keywords)

                if not new_keywords:
                    print('No additional keywords found, skipping download.')
                    with open(output.gb_file, 'w') as f:
                        f.write('')
                else:
                    new_query = hmm.build_query_additional_keywords(new_keywords, params.orig_keywords, params.organism)
                    with ThreadPoolExecutor() as pool:
                        future = pool.submit(lambda: asyncio.run(hmm.download_additional_keywords(new_query, tmp_dir)))
                        future.result()


        rule merge_new_keywords:
            input:
                fasta = f"{db_prefix}/round2/sequences.fasta",
                meta = f"{db_prefix}/round2/seq-meta.csv",
                new_kw_fasta = f"{db_prefix}/round2/new-keywords/sequences.fasta",
                new_kw_meta = f"{db_prefix}/round2/new-keywords/seq-meta.csv",
            output:
                fasta = f"{db_prefix}/round2/merged/sequences.fasta",
                meta = f"{db_prefix}/round2/merged/seq-meta.csv",
            run:
                os.makedirs(f"{db_prefix}/round2/merged", exist_ok=True)
                meta = pd.concat([pd.read_csv(input.meta), pd.read_csv(input.new_kw_meta)], ignore_index=True)
                meta.drop_duplicates(subset='accession').to_csv(output.meta, index=False)
                shell("cat {input.fasta} {input.new_kw_fasta} | seqkit rmdup -n -o {output.fasta}")


if find_no_primers:
    rule filter_no_primers:
        """
        Recovers sequences cutadapt discarded as untrimmed (no primer match) by
        salvaging them from the pre-trim source and length-filtering them.
        Handles both round 1 (round="") and, when find_missing_species is also
        enabled, round 2 (round="round2/", sourced from trim_missing's input/output).
        """
        wildcard_constraints:
            round = "|round2/"
        input:
            cutadapt_untrimmed = lambda wildcards: (
                f"{db_prefix}/primer-trim/cutadapt-info-untrimmed.tsv" if wildcards.round == ""
                else rules.trim_missing.output.info_untrimmed
            ),
            source_meta = lambda wildcards: (
                f"{db_prefix}/seq-meta.csv" if wildcards.round == "" else rules.trim_missing.input.meta
            ),
            source_fasta = lambda wildcards: (
                f"{db_prefix}/sequences.fasta" if wildcards.round == "" else rules.trim_missing.input.fasta
            ),
            r1_final_fasta = f"{db_prefix}/output/{config['db_name']}-qiime-clean.fa"
        output:
            accessions = f'{db_prefix}/{{round}}untrimmed/sequences.txt',
            fasta = f'{db_prefix}/{{round}}untrimmed/sequences.fasta',
            meta = f'{db_prefix}/{{round}}untrimmed/metadata.tsv',
            surviving_acc = f'{db_prefix}/{{round}}untrimmed/sequences-surviving.txt',
        params:
            err_margin = config['error_margin'],
            min_seq_len = config['min_seq_len'],
        run:
            os.makedirs(f"{db_prefix}/{wildcards.round}untrimmed", exist_ok=True)
            untrimmed = pd.read_csv(input.cutadapt_untrimmed, sep='\t', header=None)
            untrimmed[0].to_csv(output.accessions, header=False, index=False)
            # get max length of sequences in final DB
            max_seq_len = pd.read_csv(io.StringIO(shell("seqkit stats -T {input.r1_final_fasta}", read=True)), sep='\t')['max_len'][0]
            shell(f"seqkit grep -f {{output.accessions}} {{input.source_fasta}} | seqkit seq -m {params.min_seq_len} -M {max_seq_len + params.err_margin} -o {{output.fasta}}")
            shell("seqkit seq --name --only-id {output.fasta} -o {output.surviving_acc}")

            try:
                accs = pd.read_csv(output.surviving_acc, header=None)[0]
                meta_df = pd.read_csv(input.source_meta)
                meta_df = meta_df[meta_df['accession'].isin(accs)]
                meta_df.to_csv(output.meta, sep='\t', index=False)
            except pd.errors.EmptyDataError:
                shell("touch {output.meta}")


    module hmm_filter_no_primers:
        snakefile: "hmm-query-filter.smk"
        config:
            {**config,
             'hmm': rules.alignment_hmm_profile.output[0],
             'threshold_file': rules.alignment_score_db_seq.output.threshold_file,
             'query_fasta': rules.filter_no_primers.output.fasta.format(round=""),
             'query_meta': rules.filter_no_primers.output.meta.format(round="")}
        prefix: f"{db_prefix}/untrimmed"

    use rule * from hmm_filter_no_primers as hmm_filter_no_primers*


if find_missing_species:
    rule trim_missing:
        input:
            fasta = (
                f"{db_prefix}/round2/merged/sequences.fasta"
                if autofetch_additional_keywords
                else f"{db_prefix}/round2/sequences.fasta"
            ),
            meta = (
                f"{db_prefix}/round2/merged/seq-meta.csv"
                if autofetch_additional_keywords
                else f"{db_prefix}/round2/seq-meta.csv"
            ),
            script = "trim.py",
        output:
            fasta = f"{db_prefix}/round2/primer-trim-new/sequences.fasta",
            meta = f"{db_prefix}/round2/primer-trim-new/metadata.tsv",
            info_tmp = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-tmp.tsv",
            info = f"{db_prefix}/round2/primer-trim-new/cutadapt-info.tsv",
            info_out = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-parsed.tsv",
            info_untrimmed = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-untrimmed.tsv",
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
            shell("python filter_cutadapt_info.py {output.info_tmp} {output.info} {output.info_untrimmed}")
            shell("python {input.script} {input.meta} {output.fasta} {output.meta} {output.info} {output.info_out}")


    if find_no_primers:
        module hmm_filter_no_primers_r2:
            snakefile: "hmm-query-filter.smk"
            config:
                {**config,
                 'hmm': rules.alignment_hmm_profile.output[0],
                 'threshold_file': rules.alignment_score_db_seq.output.threshold_file,
                 'query_fasta': rules.filter_no_primers.output.fasta.format(round="round2/"),
                 'query_meta': rules.filter_no_primers.output.meta.format(round="round2/")}
            prefix: f"{db_prefix}/round2/untrimmed"

        use rule * from hmm_filter_no_primers_r2 as hmm_filter_no_primers_r2*


if find_missing_species or find_no_primers:
    rule combine_sequences:
        """
        Builds round2/primer-trim/ as the starting point for the second pass.
        Concrete output paths let Snakemake prefer this rule over the 'trim' wildcard rule
        when building the round2 source tree.
        """
        input:
            info_orig = f"{db_prefix}/primer-trim/cutadapt-info-parsed.tsv",
            # Round-1 final database: defines which primer-trim sequences form the base.
            meta_final = f"{db_prefix}/output/{config['db_name']}-meta-clean.tsv",
            fasta_final = f"{db_prefix}/output/{config['db_name']}-qiime-clean.fa",
            # Newly trimmed sequences for previously missing species (trim_missing).
            fasta_new = f"{db_prefix}/round2/primer-trim-new/sequences.fasta" if find_missing_species else [],
            meta_new = f"{db_prefix}/round2/primer-trim-new/metadata.tsv" if find_missing_species else [],
            info_new = f"{db_prefix}/round2/primer-trim-new/cutadapt-info-parsed.tsv" if find_missing_species else [],
            # Missing-species candidates recovered without primer matches (hmm_filter_no_primers_r2).
            fasta_hmm_r2 = f"{db_prefix}/round2/untrimmed/hmm/sequences-hmm-filter.fasta" if (find_missing_species and find_no_primers) else [],
            meta_hmm_r2 = f"{db_prefix}/round2/untrimmed/hmm/metadata-hmm-filter.tsv" if (find_missing_species and find_no_primers) else [],
            # Round-1 candidates recovered without primer matches (hmm_filter_no_primers).
            fasta_hmm_r1 = f"{db_prefix}/untrimmed/hmm/sequences-hmm-filter.fasta" if find_no_primers else [],
            meta_hmm_r1 = f"{db_prefix}/untrimmed/hmm/metadata-hmm-filter.tsv" if find_no_primers else [],
        output:
            fasta = f"{db_prefix}/round2/primer-trim/sequences.fasta",
            meta = f"{db_prefix}/round2/primer-trim/metadata.tsv",
            info_out = f"{db_prefix}/round2/primer-trim/cutadapt-info-parsed.tsv",
        run:
            os.makedirs(f"{db_prefix}/round2/primer-trim", exist_ok=True)

            # Round-1 final database forms the base; no overlap with any of the new sources below.
            meta_final_r1 = pd.read_csv(input.meta_final, sep='\t')
            meta_parts = [meta_final_r1]

            orig_info = pd.read_csv(input.info_orig, sep='\t')
            info_parts = [orig_info[orig_info['accession'].isin(meta_final_r1['accession'])]]

            shutil.copy2(input.fasta_final, output.fasta)

            if find_missing_species:
                new_meta = pd.read_csv(input.meta_new, sep='\t')
                meta_parts.append(new_meta)

                new_acc_file = f"{db_prefix}/round2/primer-trim-new/new_acc.txt"
                new_meta['accession'].to_csv(new_acc_file, header=False, index=False)
                shell("seqkit grep -f " + new_acc_file + " {input.fasta_new} >> {output.fasta}")

                new_info = pd.read_csv(input.info_new, sep='\t')
                info_parts.append(new_info[new_info['accession'].isin(new_meta['accession'])])

            if find_missing_species and find_no_primers:
                # No cutadapt trim coordinates for these (recovered without a primer match);
                # they're simply not used as MSA-trim references downstream (see primer.trim_msa).
                meta_parts.append(pd.read_csv(input.meta_hmm_r2, sep='\t'))
                shell("cat {input.fasta_hmm_r2} >> {output.fasta}")

            if find_no_primers:
                meta_parts.append(pd.read_csv(input.meta_hmm_r1, sep='\t'))
                shell("cat {input.fasta_hmm_r1} >> {output.fasta}")

            pd.concat(meta_parts, ignore_index=True).to_csv(output.meta, sep='\t', index=False)
            pd.concat(info_parts, ignore_index=True).to_csv(output.info_out, sep='\t', index=False)
