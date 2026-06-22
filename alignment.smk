import pandas as pd
import io
import os
import sys


db_prefix = f"{config['base_dir']}/{config['db_name']}"

sources = [db_prefix]


rule all:
    input:
        expand("{src}/untrimmed/sequences-hmm-filter.fasta", src=sources)

def _parse_hmm_tblout(path):
    cols = [
        "target_name", "target_accession", "query_name", "query_accession",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "sq_len", "strand", "e_value", "score", "bias", "description",
    ]
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, engine='python', names=cols)
    return df.sort_values('score', ascending=False).drop_duplicates(subset=['target_name'], keep='first')


rule filter_untrimmed:
    input:
        untrimmed = '{src}/primer-trim/cutadapt-info-untrimmed.tsv',
        meta = "{src}/seq-meta.csv",
        fasta = "{src}/sequences.fasta",
        fasta_res_db = f"{{src}}/output/cox2_th-qiime-clean.fa",
    output:
        accessions = '{src}/untrimmed/sequences.txt',
        fasta = '{src}/untrimmed/sequences.fasta',
        meta = '{src}/untrimmed/metadata.tsv',
        surviving_acc = '{src}/untrimmed/sequences-surviving.txt',
    params:
        err_margin = config['error_margin'],
        min_seq_len = config['min_seq_len'],
    run:
        os.makedirs(f"{wildcards.src}/untrimmed", exist_ok=True)
        untrimmed = pd.read_csv(input.untrimmed, sep='\t', header=None)
        untrimmed[0].to_csv(output.accessions, header=False, index=False)
        # get max length of sequences in final DB
        max_seq_len = pd.read_csv(io.StringIO(shell("seqkit stats -T {input.fasta_res_db}", read=True)), sep='\t')['max_len'][0]
        shell(f"seqkit grep -f {{output.accessions}} {{input.fasta}} | seqkit seq -m {params.min_seq_len} -M {max_seq_len + params.err_margin} -o {{output.fasta}}")
        shell("seqkit seq --name --only-id {output.fasta} -o {output.surviving_acc}")

        accs = pd.read_csv(output.surviving_acc, header=None)[0]
        meta_df = pd.read_csv(input.meta)
        meta_df = meta_df[meta_df['accession'].isin(accs)]
        meta_df.to_csv(output.meta, sep='\t', index=False)


rule msa:
    input: f"{{src}}/output/cox2_th-qiime-clean.fa"
    output: "{src}/untrimmed/db-msa.fasta"
    threads: workflow.cores
    shell: "muscle -super5 {input} -nt -threads {threads} -output {output}"


rule hmm_profile:
    input: rules.msa.output[0]
    output: "{src}/untrimmed/db-profile.hmm"
    threads: workflow.cores
    shell: "hmmbuild --dna --cpu {threads} {output} {input}"

rule score_db_seq:
    input:
        hmm = rules.hmm_profile.output[0],
        db_seq = f"{{src}}/output/cox2_th-qiime-clean.fa",
    output:
        tblout = "{src}/untrimmed/db-scores.tsv",
        hmm_out = "{src}/untrimmed/db-scores.hmm",
        threshold_file = "{src}/untrimmed/score-threshold.txt",
    params:
        margin = config['hmm_score_margin'],
    threads: workflow.cores
    run:
        shell("nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.db_seq}")
        df = _parse_hmm_tblout(output.tblout)
        min_score = df['score'].min()
        threshold = min_score - params.margin
        print(f"Min RefSeq score: {min_score:.1f}  →  threshold: {threshold:.1f}", file=sys.stderr)
        with open(output.threshold_file, 'w') as f:
            f.write(f"{threshold}\n")

rule score_seqs:
    input:
        hmm = rules.hmm_profile.output[0],
        query_seqs = rules.filter_untrimmed.output.fasta,
    output:
        tblout = "{src}/untrimmed/sequence-scores.tsv",
        hmm_out = "{src}/untrimmed/sequence-scores.hmm",
    threads: workflow.cores
    shell:
        "nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.query_seqs}"


rule filter_by_score:
    """Keep non-RefSeq query sequences that meet the HMM score threshold."""
    input:
        tblout = rules.score_seqs.output.tblout,
        threshold_file = rules.score_db_seq.output.threshold_file,
        fasta = rules.filter_untrimmed.output.fasta,
    output:
        accessions = "{src}/untrimmed/accessions-hmm-filter.txt",
        fasta = "{src}/untrimmed/sequences-hmm-filter.fasta",
    run:
        df = _parse_hmm_tblout(input.tblout)
        with open(input.threshold_file) as f:
            threshold = float(f.read().strip())
        good = df[df['score'].ge(threshold)]
        print(f"Sequences passing HMM score filter: {len(good)}", file=sys.stderr)
        good['target_name'].to_csv(output.accessions, header=False, index=False)
        shell("seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}")
