import sys
from src import database

# Reusable HMM-scoring component, keyed by the {src} wildcard.
# Callers must set two config keys before including/using this module:
#   config['alignment_ref_fasta']   - trusted reference sequences the HMM is built from
#   config['alignment_query_fasta'] - candidate sequences to score against that HMM
# Both are path patterns containing a literal "{src}" placeholder.


rule msa:
    input: config['alignment_ref_fasta']
    output: "{src}/hmm/db-msa.fasta"
    threads: workflow.cores
    shell: "muscle -super5 {input} -nt -threads {threads} -output {output}"


rule hmm_profile:
    input: rules.msa.output[0]
    output: "{src}/hmm/db-profile.hmm"
    threads: workflow.cores
    shell: "hmmbuild --dna --cpu {threads} {output} {input}"

rule score_db_seq:
    input:
        hmm = rules.hmm_profile.output[0],
        db_seq = config['alignment_ref_fasta'],
    output:
        tblout = "{src}/hmm/db-scores.tsv",
        hmm_out = "{src}/hmm/db-scores.hmm",
        threshold_file = "{src}/hmm/score-threshold.txt",
    params:
        margin = config['hmm_score_margin'],
    threads: workflow.cores
    run:
        shell("nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.db_seq}")
        df = database.parse_hmm_tblout(output.tblout)
        min_score = df['score'].min()
        threshold = min_score - params.margin
        print(f"Min RefSeq score: {min_score:.1f}  →  threshold: {threshold:.1f}", file=sys.stderr)
        with open(output.threshold_file, 'w') as f:
            f.write(f"{threshold}\n")

rule score_seqs:
    input:
        hmm = rules.hmm_profile.output[0],
        query_seqs = config['alignment_query_fasta'],
    output:
        tblout = "{src}/hmm/sequence-scores.tsv",
        hmm_out = "{src}/hmm/sequence-scores.hmm",
    threads: workflow.cores
    shell:
        "nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.query_seqs}"


rule filter_by_score:
    """Keep query sequences that meet the HMM score threshold."""
    input:
        tblout = rules.score_seqs.output.tblout,
        threshold_file = rules.score_db_seq.output.threshold_file,
        fasta = config['alignment_query_fasta'],
    output:
        accessions = "{src}/hmm/accessions-hmm-filter.txt",
        fasta = "{src}/hmm/sequences-hmm-filter.fasta",
    run:
        df = database.parse_hmm_tblout(input.tblout)
        with open(input.threshold_file) as f:
            threshold = float(f.read().strip())
        good = df[df['score'].ge(threshold)]
        print(f"Sequences passing HMM score filter: {len(good)}", file=sys.stderr)
        good['target_name'].to_csv(output.accessions, header=False, index=False)
        shell("seqkit grep -f {output.accessions} -o {output.fasta} {input.fasta}")
