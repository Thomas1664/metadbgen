import sys
from src.hmm import compute_score_threshold, write_threshold_file

# One-time HMM-building component: aligns a trusted reference set, builds an
# HMM profile from it, then scores that same reference set against its own
# profile to derive a score threshold. Build this once per reference set via:
#
#   module alignment:
#       snakefile: "alignment.smk"
#       config: {**config, 'ref_fasta': "path/to/reference.fasta"}
#       prefix: "path/to/output/dir"
#
#   use rule * from alignment as alignment_*
#
# The resulting profile ("path/to/output/dir/hmm/db-profile.hmm") and
# threshold ("path/to/output/dir/hmm/score-threshold.txt") can then be reused
# by hmm-query-filter.smk to score any number of different candidate sets without
# rebuilding the HMM. `prefix:` scopes this module's own rule outputs so
# multiple reference sets don't collide; ref_fasta is wrapped in local()
# because it's an ordinary caller-supplied path, not one meant to be resolved
# relative to this module's own prefix.
#
# Config keys:
#   ref_fasta - trusted reference sequences the HMM is built from


rule msa:
    input: local(config['ref_fasta'])
    output: "hmm/db-msa.fasta"
    threads: workflow.cores
    shell: "muscle -super5 {input} -nt -threads {threads} -output {output}"


rule hmm_profile:
    input: rules.msa.output[0]
    output: "hmm/db-profile.hmm"
    threads: workflow.cores
    shell: "hmmbuild --dna --cpu {threads} {output} {input}"

rule score_db_seq:
    input:
        hmm = rules.hmm_profile.output[0],
        db_seq = local(config['ref_fasta']),
    output:
        tblout = "hmm/db-scores.tsv",
        hmm_out = "hmm/db-scores.hmm",
        threshold_file = "hmm/score-threshold.txt",
    params:
        margin = config['hmm_score_margin'],
    threads: workflow.cores
    run:
        shell("nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.db_seq}")
        min_score, threshold = compute_score_threshold(output.tblout, params.margin)
        print(f"Min RefSeq score: {min_score:.1f}  →  threshold: {threshold:.1f}", file=sys.stderr)
        write_threshold_file(output.threshold_file, threshold)
