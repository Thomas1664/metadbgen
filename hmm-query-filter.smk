import os
import sys
import tempfile
from src.hmm import passing_hits, split_by_strand, read_metadata, filter_metadata_by_accession, read_threshold_file

# Reusable HMM scoring/filtering component.
#
# Scores one candidate set against an HMM profile + threshold that were
# already built once (see alignment.smk), and keeps the sequences (and their
# metadata rows) that meet the threshold. Call this once per candidate set:
#
#   module round2_hmm_filter:
#       snakefile: "hmm-query-filter.smk"
#       config:
#           {**config,
#            'hmm': "path/to/db-profile.hmm",
#            'threshold_file': "path/to/score-threshold.txt",
#            'query_fasta': "path/to/candidates.fasta",
#            'query_meta': "path/to/candidates-meta.csv"}
#       prefix: "path/to/output/dir"
#
#   use rule * from round2_hmm_filter as round2_hmm_filter_*
#
# `prefix:` scopes this module's own rule outputs to "path/to/output/dir" so
# multiple candidate sets scored against the same (or different) HMMs don't
# collide. hmm, threshold_file, query_fasta, and query_meta are wrapped in
# local() because they are ordinary caller-supplied paths, not paths meant to
# be resolved relative to this module's own prefix.
#
# Config keys:
#   hmm            - path to an already-built HMM profile
#   threshold_file - path to the matching score threshold (one float)
#   query_fasta    - candidate sequences to score against that HMM
#   query_meta     - metadata table (needs an 'accession' column) for the
#                    candidate sequences; filtered alongside the fasta


rule score_seqs:
    input:
        hmm = local(config['hmm']),
        query_seqs = local(config['query_fasta']),
    output:
        tblout = "hmm/sequence-scores.tsv",
        hmm_out = "hmm/sequence-scores.hmm",
    threads: workflow.cores
    shell:
        """
        if [ -s {input.query_seqs} ]; then
            nhmmer --cpu {threads} --tblout {output.tblout} -o {output.hmm_out} {input.hmm} {input.query_seqs}
        else
            touch {output.tblout} {output.hmm_out}
        fi
        """


rule filter_by_score:
    """Keep query sequences (and matching metadata rows) that meet the HMM score threshold.

    nhmmer scores both strands of every query independently. A hit on the
    minus strand means the query is reverse-complemented relative to the
    reference, so those survivors are flipped back with seqkit before being
    written out, to give downstream consumers a consistently-oriented fasta.
    """
    input:
        tblout = rules.score_seqs.output.tblout,
        threshold_file = local(config['threshold_file']),
        fasta = local(config['query_fasta']),
        meta = local(config['query_meta']),
    output:
        accessions = "hmm/accessions-hmm-filter.txt",
        fasta = "hmm/sequences-hmm-filter.fasta",
        meta = "hmm/metadata-hmm-filter.tsv",
    run:
        threshold = read_threshold_file(input.threshold_file)
        good = passing_hits(input.tblout, threshold)
        print(f"Sequences passing HMM score filter: {len(good)}", file=sys.stderr)
        good['target_name'].to_csv(output.accessions, header=False, index=False)

        plus, minus = split_by_strand(good)
        with tempfile.TemporaryDirectory(dir=os.path.dirname(output.fasta) or '.') as tmp:
            plus_acc = f"{tmp}/plus-accessions.txt"
            minus_acc = f"{tmp}/minus-accessions.txt"
            plus_fasta = f"{tmp}/plus.fasta"
            minus_fasta = f"{tmp}/minus.fasta"
            minus_rc_fasta = f"{tmp}/minus-rc.fasta"

            plus.to_csv(plus_acc, header=False, index=False)
            minus.to_csv(minus_acc, header=False, index=False)

            shell("seqkit grep -f {plus_acc} -o {plus_fasta} {input.fasta}")
            shell("seqkit grep -f {minus_acc} -o {minus_fasta} {input.fasta}")
            shell("seqkit seq -r -p -o {minus_rc_fasta} {minus_fasta}")
            shell("cat {plus_fasta} {minus_rc_fasta} > {output.fasta}")

        meta_df = read_metadata(input.meta)
        meta_df = filter_metadata_by_accession(meta_df, good['target_name'])
        meta_df.to_csv(output.meta, sep='\t', index=False)
