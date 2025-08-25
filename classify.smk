import glob
import os

gb_prefix = config['gb_prefix']
db_prefix = f"{config['base_dir']}/{config['db_name']}"
db_name = config['db_name']
sample_name = config['sample_name']

rule all:
    input:
        f"{db_prefix}/seq-trimmed-removed.qza",
        f"{db_prefix}/qiime/{db_name}_sequences.qza",
        f"{db_prefix}/qiime/{db_name}_tax.qza",
        f"{db_prefix}/qiime/{db_name}_trained_classifier.qza",
        f"{db_prefix}/qiime/{db_name}_taxonomy.qza",
        f"{db_prefix}/qiime/{db_name}_taxonomy.tsv"

rule qiime_import:
    input:
        f"{db_prefix}/seq-trimmed-removed.fasta"
    output:
        f"{db_prefix}/seq-trimmed-removed.qza"
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input} --output-path {output}
        """

rule import_db:
    input:
        fasta = f"{db_prefix}/{config['db_name']}-qiime.fa",
        #fasta = f"{base_dir}/taxonomy/{db_name}.fa",
        tax = f"{db_prefix}/{config['db_name']}-qiime.tsv",
    output:
        fasta = f"{db_prefix}/qiime/{db_name}_sequences.qza",
        tax = f"{db_prefix}/qiime/{db_name}_tax.qza",
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input.fasta} --output-path {output.fasta}
        qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path {input.tax} --output-path {output.tax}
        """

rule train:
    input:
        fasta = rules.import_db.output.fasta,
        tax = rules.import_db.output.tax
    output:
        f"{db_prefix}/qiime/{db_name}_trained_classifier.qza"
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads {input.fasta} --i-reference-taxonomy {input.tax} --o-classifier {output}
        """

rule assign_taxonomy:
    input:
        classifier = rules.train.output[0],
        asv = rules.qiime_import.output[0],
    output:
        f"{db_prefix}/qiime/{db_name}_taxonomy.qza"
    shell:
        "qiime feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.asv} --o-classification {output}"

rule export_taxonomy:
    input:
        taxonomy = rules.assign_taxonomy.output[0],
    output:
        taxonomy_tsv = f"{db_prefix}/qiime/{db_name}_taxonomy.tsv",
    shell:
        """
        qiime tools export --input-path {input.taxonomy} --output-path {db_prefix}/{sample_name}_taxonomy
        mv {db_prefix}/{sample_name}_taxonomy/taxonomy.tsv {output.taxonomy_tsv}
        """
