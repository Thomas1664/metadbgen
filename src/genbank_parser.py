import re
from dataclasses import dataclass
import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


@dataclass
class ParsedRecord:
    accession: str
    creation_date: str
    title: str
    authors: str
    sequence: str
    sequence_length: int
    organism_name: str
    tax_id: int
    cox2_start: str
    cox2_end: str
    cox2_strand: int


# from https://techwasti.com/python-tutorial-writing-python-data-class-to-csv
def dataclass_to_csv_row(csv_writer, row):
    csv_writer.writerow(row.__dict__.values())
def dataclass_to_csv_header(csv_writer, header):
    csv_writer.writerow(header.__annotations__.keys())


cox2_aliases = {"cox2", "coxII", "COX2", "COXII", "Cox2", "CoxII", 'cytochrome c oxidase subunit II', 'cytochrome c oxidase subunit 2', 'cytochrome oxidase subunit II', 'cytochrome oxidase subunit 2'}

def check_title(title: str):
    # TODO: check for false positives; only if cox2 pos not found in features
    unwanted_tags = ['RNA', 'ITS', 'cox1', 'coxI']
    if np.array([title.find(tag) != -1 for tag in unwanted_tags]).any():
        return np.array([title.find(cox2) != -1 for cox2 in cox2_aliases]).any()
    return True

def parse_taxid(features: list):
    for feature in features:
        if feature.type == "source" and "db_xref" in feature.qualifiers:
            for db_xref in feature.qualifiers["db_xref"]:
                if db_xref.startswith("taxon:"):
                    return db_xref.split(":")[1]

def parse_cox2_pos(features: list[SeqFeature]):
    # List of possible names for cox2 (case insensitive)
    start, end, strand = "N/A", "N/A", "N/A"
    for feature in features:
        if feature.type == "gene" and "gene" in feature.qualifiers:
            gene_name = feature.qualifiers["gene"][0]
            if any(re.fullmatch(alias, gene_name, re.IGNORECASE) for alias in cox2_aliases):
                start = feature.location.start
                end = feature.location.end
                strand = feature.location.strand
                break
    else:
        for feature in features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers["product"][0]
                if any(re.fullmatch(alias, gene_name, re.IGNORECASE) for alias in cox2_aliases):
                    start = feature.location.start
                    end = feature.location.end
                    strand = feature.location.strand
                    break
    return (start, end, strand)

def parse_record(record: SeqRecord):
    cox2_start, cox2_end, cox2_strand = parse_cox2_pos(record.features)
    return ParsedRecord(
        accession = record.id,  # Accession number
        creation_date = record.annotations.get("date", "N/A"),  # Date of creation
        authors = record.annotations.get("references", [])[0].authors if record.annotations.get("references") else "N/A",
        sequence = str(record.seq),  # Nucleotide sequence
        sequence_length = len(record.seq),  # Sequence length
        title = record.description,  # Record title
        organism_name = record.annotations.get("organism", "N/A"),  # Organism name
        tax_id = parse_taxid(record.features),
        cox2_start=cox2_start,
        cox2_end=cox2_end,
        cox2_strand=cox2_strand
    )

def parse_genbank_file(gb_path: str, csv_writer):
    seqs = SeqIO.parse(gb_path, 'genbank')
    for seq in seqs:
        record = parse_record(seq)
        dataclass_to_csv_row(csv_writer, record)
