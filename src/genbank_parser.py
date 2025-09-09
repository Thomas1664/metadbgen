from dataclasses import dataclass
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json


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
    metadata: str
    pubmed: int


# from https://techwasti.com/python-tutorial-writing-python-data-class-to-csv
def dataclass_to_csv_row(csv_writer, row):
    csv_writer.writerow(row.__dict__.values())
def dataclass_to_csv_header(csv_writer, header):
    csv_writer.writerow(header.__annotations__.keys())

def parse_taxid(features: list):
    for feature in features:
        if feature.type == "source" and "db_xref" in feature.qualifiers:
            for db_xref in feature.qualifiers["db_xref"]:
                if db_xref.startswith("taxon:"):
                    return db_xref.split(":")[1]

def parse_source(features: list):
    for feature in features:
        if feature.type == "source":
            return json.dumps(feature.qualifiers)

def parse_references(record: SeqRecord):
    pubmed_id = None
    for ref in record.annotations.get("references", []):
        if ref.pubmed_id != '':
            pubmed_id = ref.pubmed_id 
    return pubmed_id

def parse_record(record: SeqRecord):
    # some accessions have no sequence
    if not record.seq.defined:
        return None
    return ParsedRecord(
        accession = record.id,
        creation_date = record.annotations.get("date", "N/A"),
        authors = record.annotations.get("references", [])[0].authors if record.annotations.get("references") else "N/A",
        sequence = str(record.seq),
        sequence_length = len(record.seq),
        title = record.description,
        organism_name = record.annotations.get("organism", "N/A"),
        tax_id = parse_taxid(record.features),
        metadata=parse_source(record.features),
        pubmed=parse_references(record)
    )

def parse_genbank_file(gb_path: str, csv_writer):
    seqs = SeqIO.parse(gb_path, 'genbank')
    for seq in seqs:
        record = parse_record(seq)
        if record is None:
            continue
        dataclass_to_csv_row(csv_writer, record)
