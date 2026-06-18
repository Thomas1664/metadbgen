from dataclasses import dataclass
import asyncio
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import json
try:
    from .request import HttpClient
except ImportError:
    from request import HttpClient


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


async def _fetch_fasta_by_accession(accession: str) -> str | None:
    client = HttpClient()
    try:
        responses = await client.efetch_by_ids(
            db='nuccore',
            ids=[accession],
            callback=lambda text: text,
            batch_size=1,
            ret_type='fasta'
        )
    finally:
        await client._client.aclose()

    if not responses:
        return None
    fasta_response = responses[0].strip()
    if not fasta_response:
        return None

    parsed = list(SeqIO.parse(StringIO(fasta_response), 'fasta'))
    if not parsed:
        return None
    return str(parsed[0].seq)


def fetch_fasta_by_accession(accession: str) -> str | None:
    return asyncio.run(_fetch_fasta_by_accession(accession))

def parse_record(record: SeqRecord):
    accession: str = record.id

    # Some records only contain a contig annotation and no sequence body.
    if (record.seq is None or len(record.seq) == 0 or not record.seq.defined) and record.annotations.get('contig'):
        fetched_sequence = fetch_fasta_by_accession(accession)
        if fetched_sequence:
            record.seq = Seq(fetched_sequence)

    if record.seq is None or not record.seq.defined:
        return None
    return ParsedRecord(
        accession = accession,
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

def parse_genbank_file(gb_path: str, csv_writer, fasta_file):
    seqs = SeqIO.parse(gb_path, 'genbank')
    for seq in seqs:
        record = parse_record(seq)
        if record is None:
            continue
        SeqIO.write(SeqRecord(id=record.accession, seq=Seq(record.sequence), description=''), fasta_file, 'fasta')
        record.sequence = None
        dataclass_to_csv_row(csv_writer, record)
