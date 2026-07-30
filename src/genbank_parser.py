from dataclasses import dataclass
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


def _needs_fetched_sequence(record: SeqRecord) -> bool:
    """True for records that only carry a CONTIG annotation and no inline sequence."""
    is_empty = record.seq is None or len(record.seq) == 0 or not record.seq.defined
    return is_empty and bool(record.annotations.get('contig'))


def find_missing_sequence_ids(gb_path: str) -> list[str]:
    """Accessions of records in gb_path that need a separate FASTA fetch to fill their sequence."""
    return [record.id for record in SeqIO.parse(gb_path, 'genbank') if _needs_fetched_sequence(record)]


async def _fetch_missing_sequences(accessions: list[str], out_fasta_path: str, batch_size: int) -> None:
    client = HttpClient()
    try:
        with open(out_fasta_path, 'w') as fasta_file:
            await client.efetch_by_ids(
                db='nuccore',
                ids=accessions,
                callback=lambda x: fasta_file.write(x),
                batch_size=batch_size,
                ret_type='fasta',
            )
    finally:
        await client._client.aclose()


async def fetch_missing_sequences(gb_path: str, out_fasta_path: str, batch_size: int = 100) -> None:
    """Scan gb_path for records that need a fetched sequence and batch-fetch them from NCBI."""
    accessions = find_missing_sequence_ids(gb_path)
    if not accessions:
        open(out_fasta_path, 'w').close()
        return
    await _fetch_missing_sequences(accessions, out_fasta_path, batch_size)


def load_sequence_lookup(fasta_path: str) -> dict[str, str]:
    """Load a FASTA file written by fetch_missing_sequences into an {accession: sequence} dict."""
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, 'fasta')}


def parse_record(record: SeqRecord, sequence_lookup: dict[str, str] | None = None):
    accession: str = record.id

    # Some records only contain a contig annotation and no sequence body.
    if _needs_fetched_sequence(record) and sequence_lookup is not None:
        fetched_sequence = sequence_lookup.get(accession)
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

def parse_genbank_file(gb_path: str, csv_writer, fasta_file, sequence_lookup: dict[str, str] | None = None):
    seqs = SeqIO.parse(gb_path, 'genbank')
    for seq in seqs:
        record = parse_record(seq, sequence_lookup)
        if record is None:
            continue
        SeqIO.write(SeqRecord(id=record.accession, seq=Seq(record.sequence), description=''), fasta_file, 'fasta')
        record.sequence = None
        dataclass_to_csv_row(csv_writer, record)
