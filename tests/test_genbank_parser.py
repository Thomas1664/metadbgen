import os

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dotenv import load_dotenv

from src.genbank_parser import parse_record


def test_parse_record_fetches_fasta_for_empty_contig_record():
    load_dotenv()
    if not os.getenv("NCBI_API_KEY"):
        pytest.skip("NCBI_API_KEY is required for integration test")

    record = SeqRecord(Seq(None, length=0), id="NW_012157837")
    record.annotations["contig"] = "NW_012157837"
    record.annotations["date"] = "01-JAN-2000"
    record.annotations["organism"] = "test organism"

    parsed = parse_record(record)

    assert parsed is not None
    assert parsed.accession == "NW_012157837"
    assert parsed.sequence_length > 0
    assert len(parsed.sequence) == parsed.sequence_length


def test_parse_record_returns_none_for_empty_record_without_contig():
    record = SeqRecord(None, id="NO_CONTIG_001")
    record.annotations["date"] = "01-JAN-2000"
    record.annotations["organism"] = "test organism"

    parsed = parse_record(record)

    assert parsed is None


def test_parse_record_parses_normal_record():
    record = SeqRecord(Seq("ACGTACGT"), id="NORMAL_001", description="normal sequence")
    record.annotations["date"] = "01-JAN-2001"
    record.annotations["organism"] = "Testus organismus"
    record.features = [
        SeqFeature(
            FeatureLocation(0, 8),
            type="source",
            qualifiers={"db_xref": ["taxon:12345"], "country": ["Germany"]},
        )
    ]

    parsed = parse_record(record)

    assert parsed is not None
    assert parsed.accession == "NORMAL_001"
    assert parsed.sequence == "ACGTACGT"
    assert parsed.sequence_length == 8
    assert parsed.organism_name == "Testus organismus"
    assert parsed.tax_id == "12345"
