import os

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dotenv import load_dotenv

from src.genbank_parser import fetch_missing_sequences, load_sequence_lookup, parse_record


def test_parse_record_fills_sequence_from_lookup_for_empty_contig_record():
    record = SeqRecord(Seq(None, length=100), id="NW_012157837")
    record.annotations["contig"] = "NW_012157837"
    record.annotations["date"] = "01-JAN-2000"
    record.annotations["organism"] = "test organism"

    parsed = parse_record(record, sequence_lookup={"NW_012157837": "ACGTACGT"})

    assert parsed is not None
    assert parsed.accession == "NW_012157837"
    assert parsed.sequence == "ACGTACGT"
    assert parsed.sequence_length == 8


def test_parse_record_returns_none_for_empty_contig_record_without_lookup():
    record = SeqRecord(Seq(None, length=100), id="NW_012157837")
    record.annotations["contig"] = "NW_012157837"
    record.annotations["date"] = "01-JAN-2000"
    record.annotations["organism"] = "test organism"

    assert parse_record(record) is None
    assert parse_record(record, sequence_lookup={}) is None


@pytest.mark.asyncio
async def test_fetch_missing_sequences_fetches_fasta_for_empty_contig_record(tmp_path):
    load_dotenv()
    if not os.getenv("NCBI_API_KEY"):
        pytest.skip("NCBI_API_KEY is required for integration test")

    gb_path = tmp_path / "contig_only.gb"
    gb_path.write_text(
        "LOCUS       NW_012157837               1 bp    DNA     linear   CON 01-JAN-2000\n"
        "DEFINITION  test organism.\n"
        "ACCESSION   NW_012157837\n"
        "VERSION     NW_012157837.1\n"
        "SOURCE      test organism\n"
        "  ORGANISM  test organism\n"
        "CONTIG      join(NW_012157837.1:1..1)\n"
        "//\n"
    )
    out_fasta = tmp_path / "missing-sequences.fasta"

    await fetch_missing_sequences(str(gb_path), str(out_fasta))
    sequence_lookup = load_sequence_lookup(str(out_fasta))

    assert "NW_012157837.1" in sequence_lookup
    assert len(sequence_lookup["NW_012157837.1"]) > 0


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
