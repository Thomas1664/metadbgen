import xml.etree.ElementTree as ET
from lxml import etree
from dataclasses import dataclass
import csv
import pandas as pd

# from https://techwasti.com/python-tutorial-writing-python-data-class-to-csv
def dataclass_to_csv_row(csv_writer, row):
    csv_writer.writerow(row.__dict__.values())
def dataclass_to_csv_header(csv_writer, header):
    csv_writer.writerow(header.__annotations__.keys())

@dataclass
class Taxon:
    taxid: int
    name: str
    rank: str

@dataclass
class LineageInfo:
    species: Taxon
    genus: Taxon
    family: Taxon
    order: Taxon
    phylum: Taxon
    kingdom: Taxon
    domain: Taxon

def parse_taxon(taxon: ET.Element):
    """Extract TaxID, rank, scientific name, and human-readable rank from a <Taxon> element."""
    taxid = int(taxon.find("TaxId").text)
    name = str(taxon.find("ScientificName").text)
    rank = str(taxon.find("Rank").text if taxon.find("Rank") is not None else "no rank")
    return Taxon(taxid, name, rank)

def parse_xml(file_path: str, csv_writer):
    """Parse XML file and extract taxonomical information."""
    parser = etree.XMLParser(resolve_entities=False)
    tree = etree.parse(file_path, parser)
    root = tree.getroot()

    for taxon in root.findall("Taxon"):
        # Extract species-level info
        tax_data = parse_taxon(taxon)
        lineage = taxon.find('Lineage').text
        lin_arr = lineage.split('; ')
        eukayota_idx = lin_arr.index("Eukaryota")
        lin_arr = lin_arr[eukayota_idx:eukayota_idx + 7]
        csv_writer.writerow([tax_data.taxid, tax_data.name, *lin_arr])

def parse_taxonomy(input_xml: str, output_csv: str):
    with open(output_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['tax_id', 'species', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
        parse_xml(input_xml, csv_writer)

def sanitize_org_name(org_name: str):
    org_name_arr = org_name.split(' ')
    if org_name_arr[1] == 'sp.': return False
    return True

def filter_taxonomy(input_tax: str, output_tax: str):
    tax = pd.read_csv(input_tax)
    tax = tax.dropna(subset='family')
    tax = tax[~tax['family'].str.startswith('unclassified')]
    tax = tax[~tax['family'].str.startswith('environmental')]
    tax = tax.dropna(subset='genus')
    tax = tax[~tax['genus'].str.startswith('unclassified')]
    tax = tax[~tax['genus'].str.startswith('environmental')]
    tax = tax[tax['species'].apply(sanitize_org_name)]
    tax.to_csv(output_tax, sep='\t', index=False)
