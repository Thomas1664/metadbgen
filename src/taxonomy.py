import xml.etree.ElementTree as ET
from lxml import etree
from dataclasses import dataclass
import csv
import pandas as pd
from src import request


async def download_taxonomy(query: str, output_path: str):
    client = request.HttpClient()
    response = await client.esearch(query, db='taxonomy')
    await client.efetch_all(response, out_file=output_path, db='taxonomy', rettype='xml')

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

def parse_taxon(taxon: ET.Element):
    """Extract TaxID, rank, scientific name, and human-readable rank from a <Taxon> element."""
    taxid = int(taxon.find("TaxId").text)
    name = str(taxon.find("ScientificName").text)
    rank = str(taxon.find("Rank").text if taxon.find("Rank") is not None else "no rank")
    return Taxon(taxid, name, rank)

def parse_lineage(lineage_elem: ET.Element):
    wanted_ranks = {"domain", "clade", "phylum", "order", "family", "class", "genus"}
    # clade can appear multiple times
    lineage = {rank: [] if rank == "clade" else None for rank in wanted_ranks}
    # Extract Taxon entries
    for taxon in lineage_elem.findall("Taxon"):
        rank = taxon.find("Rank").text
        name = taxon.find("ScientificName").text
        if rank in wanted_ranks:
            # clade can be multiple, so keep them in a list
            if rank == "clade":
                lineage[rank].append(name)
            else:
                lineage[rank] = name
    return lineage

def parse_xml(file_path: str, csv_writer):
    """Parse XML file and extract taxonomical information."""
    parser = etree.XMLParser(resolve_entities=False)
    tree = etree.parse(file_path, parser)
    root = tree.getroot()

    for taxon in root.findall("Taxon"):
        # Extract species-level info
        tax_data = parse_taxon(taxon)
        lineage_elem = taxon.find('LineageEx')
        lineage = parse_lineage(lineage_elem)
        csv_writer.writerow([tax_data.taxid, tax_data.name, lineage['domain'], lineage['clade'][-1], lineage['phylum'], lineage['class'], lineage['order'], lineage['family'], lineage['genus']])

def parse_taxonomy(input_xml: str, output_csv: str):
    with open(output_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['tax_id', 'species', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
        parse_xml(input_xml, csv_writer)

def sanitize_org_name(org_name: str):
    org_name_arr = org_name.split(' ')
    if org_name_arr[1] == 'sp.': return False
    return True


def parse_col_tree(col_path: str, out_path: str):
    def build_lineage_dict(taxon_id):
        lineage = {}
        while taxon_id in id_to_parent and taxon_id in id_to_name:
            rank = id_to_rank[taxon_id].lower()
            lineage[rank] = id_to_name[taxon_id]
            taxon_id = id_to_parent[taxon_id]
            if taxon_id not in id_to_name and taxon_id not in id_to_parent:
                break
        return lineage

    df = pd.read_csv(col_path, sep='\t', header=0, names=[
        "id", "parent", "status", "rank", "scientificName", "authorship"
    ])

    id_to_name = dict(zip(df["id"], df["scientificName"]))
    id_to_parent = dict(zip(df["id"], df["parent"]))
    id_to_rank = dict(zip(df["id"], df["rank"]))

    lineages = []
    for _, row in df[df["rank"] == "species"].iterrows():
        lineage_dict = build_lineage_dict(row["id"])
        lineages.append(lineage_dict)

    result = pd.DataFrame(lineages)
    # Force consistent rank order
    rank_order = ["phylum", "class", "order", "family", "genus", "species"]
    result = result.reindex(columns=rank_order)
    result.to_csv(out_path, index=False)

def fill_taxonomic_blanks(row, backup, ranks):
    # Determine which ranks are missing
    missing_ranks = [r for r in ranks if pd.isna(row[r])]
    
    # For non-missing ranks, filter df2 rows where they match
    filter_mask = pd.Series(True, index=backup.index)
    for rank in ranks:
        if rank not in missing_ranks:
            filter_mask &= (backup[rank] == row[rank])
    
    # If a match exists, fill missing ranks
    if filter_mask.any():
        match = backup.loc[filter_mask].iloc[0]
        for rank in missing_ranks:
            if rank in backup.columns:
                row[rank] = match[rank]
    
    # For any remaining missing ranks, append 'incertae sedis' to the next higher rank
    ranks = ['class', 'order', 'family', 'genus']
    for i, rank in enumerate(ranks):
        if pd.isna(row[rank]):
            if rank == 'genus':
                # Extract genus from species name (first word)
                row[rank] = row['species'].split(' ')[0]
            else:
                # Look for the closest known higher rank
                higher_rank = None
                for j in range(i-1, -1, -1):
                    if pd.notna(row[ranks[j]]):
                        higher_rank = row[ranks[j]]
                        break
                if higher_rank is not None:
                    row[rank] = f"{higher_rank} incertae sedis"
                else:
                    # fallback if nothing above is known
                    row[rank] = "incertae sedis"
    return row


def _check_taxonomy_complete(meta: pd.DataFrame, tax: pd.DataFrame):
    missing_tax = meta[~meta['tax_id_sanitized'].isin(tax['tax_id'])]
    if not missing_tax.empty:
        print('Warning: No taxonomy information found for the following accessions:')
        print(missing_tax)
    return missing_tax.empty

def filter_taxonomy(meta_path: str, input_tax: str, backup_tax_path: str, output_tax: str):
    meta = pd.read_csv(meta_path, sep='\t')
    tax = pd.read_csv(input_tax)
    _check_taxonomy_complete(meta, tax)
    tax_filtered = tax[tax['tax_id'].isin(meta['tax_id_sanitized'])]
    backup_tax = pd.read_csv(backup_tax_path)
    tax_ranks = ['order', 'family', 'genus']
    tax_filled = pd.DataFrame(tax_filtered.apply(lambda row: fill_taxonomic_blanks(row, backup_tax, tax_ranks), axis=1))
    tax_filled = tax_filled[tax_filled['species'].apply(sanitize_org_name)]
    tax_filled.to_csv(output_tax, sep='\t', index=False)

def merge_taxonomy(meta_in: str, tax: str):
    meta = pd.read_csv(meta_in, sep='\t')
    tax_df = pd.read_csv(tax, sep='\t')
    meta_and_tax = meta.merge(tax_df, right_on='tax_id', left_on='tax_id_sanitized', how='inner')
    return meta_and_tax
