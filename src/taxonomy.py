from collections.abc import Iterator
from bs4 import BeautifulSoup
from pydantic import BaseModel, Field, ValidationError, model_validator
from src import request, database
import pandas as pd
import csv
import re
import time
import xmltodict


class Taxon(BaseModel):
    tax_id: int = Field(alias='TaxId')
    name: str = Field(alias='ScientificName')
    rank: str = Field(alias='Rank')

class TaxonData(BaseModel):
    tax_id: int = Field(alias='TaxId')
    name: str = Field(alias='ScientificName')
    rank: str = Field(alias='Rank')
    parent_tax_id: int = Field(alias='ParentTaxId')
    lineage: list[Taxon] = Field(alias='taxa')

    @model_validator(mode="before")
    @classmethod
    def flatten_lineage(cls, values):
        if "LineageEx" in values:
            lineage = values["LineageEx"]
            if isinstance(lineage, dict) and "Taxon" in lineage:
                values["taxa"] = lineage["Taxon"]
        return values


def parse_xml(xml_str: str):
    parsed = xmltodict.parse(xml_str)
    taxa = parsed['TaxaSet']['Taxon']
    if isinstance(taxa, dict):
        taxa = [taxa]
    try:
        result = [TaxonData(**taxon) for taxon in taxa]
    except ValidationError:
        print(taxa)
        raise
    return result


async def download_taxonomy(organism: str, output_path: str):
    client = request.HttpClient()
    # TODO: support multiple organisms?
    query = f'({organism}[Subtree]) AND species[Rank]'
    response = await client.esearch(query, db='taxonomy')
    print(f'Found {response.result_count} species in taxon {organism}')
    return await client.efetch_all(response, out_file=output_path, db='taxonomy', rettype='xml')

# from https://techwasti.com/python-tutorial-writing-python-data-class-to-csv
def dataclass_to_csv_row(csv_writer, row):
    csv_writer.writerow(row.__dict__.values())
def dataclass_to_csv_header(csv_writer, header):
    csv_writer.writerow(header.__annotations__.keys())


def _parse_lineage(taxa: list[Taxon]) -> dict[str, str | None]:
    ranks: dict[str, int | None] = {
        "domain": None,
        "kingdom": None,
        "phylum": None,
        "class": None,
        "order": None,
        "family": None,
        "genus": None,
    }
    for rank in ranks.keys():
        ranks[rank] = next((i for i, x in enumerate(taxa) if x.rank == rank), None)

    # Try find unofficial rank names by their suffix.
    # Suffixes from https://en.wikipedia.org/wiki/Taxonomic_rank#Terminations_of_names
    # See also https://www.iapt-taxon.org/nomen/pages/main/art_16.html
    n_taxa = len(taxa)
    reverse_taxa = reversed(taxa)
    def find_by_suffix(taxa: Iterator[Taxon], rank: str, suffixes: list[str]):
        if ranks[rank] is None:
            idx = next(
                (
                    i
                    for i, x in enumerate(taxa)
                    if any(x.name.endswith(suf) for suf in suffixes)
                ),
                None
            )
            if idx is not None:
                fwd_idx = n_taxa - idx - 1
                # avoid re-using names
                if fwd_idx not in ranks.values():
                    ranks[rank] = fwd_idx
    find_by_suffix(reverse_taxa, 'family', ['viridae', 'aceae', 'idae'])
    find_by_suffix(reverse_taxa, 'order', ['ales', 'ida', 'iformes'])
    find_by_suffix(reverse_taxa, 'class', ['viricetes', 'ia', 'opsida', 'phyceae', 'mycetes'])
    find_by_suffix(reverse_taxa, 'phylum', ['ota', 'ophyta'])

    # Unique for protists: Protists don't have an official kingdom.
    # Therefore, we use the lowest 'clade' rank between domain and phylum instead.
    if ranks['kingdom'] is None:
        domain_idx = ranks['domain']
        phylum_idx = ranks['phylum']
        if domain_idx is not None and phylum_idx is not None:
            if phylum_idx - domain_idx > 1:
                ranks['kingdom'] = phylum_idx - 1
    # TODO: if class is also unknown, this won't help!
    # Therefore, search lineage for phyla of all other species
    # Valid names don't contain spaces!

    return {rank: (taxa[idx].name if idx is not None else None) for rank, idx in ranks.items()}


"""
def _parse_lineage(taxa: list[Taxon]):
    wanted_ranks = {'domain', 'kingdom', 'phylum', 'order', 'family', 'class', 'genus'}
    lineage = {rank: [] if rank == "clade" else None for rank in wanted_ranks}
    for taxon in taxa:
        if taxon.rank in wanted_ranks:
            lineage[taxon.rank] = taxon.name
    return lineage
"""

def to_csv(taxa: list[TaxonData], csv_writer, header=False):
    if header:
        write_csv_header(csv_writer)
    for taxon in taxa:
        lineage = _parse_lineage(taxon.lineage)
        try:
            if taxon.rank != 'species' and lineage[taxon.rank] is None:
                lineage[taxon.rank] = taxon.name
        except KeyError:
            pass
        csv_writer.writerow([taxon.tax_id, taxon.name, lineage['domain'], lineage['kingdom'], lineage['phylum'], lineage['class'], lineage['order'], lineage['family'], lineage['genus']])

def write_csv_header(csv_writer):
    """Writes the header row. WARNING: Pay attention to the format when writing data"""
    csv_writer.writerow(['tax_id', 'species', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'])

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
                row[rank] = None
    return row


def _check_taxonomy_complete(meta: pd.DataFrame, tax: pd.DataFrame):
    missing_tax = meta[~meta['tax_id'].isin(tax['tax_id'])]
    if not missing_tax.empty:
        print('Warning: No taxonomy information found for the following accessions:')
        print(missing_tax)
    return missing_tax.empty

def filter_taxonomy(meta_path: str, input_tax: str, output_tax: str):
    meta = pd.read_csv(meta_path, sep='\t')
    tax = pd.read_csv(input_tax)
    _check_taxonomy_complete(meta, tax)
    tax_filtered = tax[tax['tax_id'].isin(meta['tax_id'])]
    tax_filled = tax_filtered.copy()
    tax_filled = tax_filled[tax_filled['species'].apply(sanitize_org_name)]
    tax_filled.to_csv(output_tax, sep='\t', index=False)

def merge_taxonomy(meta_in: str, tax: str) -> pd.DataFrame:
    meta = pd.read_csv(meta_in, sep='\t')
    tax_df = pd.read_csv(tax, sep='\t')
    meta_and_tax = meta.merge(tax_df, left_on='tax_id', right_on='tax_id', how='inner')
    meta_and_tax = meta_and_tax.rename(columns={'tax_id_x': 'tax_id'}).drop(columns=['tax_id_y'])
    return meta_and_tax


def parse_species_name(hosts: pd.Series):
    hosts = pd.Series(hosts.dropna().unique())
    hosts = hosts.str.replace(r'\(.*\)', '', regex=True)
    hosts = hosts.str.replace('  ', ' ')
    hosts = hosts.apply(database.sanitize_org_name)
    def to_genus(name: str) -> str:
        parts = name.split(' ')[:2]
        # remove non-alphabetical characters
        parts = list(map(lambda x: re.sub(r'[^a-zA-Z]', '', x), parts))
        if len(parts) == 2 and (parts[1] == 'sp' or parts[1] == 'spp'):
            return parts[0]
        return ' '.join(parts)
    hosts = hosts.apply(to_genus)
    return pd.Series(hosts.unique()).astype(str)


async def get_taxid_by_name(name: str, client: request.HttpClient) -> int | None:
    response = await client.esearch(name, 'taxonomy')
    id = None
    if response.result_count > 0:
        bs = BeautifulSoup(response.text_response, features='xml')
        id = int(bs.find('Id').string)
    return id

async def by_taxids(tax_ids: list[str]):
    if len(tax_ids) == 0:
        raise AssertionError("No tax_ids were passed to by_taxids!")
    client = request.NcbiApiClient()
    response = await client.post('taxonomy', 'dataset_report', json={
        "taxons": tax_ids
    })
    return response['reports']

async def get_taxids_by_names(names: pd.Series, client: request.HttpClient) -> dict[str, int]:
    tax_ids = {}
    for species in names:
        time.sleep(0.1)
        tax_ids[species] = await get_taxid_by_name(species, client)
    return tax_ids


async def fetch_by_tax_ids(tax_ids: list[str], client: request.HttpClient, out_file: str):
    with open(out_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        write_csv_header(csv_writer)
        def parse_output(output: str):
            taxa = parse_xml(output)
            to_csv(taxa, csv_writer)
        await client.efetch_by_ids('taxonomy', tax_ids, parse_output)

async def from_organsim_names(data: pd.Series, tax_file: str):
    org_names = parse_species_name(data)
    client = request.HttpClient()
    tax_ids = await get_taxids_by_names(org_names, client)
    tax_ids = [str(id) for id in tax_ids.values() if id is not None]
    await fetch_by_tax_ids(tax_ids, client, tax_file)

async def info_by_organism_name(organism_name: str):
    client = request.HttpClient()
    tax_id = await get_taxid_by_name(organism_name, client)
    if tax_id is None:
        raise ValueError(f'Taxon {tax_id} not found in NCBI taxonomny database')
    tax_xml = await client.efetch_by_id('taxonomy', str(tax_id))
    return parse_xml(tax_xml)[0]

def get_sisters(parent_data: dict, exclude_child: int) -> list[int]:
    sisters = []
    if len(parent_data) > 0:
        sisters = parent_data['taxonomy']['children']
        assert exclude_child in sisters
        sisters.remove(exclude_child)
    return sisters

def load_tax_qiime(path: str) -> pd.DataFrame:
    tax_df = pd.read_csv(path, sep='\t', header=None)
    tax_df[['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = tax_df[1].str.split(';', expand=True)
    return tax_df.set_index(0)
