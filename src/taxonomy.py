from bs4 import BeautifulSoup
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from src import request, database
from src.ncbi_taxonomy import NCBITaxonomyDB, TaxonData, Taxon
import pandas as pd
import csv
import re
import time


async def download_taxonomy(organism: str, db_path: str):
    # TODO: support multiple organisms?
    client = request.HttpClient()
    query_tax_id = await get_taxid_by_name(organism, client)
    if query_tax_id is None:
        raise ValueError(f'Unknown organism {organism}')
    with NCBITaxonomyDB(db_path) as db:
        tax_data = db.get_descendants_bulk(query_tax_id)
    return tax_data

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

def to_csv(taxa: Sequence[TaxonData], csv_writer, header=False):
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


def update_tax_ids(tax_db_path: str, meta_df: pd.DataFrame, index_col: str | None=None) -> pd.DataFrame:
    meta_df = meta_df.copy()
    meta_df['tax_id'] = pd.to_numeric(meta_df['tax_id'], errors='coerce')
    meta_df['tax_id'] = meta_df['tax_id'].astype(int)
    unique_tax_ids = meta_df['tax_id'].unique().tolist()
    with NCBITaxonomyDB(tax_db_path) as db:
        id_dict = {old_id: tax.tax_id for old_id, tax in db.resolve_tax_ids(unique_tax_ids).items()}
        id_df = pd.DataFrame({'old': id_dict.keys(), 'new': id_dict.values()})
        assert not id_df['new'].isna().any()
    if index_col is not None:
        meta_df.insert(loc=0, column=index_col, value=meta_df.index)
    # replace deprecated tax ids with the correct tax id
    merged = (meta_df
            .merge(id_df, how='left', left_on='tax_id', right_on='old')
            .drop(columns=['old', 'tax_id'])
            .rename(columns={'new': 'tax_id'})
    )

    with NCBITaxonomyDB(tax_db_path) as db:
        unique_tax_ids = merged['tax_id'].unique().tolist()
        species_tax = db.find_parent_at_rank(tax_ids=unique_tax_ids, rank='species')
        species_tax_ids = {old_id: tax.tax_id for old_id, tax in species_tax.items() if tax is not None}
        spec_id_df = pd.DataFrame({'old': species_tax_ids.keys(), 'new': species_tax_ids.values()})
    return (merged
            .merge(spec_id_df, how='left', left_on='tax_id', right_on='old')
            # remove sequences from taxa higher than species
            .dropna(subset=['new'])
            .drop(columns=['old', 'tax_id'])
            .rename(columns={'new': 'tax_id'})
    )


# When a species is being moved to a new genus, NCBI temporarily puts the old genus name in brackets
# for the species, e.g. [Oldgenus] species, whereas the genus is laready the new genus.
# I.e. Newgenus, [Oldgenus] species
# To avoid confusion, we replace [Oldgenus] in the species name with Newgenus:
# I.e. Newgenus, Newgenus species
def update_species_name(tax_df: pd.DataFrame) -> pd.DataFrame:
    tax_df = tax_df.copy()
    has_bracketed_genus = tax_df['species'].str.match(r'^\[.*?\]')
    tax_df.loc[has_bracketed_genus, 'species'] = tax_df.loc[has_bracketed_genus].apply(
        lambda row: re.sub(r'^\[.*?\]', row['genus'], row['species']),
        axis=1,
    )
    return tax_df


async def filter_descendants(tax_db_path: str, meta_df: pd.DataFrame, ancestor_name: str):
    client = request.HttpClient()
    ancestor_tax_id = await get_taxid_by_name(ancestor_name, client)
    if ancestor_tax_id is None:
        raise ValueError(f'Unknown organism {ancestor_name}')
    with NCBITaxonomyDB(tax_db_path) as db:
        tax_ids_in_ancestor = db.filter_descendants(meta_df['tax_id'].unique().tolist(), ancestor_tax_id)
    return meta_df[meta_df['tax_id'].isin(tax_ids_in_ancestor)]


def merge_taxonomy(meta_df: pd.DataFrame, tax_df: pd.DataFrame) -> pd.DataFrame:
    assert _check_taxonomy_complete(meta_df, tax_df)
    meta_and_tax = meta_df.merge(tax_df, left_on='tax_id', right_on='tax_id', how='inner')
    meta_and_tax['organism_name'] = meta_and_tax['species']
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

def by_taxids(tax_db_path: str, tax_ids: list[int]) -> list[TaxonData]:
    if len(tax_ids) == 0:
        raise AssertionError("No tax_ids were passed to by_taxids!")
    with NCBITaxonomyDB(tax_db_path) as db:
        return db.get_taxon_data_bulk(tax_ids)

async def get_taxids_by_names(names: pd.Series, client: request.HttpClient) -> dict[str, int]:
    tax_ids = {}
    for species in names:
        time.sleep(0.1)
        tax_ids[species] = await get_taxid_by_name(species, client)
    return tax_ids


def fetch_by_tax_ids(tax_db_path: str, tax_ids: list[int], out_file: str):
    with NCBITaxonomyDB(tax_db_path) as db:
        tax_data = db.get_taxon_data_bulk(tax_ids)
    with open(out_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        to_csv(tax_data, csv_writer, header=True)

async def from_organsim_names(tax_db_path: str, data: pd.Series, tax_file: str):
    org_names = parse_species_name(data)
    client = request.HttpClient()
    tax_ids = await get_taxids_by_names(org_names, client)
    tax_ids = [id for id in tax_ids.values() if id is not None]
    fetch_by_tax_ids(tax_db_path, tax_ids, tax_file)

def info_by_organism_name(tax_db_path: str, organism_name: str) -> TaxonData:
    with NCBITaxonomyDB(tax_db_path) as db:
        tax_id = db.find_tax_id_by_name(organism_name)
        if tax_id is None:
            raise ValueError(f'Taxon {organism_name!r} not found in local taxonomy database')
        return db.get_taxon_data(tax_id)


@dataclass(frozen=True)
class SisterTaxon:
    taxon: Taxon
    parent: Taxon

    def __eq__(self, value: object) -> bool:
        return self.taxon.tax_id == value.taxon.tax_id

    def __hash__(self) -> int:
        return self.taxon.tax_id


# drop environmental samples (child of Eukaryota) and unclassified taxa
def filter_sisters(sisters: Iterable[SisterTaxon]):
    return set(filter(lambda x: not x.taxon.name.startswith('environmental') and not x.taxon.name.startswith('unclassified') and x.taxon.name != 'other entries', sisters))


def get_sister_taxa(tax_db_path: str, target_organism: str) -> set[SisterTaxon]:
    tax_info = info_by_organism_name(tax_db_path, target_organism)

    # lineage is ordered from top to bottom, but we want to start from the bottom
    lineage = list(reversed(tax_info.lineage))
    taxon = Taxon(TaxId=tax_info.tax_id, ScientificName=tax_info.name, Rank=tax_info.rank)
    sisters: set[SisterTaxon] = set()

    # get sister tax ids
    with NCBITaxonomyDB(tax_db_path) as db:
        for parent, exclude_child in zip(lineage, [taxon, *lineage[:-1]]):
            sibling_taxa = db.get_direct_children(parent.tax_id)
            current_sisters = [c for c in sibling_taxa if c.tax_id != exclude_child.tax_id]
            sisters.update([SisterTaxon(sister, parent) for sister in current_sisters])
    return filter_sisters(sisters)


def get_sister_childs(tax_db_path: str, sisters: list[SisterTaxon]) -> set[SisterTaxon]:
    children: set[SisterTaxon] = set()
    with NCBITaxonomyDB(tax_db_path) as db:
        for item in sisters:
            child_taxa = db.get_direct_children(item.taxon.tax_id)
            print(f"Taxon {item.taxon.name}({item.taxon.tax_id}) has {len(child_taxa)} kids")
            for child in child_taxa:
                children.add(SisterTaxon(child, item.parent))
    return filter_sisters(children)


def load_tax_qiime(path: str) -> pd.DataFrame:
    tax_df = pd.read_csv(path, sep='\t', header=None)
    tax_df[['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = tax_df[1].str.split(';', expand=True)
    return tax_df.set_index(0)
