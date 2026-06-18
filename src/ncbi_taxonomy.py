import os
import sqlite3
import tarfile
from hashlib import md5
from itertools import batched
from urllib.request import urlretrieve
from pydantic import BaseModel, Field, model_validator


__all__ = ["NCBITaxonomyDB", "download_and_build_db"]

DB_VERSION = 2
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDUMP_MD5_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"

# SQLite only allows < 999 variables in a query; Use 900 as a safe margin
SQLITE_MAX_VARS = 900


class Taxon(BaseModel):
    tax_id: int = Field(alias='TaxId')
    name: str = Field(alias='ScientificName')
    rank: str = Field(alias='Rank')
    model_config = {"populate_by_name": True}

class TaxonData(BaseModel):
    tax_id: int = Field(alias='TaxId')
    name: str = Field(alias='ScientificName')
    rank: str = Field(alias='Rank')
    parent_tax_id: int = Field(alias='ParentTaxId')
    lineage: list[Taxon] = Field(alias='taxa')
    model_config = {"populate_by_name": True}

    @model_validator(mode="before")
    @classmethod
    def flatten_lineage(cls, values):
        if "LineageEx" in values:
            lineage = values["LineageEx"]
            if isinstance(lineage, dict) and "Taxon" in lineage:
                values["taxa"] = lineage["Taxon"]
        return values


def _download_taxdump(taxdump_dir: str) -> str:
    """Download taxdump.tar.gz into taxdump_dir (skipped if md5 matches). Returns local path."""
    os.makedirs(taxdump_dir, exist_ok=True)
    targz_path = os.path.join(taxdump_dir, "taxdump.tar.gz")

    md5_tmp, _ = urlretrieve(TAXDUMP_MD5_URL)
    with open(md5_tmp, "r") as f:
        expected_md5 = f.readline().split()[0]

    if os.path.exists(targz_path):
        with open(targz_path, "rb") as f:
            actual_md5 = md5(f.read()).hexdigest()
        if actual_md5 == expected_md5:
            print(f"taxdump.tar.gz is up-to-date at {targz_path}")
            return targz_path

    print(f"Downloading taxdump.tar.gz to {targz_path}...")
    urlretrieve(TAXDUMP_URL, targz_path)
    print("Download complete.")
    return targz_path


def _parse_names(
    tar: tarfile.TarFile,
) -> tuple[dict[int, str], dict[int, str], set[tuple[int, str]]]:
    """Parse names.dmp. Returns (sci_names, common_names, synonyms)."""
    SYNONYM_TYPES = {
        "synonym", "equivalent name", "genbank equivalent name",
        "anamorph", "genbank synonym", "genbank anamorph", "teleomorph",
    }
    sci_names: dict[int, str] = {}
    common_names: dict[int, str] = {}
    synonyms: set[tuple[int, str]] = set()
    seen_lower: set[tuple[int, str]] = set()

    for raw in tar.extractfile("names.dmp"):
        fields = [f.strip() for f in raw.decode().split("|")]
        taxid = int(fields[0])
        name = fields[1].strip('"')
        name_class = fields[3].lower()

        if name_class == "scientific name":
            sci_names[taxid] = name
        elif name_class == "genbank common name":
            common_names[taxid] = name
        elif name_class in SYNONYM_TYPES:
            key = (taxid, name.lower())
            if key not in seen_lower:
                seen_lower.add(key)
                synonyms.add((taxid, name))

    return sci_names, common_names, synonyms


def _parse_nodes(tar: tarfile.TarFile) -> dict[int, tuple[int, str]]:
    """Parse nodes.dmp. Returns {taxid: (parent_taxid, rank)}."""
    nodes: dict[int, tuple[int, str]] = {}
    for raw in tar.extractfile("nodes.dmp"):
        fields = raw.decode().split("|")
        taxid = int(fields[0].strip())
        parent = int(fields[1].strip())
        rank = fields[2].strip()
        nodes[taxid] = (parent, rank)
    return nodes


def _parse_merged(tar: tarfile.TarFile) -> list[tuple[int, int]]:
    """Parse merged.dmp. Returns list of (old_taxid, new_taxid)."""
    merged: list[tuple[int, int]] = []
    for raw in tar.extractfile("merged.dmp"):
        parts = [f.strip() for f in raw.decode().split("|")]
        if len(parts) >= 2 and parts[0] and parts[1]:
            merged.append((int(parts[0]), int(parts[1])))
    return merged


def _compute_tracks(parent_map: dict[int, int]) -> dict[int, str]:
    """
    Compute the lineage track string for every taxid.

    Track format: "taxid,parent_taxid,...,root_taxid"  (leaf → root),
    matching the column written by ete3's generate_table().

    Uses iterative memoisation to avoid Python recursion limits on deep trees.
    """
    tracks: dict[int, str] = {}

    for start in parent_map:
        if start in tracks:
            continue

        # Walk up until we reach a node whose track is already known
        path: list[int] = []
        current = start
        while current not in tracks:
            path.append(current)
            parent = parent_map.get(current, current)
            if parent == current:          # root node (self-referential)
                tracks[current] = str(current)
                path.pop()
                break
            current = parent

        # Fill in tracks top-down along the discovered path
        for node in reversed(path):
            parent = parent_map[node]
            tracks[node] = f"{node},{tracks[parent]}"

    return tracks


def _build_sqlite(
    db_path: str,
    sci_names: dict[int, str],
    common_names: dict[int, str],
    nodes: dict[int, tuple[int, str]],
    synonyms: set[tuple[int, str]],
    merged: list[tuple[int, int]],
    tracks: dict[int, str],
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.executescript("""
        DROP TABLE IF EXISTS stats;
        DROP TABLE IF EXISTS tax_meta;
        DROP TABLE IF EXISTS synonym;
        DROP TABLE IF EXISTS merged;

        CREATE TABLE stats (version INT PRIMARY KEY);

        CREATE TABLE tax_meta (
            taxid  INT PRIMARY KEY,
            parent INT,
            spname VARCHAR(50) COLLATE NOCASE,
            common VARCHAR(50) COLLATE NOCASE,
            rank   VARCHAR(50),
            track  TEXT
        );

        CREATE TABLE synonym (
            taxid  INT,
            spname VARCHAR(50) COLLATE NOCASE,
            PRIMARY KEY (spname, taxid)
        );

        CREATE TABLE merged (taxid_old INT, taxid_new INT);

        CREATE INDEX spname1 ON tax_meta (spname COLLATE NOCASE);
        CREATE INDEX spname2 ON synonym  (spname COLLATE NOCASE);
        CREATE INDEX idx_species_rank ON tax_meta (rank);
    """)

    conn.execute("INSERT INTO stats (version) VALUES (?)", (DB_VERSION,))

    print("Inserting species rows...")
    species_rows = [
        (
            tid,
            parent,
            sci_names.get(tid, ""),
            common_names.get(tid, ""),
            rank,
            tracks.get(tid, str(tid)),
        )
        for tid, (parent, rank) in nodes.items()
    ]
    conn.executemany(
        "INSERT INTO tax_meta (taxid, parent, spname, common, rank, track) VALUES (?,?,?,?,?,?)",
        species_rows,
    )

    print("Inserting synonym rows...")
    conn.executemany("INSERT OR IGNORE INTO synonym (taxid, spname) VALUES (?,?)", synonyms)

    print("Inserting merged rows...")
    conn.executemany("INSERT INTO merged (taxid_old, taxid_new) VALUES (?,?)", merged)

    conn.commit()
    conn.close()


def download_and_build_db(
    db_path: str,
    taxdump_dir: str,
    force_update: bool = False,
) -> None:
    """
    Download taxdump.tar.gz from NCBI into taxdump_dir, parse it, and write
    the taxonomy SQLite database to db_path.

    Skips downloading (and rebuilding) if the existing database is already at
    the current schema version, unless force_update=True.

    Args:
        db_path:      Path where the SQLite file will be written.
        taxdump_dir:  Directory where taxdump.tar.gz is downloaded/cached.
        force_update: Re-download and rebuild even if the DB appears current.
    """
    if not force_update and os.path.exists(db_path):
        try:
            conn = sqlite3.connect(db_path)
            row = conn.execute("SELECT version FROM stats").fetchone()
            conn.close()
            if row and row[0] == DB_VERSION:
                print(f"Database at {db_path} is already up-to-date.")
                return
        except sqlite3.OperationalError:
            pass

    targz_path = _download_taxdump(taxdump_dir)

    print("Parsing taxdump archive...")
    with tarfile.open(targz_path, "r:gz") as tar:
        sci_names, common_names, synonyms = _parse_names(tar)
        nodes = _parse_nodes(tar)
        merged = _parse_merged(tar)

    print(f"  {len(sci_names):,} names  |  {len(nodes):,} nodes  |  {len(synonyms):,} synonyms")

    print("Computing lineage tracks...")
    parent_map = {tid: parent for tid, (parent, _) in nodes.items()}
    tracks = _compute_tracks(parent_map)

    print(f"Writing database to {db_path}...")
    _build_sqlite(db_path, sci_names, common_names, nodes, synonyms, merged, tracks)
    print("Done.")


class NCBITaxonomyDB:
    def __init__(self, db_path: str):
        """
        Open an existing taxonomy SQLite database for querying.

        The database must be built beforehand with download_and_build_db().

        Args:
            db_path: Path to the SQLite taxonomy database file.

        Raises:
            FileNotFoundError: If db_path does not exist.
        """
        if not os.path.exists(db_path):
            raise FileNotFoundError(
                f"Taxonomy database not found at '{db_path}'. "
                "Build it first with download_and_build_db()."
            )
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA cache_size=-65536")
        self.conn.commit()

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def _parse_track(self, track: str) -> list[int]:
        """Return lineage as [root_taxid, ..., parent_taxid, taxid] (root first)."""
        return [int(x) for x in reversed(track.split(","))]

    def _query_species(self, ids: list[int]) -> list[sqlite3.Row]:
        rows: list[sqlite3.Row] = []
        for chunk in batched(ids, SQLITE_MAX_VARS):
            placeholders = ",".join("?" * len(chunk))
            rows.extend(
                self.conn.execute(
                    f"SELECT taxid, parent, spname, rank, track FROM tax_meta "
                    f"WHERE taxid IN ({placeholders})",
                    list(chunk),
                ).fetchall()
            )
        return rows


    def resolve_tax_ids(self, tax_ids: list[int]) -> dict[int, Taxon]:
        """
        Map a list of (possibly deprecated) tax IDs to their current equivalents.

        Returns a dict keyed by the original input IDs. IDs that are already
        current map to themselves. Raises ValueError for any ID not found after
        resolution.
        """
        merged: dict[int, int] = {}
        for chunk in batched(tax_ids, SQLITE_MAX_VARS):
            placeholders = ",".join("?" * len(chunk))
            rows = self.conn.execute(
                f"SELECT taxid_old, taxid_new FROM merged WHERE taxid_old IN ({placeholders})",
                list(chunk),
            ).fetchall()
            merged.update({r["taxid_old"]: r["taxid_new"] for r in rows})

        def follow(tid: int) -> int:
            seen: set[int] = set()
            while tid in merged and tid not in seen:
                seen.add(tid)
                tid = merged[tid]
            return tid

        resolved: dict[int, int] = {tid: follow(tid) for tid in tax_ids}
        current_ids = list(set(resolved.values()))
        species_rows: dict[int, sqlite3.Row] = {
            row["taxid"]: row for row in self._query_species(current_ids)
        }

        missing = [current_id for current_id in resolved.values() if current_id not in species_rows]
        if missing:
            raise ValueError(f"Resolved tax IDs not found in database: {missing}")

        return {
            original_id: Taxon.model_validate({
                "TaxId": current_id,
                "ScientificName": species_rows[current_id]["spname"],
                "Rank": species_rows[current_id]["rank"],
            })
            for original_id, current_id in resolved.items()
        }


    def get_taxon_data(self, tax_id: int) -> TaxonData:
        """Single taxon lookup. Expects a current (non-deprecated) tax ID."""
        return self.get_taxon_data_bulk([tax_id])[0]

    def get_taxon_data_bulk(self, tax_ids: list[int]) -> list[TaxonData]:
        """
        Bulk fetch taxonomy data including full lineage for each taxid.

        Expects current (non-deprecated) tax IDs. Call resolve_tax_ids() first
        if your input may contain deprecated IDs.

        Returns results in the same order as the input list.
        """
        primary_rows: dict[int, sqlite3.Row] = {
            row["taxid"]: row for row in self._query_species(tax_ids)
        }

        all_lineage_ids: set[int] = set()
        parsed_tracks: dict[int, list[int]] = {}

        for tax_id in tax_ids:
            if tax_id not in primary_rows:
                continue
            track_ids = self._parse_track(primary_rows[tax_id]["track"])
            parsed_tracks[tax_id] = track_ids
            all_lineage_ids.update(track_ids[:-1])  # exclude the root placeholder

        lineage_ids_list = list(all_lineage_ids - set(tax_ids))
        lineage_lookup: dict[int, sqlite3.Row] = dict(primary_rows)
        if lineage_ids_list:
            lineage_lookup.update({
                row["taxid"]: row
                for row in self._query_species(lineage_ids_list)
            })

        results: list[TaxonData] = []
        for tax_id in tax_ids:
            if tax_id not in primary_rows:
                raise ValueError(
                    f"taxid {tax_id} not found — did you forget resolve_tax_ids()?"
                )

            row = primary_rows[tax_id]
            track_ids = parsed_tracks[tax_id]

            lineage_taxa = [
                Taxon.model_validate({
                    "TaxId": tid,
                    "ScientificName": lineage_lookup[tid]["spname"] if tid in lineage_lookup else "unknown",
                    "Rank": lineage_lookup[tid]["rank"] if tid in lineage_lookup else "no rank",
                })
                for tid in track_ids[:-1]
            ]

            results.append(TaxonData.model_validate({
                "TaxId": tax_id,
                "ScientificName": row["spname"],
                "Rank": row["rank"],
                "ParentTaxId": row["parent"],
                "taxa": [t.model_dump(by_alias=True) for t in lineage_taxa],
            }))

        return results

    def get_descendants_bulk(self, ancestor_tax_id: int, rank: str = "species") -> list[TaxonData]:
        """
        Return TaxonData for all taxa of the given rank that descend from
        ancestor_tax_id (inclusive of ancestor itself if it matches the rank).

        Uses the pre-computed track column — no recursive queries needed.
        """
        rows = self.conn.execute(
            """
            SELECT taxid FROM tax_meta
            WHERE rank = ?
            AND INSTR(',' || track || ',', ',' || ? || ',') > 0
            """,
            (rank, str(ancestor_tax_id)),
        ).fetchall()
        descendant_ids = [r["taxid"] for r in rows]
        if not descendant_ids:
            return []
        return self.get_taxon_data_bulk(descendant_ids)

    def get_direct_children(self, tax_id: int) -> list[Taxon]:
        """Return immediate children of tax_id (one level only)."""
        rows = self.conn.execute(
            "SELECT taxid, spname, rank FROM tax_meta WHERE parent = ? AND taxid != ?",
            (tax_id, tax_id),
        ).fetchall()
        return [
            Taxon.model_validate({"TaxId": r["taxid"], "ScientificName": r["spname"], "Rank": r["rank"]})
            for r in rows
        ]

    def find_tax_id_by_name(self, name: str) -> int | None:
        """Look up a tax ID by scientific name or synonym (case-insensitive)."""
        row = self.conn.execute(
            "SELECT taxid FROM tax_meta WHERE spname = ? LIMIT 1", (name,)
        ).fetchone()
        if row:
            return row["taxid"]
        row = self.conn.execute(
            "SELECT taxid FROM synonym WHERE spname = ? LIMIT 1", (name,)
        ).fetchone()
        return row["taxid"] if row else None

    def filter_descendants(self, tax_ids: list[int], ancestor_tax_id: int) -> list[int]:
        """
        Return the subset of tax_ids that are descendants of ancestor_tax_id
        (inclusive — ancestor_tax_id itself is returned if present in tax_ids).

        Uses the pre-computed track column: a taxid is a descendant if
        ancestor_tax_id appears in its lineage string. Unknown tax IDs are
        silently excluded. Output order is not guaranteed to match input order.

        Args:
            tax_ids:         Candidate tax IDs to test.
            ancestor_tax_id: The higher-rank ancestor to test membership against.

        Returns:
            List of tax IDs from tax_ids that are descendants of ancestor_tax_id.
        """
        ancestor_str = str(ancestor_tax_id)
        matched: list[int] = []
        for chunk in batched(tax_ids, SQLITE_MAX_VARS):
            placeholders = ",".join("?" * len(chunk))
            rows = self.conn.execute(
                f"SELECT taxid FROM tax_meta "
                f"WHERE taxid IN ({placeholders}) "
                f"AND INSTR(',' || track || ',', ',' || ? || ',') > 0",
                list(chunk) + [ancestor_str],
            ).fetchall()
            matched.extend(r["taxid"] for r in rows)
        return matched

    def find_parent_at_rank(self, tax_ids: list[int], rank: str) -> dict[int, Taxon | None]:
        """
        For each tax ID, return the nearest ancestor (including itself) that has
        the given NCBI rank (e.g. "species", "genus", "family").

        Returns a dict keyed by the original input IDs. The value is None when
        no ancestor with that rank exists in the lineage, or when the tax ID is
        not found in the database.

        All queries are batched — the full set of lineage IDs is resolved in
        two SQL round-trips regardless of how many tax IDs are provided.
        """
        # 1. Fetch full rows (incl. rank) for all input IDs
        input_data: dict[int, sqlite3.Row] = {
            row["taxid"]: row for row in self._query_species(tax_ids)
        }

        # 2. Short-circuit: inputs already at the target rank need no ancestor lookup
        result: dict[int, Taxon | None] = {}
        needs_lookup: list[int] = []
        for tid in tax_ids:
            if tid not in input_data.keys():
                result[tid] = None
            elif input_data[tid]["rank"] == rank:
                r = input_data[tid]
                result[tid] = Taxon.model_validate({
                    "TaxId": r["taxid"],
                    "ScientificName": r["spname"],
                    "Rank": r["rank"],
                })
            else:
                needs_lookup.append(tid)

        if len(needs_lookup) == 0:
            return result

        # 3. Collect ancestor IDs for the remaining inputs, excluding already-fetched taxids
        lineages: dict[int, list[int]] = {}
        ancestor_ids: set[int] = set()
        for tid in needs_lookup:
            ids = self._parse_track(input_data[tid]["track"])  # [root, ..., parent, tid]
            lineages[tid] = ids
            ancestor_ids.update(ids)
        ancestor_ids -= set(input_data)

        # 4. Fetch ancestor rows via _query_species (primary-key lookup), filter by rank in Python
        ancestor_data: dict[int, sqlite3.Row] = {
            row["taxid"]: row for row in self._query_species(list(ancestor_ids))
        }
        ranked: dict[int, sqlite3.Row] = {
            tid: row
            for tid, row in (*input_data.items(), *ancestor_data.items())
            if row["rank"] == rank
        }

        # 5. For each remaining input walk leaf → root, return deepest ranked match
        for tid in needs_lookup:
            match = None
            for ancestor in reversed(lineages[tid]):  # reversed = leaf → root
                if ancestor in ranked:
                    r = ranked[ancestor]
                    match = Taxon.model_validate({
                        "TaxId": r["taxid"],
                        "ScientificName": r["spname"],
                        "Rank": r["rank"],
                    })
                    break
            result[tid] = match

        return result
