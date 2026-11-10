"""The CheckAMG-derived AVG reference tables, read from the phu-db bundle.

The bundle's ``avg`` database holds the upstream CheckAMG tables under stable
names. load_avg_reference() turns them into the rows the AVG workflow uses:
Pfam and KOfam entries only, accessions normalised, one row per identity.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

SCORES_FILE = "v_scores.tsv"
_POSITIVE_SOURCES = {
    "positive_amg.tsv": "putative_amg",
    "positive_apg.tsv": "putative_apg",
    "positive_areg.tsv": "putative_areg",
}
_FILTER_SOURCES = {
    "filters_amg.tsv": "putative_amg",
    "filters_apg.tsv": "putative_apg",
    "filters_areg.tsv": "putative_areg",
}
_FILTER_COLUMNS = (
    "filter_essential",
    "filter_glucan",
    "filter_nucleotide",
    "filter_methyl",
    "filter_lipid",
)


@dataclass(frozen=True)
class AvgReference:
    """Normalised reference rows, each sorted by identity."""

    scores: list[dict[str, str]]
    positive: list[dict[str, str]]
    filters: list[dict[str, str]]


def _normalize_database(value: str) -> str | None:
    database = value.strip().lower()
    if database == "pfam":
        return "pfam"
    if database in {"kegg", "ko", "kofam"}:
        return "kofam"
    return None


def normalize_accession(database: str, value: str) -> str:
    accession = value.strip().upper()
    if database == "pfam":
        accession = accession.split(".", 1)[0]
        if not accession.startswith("PF") or len(accession) != 7 or not accession[2:].isdigit():
            raise ValueError(f"Invalid Pfam accession: {value}")
    elif database == "kofam":
        if not accession.startswith("K") or len(accession) != 6 or not accession[1:].isdigit():
            raise ValueError(f"Invalid KOfam accession: {value}")
    else:
        raise ValueError(f"Unsupported AVG database: {database}")
    return accession


def _read_rows(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"AVG reference table has no header: {path.name}")
        required = {"id", "db", "name"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"AVG reference table {path.name} is missing columns: {', '.join(sorted(missing))}"
            )
        yield from reader


def _selected_row(row: dict[str, str]) -> tuple[str, str] | None:
    database = _normalize_database(row.get("db", ""))
    if database is None:
        return None
    return database, normalize_accession(database, row.get("id", ""))


def load_avg_reference(root: Path) -> AvgReference:
    """Read and normalise the ``avg`` database in ``root`` (see resolve_db)."""
    root = Path(root)
    positive_rows: dict[tuple[str, str, str], dict[str, str]] = {}
    filter_rows: dict[tuple[str, str, str, str], dict[str, str]] = {}
    score_rows: dict[tuple[str, str], dict[str, str]] = {}

    for source_name, classification in _POSITIVE_SOURCES.items():
        for row in _read_rows(root / source_name):
            key = _selected_row(row)
            if key is None:
                continue
            database, accession = key
            normalized = {
                "database": database,
                "accession": accession,
                "original_accession": row["id"].strip(),
                "proposed_class": classification,
                "name": row.get("name", "").strip(),
                "metabolic_paths": row.get("metabolic_paths", "").strip(),
                "total_paths": row.get("total_paths", "").strip(),
                "path_source": row.get("path_source", "").strip(),
                "metabolic_ratio": row.get("metabolic_ratio", "").strip(),
                "vl_score": row.get("VL-score", "").strip(),
                "amg_weight": row.get("AMG_weight", "").strip(),
                "amg_level": row.get("amg_level", "").strip(),
            }
            identity = (database, accession, classification)
            previous = positive_rows.get(identity)
            if previous is not None and previous != normalized:
                raise ValueError(f"Conflicting AVG positive rows for {identity}")
            positive_rows[identity] = normalized

    for source_name, classification in _FILTER_SOURCES.items():
        for row in _read_rows(root / source_name):
            key = _selected_row(row)
            if key is None:
                continue
            database, accession = key
            for category in _FILTER_COLUMNS:
                if row.get(category, "").strip().lower() != "true":
                    continue
                normalized = {
                    "database": database,
                    "accession": accession,
                    "original_accession": row["id"].strip(),
                    "proposed_class": classification,
                    "name": row.get("name", "").strip(),
                    "filter_category": category,
                }
                identity = (database, accession, classification, category)
                previous = filter_rows.get(identity)
                if previous is not None and previous != normalized:
                    raise ValueError(f"Conflicting AVG filter rows for {identity}")
                filter_rows[identity] = normalized

    for row in _read_rows(root / SCORES_FILE):
        key = _selected_row(row)
        if key is None:
            continue
        database, accession = key
        normalized = {
            "database": database,
            "accession": accession,
            "original_accession": row["id"].strip(),
            "name": row.get("name", "").strip(),
            "v_score": row.get("V-score", "").strip(),
            "vl_score": row.get("VL-score", "").strip(),
        }
        previous = score_rows.get(key)
        if previous is not None and previous != normalized:
            raise ValueError(f"Conflicting AVG score rows for {key}")
        score_rows[key] = normalized

    return AvgReference(
        scores=[score_rows[key] for key in sorted(score_rows)],
        positive=[positive_rows[key] for key in sorted(positive_rows)],
        filters=[filter_rows[key] for key in sorted(filter_rows)],
    )
