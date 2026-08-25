from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

AVG_RELEASE = "v1.1.1"
AVG_NAME = "avg"
AVG_SOURCE_REPOSITORY = "https://github.com/AnantharamanLab/CheckAMG"
AVG_SOURCE_ROOT = (
    "https://raw.githubusercontent.com/AnantharamanLab/CheckAMG/"
    f"{AVG_RELEASE}/CheckAMG/files"
)
AVG_SOURCE_FILES = (
    "AMGs.tsv",
    "APGs.tsv",
    "AReGs.tsv",
    "AMG_filters.tsv",
    "APG_filters.tsv",
    "AReG_filters.tsv",
    "vscores.tsv",
)
NORMALIZED_FILES = ("v_scores.tsv", "avg_positive.tsv", "avg_filters.tsv")
_ALLOWED_DATABASES = {"pfam", "kofam"}
_POSITIVE_SOURCES = {
    "AMGs.tsv": "putative_amg",
    "APGs.tsv": "putative_apg",
    "AReGs.tsv": "putative_areg",
}
_FILTER_SOURCES = {
    "AMG_filters.tsv": "putative_amg",
    "APG_filters.tsv": "putative_apg",
    "AReG_filters.tsv": "putative_areg",
}
_FILTER_COLUMNS = (
    "filter_essential",
    "filter_glucan",
    "filter_nucleotide",
    "filter_methyl",
    "filter_lipid",
)


def _db_root() -> Path:
    if db_env := os.environ.get("PHU_DB_FOLDER"):
        return Path(db_env)
    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


def _avg_root() -> Path:
    return _db_root() / AVG_NAME


def _atomic_write_bytes(destination: Path, content: bytes) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as handle:
        temporary_path = Path(handle.name)
        handle.write(content)
    temporary_path.replace(destination)


def _atomic_write_text(destination: Path, content: str) -> None:
    _atomic_write_bytes(destination, content.encode("utf-8"))


def _download(url: str, destination: Path) -> None:
    with urlopen(url) as response:
        _atomic_write_bytes(destination, response.read())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> int:
    from io import StringIO

    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    count = 0
    for row in rows:
        writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fieldnames})
        count += 1
    _atomic_write_text(path, buffer.getvalue())
    return count


def _compile_tables(root: Path) -> dict[str, int]:
    positive_rows: dict[tuple[str, str, str], dict[str, object]] = {}
    filter_rows: dict[tuple[str, str, str, str], dict[str, object]] = {}
    score_rows: dict[tuple[str, str], dict[str, object]] = {}

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

    for row in _read_rows(root / "vscores.tsv"):
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

    counts = {
        "positive": _write_tsv(
            root / "avg_positive.tsv",
            list(next(iter(positive_rows.values()), {"database": ""}).keys()),
            sorted(positive_rows.values(), key=lambda row: (row["database"], row["accession"], row["proposed_class"])),
        ),
        "filters": _write_tsv(
            root / "avg_filters.tsv",
            list(next(iter(filter_rows.values()), {"database": ""}).keys()),
            sorted(filter_rows.values(), key=lambda row: (row["database"], row["accession"], row["proposed_class"], row["filter_category"])),
        ),
        "scores": _write_tsv(
            root / "v_scores.tsv",
            list(next(iter(score_rows.values()), {"database": ""}).keys()),
            sorted(score_rows.values(), key=lambda row: (row["database"], row["accession"])),
        ),
    }
    return counts


def _manifest(root: Path, counts: dict[str, int]) -> dict[str, object]:
    files = {}
    for name in (*AVG_SOURCE_FILES, *NORMALIZED_FILES):
        path = root / name
        files[name] = {"path": str(path), "size": path.stat().st_size, "sha256": _sha256(path)}
    return {
        "name": AVG_NAME,
        "release": AVG_RELEASE,
        "source_repository": AVG_SOURCE_REPOSITORY,
        "source_urls": {name: f"{AVG_SOURCE_ROOT}/{name}" for name in AVG_SOURCE_FILES},
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "license": "GPL-3.0",
        "citation": "See the upstream release and phu documentation for attribution.",
        "files": files,
        "normalized_row_counts": counts,
    }


def prepare_avg_database(*, force_refresh: bool = False) -> dict[str, object]:
    root = _avg_root()
    root.mkdir(parents=True, exist_ok=True)
    for name in AVG_SOURCE_FILES:
        path = root / name
        if force_refresh or not path.exists():
            _download(f"{AVG_SOURCE_ROOT}/{name}", path)
    counts = _compile_tables(root)
    manifest = _manifest(root, counts)
    _atomic_write_text(root / "manifest.json", json.dumps(manifest, indent=2) + "\n")
    return manifest


def ensure_avg_database(*, force_refresh: bool = False) -> dict[str, object]:
    root = _avg_root()
    if force_refresh or not (root / "manifest.json").exists():
        return prepare_avg_database(force_refresh=force_refresh)
    try:
        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        for name in (*AVG_SOURCE_FILES, *NORMALIZED_FILES):
            if not (root / name).exists() or _sha256(root / name) != manifest["files"][name]["sha256"]:
                return prepare_avg_database(force_refresh=True)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return prepare_avg_database(force_refresh=True)
    return manifest


def refresh_avg_database() -> dict[str, object]:
    return prepare_avg_database(force_refresh=True)


def remove_avg_database() -> bool:
    root = _avg_root()
    if not root.exists():
        return False
    for path in sorted(root.iterdir(), reverse=True):
        if path.is_file() or path.is_symlink():
            path.unlink()
    root.rmdir()
    return True


def get_avg_database_status() -> dict[str, object]:
    root = _avg_root()
    manifest_path = root / "manifest.json"
    manifest: dict[str, object] = {}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            manifest = {}
    ready = bool(manifest) and all((root / name).exists() for name in NORMALIZED_FILES)
    counts = manifest.get("normalized_row_counts", {}) if isinstance(manifest, dict) else {}
    return {
        "name": AVG_NAME,
        "release": manifest.get("release"),
        "root": str(root),
        "manifest_path": str(manifest_path),
        "manifest_exists": manifest_path.exists(),
        "downloaded": all((root / name).exists() for name in AVG_SOURCE_FILES),
        "indexed": ready,
        "ready": ready,
        "normalized_row_counts": counts,
    }
