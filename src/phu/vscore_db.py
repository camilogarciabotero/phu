from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

V_SCORE_URL = (
    "https://raw.githubusercontent.com/AnantharamanLab/V-Score-Search/"
    "refs/heads/main/Software/VScoreDataNormalized.csv"
)
V_SCORE_NAME = "vscore"
V_SCORE_REQUIRED_COLUMNS = {
    "Accession",
    "Protein Function",
    "V-Score",
    "Log10[Hit Number]",
    "Database Origin",
}
V_SCORE_OPTIONAL_SCORE_COLUMNS = {
    "Normalized VL-score",
    "Normalized V-score",
}


@dataclass(frozen=True)
class VScoreRecord:
    accession: str
    protein_function: str
    v_score: float
    vl_score: float
    log10_hit_number: float
    database_origin: str


def _db_root() -> Path:
    if db_env := os.environ.get("PHU_DB_FOLDER"):
        return Path(db_env)
    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


def _vscore_root() -> Path:
    return _db_root() / V_SCORE_NAME


def _vscore_csv_path() -> Path:
    return _vscore_root() / "VScoreData.csv"


def _manifest_path() -> Path:
    return _vscore_root() / "manifest.json"


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as tmp:
        temporary_path = Path(tmp.name)
        with urlopen(url) as response:
            shutil.copyfileobj(response, tmp, length=1024 * 1024)
    temporary_path.replace(destination)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_vl_score_column(fieldnames: list[str] | None) -> str:
    if not fieldnames:
        raise ValueError("Incompatible V-score CSV: no header row found")
    for candidate in ("Normalized VL-score", "Normalized V-score"):
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        "Incompatible V-score CSV: missing required VL-score column "
        "(expected 'Normalized VL-score' or legacy 'Normalized V-score')"
    )


def parse_vscore_csv(path: Path) -> dict[str, VScoreRecord]:
    """Parse and validate a V-score CSV, keyed by normalized accession."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = V_SCORE_REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                "Incompatible V-score CSV: missing required columns: "
                + ", ".join(sorted(missing))
            )

        vl_score_column = _resolve_vl_score_column(reader.fieldnames)

        records: dict[str, VScoreRecord] = {}
        for row_number, row in enumerate(reader, start=2):
            accession = (row["Accession"] or "").strip().upper()
            if not accession:
                raise ValueError(f"Invalid V-score CSV row {row_number}: empty Accession")
            try:
                record = VScoreRecord(
                    accession=accession,
                    protein_function=(row["Protein Function"] or "").strip(),
                    v_score=float(row["V-Score"]),
                    vl_score=float(row[vl_score_column]),
                    log10_hit_number=float(row["Log10[Hit Number]"]),
                    database_origin=(row["Database Origin"] or "").strip(),
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid V-score CSV row {row_number}: numeric fields are malformed"
                ) from exc
            records[accession] = record
    return records


def ensure_vscore_database(force_refresh: bool = False) -> dict[str, str]:
    """Ensure the pinned V-score CSV is present and schema-compatible locally."""
    csv_path = _vscore_csv_path()
    if force_refresh or not csv_path.exists():
        _download(V_SCORE_URL, csv_path)
    else:
        try:
            parse_vscore_csv(csv_path)
        except ValueError:
            _download(V_SCORE_URL, csv_path)
            parse_vscore_csv(csv_path)
    records = parse_vscore_csv(csv_path)
    manifest = {
        "name": V_SCORE_NAME,
        "source_url": V_SCORE_URL,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "csv_path": str(csv_path),
        "sha256": _sha256(csv_path),
        "record_count": len(records),
    }
    _manifest_path().write_text(json.dumps(manifest, indent=2) + "\n")
    return {key: str(value) for key, value in manifest.items()}


def get_vscore_map() -> dict[str, VScoreRecord]:
    """Load the local V-score table, refreshing stale/invalid files automatically."""
    csv_path = _vscore_csv_path()
    if not csv_path.exists():
        raise FileNotFoundError(
            f"V-score database not found: {csv_path}. Run `phu dbs prepare vscore` first."
        )
    try:
        return parse_vscore_csv(csv_path)
    except ValueError:
        ensure_vscore_database(force_refresh=True)
        return parse_vscore_csv(csv_path)


def get_vscore_database_status() -> dict[str, object]:
    csv_path = _vscore_csv_path()
    manifest = _manifest_path()
    record_count = 0
    if csv_path.exists():
        try:
            record_count = len(parse_vscore_csv(csv_path))
        except ValueError:
            record_count = 0
    return {
        "name": V_SCORE_NAME,
        "root": str(_vscore_root()),
        "csv_path": str(csv_path),
        "manifest_path": str(manifest),
        "downloaded": csv_path.exists(),
        "manifest_exists": manifest.exists(),
        "record_count": record_count,
    }


def refresh_vscore_database() -> dict[str, object]:
    ensure_vscore_database(force_refresh=False)
    status = get_vscore_database_status()
    status["refreshed"] = True
    return status


def remove_vscore_database() -> bool:
    root = _vscore_root()
    if not root.exists():
        return False
    shutil.rmtree(root)
    return True