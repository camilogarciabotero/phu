"""The PHROG database of the phu-db bundle: HMMER models and annotations.

PHROGs are published as alignments and HH-suite profiles. The bundle carries
HMMER models built from the official alignments, named ``phrog_<number>``, and
the upstream annotation table, which lists every PHROG including the many
without a model.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from phu.db.bundle import derived_dir, derived_files_lock, resolve_db, write_atomically

PHROG_NAME = "PHROG"
PHROG_ID_PREFIX = "phrog_"
ANNOTATIONS_FILE = "annotations.tsv"
_UNANNOTATED = {"", "NA"}


@dataclass(frozen=True)
class PhrogAnnotation:
    phrog_id: str
    annotation: str
    category: str


def is_phrog_id(token: str) -> bool:
    """True for ``phrog_216`` in any letter case."""
    token = token.strip().lower()
    return token.startswith(PHROG_ID_PREFIX) and token[len(PHROG_ID_PREFIX) :].isdigit()


def normalize_phrog_id(token: str) -> str:
    """Normalise a PHROG identifier to ``phrog_<number>`` without leading zeros."""
    if not is_phrog_id(token):
        raise ValueError(f"Invalid PHROG identifier: {token}")
    return f"{PHROG_ID_PREFIX}{int(token.strip()[len(PHROG_ID_PREFIX):])}"


def ensure_phrog_database() -> dict[str, str]:
    """Locate PHROG in the phu-db bundle, installing the bundle on first use."""
    root = resolve_db("phrog")
    return {
        "name": PHROG_NAME,
        "hmm_path": str(root / "profiles.hmm"),
        "annotations_path": str(root / ANNOTATIONS_FILE),
    }


def read_phrog_annotations(path: Path) -> dict[str, PhrogAnnotation]:
    """Read the annotation table; an ``NA`` annotation becomes an empty string."""
    annotations: dict[str, PhrogAnnotation] = {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"phrog", "annot", "category"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"PHROG annotation table is missing columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            number = (row["phrog"] or "").strip()
            if not number.isdigit():
                raise ValueError(f"Invalid PHROG number in annotation table: {number!r}")
            phrog_id = f"{PHROG_ID_PREFIX}{int(number)}"
            annotation = (row["annot"] or "").strip()
            annotations[phrog_id] = PhrogAnnotation(
                phrog_id=phrog_id,
                annotation="" if annotation in _UNANNOTATED else annotation,
                category=(row["category"] or "").strip(),
            )
    return annotations


@lru_cache(maxsize=2)
def _cached_annotations(path: str, size: int, mtime_ns: int) -> dict[str, PhrogAnnotation]:
    return read_phrog_annotations(Path(path))


def get_all_phrog_annotations() -> dict[str, PhrogAnnotation]:
    """Annotations of every PHROG in the installed bundle, keyed by identifier."""
    path = Path(ensure_phrog_database()["annotations_path"])
    stat = path.stat()
    return _cached_annotations(str(path), stat.st_size, stat.st_mtime_ns)


# Extraction of single models, for `phu screen phrog_<n>`.

OFFSETS_NAME = "offsets.json"
_OFFSETS_SCHEMA_VERSION = 1


def _offsets_path(hmm_path: Path) -> Path:
    return derived_dir(hmm_path.parent) / OFFSETS_NAME


def _source_stamp(hmm_path: Path) -> dict[str, object]:
    stat = hmm_path.stat()
    return {
        "schema_version": _OFFSETS_SCHEMA_VERSION,
        "source_path": str(hmm_path),
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
    }


def _read_offsets(hmm_path: Path) -> dict[str, list[int]] | None:
    """The stored byte range of every model, or None if absent or out of date."""
    try:
        data = json.loads(_offsets_path(hmm_path).read_text())
    except (OSError, ValueError):
        return None
    stamp = _source_stamp(hmm_path)
    if not isinstance(data, dict) or any(data.get(key) != value for key, value in stamp.items()):
        return None
    offsets = data.get("offsets")
    return offsets if isinstance(offsets, dict) else None


def _scan_offsets(hmm_path: Path) -> dict[str, list[int]]:
    """Byte range of every model in a HMMER3 text file, keyed by PHROG identifier."""
    offsets: dict[str, list[int]] = {}
    with hmm_path.open("rb") as handle:
        start = 0
        name: str | None = None
        position = 0
        for line in handle:
            position += len(line)
            if line.startswith(b"NAME"):
                token = line.split(None, 1)[1].decode(errors="replace").strip()
                name = normalize_phrog_id(token) if is_phrog_id(token) else None
            elif line.strip() == b"//":
                if name is not None and name not in offsets:
                    offsets[name] = [start, position]
                start = position
                name = None
    return offsets


def _offsets(hmm_path: Path) -> dict[str, list[int]]:
    offsets = _read_offsets(hmm_path)
    if offsets is not None:
        return offsets
    place = derived_dir(hmm_path.parent)
    place.mkdir(parents=True, exist_ok=True)
    with derived_files_lock(place):
        # Another job may have indexed the file while this one waited.
        offsets = _read_offsets(hmm_path)
        if offsets is None:
            offsets = _scan_offsets(hmm_path)
            write_atomically(
                _offsets_path(hmm_path),
                json.dumps({**_source_stamp(hmm_path), "offsets": offsets}).encode(),
            )
    return offsets


def extract_phrog_models(
    requested_ids: Iterable[str], output_dir: Path
) -> tuple[list[Path], list[str]]:
    """Write each requested PHROG model to ``output_dir/<id>.hmm``.

    Returns the written paths, in request order, and the identifiers that have
    no model: most PHROGs above 38,880 are singletons without an alignment.
    """
    hmm_path = Path(ensure_phrog_database()["hmm_path"])
    offsets = _offsets(hmm_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    missing: list[str] = []
    with hmm_path.open("rb") as handle:
        for phrog_id in dict.fromkeys(normalize_phrog_id(token) for token in requested_ids):
            span = offsets.get(phrog_id)
            if span is None:
                missing.append(phrog_id)
                continue
            handle.seek(span[0])
            path = output_dir / f"{phrog_id}.hmm"
            path.write_bytes(handle.read(span[1] - span[0]))
            paths.append(path)
    return paths, missing
