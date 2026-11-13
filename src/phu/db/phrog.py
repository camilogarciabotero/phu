"""The PHROG database of the phu-db bundle: HMMER models and annotations.

PHROGs are published as alignments and HH-suite profiles. The bundle carries
HMMER models built from the official alignments, named ``phrog_<number>``, and
the upstream annotation table, which lists every PHROG including the many
without a model.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from phu.db.bundle import resolve_db

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
