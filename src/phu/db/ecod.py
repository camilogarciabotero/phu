"""The ECOD database of the phu-db bundle: one HMMER model per F-group.

ECOD (Evolutionary Classification of protein Domains) groups domains into X
(possible homology), H (homology), T (topology) and F (family) groups, with
dotted ids: F-group ``2004.1.1.5`` lies in T-group ``2004.1.1``, H-group
``2004.1`` and X-group ``2004``. The bundle carries one model per F-group, built
from ECOD's F70 representatives, and ``families.tsv`` with each family's names.

Users name ECOD groups as ``ecod:<id>`` at any level; a T, H or X id stands for
every F-group inside it. Model names inside ``profiles.hmm`` are bare F-group
ids; extracted models are renamed to ``ecod:<f_id>`` so outputs carry the same
ids users type.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from phu.db.bundle import resolve_db
from phu.db.phrog import _offsets

ECOD_NAME = "ECOD"
ECOD_PREFIX = "ecod:"
FAMILIES_FILE = "families.tsv"
LEVELS = {1: "X", 2: "H", 3: "T", 4: "F"}
_ID = re.compile(r"\d+(?:\.\d+){0,3}")
_FAMILY_COLUMNS = (
    "f_id", "f_name", "t_id", "t_name", "h_id", "h_name", "x_id", "x_name",
    "architecture", "pfam_acc", "n_domains",
)  # fmt: skip


@dataclass(frozen=True)
class EcodFamily:
    f_id: str
    f_name: str
    t_id: str
    t_name: str
    h_id: str
    h_name: str
    x_id: str
    x_name: str
    architecture: str
    pfam_acc: str
    n_domains: int

    @property
    def ecod_id(self) -> str:
        return f"{ECOD_PREFIX}{self.f_id}"


def is_ecod_id(token: str) -> bool:
    """True for ``ecod:2004``, ``ecod:2004.1``, ``ecod:2004.1.1`` or ``ecod:2004.1.1.5``."""
    token = token.strip()
    return token[: len(ECOD_PREFIX)].lower() == ECOD_PREFIX and bool(
        _ID.fullmatch(token[len(ECOD_PREFIX) :])
    )


def normalize_ecod_id(token: str) -> str:
    """``ECOD:02004.1`` -> ``ecod:2004.1``."""
    if not is_ecod_id(token):
        raise ValueError(
            f"Invalid ECOD identifier: {token} (expected ecod:<X>[.<H>[.<T>[.<F>]]])"
        )
    parts = token.strip()[len(ECOD_PREFIX) :].split(".")
    return ECOD_PREFIX + ".".join(str(int(part)) for part in parts)


def ecod_level(ecod_id: str) -> str:
    """The ECOD level of a normalised id: X, H, T or F."""
    return LEVELS[ecod_id[len(ECOD_PREFIX) :].count(".") + 1]


def ensure_ecod_database() -> dict[str, str]:
    """Locate ECOD in the phu-db bundle, installing the bundle on first use."""
    root = resolve_db("ecod")
    return {
        "name": ECOD_NAME,
        "hmm_path": str(root / "profiles.hmm"),
        "families_path": str(root / FAMILIES_FILE),
    }


def read_ecod_families(path: Path) -> dict[str, EcodFamily]:
    """``families.tsv``, keyed by F-group id."""
    families: dict[str, EcodFamily] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = set(_FAMILY_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"ECOD families table is missing columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            values = {column: (row[column] or "").strip() for column in _FAMILY_COLUMNS}
            values["n_domains"] = int(values["n_domains"] or 0)
            families[values["f_id"]] = EcodFamily(**values)
    return families


@lru_cache(maxsize=2)
def _cached_families(path: str, size: int, mtime_ns: int) -> dict[str, EcodFamily]:
    return read_ecod_families(Path(path))


def get_ecod_families() -> dict[str, EcodFamily]:
    """Every F-group of the installed bundle, keyed by bare F-group id."""
    path = Path(ensure_ecod_database()["families_path"])
    stat = path.stat()
    return _cached_families(str(path), stat.st_size, stat.st_mtime_ns)


def expand_ecod_ids(
    tokens: Iterable[str], families: dict[str, EcodFamily] | None = None
) -> tuple[dict[str, list[str]], list[str]]:
    """Map each requested id to the F-groups it stands for.

    Returns ``{normalised id: [ecod:<f_id>, ...]}`` in request order (F-groups
    in id order) and the ids that match no F-group.
    """
    families = get_ecod_families() if families is None else families
    expanded: dict[str, list[str]] = {}
    missing: list[str] = []
    for ecod_id in dict.fromkeys(normalize_ecod_id(token) for token in tokens):
        bare = ecod_id[len(ECOD_PREFIX) :]
        if ecod_level(ecod_id) == "F":
            members = [bare] if bare in families else []
        else:
            members = [f_id for f_id in families if f_id.startswith(bare + ".")]
        if not members:
            missing.append(ecod_id)
            continue
        expanded[ecod_id] = [
            ECOD_PREFIX + f_id
            for f_id in sorted(members, key=lambda f: tuple(int(p) for p in f.split(".")))
        ]
    return expanded, missing


def _scan_offsets(hmm_path: Path) -> dict[str, list[int]]:
    """Byte range of every model, keyed by bare F-group id."""
    offsets: dict[str, list[int]] = {}
    with hmm_path.open("rb") as handle:
        start = position = 0
        name: str | None = None
        for line in handle:
            position += len(line)
            if line.startswith(b"NAME"):
                name = line.split(None, 1)[1].decode(errors="replace").strip()
            elif line.strip() == b"//":
                if name and name not in offsets:
                    offsets[name] = [start, position]
                start = position
                name = None
    return offsets


def extract_ecod_models(
    ecod_ids: Iterable[str], output_dir: Path
) -> tuple[list[Path], list[str]]:
    """Write each F-group model to ``output_dir/ecod_<f_id>.hmm``, named ``ecod:<f_id>``.

    ``ecod_ids`` are F-level ids (``ecod:<f_id>``, as from ``expand_ecod_ids``).
    Returns the written paths in request order and the ids without a model.
    """
    hmm_path = Path(ensure_ecod_database()["hmm_path"])
    offsets = _offsets(hmm_path, _scan_offsets)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    missing: list[str] = []
    with hmm_path.open("rb") as handle:
        for ecod_id in dict.fromkeys(normalize_ecod_id(token) for token in ecod_ids):
            bare = ecod_id[len(ECOD_PREFIX) :]
            span = offsets.get(bare)
            if span is None:
                missing.append(ecod_id)
                continue
            handle.seek(span[0])
            model = handle.read(span[1] - span[0])
            model = re.sub(
                rb"^NAME( +)" + re.escape(bare.encode()) + rb"$",
                b"NAME\\1" + ecod_id.encode(),
                model,
                count=1,
                flags=re.MULTILINE,
            )
            path = output_dir / f"ecod_{bare}.hmm"
            path.write_bytes(model)
            paths.append(path)
    return paths, missing
