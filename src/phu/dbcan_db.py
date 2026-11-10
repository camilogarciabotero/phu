from __future__ import annotations

import functools
import json
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .db_bundle import resolve_db

DBCAN_NAME = "dbcan"
DBCAN_SCHEMA_VERSION = 3
PUL_RULES_HEADER = ("pul_id", "substrate", "cazymes_predicted_dbcan")
DBCAN_FAMILY_RE = re.compile(r"^[A-Z]{2,}\d+(?:_\d+)?$")
DBCAN_PUL_RE = re.compile(r"^PUL\d+$")
DBCAN_CANONICAL_RE = re.compile(r"^(?:AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?$")


@dataclass(frozen=True)
class PULRule:
    pul_id: str
    substrate: str
    raw_rule: str
    families: tuple[str, ...]
    unresolved_tokens: tuple[str, ...]

    @property
    def resolved(self) -> bool:
        return bool(self.families) and not self.unresolved_tokens


def is_dbcan_id(token: str) -> bool:
    value = token.strip().upper()
    value = value.removesuffix(".HMM")
    return bool(DBCAN_FAMILY_RE.fullmatch(value))


def normalize_dbcan_id(token: str) -> str:
    value = token.strip().upper()
    value = value.removesuffix(".HMM")
    if not is_dbcan_id(value):
        raise ValueError(f"Invalid dbCAN family: {token}")
    return value


def is_canonical_dbcan_id(token: str) -> bool:
    """Return whether a model belongs to the canonical CAZyme classes."""
    value = token.strip().upper().removesuffix(".HMM")
    return bool(DBCAN_CANONICAL_RE.fullmatch(value))


def is_dbcan_pul_id(token: str) -> bool:
    return bool(DBCAN_PUL_RE.fullmatch(token.strip().upper()))


def normalize_dbcan_pul_id(token: str) -> str:
    value = token.strip().upper()
    if not is_dbcan_pul_id(value):
        raise ValueError(f"Invalid dbCAN PUL identifier: {token}")
    return value


@functools.lru_cache(maxsize=1)
def _dbcan_root() -> Path:
    """dbCAN in the phu-db bundle, installed on first use; checked once per process."""
    return resolve_db(DBCAN_NAME)


def _hmm_path() -> Path:
    return _dbcan_root() / "profiles.hmm"


def _pul_rules_path() -> Path:
    return _dbcan_root() / "pul_rules.tsv"


def _offsets_path() -> Path:
    return _dbcan_root() / "family_offsets.json"


def _atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=path.parent) as handle:
        temporary_path = Path(handle.name)
        handle.write(content)
    temporary_path.replace(path)


def _models_dir() -> Path:
    return _dbcan_root() / "models"


def _clear_sparse_models_cache() -> None:
    models = _models_dir()
    if not models.exists():
        return
    for path in models.glob("*.hmm"):
        path.unlink()


def _model_from_name(name: str) -> str | None:
    value = name.strip().split()[0].upper().removesuffix(".HMM")
    return value if value and re.fullmatch(r"[A-Z][A-Z0-9_]*", value) else None


def build_family_offsets(hmm_path: Path) -> dict[str, list[int]]:
    offsets: dict[str, list[int]] = {}
    block_start = 0
    family: str | None = None
    with hmm_path.open("rb") as handle:
        while raw_line := handle.readline():
            if raw_line.startswith(b"NAME"):
                parts = raw_line.decode("utf-8", errors="replace").split()
                family = _model_from_name(parts[1]) if len(parts) > 1 else None
            if raw_line.strip() == b"//":
                block_end = handle.tell()
                if family and family not in offsets:
                    offsets[family] = [block_start, block_end]
                block_start = block_end
                family = None
    return offsets


def _write_json(path: Path, value: object) -> None:
    _atomic_write(path, (json.dumps(value, indent=2) + "\n").encode())


def parse_pul_rules(path: Path, family_ids: set[str]) -> dict[str, PULRule]:
    """Read the bundle's pul_rules.tsv; ``family_ids`` are the canonical models.

    A rule's tokens are split on ``,`` and ``|``. Those naming a canonical
    model become its families; the rest stay as unresolved tokens. A PUL listed
    twice keeps its first position and its last rule, as the workbook reader did.
    """
    rules: dict[str, PULRule] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != PUL_RULES_HEADER:
            raise ValueError(f"unexpected PUL rules header in {Path(path).name}: {header}")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            fields += [""] * (len(PUL_RULES_HEADER) - len(fields))
            pul_id = fields[0].strip().upper()
            if not is_dbcan_pul_id(pul_id):
                continue
            raw_rule = fields[2].strip()
            tokens = tuple(
                token.strip().upper() for token in re.split(r"[,|]", raw_rule) if token.strip()
            )
            families = tuple(dict.fromkeys(token for token in tokens if token in family_ids))
            unresolved = tuple(dict.fromkeys(token for token in tokens if token not in family_ids))
            rules[pul_id] = PULRule(pul_id, fields[1].strip(), raw_rule, families, unresolved)
    return rules


def _load_offsets() -> dict[str, list[int]]:
    """Model name -> byte span in profiles.hmm, from an index kept beside it.

    The index is rebuilt, and the sparse model cache cleared, whenever it does
    not describe the profiles.hmm that is there now.
    """
    hmm_path = _hmm_path()
    stat = hmm_path.stat()
    source = {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    try:
        index = json.loads(_offsets_path().read_text())
        if (
            index["schema_version"] == DBCAN_SCHEMA_VERSION
            and index["source"] == source
            and isinstance(index["offsets"], dict)
        ):
            return index["offsets"]
    except (OSError, KeyError, TypeError, ValueError):
        pass
    _clear_sparse_models_cache()
    offsets = build_family_offsets(hmm_path)
    index = {"schema_version": DBCAN_SCHEMA_VERSION, "source": source, "offsets": offsets}
    _atomic_write(_offsets_path(), (json.dumps(index, indent=2) + "\n").encode())
    return offsets


def ensure_dbcan_database() -> dict[str, str]:
    """Locate dbCAN in the phu-db bundle, installing the bundle on first use."""
    return {
        "name": DBCAN_NAME,
        "hmm_path": str(_hmm_path()),
        "pul_rules_path": str(_pul_rules_path()),
    }


def extract_dbcan_models(requested_ids: list[str], output_dir: Path) -> tuple[list[Path], list[str]]:
    offsets = _load_offsets()
    output_dir.mkdir(parents=True, exist_ok=True)
    models = _models_dir()
    models.mkdir(parents=True, exist_ok=True)
    normalized = list(dict.fromkeys(normalize_dbcan_id(value) for value in requested_ids))
    paths: list[Path] = []
    missing: list[str] = []
    with _hmm_path().open("rb") as handle:
        for family in normalized:
            cache_path = models / f"{family}.hmm"
            if not cache_path.exists():
                span = offsets.get(family)
                if not isinstance(span, list) or len(span) != 2:
                    missing.append(family)
                    continue
                handle.seek(span[0])
                cache_path.write_bytes(handle.read(span[1] - span[0]))
            output_path = output_dir / cache_path.name
            output_path.write_bytes(cache_path.read_bytes())
            paths.append(output_path)
    return paths, missing


def get_dbcan_model_inventory() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return indexed canonical models and excluded ancillary model names."""
    models = tuple(_load_offsets())
    canonical = tuple(model for model in models if is_canonical_dbcan_id(model))
    excluded = tuple(model for model in models if not is_canonical_dbcan_id(model))
    return canonical, excluded


def get_dbcan_pul_rules() -> dict[str, PULRule]:
    """Every PUL rule in the bundle, resolved against the canonical models."""
    canonical, _ = get_dbcan_model_inventory()
    return parse_pul_rules(_pul_rules_path(), set(canonical))


def get_dbcan_pul(pul_id: str) -> PULRule | None:
    return get_dbcan_pul_rules().get(normalize_dbcan_pul_id(pul_id))
