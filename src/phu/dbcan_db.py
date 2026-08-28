from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
from xml.etree import ElementTree

DBCAN_NAME = "dbcan"
DBCAN_SCHEMA_VERSION = 2
DBCAN_HMM_URL = "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/V15/dbCAN-HMMdb-V15.txt"
DBCAN_PUL_URL = "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/dbCAN-PUL.xlsx"
DBCAN_SUBSTRATE_URL = "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/dbCAN-PUL_12-12-2023.txt"
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


def _db_root() -> Path:
    if db_env := os.environ.get("PHU_DB_FOLDER"):
        return Path(db_env)
    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


def _dbcan_root() -> Path:
    return _db_root() / DBCAN_NAME


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=path.parent) as handle:
        temporary_path = Path(handle.name)
        handle.write(content)
    temporary_path.replace(path)


def _download(url: str, destination: Path) -> None:
    with urlopen(url) as response:
        _atomic_write(destination, response.read())


def _hmm_path() -> Path:
    return _dbcan_root() / "dbCAN-HMMdb-V15.txt"


def _pul_path() -> Path:
    return _dbcan_root() / "dbCAN-PUL.xlsx"


def _manifest_path() -> Path:
    return _dbcan_root() / "manifest.json"


def _offsets_path() -> Path:
    return _dbcan_root() / "family_offsets.json"


def _pul_index_path() -> Path:
    return _dbcan_root() / "pul_rules.json"


def _models_dir() -> Path:
    return _dbcan_root() / "models"


def _clear_sparse_models_cache() -> None:
    models = _models_dir()
    if not models.exists():
        return
    for path in models.glob("*.hmm"):
        path.unlink()


def _family_from_name(name: str) -> str | None:
    value = name.strip().split()[0].removesuffix(".hmm")
    return normalize_dbcan_id(value) if is_dbcan_id(value) else None


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


def _xlsx_values(path: Path) -> list[list[str]]:
    namespace = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in root.findall(f"{{{namespace}}}si")]
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in rels
        }
        sheet = workbook.find(f".//{{{namespace}}}sheet")
        if sheet is None:
            raise ValueError("PUL workbook has no worksheets")
        target = rel_map[sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]]
        sheet_path = "xl/" + target.removeprefix("/") if not target.startswith("xl/") else target
        root = ElementTree.fromstring(archive.read(sheet_path))
        rows: list[list[str]] = []
        for row in root.findall(f".//{{{namespace}}}row"):
            values: dict[int, str] = {}
            for cell in row.findall(f"{{{namespace}}}c"):
                ref = cell.attrib.get("r", "")
                match = re.search(r"[A-Z]+", ref)
                if not match:
                    continue
                column = 0
                for character in match.group():
                    column = column * 26 + ord(character) - 64
                column -= 1
                value = cell.find(f"{{{namespace}}}v")
                text = "" if value is None else value.text or ""
                if cell.attrib.get("t") == "s" and text:
                    text = shared[int(text)]
                elif cell.attrib.get("t") == "inlineStr":
                    inline = cell.find(f"{{{namespace}}}is")
                    text = "" if inline is None else "".join(inline.itertext())
                values[column] = text
            if values:
                rows.append([values.get(i, "") for i in range(max(values) + 1)])
        return rows


def parse_pul_workbook(path: Path, family_ids: set[str]) -> dict[str, PULRule]:
    rows = _xlsx_values(path)
    if not rows:
        raise ValueError("PUL workbook is empty")
    headers = {value.strip().lower(): index for index, value in enumerate(rows[0])}
    pul_column = next(
        (
            headers[name]
            for name in ("id", "pul", "pul_id", "pul id")
            if name in headers
        ),
        None,
    )
    rule_column = headers.get("cazymes_predicted_dbcan")
    if pul_column is None:
        raise ValueError("PUL workbook is missing columns: pul")
    if rule_column is None:
        raise ValueError("PUL workbook is missing columns: cazymes_predicted_dbcan")
    substrate_column = next((index for name, index in headers.items() if "substrate" in name), None)
    rules: dict[str, PULRule] = {}
    for row in rows[1:]:
        if pul_column >= len(row):
            continue
        pul_id = row[pul_column].strip().upper()
        if not is_dbcan_pul_id(pul_id):
            continue
        raw_rule = row[rule_column].strip() if rule_column < len(row) else ""
        tokens = tuple(
            token.strip().upper()
            for token in re.split(r"[,|]", raw_rule)
            if token.strip()
        )
        families = tuple(dict.fromkeys(token for token in tokens if token in family_ids))
        unresolved = tuple(dict.fromkeys(token for token in tokens if token not in family_ids))
        substrate = row[substrate_column].strip() if substrate_column is not None and substrate_column < len(row) else ""
        rules[pul_id] = PULRule(pul_id, substrate, raw_rule, families, unresolved)
    return rules


def _write_json(path: Path, value: object) -> None:
    _atomic_write(path, (json.dumps(value, indent=2) + "\n").encode())


def prepare_dbcan_database(*, force_refresh: bool = False) -> dict[str, object]:
    root = _dbcan_root()
    root.mkdir(parents=True, exist_ok=True)
    hmm_path, pul_path = _hmm_path(), _pul_path()
    if force_refresh:
        _clear_sparse_models_cache()
    if force_refresh or not hmm_path.exists():
        _download(DBCAN_HMM_URL, hmm_path)
    if force_refresh or not pul_path.exists():
        _download(DBCAN_PUL_URL, pul_path)
    offsets = build_family_offsets(hmm_path)
    _write_json(_offsets_path(), {"schema_version": DBCAN_SCHEMA_VERSION, "offsets": offsets})
    canonical_ids = {model for model in offsets if is_canonical_dbcan_id(model)}
    rules = parse_pul_workbook(pul_path, canonical_ids)
    _write_json(_pul_index_path(), {"schema_version": DBCAN_SCHEMA_VERSION, "rules": {key: value.__dict__ for key, value in rules.items()}})
    manifest = {
        "name": DBCAN_NAME,
        "schema_version": DBCAN_SCHEMA_VERSION,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "source_urls": {"hmm": DBCAN_HMM_URL, "pul": DBCAN_PUL_URL, "substrate_metadata": DBCAN_SUBSTRATE_URL},
        "files": {
            "hmm": {"path": str(hmm_path), "sha256": _sha256(hmm_path)},
            "pul": {"path": str(pul_path), "sha256": _sha256(pul_path)},
            "family_offsets": {"path": str(_offsets_path()), "sha256": _sha256(_offsets_path())},
            "pul_rules": {"path": str(_pul_index_path()), "sha256": _sha256(_pul_index_path())},
        },
        "family_count": len(offsets),
        "pul_count": len(rules),
    }
    _write_json(_manifest_path(), manifest)
    return manifest


def ensure_dbcan_database(*, force_refresh: bool = False) -> dict[str, object]:
    if force_refresh or not _manifest_path().exists():
        return prepare_dbcan_database(force_refresh=force_refresh)
    try:
        manifest = json.loads(_manifest_path().read_text())
        if manifest.get("schema_version") != DBCAN_SCHEMA_VERSION:
            _clear_sparse_models_cache()
            return prepare_dbcan_database(force_refresh=False)
        for name in ("hmm", "pul", "family_offsets", "pul_rules"):
            item = manifest["files"][name]
            path = Path(item["path"])
            if not path.exists() or _sha256(path) != item["sha256"]:
                _clear_sparse_models_cache()
                return prepare_dbcan_database(force_refresh=True)
        return manifest
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return prepare_dbcan_database(force_refresh=True)


def refresh_dbcan_database() -> dict[str, object]:
    return prepare_dbcan_database(force_refresh=True)


def remove_dbcan_database() -> bool:
    root = _dbcan_root()
    if not root.exists():
        return False
    shutil.rmtree(root)
    return True


def get_dbcan_database_status() -> dict[str, object]:
    root = _dbcan_root()
    manifest = _manifest_path()
    data: dict[str, object] = {}
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text())
        except (OSError, json.JSONDecodeError):
            data = {}
    offsets = _offsets_path()
    try:
        family_count = len(json.loads(offsets.read_text()).get("offsets", {}))
    except (OSError, json.JSONDecodeError, AttributeError):
        family_count = 0
    models = _models_dir()
    return {
        "name": DBCAN_NAME,
        "root": str(root),
        "hmm_path": str(_hmm_path()),
        "pul_path": str(_pul_path()),
        "offsets_path": str(offsets),
        "pul_index_path": str(_pul_index_path()),
        "manifest_path": str(manifest),
        "downloaded": _hmm_path().exists() and _pul_path().exists(),
        "indexed": offsets.exists() and _pul_index_path().exists(),
        "manifest_exists": manifest.exists(),
        "model_count": family_count,
        "pul_count": data.get("pul_count", 0),
        "sparse_cached_count": len(list(models.glob("*.hmm"))) if models.exists() else 0,
    }


def extract_dbcan_models(requested_ids: list[str], output_dir: Path) -> tuple[list[Path], list[str]]:
    ensure_dbcan_database()
    output_dir.mkdir(parents=True, exist_ok=True)
    offsets = json.loads(_offsets_path().read_text())["offsets"]
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


def get_dbcan_pul(pul_id: str) -> PULRule | None:
    ensure_dbcan_database()
    payload = json.loads(_pul_index_path().read_text()).get("rules", {}).get(normalize_dbcan_pul_id(pul_id))
    return None if payload is None else PULRule(**payload)


def get_dbcan_pul_rules() -> dict[str, PULRule]:
    """Load every PUL rule from the installed dbCAN rule index."""
    ensure_dbcan_database()
    payload = json.loads(_pul_index_path().read_text()).get("rules", {})
    return {
        pul_id: PULRule(**rule)
        for pul_id, rule in payload.items()
        if isinstance(rule, dict)
    }


def get_dbcan_model_inventory() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return indexed canonical models and excluded ancillary model names."""
    ensure_dbcan_database()
    offsets = json.loads(_offsets_path().read_text()).get("offsets", {})
    models = tuple(offsets) if isinstance(offsets, dict) else ()
    canonical = tuple(model for model in models if is_canonical_dbcan_id(model))
    excluded = tuple(model for model in models if not is_canonical_dbcan_id(model))
    return canonical, excluded