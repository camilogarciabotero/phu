from __future__ import annotations

import json
import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import click

from phu.core.progress import run_click_task
from phu.db.bundle import derived_dir, derived_files_lock, resolve_db_cached, write_atomically

logger = logging.getLogger(__name__)

KOFAM_ID_PREFIX = "K"
KOFAM_NAME = "kofam"
KOFAM_OFFSETS_SCHEMA_VERSION = 1


def is_kofam_id(token: str) -> bool:
    """Return True when token looks like a KO identifier (e.g. K00001)."""
    token = token.strip().upper()
    return len(token) == 6 and token.startswith(KOFAM_ID_PREFIX) and token[1:].isdigit()


def normalize_kofam_id(token: str) -> str:
    """Normalize KO identifier to upper-case K##### format."""
    token = token.strip().upper()
    if not is_kofam_id(token):
        raise ValueError(f"Invalid KO identifier: {token}")
    return token


def _kofam_root() -> Path:
    """KOfam in the phu-db bundle, installed on first use; checked once per process."""
    return resolve_db_cached("kofam")


def ensure_kofam_database() -> dict[str, str]:
    """Locate KOfam in the phu-db bundle, installing the bundle on first use."""
    return {
        "name": KOFAM_NAME,
        "hmm_path": str(_kofam_hmm_path()),
        "ko_list_path": str(_kofam_ko_list_path()),
    }


def _kofam_ko_list_path() -> Path:
    return _kofam_root() / "ko_list.tsv"


def _kofam_hmm_path() -> Path:
    return _kofam_root() / "profiles.hmm"


def _kofam_models_dir() -> Path:
    return derived_dir(_kofam_root()) / "models"


def _kofam_metadata_index_path() -> Path:
    return derived_dir(_kofam_root()) / "ko_metadata.json"


def _kofam_offsets_index_path() -> Path:
    return derived_dir(_kofam_root()) / "offsets.json"


def _write_manifest_atomically(
    manifest_path: Path, metadata: dict[str, object]
) -> None:
    """Write manifest JSON atomically to avoid corruption on interruption."""
    write_atomically(manifest_path, json.dumps(metadata, indent=2).encode())


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


@dataclass(frozen=True)
class KOFamMetadata:
    ko_id: str
    threshold: Optional[float]
    score_type: str
    profile_type: str
    definition: str


def _parse_float(token: str) -> Optional[float]:
    token = token.strip()
    if token in {"", "-", "NA", "N/A"}:
        return None
    return float(token)


def _parse_ko_list(ko_list_path: Path) -> dict[str, KOFamMetadata]:
    """Parse ko_list into KO metadata map."""
    metadata: dict[str, KOFamMetadata] = {}

    with ko_list_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 12:
                continue

            ko_id = parts[0].strip().upper()
            if not is_kofam_id(ko_id):
                continue

            threshold = _parse_float(parts[1])
            score_type = parts[2].strip().lower()
            profile_type = parts[3].strip().lower()
            definition = parts[11].strip()

            metadata[ko_id] = KOFamMetadata(
                ko_id=ko_id,
                threshold=threshold,
                score_type=score_type if score_type in {"full", "domain"} else "full",
                profile_type=profile_type,
                definition=definition,
            )

    return metadata


def _build_kofam_metadata_index(ko_list_path: Path) -> dict[str, KOFamMetadata]:
    metadata_map = _parse_ko_list(ko_list_path)
    serializable = {
        ko_id: {
            "ko_id": item.ko_id,
            "threshold": item.threshold,
            "score_type": item.score_type,
            "profile_type": item.profile_type,
            "definition": item.definition,
        }
        for ko_id, item in metadata_map.items()
    }

    _write_manifest_atomically(_kofam_metadata_index_path(), {"metadata": serializable})
    return metadata_map


def _load_kofam_metadata() -> dict[str, KOFamMetadata]:
    index_path = _kofam_metadata_index_path()
    ko_list_path = _kofam_ko_list_path()

    if not index_path.exists():
        if not ko_list_path.exists():
            raise FileNotFoundError(f"KOFam metadata not found: {ko_list_path}")
        return _build_kofam_metadata_index(ko_list_path)

    data = _read_json(index_path)
    metadata_obj = data.get("metadata")
    if not isinstance(metadata_obj, dict):
        return _build_kofam_metadata_index(ko_list_path)

    parsed: dict[str, KOFamMetadata] = {}
    for ko_id, payload in metadata_obj.items():
        if not isinstance(payload, dict):
            continue
        if not is_kofam_id(ko_id):
            continue
        threshold_obj = payload.get("threshold")
        threshold = threshold_obj if isinstance(threshold_obj, (int, float)) else None
        parsed[ko_id] = KOFamMetadata(
            ko_id=ko_id,
            threshold=float(threshold) if threshold is not None else None,
            score_type=str(payload.get("score_type", "full")),
            profile_type=str(payload.get("profile_type", "")),
            definition=str(payload.get("definition", "")),
        )

    if not parsed:
        return _build_kofam_metadata_index(ko_list_path)

    return parsed


def _clear_sparse_models_cache() -> None:
    models_dir = _kofam_models_dir()
    if not models_dir.exists():
        return
    for model_file in models_dir.glob("*.hmm"):
        model_file.unlink()


def _is_offsets_index_valid(hmm_db_path: Path) -> bool:
    offsets_index = _kofam_offsets_index_path()
    if not offsets_index.exists() or not hmm_db_path.exists():
        return False

    try:
        data = _read_json(offsets_index)
    except (json.JSONDecodeError, OSError, ValueError):
        return False

    required_keys = {
        "schema_version",
        "source_path",
        "source_size",
        "source_mtime_ns",
        "offsets",
        "complete",
    }
    if not required_keys.issubset(data.keys()):
        return False

    offsets = data.get("offsets")
    if not isinstance(offsets, dict):
        return False

    stat = hmm_db_path.stat()
    return (
        data["schema_version"] == KOFAM_OFFSETS_SCHEMA_VERSION
        and data["source_path"] == str(hmm_db_path)
        and data["source_size"] == stat.st_size
        and data["source_mtime_ns"] == stat.st_mtime_ns
        and data["complete"] is True
    )


def _build_offsets_index(hmm_db_path: Path) -> None:
    """Index KO model byte offsets by identifier for fast sparse extraction."""
    if not hmm_db_path.exists():
        raise FileNotFoundError(f"KOFam HMM database not found: {hmm_db_path}")

    source_stat = hmm_db_path.stat()
    offsets: dict[str, list[int]] = {}

    def _scan_offsets(update_progress: Optional[Callable[[int], None]] = None) -> None:
        block_start = 0
        block_ko: Optional[str] = None
        pending_progress = 0
        progress_chunk = 1024 * 1024

        with hmm_db_path.open("rb") as handle:
            block_start = handle.tell()
            while True:
                raw_line = handle.readline()
                if not raw_line:
                    break

                if update_progress is not None:
                    pending_progress += len(raw_line)
                    if pending_progress >= progress_chunk:
                        update_progress(pending_progress)
                        pending_progress = 0

                if raw_line.startswith(b"ACC"):
                    parts = raw_line.split()
                    if len(parts) >= 2:
                        token = (
                            parts[1].decode("utf-8", errors="replace").strip().upper()
                        )
                        if is_kofam_id(token):
                            block_ko = token
                elif raw_line.startswith(b"NAME") and block_ko is None:
                    parts = raw_line.split()
                    if len(parts) >= 2:
                        token = (
                            parts[1].decode("utf-8", errors="replace").strip().upper()
                        )
                        if is_kofam_id(token):
                            block_ko = token

                if raw_line.strip() == b"//":
                    block_end = handle.tell()
                    if block_ko is not None and block_ko not in offsets:
                        offsets[block_ko] = [block_start, block_end]
                    block_start = handle.tell()
                    block_ko = None

        if update_progress is not None and pending_progress:
            update_progress(pending_progress)

    if source_stat.st_size > 0:
        with click.progressbar(
            length=source_stat.st_size,
            label="Indexing KOFam HMMs",
            show_pos=True,
            show_percent=True,
            show_eta=True,
        ) as bar:
            _scan_offsets(bar.update)
    else:
        run_click_task("Indexing KOFam HMMs", _scan_offsets)

    metadata: dict[str, object] = {
        "schema_version": KOFAM_OFFSETS_SCHEMA_VERSION,
        "name": KOFAM_NAME,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "source_path": str(hmm_db_path),
        "source_size": source_stat.st_size,
        "source_mtime_ns": source_stat.st_mtime_ns,
        "model_count": len(offsets),
        "offsets": offsets,
        "complete": True,
    }
    _write_manifest_atomically(_kofam_offsets_index_path(), metadata)


def extract_kofam_models(
    requested_ids: Iterable[str],
    output_dir: Path,
) -> tuple[list[Path], list[str]]:
    """Extract requested KO models into output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    hmm_db_path = _kofam_hmm_path()
    if not hmm_db_path.exists():
        raise FileNotFoundError(f"KOFam HMM database not found: {hmm_db_path}")

    # Jobs sharing the install do this one at a time; one that waited finds the
    # index valid and leaves the other's cache alone.
    if not _is_offsets_index_valid(hmm_db_path):
        with derived_files_lock(derived_dir(hmm_db_path.parent)):
            if not _is_offsets_index_valid(hmm_db_path):
                _clear_sparse_models_cache()
                _build_offsets_index(hmm_db_path)

    models_dir = _kofam_models_dir()
    models_dir.mkdir(parents=True, exist_ok=True)

    normalized_ids = list(
        dict.fromkeys(normalize_kofam_id(token) for token in requested_ids)
    )

    out_paths_by_id: dict[str, Path] = {}
    missing: list[str] = []

    needs_index_read = False
    for ko_id in normalized_ids:
        model_path = models_dir / f"{ko_id}.hmm"
        if not model_path.exists():
            needs_index_read = True
            continue
        out_path = output_dir / f"{ko_id}.hmm"
        out_path.write_bytes(model_path.read_bytes())
        out_paths_by_id[ko_id] = out_path

    offsets_map: dict[str, object] = {}
    if needs_index_read:
        index_data = _read_json(_kofam_offsets_index_path())
        offsets_obj = index_data.get("offsets")
        if isinstance(offsets_obj, dict):
            offsets_map = offsets_obj

        with hmm_db_path.open("rb") as handle:
            for ko_id in normalized_ids:
                if ko_id in out_paths_by_id:
                    continue

                span_obj = offsets_map.get(ko_id)
                if not isinstance(span_obj, list) or len(span_obj) != 2:
                    missing.append(ko_id)
                    continue

                start, end = span_obj
                if (
                    not isinstance(start, int)
                    or not isinstance(end, int)
                    or end <= start
                ):
                    missing.append(ko_id)
                    continue

                handle.seek(start)
                model_blob = handle.read(end - start)

                cache_path = models_dir / f"{ko_id}.hmm"
                write_atomically(cache_path, model_blob)

                out_path = output_dir / f"{ko_id}.hmm"
                out_path.write_bytes(model_blob)
                out_paths_by_id[ko_id] = out_path

    extracted_paths = [
        out_paths_by_id[ko_id] for ko_id in normalized_ids if ko_id in out_paths_by_id
    ]
    return extracted_paths, missing


def get_kofam_metadata(ko_id: str) -> Optional[KOFamMetadata]:
    """Get metadata for one KO identifier."""
    ko_id = normalize_kofam_id(ko_id)
    metadata = _load_kofam_metadata()
    return metadata.get(ko_id)


def get_kofam_metadata_map(ko_ids: Iterable[str]) -> dict[str, KOFamMetadata]:
    """Get metadata map for KO identifiers (missing entries are omitted)."""
    metadata = _load_kofam_metadata()
    selected: dict[str, KOFamMetadata] = {}
    for token in ko_ids:
        ko_id = normalize_kofam_id(token)
        if ko_id in metadata:
            selected[ko_id] = metadata[ko_id]
    return selected


def get_all_kofam_metadata() -> dict[str, KOFamMetadata]:
    """Return full KO metadata map loaded from the local metadata index."""
    return _load_kofam_metadata()
