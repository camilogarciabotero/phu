from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple
from urllib.request import urlopen

import click

from ._click import run_click_task

PFAM_HMM_GZ_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_ID_PREFIX = "PF"
PFAM_NAME = "Pfam-A"
PFAM_OFFSETS_SCHEMA_VERSION = 1

def is_pfam_id(token: str) -> bool:
    """Return True if token looks like a PFAM accession (with optional version)."""
    token = token.strip().upper()
    if not token.startswith(PFAM_ID_PREFIX):
        return False
    body = token[2:]
    if "." in body:
        acc, version = body.split(".", 1)
        return len(acc) == 5 and acc.isdigit() and version.isdigit() and len(version) > 0
    return len(body) == 5 and body.isdigit()


def normalize_pfam_id(token: str) -> str:
    """Normalize PFAM accession to unversioned upper-case format (e.g. PF00001)."""
    token = token.strip().upper()
    if not is_pfam_id(token):
        raise ValueError(f"Invalid PFAM accession: {token}")
    return token.split(".", 1)[0]


def _db_root() -> Path:
    """Resolve PHU DB root directory from env var or XDG-like defaults."""
    if db_env := os.environ.get("PHU_DB_FOLDER"):
        return Path(db_env)

    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


def _pfam_root() -> Path:
    return _db_root() / "pfam"


def _pfam_hmm_gz_path() -> Path:
    return _pfam_root() / "Pfam-A.hmm.gz"


def _pfam_hmm_path() -> Path:
    return _pfam_root() / "Pfam-A.hmm"


def _pfam_manifest_path() -> Path:
    return _pfam_root() / "manifest.json"


def _pfam_models_dir(hmm_db_path: Path) -> Path:
    return hmm_db_path.parent / "models"


def _pfam_offsets_index_path(hmm_db_path: Path) -> Path:
    return hmm_db_path.parent / "offsets.json"


def _stream_download_to_path(url: str, destination: Path) -> None:
    """Download URL content to destination atomically via a temporary file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as tmp:
        tmp_path = Path(tmp.name)
        with urlopen(url) as response:
            content_length = response.info().get("Content-Length")
            total_size = int(content_length) if content_length and content_length.isdigit() else None

            if total_size is not None and total_size > 0:
                with click.progressbar(
                    length=total_size,
                    label="Downloading Pfam-A",
                    show_pos=True,
                    show_percent=True,
                    show_eta=True,
                ) as bar:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        tmp.write(chunk)
                        bar.update(len(chunk))
            else:
                def _copy_stream() -> None:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        tmp.write(chunk)

                run_click_task("Downloading Pfam-A", _copy_stream)

    tmp_path.replace(destination)


def _decompress_gzip_to_path(gz_path: Path, out_path: Path) -> None:
    """Decompress gz_path into out_path atomically."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=out_path.parent) as tmp:
        tmp_path = Path(tmp.name)
        with gzip.open(gz_path, "rb") as src:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                tmp.write(chunk)
    tmp_path.replace(out_path)


def _sha256(path: Path) -> str:
    """Compute SHA256 for a file path."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _write_manifest_atomically(manifest_path: Path, metadata: Dict[str, object]) -> None:
    """Write manifest JSON atomically to avoid corruption on interruption."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=manifest_path.parent, suffix=".json"
    ) as tmp:
        tmp_path = Path(tmp.name)
        json.dump(metadata, tmp, indent=2)
    tmp_path.replace(manifest_path)


def _read_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def _scan_hmm_blocks(
    hmm_db_path: Path,
    *,
    label: str,
    on_block: Callable[[str | None, List[str]], None],
) -> None:
    """Scan an HMM database and invoke on_block(accession, lines) for each model block."""
    total_size = hmm_db_path.stat().st_size if hmm_db_path.exists() else None
    block_lines: List[str] = []
    block_acc: str | None = None

    def _flush_block() -> None:
        nonlocal block_lines, block_acc
        on_block(block_acc, block_lines)
        block_lines = []
        block_acc = None

    def _scan_stream(update_progress: Callable[[int], None] | None = None) -> None:
        nonlocal block_acc
        pending_progress = 0
        progress_chunk = 1024 * 1024
        with hmm_db_path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                block_lines.append(raw_line)
                if update_progress is not None:
                    pending_progress += len(raw_line)
                    if pending_progress >= progress_chunk:
                        update_progress(pending_progress)
                        pending_progress = 0

                if line.startswith("ACC"):
                    parts = line.split()
                    if len(parts) >= 2:
                        acc_raw = parts[1].strip()
                        try:
                            block_acc = normalize_pfam_id(acc_raw)
                        except ValueError:
                            block_acc = None

                if line.strip() == "//":
                    _flush_block()

        if update_progress is not None and pending_progress:
            update_progress(pending_progress)

    if total_size is not None and total_size > 0:
        with click.progressbar(
            length=total_size,
            label=label,
            show_pos=True,
            show_percent=True,
            show_eta=True,
        ) as bar:
            _scan_stream(bar.update)
    else:
        run_click_task(label, _scan_stream)


def _is_offsets_index_valid(hmm_db_path: Path) -> bool:
    offsets_index = _pfam_offsets_index_path(hmm_db_path)
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
        data["schema_version"] == PFAM_OFFSETS_SCHEMA_VERSION
        and data["source_path"] == str(hmm_db_path)
        and data["source_size"] == stat.st_size
        and data["source_mtime_ns"] == stat.st_mtime_ns
        and data["complete"] is True
    )


def _build_offsets_index(hmm_db_path: Path) -> None:
    """Index PFAM model byte offsets by accession for fast sparse extraction."""
    if not hmm_db_path.exists():
        raise FileNotFoundError(f"PFAM HMM database not found: {hmm_db_path}")

    source_stat = hmm_db_path.stat()
    offsets: Dict[str, List[int]] = {}

    def _scan_offsets(update_progress: Callable[[int], None] | None = None) -> None:
        block_start = 0
        block_acc: str | None = None
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
                        acc_raw = parts[1].decode("utf-8", errors="replace").strip()
                        try:
                            block_acc = normalize_pfam_id(acc_raw)
                        except ValueError:
                            block_acc = None

                if raw_line.strip() == b"//":
                    block_end = handle.tell()
                    if block_acc is not None and block_acc not in offsets:
                        offsets[block_acc] = [block_start, block_end]
                    block_start = handle.tell()
                    block_acc = None

        if update_progress is not None and pending_progress:
            update_progress(pending_progress)

    if source_stat.st_size > 0:
        with click.progressbar(
            length=source_stat.st_size,
            label="Indexing Pfam-A",
            show_pos=True,
            show_percent=True,
            show_eta=True,
        ) as bar:
            _scan_offsets(bar.update)
    else:
        run_click_task("Indexing Pfam-A", _scan_offsets)

    metadata: Dict[str, object] = {
        "schema_version": PFAM_OFFSETS_SCHEMA_VERSION,
        "name": PFAM_NAME,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "source_path": str(hmm_db_path),
        "source_size": source_stat.st_size,
        "source_mtime_ns": source_stat.st_mtime_ns,
        "model_count": len(offsets),
        "offsets": offsets,
        "complete": True,
    }
    _write_manifest_atomically(_pfam_offsets_index_path(hmm_db_path), metadata)


def _extract_from_split_cache(
    hmm_db_path: Path,
    requested_ids: Iterable[str],
    output_dir: Path,
) -> Tuple[List[Path], List[str]]:
    """Resolve requested PFAM IDs via offset index and sparse local model cache."""
    models_dir = _pfam_models_dir(hmm_db_path)
    models_dir.mkdir(parents=True, exist_ok=True)
    normalized_ordered = list(dict.fromkeys(normalize_pfam_id(token) for token in requested_ids))

    out_paths_by_id: Dict[str, Path] = {}
    missing: List[str] = []

    needs_index = False
    for pfam_id in normalized_ordered:
        model_path = models_dir / f"{pfam_id}.hmm"
        if not model_path.exists():
            needs_index = True
            continue
        out_path = output_dir / f"{pfam_id}.hmm"
        out_path.write_bytes(model_path.read_bytes())
        out_paths_by_id[pfam_id] = out_path

    offsets: Dict[str, object] = {}
    if needs_index:
        if not _is_offsets_index_valid(hmm_db_path):
            _build_offsets_index(hmm_db_path)

        index_data = _read_json(_pfam_offsets_index_path(hmm_db_path))
        offsets_obj = index_data.get("offsets")
        if isinstance(offsets_obj, dict):
            offsets = offsets_obj

        with hmm_db_path.open("rb") as handle:
            for pfam_id in normalized_ordered:
                if pfam_id in out_paths_by_id:
                    continue
                span_obj = offsets.get(pfam_id)
                if not isinstance(span_obj, list) or len(span_obj) != 2:
                    missing.append(pfam_id)
                    continue

                start, end = span_obj
                if not isinstance(start, int) or not isinstance(end, int) or end <= start:
                    missing.append(pfam_id)
                    continue

                handle.seek(start)
                model_blob = handle.read(end - start)
                cache_path = models_dir / f"{pfam_id}.hmm"
                cache_path.write_bytes(model_blob)

                out_path = output_dir / f"{pfam_id}.hmm"
                out_path.write_bytes(model_blob)
                out_paths_by_id[pfam_id] = out_path

    out_paths = [out_paths_by_id[pfam_id] for pfam_id in normalized_ordered if pfam_id in out_paths_by_id]
    return out_paths, missing


def ensure_pfam_database(force_refresh: bool = False) -> Dict[str, str]:
    """
    Ensure local Pfam-A database exists and return metadata.

    On first use (or force_refresh), downloads and decompresses latest Pfam-A.
    """
    pfam_root = _pfam_root()
    pfam_root.mkdir(parents=True, exist_ok=True)

    hmm_gz = _pfam_hmm_gz_path()
    hmm = _pfam_hmm_path()
    manifest = _pfam_manifest_path()

    needs_fetch = force_refresh or not hmm.exists()
    if needs_fetch:
        _stream_download_to_path(PFAM_HMM_GZ_URL, hmm_gz)
        _decompress_gzip_to_path(hmm_gz, hmm)

        metadata = {
            "name": PFAM_NAME,
            "source_url": PFAM_HMM_GZ_URL,
            "downloaded_at": datetime.now(timezone.utc).isoformat(),
            "hmm_path": str(hmm),
            "hmm_gz_path": str(hmm_gz),
            "hmm_sha256": _sha256(hmm),
            "hmm_gz_sha256": _sha256(hmm_gz),
        }
        _write_manifest_atomically(manifest, metadata)

    if manifest.exists():
        data = json.loads(manifest.read_text())
    else:
        data = {
            "name": PFAM_NAME,
            "source_url": PFAM_HMM_GZ_URL,
            "downloaded_at": "unknown",
            "hmm_path": str(hmm),
        }

    data["hmm_path"] = str(hmm)
    return data


def prepare_pfam_database(
    download: bool = True,
    index: bool = True,
    force_refresh: bool = False,
) -> Dict[str, object]:
    """Prepare the local Pfam database by downloading and/or indexing it."""
    hmm = _pfam_hmm_path()

    if download or force_refresh:
        pfam_meta = ensure_pfam_database(force_refresh=force_refresh)
    else:
        if not hmm.exists():
            raise FileNotFoundError(
                f"PFAM HMM database not found: {hmm}. Run `phu prepare-dbs` with download enabled first."
            )
        pfam_meta = ensure_pfam_database(force_refresh=False)

    result: Dict[str, object] = dict(pfam_meta)
    result["downloaded"] = bool(download or force_refresh)

    if index:
        if not _is_offsets_index_valid(hmm):
            _build_offsets_index(hmm)
        result["offsets_path"] = str(_pfam_offsets_index_path(hmm))
        result["indexed"] = True
    else:
        result["indexed"] = False

    return result


def refresh_pfam_database() -> Dict[str, object]:
    """Refresh PFAM metadata/index integrity and repair incomplete local state."""
    hmm = _pfam_hmm_path()
    hmm_gz = _pfam_hmm_gz_path()
    manifest = _pfam_manifest_path()

    if not hmm.exists():
        result = prepare_pfam_database(download=True, index=True, force_refresh=False)
        result["refreshed"] = True
        return result

    if not manifest.exists():
        metadata: Dict[str, object] = {
            "name": PFAM_NAME,
            "source_url": PFAM_HMM_GZ_URL,
            "downloaded_at": datetime.now(timezone.utc).isoformat(),
            "hmm_path": str(hmm),
            "hmm_sha256": _sha256(hmm),
        }
        if hmm_gz.exists():
            metadata["hmm_gz_path"] = str(hmm_gz)
            metadata["hmm_gz_sha256"] = _sha256(hmm_gz)
        _write_manifest_atomically(manifest, metadata)

    manifest_data = _read_json(manifest)
    expected_sha = manifest_data.get("hmm_sha256")
    if isinstance(expected_sha, str):
        current_sha = _sha256(hmm)
        if current_sha != expected_sha:
            result = prepare_pfam_database(download=True, index=True, force_refresh=True)
            result["refreshed"] = True
            return result

    result = prepare_pfam_database(download=False, index=True, force_refresh=False)
    result["refreshed"] = True
    return result


def remove_pfam_database() -> bool:
    """Remove the full local PFAM database directory."""
    pfam_root = _pfam_root()
    if not pfam_root.exists():
        return False
    shutil.rmtree(pfam_root)
    return True


def get_pfam_database_status() -> Dict[str, object]:
    """Return local PFAM database status metadata."""
    pfam_root = _pfam_root()
    hmm = _pfam_hmm_path()
    manifest = _pfam_manifest_path()
    offsets = _pfam_offsets_index_path(hmm)
    models_dir = _pfam_models_dir(hmm)

    index_valid = _is_offsets_index_valid(hmm)
    model_count = 0
    if offsets.exists():
        try:
            index_data = _read_json(offsets)
            offset_map = index_data.get("offsets")
            if isinstance(offset_map, dict):
                model_count = len(offset_map)
        except (json.JSONDecodeError, OSError, ValueError):
            model_count = 0

    sparse_cached_count = 0
    if models_dir.exists():
        sparse_cached_count = len(list(models_dir.glob("*.hmm")))

    return {
        "name": "pfam",
        "root": str(pfam_root),
        "hmm_path": str(hmm),
        "manifest_path": str(manifest),
        "offsets_path": str(offsets),
        "models_dir": str(models_dir),
        "downloaded": hmm.exists(),
        "manifest_exists": manifest.exists(),
        "indexed": index_valid,
        "model_count": model_count,
        "sparse_cached_count": sparse_cached_count,
    }


def extract_pfam_models(
    hmm_db_path: Path,
    requested_ids: Iterable[str],
    output_dir: Path,
) -> Tuple[List[Path], List[str]]:
    """
    Extract requested PFAM models from a Pfam-A HMM database into individual files.

    Returns:
        (list_of_extracted_paths, list_of_missing_normalized_ids)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    return _extract_from_split_cache(
        hmm_db_path=hmm_db_path,
        requested_ids=requested_ids,
        output_dir=output_dir,
    )
