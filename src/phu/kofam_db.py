from __future__ import annotations

import gzip
import hashlib
import json
import logging
import lzma
import os
import shutil
import tempfile
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import click
import typer

from ._click import run_click_task

logger = logging.getLogger(__name__)

KOFAM_BASE_URL = "https://zenodo.org/records/19503464/files/"
KOFAM_KO_LIST_GZ_URL = KOFAM_BASE_URL + "ko_list.gz"
KOFAM_HMM_XZ_URL = KOFAM_BASE_URL + "kofam.hmm.xz"
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


def _db_root() -> Path:
    if db_env := os.environ.get("PHU_DB_FOLDER"):
        return Path(db_env)

    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


def _kofam_root() -> Path:
    return _db_root() / "kofam"


def _kofam_ko_list_gz_path() -> Path:
    return _kofam_root() / "ko_list.gz"


def _kofam_ko_list_path() -> Path:
    return _kofam_root() / "ko_list"


def _kofam_hmm_xz_path() -> Path:
    return _kofam_root() / "kofam.hmm.xz"


def _kofam_hmm_path() -> Path:
    return _kofam_root() / "kofam.hmm"


def _kofam_models_dir() -> Path:
    return _kofam_root() / "models"


def _kofam_manifest_path() -> Path:
    return _kofam_root() / "manifest.json"


def _kofam_metadata_index_path() -> Path:
    return _kofam_root() / "ko_metadata.json"


def _kofam_offsets_index_path() -> Path:
    return _kofam_root() / "offsets.json"


def _get_file_size(url: str) -> Optional[int]:
    """Get remote file size from Content-Length when available."""
    try:
        with urlopen(url) as response:
            content_length = response.info().get("Content-Length")
            if content_length and content_length.isdigit():
                return int(content_length)
    except (HTTPError, URLError, OSError, ValueError) as exc:
        logger.debug("Could not determine remote file size for %s: %s", url, exc)
    return None


def _stream_download_to_path(url: str, destination: Path, label: str) -> None:
    """Download URL content to destination atomically via temporary file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as tmp:
        tmp_path = Path(tmp.name)
        with urlopen(url) as response:
            content_length = response.info().get("Content-Length")
            total_size = (
                int(content_length)
                if content_length and content_length.isdigit()
                else None
            )

            if total_size is not None and total_size > 0:
                with click.progressbar(
                    length=total_size,
                    label=label,
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

                run_click_task(label, _copy_stream)

    tmp_path.replace(destination)


def _download_range(
    url: str,
    start: int,
    end: int,
    destination: Path,
    progress_callback: Optional[Callable[[int], None]] = None,
) -> None:
    """Download one byte-range chunk into a temporary part file."""
    req = Request(url)
    req.add_header("Range", f"bytes={start}-{end}")

    with urlopen(req) as response:
        status = getattr(response, "status", response.getcode())
        content_range = response.headers.get("Content-Range")
        if status != 206 or not content_range:
            raise RuntimeError(
                "Server did not honor byte-range request; falling back to single-stream download"
            )

        chunk_file = destination.parent / f"{destination.name}.part.{start}-{end}"
        with chunk_file.open("wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                if progress_callback is not None:
                    progress_callback(len(chunk))


def _range_download_supported(url: str) -> bool:
    """Probe whether the server supports partial-content downloads."""
    req = Request(url)
    req.add_header("Range", "bytes=0-0")

    try:
        with urlopen(req) as response:
            status = getattr(response, "status", response.getcode())
            content_range = response.headers.get("Content-Range", "")
            return status == 206 and content_range.startswith("bytes 0-0/")
    except (HTTPError, URLError, OSError, ValueError):
        return False


def _download_parallel_chunked(
    url: str, destination: Path, label: str, num_chunks: int = 4
) -> None:
    """Download file with HTTP range requests (pget-like parallel mode)."""
    destination.parent.mkdir(parents=True, exist_ok=True)

    file_size = _get_file_size(url)
    if file_size is None or file_size < 1024 * 1024:
        typer.secho(f"Downloading {label} (single connection)...", fg="cyan")
        _stream_download_to_path(url, destination, label)
        return

    if not _range_download_supported(url):
        typer.secho(f"Downloading {label} (single connection)...", fg="cyan")
        _stream_download_to_path(url, destination, label)
        return

    typer.secho(
        f"Downloading {label} ({num_chunks} parallel connections)...", fg="cyan"
    )
    chunk_size = file_size // num_chunks

    with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as tmp:
        tmp_path = Path(tmp.name)

    ranges: list[tuple[int, int]] = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = file_size - 1 if i == num_chunks - 1 else start + chunk_size - 1
        ranges.append((start, end))

    progress_lock = Lock()
    downloaded = [0]

    def _progress_callback(size: int) -> None:
        with progress_lock:
            downloaded[0] += size
            pct = int(100 * downloaded[0] / file_size)
            typer.secho(
                f"\r{label}: {pct}% ({downloaded[0] // 1024 // 1024}MB / {file_size // 1024 // 1024}MB)",
                nl=False,
            )

    with ThreadPoolExecutor(max_workers=num_chunks) as executor:
        futures = [
            executor.submit(
                _download_range, url, start, end, destination, _progress_callback
            )
            for start, end in ranges
        ]
        for future in as_completed(futures):
            future.result()

    typer.secho("")

    with tmp_path.open("wb") as out:
        for start, end in ranges:
            chunk_file = destination.parent / f"{destination.name}.part.{start}-{end}"
            with chunk_file.open("rb") as src:
                shutil.copyfileobj(src, out, length=1024 * 1024)
            chunk_file.unlink()

    tmp_path.replace(destination)


def _decompress_gzip_to_path(gz_path: Path, out_path: Path) -> None:
    """Decompress a gzip file atomically."""
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


def _decompress_xz_to_path(xz_path: Path, out_path: Path) -> None:
    """Decompress an xz file atomically."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=out_path.parent) as tmp:
        tmp_path = Path(tmp.name)
        with lzma.open(xz_path, "rb") as src:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                tmp.write(chunk)
    tmp_path.replace(out_path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _write_manifest_atomically(
    manifest_path: Path, metadata: dict[str, object]
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=manifest_path.parent, suffix=".json"
    ) as tmp:
        tmp_path = Path(tmp.name)
        json.dump(metadata, tmp, indent=2)
    tmp_path.replace(manifest_path)


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
            raise FileNotFoundError(
                f"KOFam metadata not found: {ko_list_path}. Run `phu dbs prepare kofam` first."
            )
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


def ensure_kofam_database(force_refresh: bool = False) -> dict[str, str]:
    """Ensure local KOFam database exists and return metadata."""
    root = _kofam_root()
    root.mkdir(parents=True, exist_ok=True)

    ko_list_gz = _kofam_ko_list_gz_path()
    ko_list = _kofam_ko_list_path()
    hmm_xz = _kofam_hmm_xz_path()
    hmm = _kofam_hmm_path()
    manifest = _kofam_manifest_path()

    needs_fetch = force_refresh or not ko_list.exists() or not hmm.exists()

    if needs_fetch:
        typer.secho("Downloading KOFam database files...", fg="cyan")

        _stream_download_to_path(KOFAM_KO_LIST_GZ_URL, ko_list_gz, "ko_list.gz")
        _download_parallel_chunked(
            KOFAM_HMM_XZ_URL, hmm_xz, "kofam.hmm.xz", num_chunks=4
        )

        typer.secho("Decompressing KOFam files...", fg="cyan")
        _decompress_gzip_to_path(ko_list_gz, ko_list)
        _decompress_xz_to_path(hmm_xz, hmm)

        ko_metadata = _build_kofam_metadata_index(ko_list)
        _clear_sparse_models_cache()
        _build_offsets_index(hmm)

        metadata = {
            "name": KOFAM_NAME,
            "source_url": KOFAM_BASE_URL,
            "downloaded_at": datetime.now(timezone.utc).isoformat(),
            "ko_list_path": str(ko_list),
            "ko_list_gz_path": str(ko_list_gz),
            "hmm_path": str(hmm),
            "hmm_xz_path": str(hmm_xz),
            "ko_list_sha256": _sha256(ko_list),
            "ko_list_gz_sha256": _sha256(ko_list_gz),
            "hmm_sha256": _sha256(hmm),
            "hmm_xz_sha256": _sha256(hmm_xz),
            "model_count": len(ko_metadata),
            "ko_metadata_count": len(ko_metadata),
        }
        _write_manifest_atomically(manifest, metadata)

    if manifest.exists():
        data = _read_json(manifest)
    else:
        data = {
            "name": KOFAM_NAME,
            "source_url": KOFAM_BASE_URL,
            "downloaded_at": "unknown",
            "hmm_path": str(hmm),
            "ko_list_path": str(ko_list),
        }

    data["hmm_path"] = str(hmm)
    data["ko_list_path"] = str(ko_list)
    return {k: str(v) for k, v in data.items()}


def prepare_kofam_database(force_refresh: bool = False) -> dict[str, object]:
    result = ensure_kofam_database(force_refresh=force_refresh)
    hmm = _kofam_hmm_path()
    if not _is_offsets_index_valid(hmm):
        _clear_sparse_models_cache()
        _build_offsets_index(hmm)
    result["offsets_path"] = str(_kofam_offsets_index_path())
    result["prepared"] = True
    return result


def refresh_kofam_database() -> dict[str, object]:
    ko_list = _kofam_ko_list_path()
    hmm = _kofam_hmm_path()
    if not ko_list.exists() or not hmm.exists():
        return prepare_kofam_database(force_refresh=False)

    _build_kofam_metadata_index(ko_list)
    if not _is_offsets_index_valid(hmm):
        _clear_sparse_models_cache()
        _build_offsets_index(hmm)

    status = get_kofam_database_status()
    status["refreshed"] = True
    return status


def remove_kofam_database() -> bool:
    root = _kofam_root()
    if not root.exists():
        return False
    shutil.rmtree(root)
    return True


def get_kofam_database_status() -> dict[str, object]:
    root = _kofam_root()
    manifest = _kofam_manifest_path()
    ko_list = _kofam_ko_list_path()
    hmm = _kofam_hmm_path()
    hmm_xz = _kofam_hmm_xz_path()
    metadata_index = _kofam_metadata_index_path()
    offsets = _kofam_offsets_index_path()
    models_dir = _kofam_models_dir()

    model_count = 0
    if offsets.exists():
        try:
            index_data = _read_json(offsets)
            offset_map = index_data.get("offsets")
            if isinstance(offset_map, dict):
                model_count = len(offset_map)
        except (json.JSONDecodeError, OSError, ValueError):
            model_count = 0

    sparse_cached_count = (
        len(list(models_dir.glob("*.hmm"))) if models_dir.exists() else 0
    )

    return {
        "name": KOFAM_NAME,
        "root": str(root),
        "hmm_path": str(hmm),
        "hmm_xz_path": str(hmm_xz),
        "ko_list_path": str(ko_list),
        "manifest_path": str(manifest),
        "metadata_index_path": str(metadata_index),
        "offsets_path": str(offsets),
        "downloaded": ko_list.exists() and hmm.exists(),
        "manifest_exists": manifest.exists(),
        "indexed": _is_offsets_index_valid(hmm),
        "model_count": model_count,
        "sparse_cached_count": sparse_cached_count,
    }


def extract_kofam_models(
    requested_ids: Iterable[str],
    output_dir: Path,
) -> tuple[list[Path], list[str]]:
    """Extract requested KO models into output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    hmm_db_path = _kofam_hmm_path()
    if not hmm_db_path.exists():
        raise FileNotFoundError(f"KOFam HMM database not found: {hmm_db_path}")

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
                cache_path.write_bytes(model_blob)

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
