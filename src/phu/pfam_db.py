from __future__ import annotations

import json
from collections.abc import Callable, Iterable
from datetime import datetime, timezone
from pathlib import Path

import click

from ._click import run_click_task
from .db_bundle import derived_dir, derived_files_lock, resolve_db, write_atomically

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
        return (
            len(acc) == 5 and acc.isdigit() and version.isdigit() and len(version) > 0
        )
    return len(body) == 5 and body.isdigit()


def normalize_pfam_id(token: str) -> str:
    """Normalize PFAM accession to unversioned upper-case format (e.g. PF00001)."""
    token = token.strip().upper()
    if not is_pfam_id(token):
        raise ValueError(f"Invalid PFAM accession: {token}")
    return token.split(".", 1)[0]


def ensure_pfam_database() -> dict[str, str]:
    """Locate Pfam in the phu-db bundle, installing the bundle on first use."""
    return {"name": PFAM_NAME, "hmm_path": str(resolve_db("pfam") / "profiles.hmm")}


def _pfam_models_dir(hmm_db_path: Path) -> Path:
    return derived_dir(hmm_db_path.parent) / "models"


def _pfam_offsets_index_path(hmm_db_path: Path) -> Path:
    return derived_dir(hmm_db_path.parent) / "offsets.json"


def _clear_sparse_models_cache(hmm_db_path: Path) -> None:
    """Remove cached per-model HMM files to avoid serving stale entries."""
    models_dir = _pfam_models_dir(hmm_db_path)
    if not models_dir.exists():
        return
    for model_file in models_dir.glob("*.hmm"):
        model_file.unlink()


def _write_manifest_atomically(
    manifest_path: Path, metadata: dict[str, object]
) -> None:
    """Write manifest JSON atomically to avoid corruption on interruption."""
    write_atomically(manifest_path, json.dumps(metadata, indent=2).encode())


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


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
    offsets: dict[str, list[int]] = {}

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

    metadata: dict[str, object] = {
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
) -> tuple[list[Path], list[str]]:
    """Resolve requested PFAM IDs via offset index and sparse local model cache."""
    # Validate index against source DB before serving any sparse cached models.
    # If source changed, clear stale per-model cache and rebuild offsets.
    # Jobs sharing the install do this one at a time; one that waited finds the
    # index valid and leaves the other's cache alone.
    if not _is_offsets_index_valid(hmm_db_path):
        with derived_files_lock(derived_dir(hmm_db_path.parent)):
            if not _is_offsets_index_valid(hmm_db_path):
                _clear_sparse_models_cache(hmm_db_path)
                _build_offsets_index(hmm_db_path)

    models_dir = _pfam_models_dir(hmm_db_path)
    models_dir.mkdir(parents=True, exist_ok=True)
    normalized_ordered = list(
        dict.fromkeys(normalize_pfam_id(token) for token in requested_ids)
    )

    out_paths_by_id: dict[str, Path] = {}
    missing: list[str] = []

    needs_index = False
    for pfam_id in normalized_ordered:
        model_path = models_dir / f"{pfam_id}.hmm"
        if not model_path.exists():
            needs_index = True
            continue
        out_path = output_dir / f"{pfam_id}.hmm"
        out_path.write_bytes(model_path.read_bytes())
        out_paths_by_id[pfam_id] = out_path

    offsets: dict[str, object] = {}
    if needs_index:
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
                if (
                    not isinstance(start, int)
                    or not isinstance(end, int)
                    or end <= start
                ):
                    missing.append(pfam_id)
                    continue

                handle.seek(start)
                model_blob = handle.read(end - start)
                cache_path = models_dir / f"{pfam_id}.hmm"
                write_atomically(cache_path, model_blob)

                out_path = output_dir / f"{pfam_id}.hmm"
                out_path.write_bytes(model_blob)
                out_paths_by_id[pfam_id] = out_path

    out_paths = [
        out_paths_by_id[pfam_id]
        for pfam_id in normalized_ordered
        if pfam_id in out_paths_by_id
    ]
    return out_paths, missing


def extract_pfam_models(
    hmm_db_path: Path,
    requested_ids: Iterable[str],
    output_dir: Path,
) -> tuple[list[Path], list[str]]:
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
