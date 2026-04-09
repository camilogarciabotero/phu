from __future__ import annotations

import gzip
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple
from urllib.request import urlopen

import click

from ._click import run_click_task

PFAM_HMM_GZ_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_ID_PREFIX = "PF"
PFAM_NAME = "Pfam-A"

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
        manifest.write_text(json.dumps(metadata, indent=2))

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
    normalized_ids: Set[str] = {normalize_pfam_id(token) for token in requested_ids}
    output_dir.mkdir(parents=True, exist_ok=True)

    out_paths: List[Path] = []
    found: Set[str] = set()
    total_size = hmm_db_path.stat().st_size if hmm_db_path.exists() else None

    block_lines: List[str] = []
    block_acc: str | None = None

    def _flush_block() -> None:
        nonlocal block_lines, block_acc
        if block_acc is None or block_acc not in normalized_ids:
            block_lines = []
            block_acc = None
            return

        out_path = output_dir / f"{block_acc}.hmm"
        out_path.write_text("".join(block_lines))
        out_paths.append(out_path)
        found.add(block_acc)
        block_lines = []
        block_acc = None

    if total_size is not None and total_size > 0:
        with click.progressbar(
            length=total_size,
            label="Scanning Pfam-A",
            show_pos=True,
            show_percent=True,
            show_eta=True,
        ) as bar:
            pending_progress = 0
            progress_chunk = 1024 * 1024
            with hmm_db_path.open("r", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    line = raw_line.rstrip("\n")
                    block_lines.append(raw_line)
                    pending_progress += len(raw_line)
                    if pending_progress >= progress_chunk:
                        bar.update(pending_progress)
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

            if pending_progress:
                bar.update(pending_progress)
    else:
        def _scan_file() -> None:
            nonlocal block_acc, block_lines
            with hmm_db_path.open("r", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    line = raw_line.rstrip("\n")
                    block_lines.append(raw_line)

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

        run_click_task("Scanning Pfam-A", _scan_file)

    missing = sorted(normalized_ids - found)
    return out_paths, missing
