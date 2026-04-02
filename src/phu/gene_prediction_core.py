"""
Gene prediction caching core module.

Provides transparent, mmseqs-style caching for protein predictions.
Cache is implicit and environment-controlled (PHU_CACHE_DIR, PHU_CACHE).

Key features:
- Deterministic cache key from contigs + prediction params
- Crash-safe atomic operations with .lock and .partial handling
- XDG-compliant cache directory defaults
- Full backward compatibility with existing screen/jack interfaces
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Tuple

from ._exec import CmdNotFound


@dataclass
class PredictionInputs:
    """Deterministic inputs to gene prediction (used for cache key)."""

    input_contigs: Path
    mode: str = "meta"
    min_gene_len: int = 90
    min_protein_len_aa: int = 30
    translation_table: int = 11

    def __post_init__(self) -> None:
        """Validate inputs."""
        self.input_contigs = Path(self.input_contigs).resolve()
        if not self.input_contigs.exists():
            raise FileNotFoundError(f"Input contigs not found: {self.input_contigs}")
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        if self.min_gene_len < 1:
            raise ValueError("min_gene_len must be >= 1")
        if self.min_protein_len_aa < 1:
            raise ValueError("min_protein_len_aa must be >= 1")


@dataclass
class CacheArtifact:
    """Result of prediction: where to find proteins + metadata."""

    proteins_path: Path
    protein_count: int
    cache_hit: bool
    cache_key: str
    cache_dir: Optional[Path] = None


def compute_cache_key(inputs: PredictionInputs) -> str:
    """
    Generate stable, deterministic cache key from prediction inputs.

    Key includes file identity (mtime as quick invalidation) and all prediction params.
    Seeds/HMMs/markers are explicitly excluded to allow reuse across different searches.
    """
    stat = inputs.input_contigs.stat()

    seed = (
        str(inputs.input_contigs)
        + str(stat.st_mtime_ns)  # Quick invalidation on file change
        + inputs.mode
        + str(inputs.min_gene_len)
        + str(inputs.min_protein_len_aa)
        + str(inputs.translation_table)
    )
    return hashlib.sha256(seed.encode()).hexdigest()[:16]


def get_cache_dir() -> Path:
    """Resolve cache directory from env or defaults (XDG-compliant)."""
    if cache_env := os.environ.get("PHU_CACHE_DIR"):
        return Path(cache_env)

    # XDG-compliant fallback
    if xdg := os.environ.get("XDG_CACHE_HOME"):
        return Path(xdg) / "phu"
    return Path.home() / ".cache" / "phu"


def clean_prediction_cache() -> Tuple[Path, bool]:
    """Remove the full prediction cache directory and report whether it existed."""
    cache_dir = get_cache_dir()
    existed = cache_dir.exists()
    if existed:
        shutil.rmtree(cache_dir)
    return cache_dir, existed


def _acquire_lock(lock_file: Path) -> Optional:
    """
    Acquire exclusive lock on cache.

    Returns file handle (must be kept open during critical section).
    On Windows or if fcntl unavailable, returns None (no-op lock).
    """
    lock_file.parent.mkdir(parents=True, exist_ok=True)
    fh = open(lock_file, "a")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        return fh
    except (AttributeError, OSError):
        # fcntl not available (Windows) or locking failed
        # Fall back to no-op; assume single-threaded or separate processes won't collide
        return fh


def _release_lock(lock_fh: Optional) -> None:
    """Release lock and close file handle."""
    if lock_fh is None:
        return
    try:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
    except (AttributeError, OSError):
        pass
    lock_fh.close()


def get_or_predict_proteins(
    inputs: PredictionInputs,
    use_cache: bool = True,
    threads: int = 1,
) -> CacheArtifact:
    """
    Predict proteins, reusing cache if possible.

    If use_cache=False, always predict fresh and return tempfile path (no caching).
    If use_cache=True:
      - Check for existing cache entry
      - If hit, return path to cached proteins
      - If miss, predict to temp location, atomically promote to cache, return path

    Crash safety: .partial/ directories are cleaned up on next run.

    Args:
        inputs: PredictionInputs with contigs path and prediction parameters
        use_cache: Whether to use cache (default True; disable with PHU_CACHE=off)
        threads: Number of threads for gene prediction

    Returns:
        CacheArtifact with proteins_path, protein_count, cache_hit status, and cache metadata
    """
    # Import here to avoid circular dependency
    from .screen import _predict_proteins_pyrodigal

    inputs.__post_init__()  # Validate

    if not use_cache:
        # No caching: predict to tempfile and return immediately
        temp_dir = Path(tempfile.mkdtemp(prefix="phu_pred_"))
        temp_prot = temp_dir / "proteins.faa"

        n_prot = _predict_proteins_pyrodigal(
            inputs.input_contigs,
            temp_prot,
            mode=inputs.mode,
            min_len=inputs.min_gene_len,
            min_protein_len_aa=inputs.min_protein_len_aa,
            translation_table=inputs.translation_table,
            threads=threads,
        )

        return CacheArtifact(
            proteins_path=temp_prot,
            protein_count=n_prot,
            cache_hit=False,
            cache_key="",
            cache_dir=None,
        )

    # Cache-aware path
    cache_dir = get_cache_dir()
    cache_root = cache_dir / "v1"
    cache_key = compute_cache_key(inputs)
    cache_subdir = cache_root / cache_key
    cache_proteins = cache_subdir / "proteins.faa"
    cache_manifest = cache_subdir / "manifest.json"
    partial_dir = cache_root / f"{cache_key}.partial"
    lock_file = cache_subdir / ".lock"

    # Ensure cache root exists
    cache_root.mkdir(parents=True, exist_ok=True)

    # Acquire exclusive lock for this cache key
    lock_fh = _acquire_lock(lock_file)

    try:
        # Check for cache hit (must recheck after lock acquisition in case another process won)
        if cache_proteins.exists() and cache_manifest.exists():
            try:
                manifest = json.loads(cache_manifest.read_text())
                n_prot = manifest.get("protein_count", 0)
                return CacheArtifact(
                    proteins_path=cache_proteins,
                    protein_count=n_prot,
                    cache_hit=True,
                    cache_key=cache_key,
                    cache_dir=cache_subdir,
                )
            except (json.JSONDecodeError, KeyError):
                # Corrupted manifest; treat as miss and rebuild
                pass

        # Cache miss or incomplete: rebuild
        # Clean up any stale partial from previous interrupted run
        if partial_dir.exists():
            shutil.rmtree(partial_dir, ignore_errors=True)

        partial_dir.mkdir(parents=True, exist_ok=True)
        temp_prot = partial_dir / "proteins.faa"

        n_prot = _predict_proteins_pyrodigal(
            inputs.input_contigs,
            temp_prot,
            mode=inputs.mode,
            min_len=inputs.min_gene_len,
            min_protein_len_aa=inputs.min_protein_len_aa,
            translation_table=inputs.translation_table,
            threads=threads,
        )

        # Atomically promote partial to cache
        cache_subdir.mkdir(parents=True, exist_ok=True)
        temp_prot.rename(cache_proteins)

        # Write manifest with prediction metadata
        manifest_data = {
            "input_contigs": str(inputs.input_contigs),
            "mode": inputs.mode,
            "min_gene_len": inputs.min_gene_len,
            "min_protein_len_aa": inputs.min_protein_len_aa,
            "translation_table": inputs.translation_table,
            "protein_count": n_prot,
            "cache_key": cache_key,
        }
        cache_manifest.write_text(json.dumps(manifest_data, indent=2))

        return CacheArtifact(
            proteins_path=cache_proteins,
            protein_count=n_prot,
            cache_hit=False,
            cache_key=cache_key,
            cache_dir=cache_subdir,
        )

    finally:
        _release_lock(lock_fh)


def write_prediction_metadata(
    path: Path,
    cache_hit: bool,
    cache_key: str,
    cache_dir: Optional[Path],
) -> None:
    """
    Write optional prediction metadata JSON to output folder.

    Useful for debugging and reproducibility tracking.
    """
    meta = {
        "cache_hit": cache_hit,
        "cache_key": cache_key,
        "cache_dir": str(cache_dir) if cache_dir else None,
    }
    path.write_text(json.dumps(meta, indent=2))
