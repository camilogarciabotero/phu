from __future__ import annotations

import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


def get_cache_dir() -> Path:
    """Resolve the shared PHU cache directory."""
    if cache_env := os.environ.get("PHU_CACHE_DIR"):
        return Path(cache_env)
    if xdg := os.environ.get("XDG_CACHE_HOME"):
        return Path(xdg) / "phu"
    return Path.home() / ".cache" / "phu"


@dataclass
class AvgerInputs:
    """Inputs that should determine avger cache validity."""

    input_contigs: Path
    mode: str = "meta"
    threads: int = 1
    max_evalue: float = 1e-5
    min_gene_len: int = 90
    min_protein_len_aa: int = 30
    translation_table: int = 11
    keep_all_hits: bool = False
    use_vscore: bool = True
    require_flank_support: bool = False

    def __post_init__(self) -> None:
        self.input_contigs = Path(self.input_contigs).resolve()
        if not self.input_contigs.exists():
            raise FileNotFoundError(f"Input contigs not found: {self.input_contigs}")
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.max_evalue < 0:
            raise ValueError("max_evalue must be >= 0")
        if self.min_gene_len < 1:
            raise ValueError("min_gene_len must be >= 1")
        if self.min_protein_len_aa < 1:
            raise ValueError("min_protein_len_aa must be >= 1")


def compute_avger_cache_key(inputs: AvgerInputs) -> str:
    """Serialize avger inputs into a stable cache key."""
    stat = inputs.input_contigs.stat()
    seed = (
        str(inputs.input_contigs)
        + str(stat.st_mtime_ns)
        + inputs.mode
        + str(inputs.threads)
        + str(inputs.max_evalue)
        + str(inputs.min_gene_len)
        + str(inputs.min_protein_len_aa)
        + str(inputs.translation_table)
        + str(bool(inputs.keep_all_hits))
        + str(bool(inputs.use_vscore))
        + str(bool(inputs.require_flank_support))
    )
    return hashlib.sha256(seed.encode()).hexdigest()[:16]


def _avger_cache_root() -> Path:
    return get_cache_dir() / "avger" / "v1"


def _cached_output_paths(cache_key: str) -> tuple[Path, Path]:
    root = _avger_cache_root() / cache_key
    return root / "best_hits.tsv", root / "manifest.json"


def get_or_run_avger(
    inputs: AvgerInputs,
    output_folder: Path,
    run_func: Callable[[AvgerInputs, Path], Path],
) -> tuple[Path, bool]:
    """Run avger once per deterministic input and cache the final best_hits TSV."""
    inputs.__post_init__()
    cache_key = compute_avger_cache_key(inputs)
    cache_best_hits, cache_manifest = _cached_output_paths(cache_key)
    cache_best_hits.parent.mkdir(parents=True, exist_ok=True)

    if cache_best_hits.exists() and cache_manifest.exists():
        output_folder.mkdir(parents=True, exist_ok=True)
        target_best = output_folder / "best_hits.tsv"
        if target_best != cache_best_hits:
            shutil.copy2(cache_best_hits, target_best)
        return target_best, True

    output_folder.mkdir(parents=True, exist_ok=True)
    generated = run_func(inputs, output_folder)

    if generated.exists():
        shutil.copy2(generated, cache_best_hits)

    manifest = {
        "input_contigs": str(inputs.input_contigs),
        "mode": inputs.mode,
        "threads": inputs.threads,
        "max_evalue": inputs.max_evalue,
        "min_gene_len": inputs.min_gene_len,
        "min_protein_len_aa": inputs.min_protein_len_aa,
        "translation_table": inputs.translation_table,
        "keep_all_hits": inputs.keep_all_hits,
        "use_vscore": inputs.use_vscore,
        "require_flank_support": inputs.require_flank_support,
        "cache_key": cache_key,
    }
    cache_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    return generated, False
