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

try:
    import fcntl
except ModuleNotFoundError:
    fcntl = None
import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from pyrodigal_gv import ViralGeneFinder


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


@dataclass(frozen=True)
class PredictedGene:
    """Typed representation of one predicted coding sequence."""

    contig_id: str
    gene_id: str
    start: int
    end: int
    strand: int
    ordinal: int
    nucleotide_sequence: str
    amino_acid_sequence: str


@dataclass
class CacheArtifact:
    """Result of prediction: where to find proteins + metadata."""

    proteins_path: Path
    protein_count: int
    cache_hit: bool
    cache_key: str
    cache_dir: Optional[Path] = None
    temp_dir: Optional[Path] = None
    genes: Optional[list[PredictedGene]] = None


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


def clean_prediction_cache() -> tuple[Path, bool]:
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


def _read_fasta_python(fp: Path) -> Iterable[tuple[str, str]]:
    """Read FASTA records using Python text IO (supports .gz)."""
    import gzip

    opener = gzip.open if fp.suffix == ".gz" else open
    with opener(fp, "rt") as handle:
        seq_id: Optional[str] = None
        seq_chunks: list[str] = []

        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_chunks)
                # FASTA header id is first token after '>'
                seq_id = line[1:].split(None, 1)[0]
                seq_chunks = []
                continue

            if seq_id is None:
                raise ValueError(
                    f"Invalid FASTA format in {fp}: sequence before header"
                )

            seq_chunks.append(line)

        if seq_id is not None:
            yield seq_id, "".join(seq_chunks)


def predict_genes_pyrodigal(
    inputs: PredictionInputs,
    *,
    threads: int = 1,
) -> list[PredictedGene]:
    """Predict genes for all contigs and return typed gene records."""

    def _predict_for_contig(contig_id: str, seq: str) -> list[PredictedGene]:
        max_overlap = max(0, inputs.min_gene_len - 1)
        finder = ViralGeneFinder(
            meta=(inputs.mode == "meta"),
            min_gene=inputs.min_gene_len,
            max_overlap=max_overlap,
        )
        if inputs.mode == "single":
            train_seq = seq
            if len(seq) < 100000:
                # pyrodigal-gv single-mode requires long training context.
                train_seq = seq + ("A" * (100000 - len(seq)))
            finder.train(train_seq, translation_table=inputs.translation_table)

        predicted: list[PredictedGene] = []
        for ordinal, gene in enumerate(finder.find_genes(seq), start=1):
            aa = gene.translate(translation_table=inputs.translation_table)
            if not aa or len(aa) < inputs.min_protein_len_aa:
                continue

            gene_id = f"{contig_id}|gene{ordinal}"
            predicted.append(
                PredictedGene(
                    contig_id=contig_id,
                    gene_id=gene_id,
                    start=int(gene.begin),
                    end=int(gene.end),
                    strand=int(gene.strand),
                    ordinal=ordinal,
                    nucleotide_sequence=str(gene.sequence()),
                    amino_acid_sequence=str(aa),
                )
            )

        return predicted

    contigs = list(_read_fasta_python(inputs.input_contigs))
    if not contigs:
        return []

    if threads > 1:
        from multiprocessing.pool import ThreadPool

        with ThreadPool(processes=threads) as pool:
            nested = pool.starmap(_predict_for_contig, contigs)
    else:
        nested = [_predict_for_contig(contig_id, seq) for contig_id, seq in contigs]

    out: list[PredictedGene] = []
    for genes in nested:
        out.extend(genes)
    return out


def write_predicted_proteins_fasta(genes: Iterable[PredictedGene], output_path: Path) -> int:
    """Write predicted proteins to FASTA and return record count."""
    count = 0
    with output_path.open("w") as out:
        for gene in genes:
            out.write(f">{gene.gene_id}\n{gene.amino_acid_sequence}\n")
            count += 1
    return count


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
    inputs.__post_init__()  # Validate

    if not use_cache:
        # No caching: predict to tempfile and return immediately
        temp_dir = Path(tempfile.mkdtemp(prefix="phu_pred_"))
        temp_prot = temp_dir / "proteins.faa"

        genes = predict_genes_pyrodigal(inputs, threads=threads)
        n_prot = write_predicted_proteins_fasta(genes, temp_prot)

        return CacheArtifact(
            proteins_path=temp_prot,
            protein_count=n_prot,
            cache_hit=False,
            cache_key="",
            cache_dir=None,
            temp_dir=temp_dir,
            genes=genes,
        )

    # Cache-aware path
    cache_dir = get_cache_dir()
    cache_root = cache_dir / "v1"
    cache_key = compute_cache_key(inputs)
    cache_subdir = cache_root / cache_key
    cache_proteins = cache_subdir / "proteins.faa"
    cache_manifest = cache_subdir / "manifest.json"
    cache_genes = cache_subdir / "genes.json"
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
                    genes=[PredictedGene(**item) for item in json.loads(cache_genes.read_text())]
                    if cache_genes.exists()
                    else None,
                )
            except (json.JSONDecodeError, KeyError):
                # Corrupted manifest; treat as miss and rebuild
                cache_manifest.unlink(missing_ok=True)

        # Cache miss or incomplete: rebuild
        # Clean up any stale partial from previous interrupted run
        if partial_dir.exists():
            shutil.rmtree(partial_dir, ignore_errors=True)

        partial_dir.mkdir(parents=True, exist_ok=True)
        temp_prot = partial_dir / "proteins.faa"

        genes = predict_genes_pyrodigal(inputs, threads=threads)
        n_prot = write_predicted_proteins_fasta(genes, temp_prot)

        # Atomically promote partial to cache
        cache_subdir.mkdir(parents=True, exist_ok=True)
        temp_prot.replace(cache_proteins)
        cache_genes.write_text(json.dumps([gene.__dict__ for gene in genes], indent=2))

        # Clean up partial directory after successful promotion
        shutil.rmtree(partial_dir, ignore_errors=True)

        # Write manifest with prediction metadata
        # Write manifest atomically to avoid partial writes on crash
        manifest_data = {
            "input_contigs": str(inputs.input_contigs),
            "mode": inputs.mode,
            "min_gene_len": inputs.min_gene_len,
            "min_protein_len_aa": inputs.min_protein_len_aa,
            "translation_table": inputs.translation_table,
            "protein_count": n_prot,
            "cache_key": cache_key,
        }
        tmp_manifest = cache_subdir / ".manifest.tmp"
        tmp_manifest.write_text(json.dumps(manifest_data, indent=2))
        tmp_manifest.replace(cache_manifest)

        return CacheArtifact(
            proteins_path=cache_proteins,
            protein_count=n_prot,
            cache_hit=False,
            cache_key=cache_key,
            cache_dir=cache_subdir,
            genes=genes,
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
