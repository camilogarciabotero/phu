from __future__ import annotations

import os
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyhmmer.easel
import pyhmmer.hmmer

from ._exec import _executable
from .screen import Hit, _read_fasta, _seqkit_extract
from .gene_prediction_core import (
    PredictionInputs,
    get_or_predict_proteins,
    write_prediction_metadata,
)


@dataclass
class JackConfig:
    """Configuration for seed-based iterative jackhmmer screening."""

    input_contigs: Path
    seed_marker: Path
    outdir: Path = Path("phu-jack")
    mode: str = "meta"
    threads: int = 1
    iterations: int = 5
    inc_evalue: float = 1e-3
    max_evalue: Optional[float] = 1e-5
    top_per_contig: int = 1
    min_gene_len: int = 90
    translation_table: int = 11
    min_protein_len_aa: int = 30
    keep_proteins: bool = False
    save_hmm: bool = False
    combine_mode: str = "any"
    min_seed_hits: int = 1

    def __post_init__(self) -> None:
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.iterations < 1:
            raise ValueError("iterations must be >= 1")
        if self.top_per_contig < 1:
            raise ValueError("top_per_contig must be >= 1")
        if self.combine_mode not in {"any", "all", "threshold"}:
            raise ValueError("combine_mode must be 'any', 'all', or 'threshold'")
        if self.min_seed_hits < 1:
            raise ValueError("min_seed_hits must be >= 1")
        if self.min_protein_len_aa < 1:
            raise ValueError("min_protein_len_aa must be >= 1")


def _read_seed_queries(seed_marker: Path):
    """Read one or more seed proteins from FASTA and digitize them."""
    records = list(_read_fasta(seed_marker))
    if not records:
        raise ValueError(f"Seed marker file is empty: {seed_marker}")

    seen = set()
    out = []
    alphabet = pyhmmer.easel.Alphabet.amino()
    for seq_id, seq in records:
        if seq_id in seen:
            raise ValueError(
                f"Duplicate seed sequence ID found: {seq_id}. "
                "Seed IDs must be unique."
            )
        seen.add(seq_id)
        text_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=seq)
        out.append((seq_id, text_seq.digitize(alphabet), alphabet))

    return out


def _run_jackhmmer(
    query,
    alphabet,
    seed_id: str,
    proteins_fa: Path,
    iterations: int,
    inc_evalue: float,
    max_evalue: Optional[float],
    threads: int,
    hmm_output_path: Optional[Path] = None,
) -> Tuple[List[Hit], List[Dict[str, object]]]:
    """Run pyhmmer.hmmer.jackhmmer and return final hits plus iteration summary."""
    with pyhmmer.easel.SequenceFile(
        str(proteins_fa), digital=True, alphabet=alphabet
    ) as seq_file:
        targets = seq_file.read_block()

    results = list(
        pyhmmer.hmmer.jackhmmer(
            query,
            targets,
            max_iterations=iterations,
            checkpoints=True,
            cpus=threads,
            incE=inc_evalue,
            incdomE=inc_evalue,
        )
    )

    if not results:
        return [], []
        
    # jackhmmer returns an iterator over results per query sequence.
    # We provided exactly one query, so there is exactly one result item.
    # Since checkpoints=True, that item is a list of IterationResult objects.
    iterations_out = results[0]
    
    if not iterations_out:
        return [], []

    def _iter_hits(iteration_obj):
        """Return hits for either pyhmmer iteration objects or plain list outputs."""
        return getattr(iteration_obj, "hits", iteration_obj)

    def _iter_meta(iteration_obj, idx: int) -> Tuple[int, bool]:
        """Return (iteration_index, converged) with sensible fallbacks."""
        iteration = getattr(iteration_obj, "iteration", idx)
        converged = getattr(iteration_obj, "converged", False)
        return int(iteration), bool(converged)

    summary: List[Dict[str, object]] = []
    for idx, it in enumerate(iterations_out, start=1):
        hits_in_iter = _iter_hits(it)
        iter_index, iter_converged = _iter_meta(it, idx)
        n_hits = 0
        n_included = 0
        for hit in hits_in_iter:
            n_hits += 1
            if hit.included:
                n_included += 1
        summary.append(
            {
                "iteration": iter_index,
                "n_hits": n_hits,
                "n_included": n_included,
                "converged": iter_converged,
            }
        )

    final = iterations_out[-1]
    if hmm_output_path is not None:
        final_hmm = getattr(final, "hmm", None)
        if final_hmm is not None:
            with hmm_output_path.open("wb") as out_hmm:
                final_hmm.write(out_hmm)

    final_hits = _iter_hits(final)
    kept: List[Hit] = []
    gene_pattern = re.compile(r"\|gene\d+$")

    for hit in final_hits:
        if not hit.included:
            continue
        if max_evalue is not None and hit.evalue > max_evalue:
            continue

        prot_id = hit.name if isinstance(hit.name, str) else hit.name.decode()
        match = gene_pattern.search(prot_id)
        if match:
            contig = prot_id[: match.start()]
        elif "|" in prot_id:
            contig = prot_id.rsplit("|", 1)[0]
        else:
            contig = prot_id

        kept.append(
            Hit(
                contig=contig,
                prot_id=prot_id,
                model=seed_id,
                bitscore=hit.score,
                evalue=hit.evalue,
            )
        )

    return kept, summary


def _choose_top_hits_per_contig(
    hits: List[Hit],
    top_per_contig: int,
    combine_mode: str,
    min_seed_hits: int,
    total_seeds: int,
) -> Tuple[List[Hit], List[str]]:
    """Combine seed hits and keep top hits per contig based on combine mode."""
    per_contig: Dict[str, List[Hit]] = defaultdict(list)
    for hit in hits:
        per_contig[hit.contig].append(hit)

    kept_hits: List[Hit] = []
    kept_contigs: List[str] = []

    for contig, contig_hits in per_contig.items():
        seed_ids = set(hit.model for hit in contig_hits)

        if combine_mode == "any":
            contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
            kept_hits.extend(contig_hits[:top_per_contig])
            kept_contigs.append(contig)
            continue

        if combine_mode == "all":
            if len(seed_ids) != total_seeds:
                continue
            hits_per_seed: Dict[str, List[Hit]] = defaultdict(list)
            for hit in contig_hits:
                hits_per_seed[hit.model].append(hit)
            for seed_hits in hits_per_seed.values():
                seed_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                kept_hits.extend(seed_hits[:1])
            kept_contigs.append(contig)
            continue

        # threshold mode
        if len(seed_ids) >= min_seed_hits:
            contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
            kept_hits.extend(contig_hits[:top_per_contig])
            kept_contigs.append(contig)

    return kept_hits, kept_contigs


def _write_iteration_summary(path: Path, rows: List[Dict[str, object]]) -> None:
    path.write_text("seed_id\titeration\tn_hits\tn_included\tconverged\n")
    if not rows:
        return
    with path.open("a") as out:
        for row in rows:
            out.write(
                f"{row['seed_id']}\t{row['iteration']}\t{row['n_hits']}\t{row['n_included']}\t{row['converged']}\n"
            )


def _write_hits_table(path: Path, hits: List[Hit]) -> None:
    path.write_text("contig\tprotein\tseed_id\tbitscore\tevalue\n")
    if not hits:
        return
    with path.open("a") as out:
        for hit in hits:
            out.write(
                f"{hit.contig}\t{hit.prot_id}\t{hit.model}\t{hit.bitscore:.4f}\t{hit.evalue:.6g}\n"
            )


def _safe_seed_filename(seed_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", seed_id)


def _jack(cfg: JackConfig) -> None:
    """Run seed-based iterative jackhmmer screening on predicted proteins."""
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")
    if not cfg.seed_marker.exists():
        raise FileNotFoundError(f"Seed marker file not found: {cfg.seed_marker}")

    seqkit_bin = _executable(["seqkit"])
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    proteins_fa = cfg.outdir / "proteins.faa"
    kept_ids = cfg.outdir / "kept_contigs.txt"
    out_contigs = cfg.outdir / "screened_contigs.fasta"
    iter_tsv = cfg.outdir / "jackhmmer_iterations.tsv"
    hits_tsv = cfg.outdir / "jackhmmer_hits.tsv"
    final_hmm_path = cfg.outdir / "last_iteration.hmm"
    hmm_dir = cfg.outdir / "last_iteration_hmms"

    if cfg.save_hmm:
        if final_hmm_path.exists():
            final_hmm_path.unlink()
        if hmm_dir.exists():
            for p in hmm_dir.glob("*.hmm"):
                p.unlink()

    # Use cache-aware protein prediction
    pred_inputs = PredictionInputs(
        input_contigs=cfg.input_contigs,
        mode=cfg.mode,
        min_gene_len=cfg.min_gene_len,
        translation_table=cfg.translation_table,
        min_protein_len_aa=cfg.min_protein_len_aa,
    )
    cache_enabled = os.environ.get("PHU_CACHE", "on") != "off"
    cache_artifact = get_or_predict_proteins(
        pred_inputs,
        use_cache=cache_enabled,
        threads=cfg.threads,
    )

    print(
        f"Predicting proteins with pyrodigal…"
        + (" [cache hit]" if cache_artifact.cache_hit else "")
    )
    print(f"  Proteins predicted: {cache_artifact.protein_count}")

    n_prot = cache_artifact.protein_count
    proteins_fa_actual = cache_artifact.proteins_path

    if n_prot == 0:
        out_contigs.write_text("")
        kept_ids.write_text("")
        _write_iteration_summary(iter_tsv, [])
        _write_hits_table(hits_tsv, [])
        print("No proteins predicted. Exiting with empty outputs.")
        return

    seeds = _read_seed_queries(cfg.seed_marker)
    n_seeds = len(seeds)
    if cfg.combine_mode == "threshold" and cfg.min_seed_hits > n_seeds:
        raise ValueError(
            f"min_seed_hits ({cfg.min_seed_hits}) cannot exceed total seeds ({n_seeds})"
        )

    print(
        f"Running iterative pyhmmer.hmmer.jackhmmer for {n_seeds} seed(s) "
        f"(combine_mode={cfg.combine_mode})…"
    )

    all_hits: List[Hit] = []
    iter_rows: List[Dict[str, object]] = []

    for seed_id, query, alphabet in seeds:
        print(f"  Seed: {seed_id}")
        hmm_output_path: Optional[Path] = None
        if cfg.save_hmm:
            if n_seeds == 1:
                hmm_output_path = final_hmm_path
            else:
                hmm_dir.mkdir(parents=True, exist_ok=True)
                hmm_output_path = hmm_dir / f"{_safe_seed_filename(seed_id)}.hmm"

        seed_hits, seed_iter_rows = _run_jackhmmer(
            query=query,
            alphabet=alphabet,
            seed_id=seed_id,
            proteins_fa=proteins_fa_actual,
            iterations=cfg.iterations,
            inc_evalue=cfg.inc_evalue,
            max_evalue=cfg.max_evalue,
            threads=cfg.threads,
            hmm_output_path=hmm_output_path,
        )
        print(f"    Included final hits: {len(seed_hits)}")
        all_hits.extend(seed_hits)
        for row in seed_iter_rows:
            row["seed_id"] = seed_id
        iter_rows.extend(seed_iter_rows)

    kept_hits, contig_ids = _choose_top_hits_per_contig(
        all_hits,
        top_per_contig=cfg.top_per_contig,
        combine_mode=cfg.combine_mode,
        min_seed_hits=cfg.min_seed_hits,
        total_seeds=n_seeds,
    )

    kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))
    _write_iteration_summary(iter_tsv, iter_rows)
    _write_hits_table(hits_tsv, kept_hits)

    print(f"Extracting {len(contig_ids)} contig(s) with seqkit…")
    _seqkit_extract(cfg.input_contigs, contig_ids, out_contigs, seqkit_bin)

    # Output handling: copy proteins to output folder if requested
    if cfg.keep_proteins:
        shutil.copy(proteins_fa_actual, proteins_fa)
        write_prediction_metadata(
            cfg.outdir / ".phu_prediction_metadata.json",
            cache_hit=cache_artifact.cache_hit,
            cache_key=cache_artifact.cache_key,
            cache_dir=cache_artifact.cache_dir,
        )

    print(f"Done. Output FASTA: {out_contigs}")
    files_msg = "Also wrote: kept_contigs.txt, jackhmmer_hits.tsv, jackhmmer_iterations.tsv"
    if cfg.keep_proteins:
        if cache_artifact.cache_hit:
            files_msg += ", proteins.faa (cached)"
        else:
            files_msg += ", proteins.faa"
    if cfg.save_hmm and final_hmm_path.exists():
        files_msg += ", last_iteration.hmm"
    elif cfg.save_hmm and hmm_dir.exists():
        files_msg += ", last_iteration_hmms/*.hmm"
    print(f"{files_msg}.")
