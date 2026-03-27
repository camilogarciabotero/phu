from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ._exec import _executable
from .screen import Hit, _predict_proteins_pyrodigal, _read_fasta, _seqkit_extract


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
    keep_proteins: bool = False

    def __post_init__(self) -> None:
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.iterations < 1:
            raise ValueError("iterations must be >= 1")
        if self.top_per_contig < 1:
            raise ValueError("top_per_contig must be >= 1")


def _read_single_seed_query(seed_marker: Path):
    """
    Read exactly one protein seed sequence from FASTA and digitize it.

    This strict v1 contract avoids ambiguity around multi-seed behavior.
    """
    import pyhmmer.easel

    records = list(_read_fasta(seed_marker))
    if not records:
        raise ValueError(f"Seed marker file is empty: {seed_marker}")
    if len(records) != 1:
        raise ValueError(
            "Seed marker file must contain exactly one sequence in this version. "
            "Provide a FASTA with a single protein sequence."
        )

    seq_id, seq = records[0]
    alphabet = pyhmmer.easel.Alphabet.amino()
    text_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=seq)
    return text_seq.digitize(alphabet), alphabet


def _run_jackhmmer(
    query,
    alphabet,
    proteins_fa: Path,
    iterations: int,
    inc_evalue: float,
    max_evalue: Optional[float],
    threads: int,
) -> Tuple[List[Hit], List[Dict[str, object]]]:
    """Run pyhmmer.hmmer.jackhmmer and return final hits plus iteration summary."""
    import pyhmmer.easel
    import pyhmmer.hmmer

    with pyhmmer.easel.SequenceFile(
        str(proteins_fa), digital=True, alphabet=alphabet
    ) as seq_file:
        targets = seq_file.read_block()

    iterations_out = list(
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

    if not iterations_out:
        return [], []

    summary: List[Dict[str, object]] = []
    for it in iterations_out:
        n_hits = 0
        n_included = 0
        for hit in it.hits:
            n_hits += 1
            if hit.included:
                n_included += 1
        summary.append(
            {
                "iteration": int(it.iteration),
                "n_hits": n_hits,
                "n_included": n_included,
                "converged": bool(it.converged),
            }
        )

    final = iterations_out[-1]
    kept: List[Hit] = []
    gene_pattern = re.compile(r"\|gene\d+$")

    for hit in final.hits:
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
                model="jackhmmer",
                bitscore=hit.score,
                evalue=hit.evalue,
            )
        )

    return kept, summary


def _choose_top_hits_per_contig(
    hits: List[Hit], top_per_contig: int
) -> Tuple[List[Hit], List[str]]:
    """Keep top-N hits per contig by bitscore/evalue."""
    per_contig: Dict[str, List[Hit]] = defaultdict(list)
    for hit in hits:
        per_contig[hit.contig].append(hit)

    kept_hits: List[Hit] = []
    kept_contigs: List[str] = []
    for contig, contig_hits in per_contig.items():
        contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
        kept_hits.extend(contig_hits[:top_per_contig])
        kept_contigs.append(contig)

    return kept_hits, kept_contigs


def _write_iteration_summary(path: Path, rows: List[Dict[str, object]]) -> None:
    path.write_text("iteration\tn_hits\tn_included\tconverged\n")
    if not rows:
        return
    with path.open("a") as out:
        for row in rows:
            out.write(
                f"{row['iteration']}\t{row['n_hits']}\t{row['n_included']}\t{row['converged']}\n"
            )


def _write_hits_table(path: Path, hits: List[Hit]) -> None:
    path.write_text("contig\tprotein\tbitscore\tevalue\n")
    if not hits:
        return
    with path.open("a") as out:
        for hit in hits:
            out.write(
                f"{hit.contig}\t{hit.prot_id}\t{hit.bitscore:.4f}\t{hit.evalue:.6g}\n"
            )


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

    print("Predicting proteins with pyrodigal…")
    n_prot = _predict_proteins_pyrodigal(
        cfg.input_contigs,
        proteins_fa,
        mode=cfg.mode,
        min_len=cfg.min_gene_len,
        translation_table=cfg.translation_table,
        threads=cfg.threads,
    )
    print(f"  Proteins predicted: {n_prot}")

    if n_prot == 0:
        out_contigs.write_text("")
        kept_ids.write_text("")
        _write_iteration_summary(iter_tsv, [])
        _write_hits_table(hits_tsv, [])
        if not cfg.keep_proteins and proteins_fa.exists():
            proteins_fa.unlink()
        print("No proteins predicted. Exiting with empty outputs.")
        return

    query, alphabet = _read_single_seed_query(cfg.seed_marker)

    print("Running iterative pyhmmer.hmmer.jackhmmer…")
    final_hits, iter_rows = _run_jackhmmer(
        query=query,
        alphabet=alphabet,
        proteins_fa=proteins_fa,
        iterations=cfg.iterations,
        inc_evalue=cfg.inc_evalue,
        max_evalue=cfg.max_evalue,
        threads=cfg.threads,
    )
    print(f"  Included final hits: {len(final_hits)}")

    kept_hits, contig_ids = _choose_top_hits_per_contig(final_hits, cfg.top_per_contig)

    kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))
    _write_iteration_summary(iter_tsv, iter_rows)
    _write_hits_table(hits_tsv, kept_hits)

    print(f"Extracting {len(contig_ids)} contig(s) with seqkit…")
    _seqkit_extract(cfg.input_contigs, contig_ids, out_contigs, seqkit_bin)

    if not cfg.keep_proteins and proteins_fa.exists():
        proteins_fa.unlink()

    print(f"Done. Output FASTA: {out_contigs}")
    print(
        "Also wrote: kept_contigs.txt, jackhmmer_hits.tsv, "
        "jackhmmer_iterations.tsv."
    )
