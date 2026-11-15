from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Optional

import pyhmmer.easel
import pyhmmer.hmmer
import pyhmmer.plan7

from phu import __version__
from phu.db.bundle import bundle_provenance
from phu.db.phrog import extract_phrog_models, phrogs_with_models
from phu.db.phrog_lookup import Resolved, load_aliases, resolve

from phu.core.log import log_call
from phu.core.gene_prediction import (
    PredictionInputs,
    get_or_predict_proteins,
    write_prediction_metadata,
)
from phu.commands.screen import (
    MATCHED_CONTIG_IDS,
    MATCHED_CONTIGS,
    PREDICTED_PROTEINS,
    PREDICTION_METADATA,
    Hit,
    _read_fasta,
    _extract_contigs,
)


@dataclass
class JackConfig:
    """Configuration for seed-based iterative jackhmmer screening."""

    input_contigs: Path
    seed_marker: Optional[Path] = None
    outdir: Path = Path("phu-jack")
    mode: str = "meta"
    threads: int = 1
    iterations: int = 5
    inc_evalue: float = 1e-3
    max_evalue: Optional[float] = 1e-5
    top_per_contig: int = 1
    min_gene_len: int = 90
    translation_table: int | None = None
    min_protein_len_aa: int = 30
    keep_proteins: bool = False
    save_hmm: bool = False
    combine_mode: str = "any"
    min_seed_hits: int = 1
    # PHROG families used as seeds: queries as for `phu phrogs resolve`
    # (phrog_<n>, alias:<name>, category:<name>, annot:<text>).
    phrog_queries: tuple[str, ...] = ()
    max_phrog_seeds: int = 50

    def __post_init__(self) -> None:
        self.phrog_queries = tuple(self.phrog_queries)
        if self.seed_marker is None and not self.phrog_queries:
            raise ValueError("give a seed FASTA, --phrog queries, or both")
        if self.max_phrog_seeds < 1:
            raise ValueError("max_phrog_seeds must be >= 1")
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
                f"Duplicate seed sequence ID found: {seq_id}. Seed IDs must be unique."
            )
        seen.add(seq_id)
        text_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=seq)
        out.append((seq_id, text_seq.digitize(alphabet), alphabet))

    return out


PhrogSeeds = tuple[list[Resolved], list[str]]


def resolve_phrog_seeds(queries: tuple[str, ...], cap: int) -> PhrogSeeds:
    """PHROGs selected by ``queries``, refused when more than ``cap``.

    Each PHROG seed is one iterative search, and one category can hold over
    a thousand PHROGs, so a large selection stops before any work is done.
    """
    selected, notes = resolve(queries, phrogs_with_models(), load_aliases())
    if len(selected) > cap:
        raise ValueError(
            f"--phrog selects {len(selected)} PHROGs, more than --max-phrog-seeds "
            f"{cap}; each is one iterative search. Narrow the query (see "
            "`phu phrogs resolve ...`) or raise --max-phrog-seeds."
        )
    return selected, notes


def _load_phrog_hmms(phrog_ids: list[str]) -> list[tuple[str, object]]:
    """The PHROG models for ``phrog_ids``, in that order."""
    with tempfile.TemporaryDirectory(prefix="phu-jack-phrog-") as tmp:
        paths, missing = extract_phrog_models(phrog_ids, Path(tmp))
        if missing:
            raise ValueError(f"PHROG models not found: {', '.join(missing)}")
        models = []
        for phrog_id, path in zip(phrog_ids, paths):
            with pyhmmer.plan7.HMMFile(path) as hmm_file:
                models.append((phrog_id, hmm_file.read()))
    return models


def _read_targets(proteins_fa: Path, alphabet) -> "pyhmmer.easel.DigitalSequenceBlock":
    """Load the predicted proteins once; every seed searches the same block."""
    with pyhmmer.easel.SequenceFile(
        str(proteins_fa), digital=True, alphabet=alphabet
    ) as seq_file:
        return seq_file.read_block()


def _run_jackhmmer(
    query,
    targets,
    seed_id: str,
    iterations: int,
    inc_evalue: float,
    max_evalue: Optional[float],
    threads: int,
    hmm_output_path: Optional[Path] = None,
) -> tuple[list[Hit], list[dict[str, object]]]:
    """Iterate from one seed and return final hits plus an iteration summary.

    A protein seed runs pyhmmer.hmmer.jackhmmer. A PHROG seed is a profile:
    the first round searches with the PHROG model itself, later rounds with
    models rebuilt from the included hits (pyhmmer Pipeline.iterate_hmm).
    """
    if isinstance(query, pyhmmer.plan7.HMM):
        iterations_out = _iterate_from_hmm(
            query, targets, seed_id, iterations, inc_evalue, max_evalue
        )
    else:
        log_call(
            "pyhmmer.hmmer.jackhmmer",
            seed=seed_id,
            queries=1,
            sequences=len(targets),
            max_iterations=iterations,
            incE=inc_evalue,
            incdomE=inc_evalue,
            cpus=threads,
            checkpoints=True,
            max_evalue=max_evalue,
        )
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
        # One query sequence gives one result: with checkpoints=True, the
        # list of its IterationResult objects.
        iterations_out = results[0] if results else []

    if not iterations_out:
        return [], []

    def _iter_hits(iteration_obj):
        """Return hits for either pyhmmer iteration objects or plain list outputs."""
        return getattr(iteration_obj, "hits", iteration_obj)

    def _iter_meta(iteration_obj, idx: int) -> tuple[int, bool]:
        """Return (iteration_index, converged) with sensible fallbacks."""
        iteration = getattr(iteration_obj, "iteration", idx)
        converged = getattr(iteration_obj, "converged", False)
        return int(iteration), bool(converged)

    summary: list[dict[str, object]] = []
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
    kept: list[Hit] = []
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


def _iterate_from_hmm(
    hmm, targets, seed_id: str, iterations: int, inc_evalue: float, max_evalue
) -> list:
    """Rounds of an iterative search started from a profile HMM."""
    log_call(
        "pyhmmer.plan7.Pipeline.iterate_hmm",
        seed=seed_id,
        sequences=len(targets),
        max_iterations=iterations,
        incE=inc_evalue,
        incdomE=inc_evalue,
        max_evalue=max_evalue,
    )
    alphabet = hmm.alphabet
    pipeline = pyhmmer.plan7.Pipeline(
        alphabet,
        background=pyhmmer.plan7.Background(alphabet),
        incE=inc_evalue,
        incdomE=inc_evalue,
    )
    # iterate_hmm only accepts a builder with "hand" architecture.
    builder = pyhmmer.plan7.Builder(alphabet, architecture="hand")
    rounds = []
    search = pipeline.iterate_hmm(hmm, targets, builder=builder)
    while len(rounds) < iterations:
        try:
            result = next(search)
        except StopIteration:
            break
        except ValueError as error:
            # pyhmmer cannot rebuild a model when a round includes nothing,
            # and says so before returning that round. Report it as jackhmmer
            # does for a sequence seed: a last, converged round with no
            # included hits.
            if "No included domains" not in str(error):
                raise
            fresh = pyhmmer.plan7.Pipeline(
                alphabet,
                background=pyhmmer.plan7.Background(alphabet),
                incE=inc_evalue,
                incdomE=inc_evalue,
            )
            last = rounds[-1].hmm if rounds else hmm
            rounds.append(
                SimpleNamespace(
                    hits=fresh.search_hmm(last, targets),
                    iteration=len(rounds) + 1,
                    converged=True,
                    hmm=last,
                )
            )
            break
        rounds.append(result)
        if result.converged:
            break
    return rounds


def _choose_top_hits_per_contig(
    hits: list[Hit],
    top_per_contig: int,
    combine_mode: str,
    min_seed_hits: int,
    total_seeds: int,
) -> tuple[list[Hit], list[str]]:
    """Combine seed hits and keep top hits per contig based on combine mode."""
    per_contig: dict[str, list[Hit]] = defaultdict(list)
    for hit in hits:
        per_contig[hit.contig].append(hit)

    kept_hits: list[Hit] = []
    kept_contigs: list[str] = []

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
            hits_per_seed: dict[str, list[Hit]] = defaultdict(list)
            for hit in contig_hits:
                hits_per_seed[hit.model].append(hit)
            for seed_hits in hits_per_seed.values():
                seed_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                kept_hits.extend(seed_hits[:top_per_contig])
            kept_contigs.append(contig)
            continue

        # threshold mode
        if len(seed_ids) >= min_seed_hits:
            contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
            kept_hits.extend(contig_hits[:top_per_contig])
            kept_contigs.append(contig)

    return kept_hits, kept_contigs


def _write_iteration_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.write_text("seed_id\tseed_source\titeration\tn_hits\tn_included\tconverged\n")
    if not rows:
        return
    with path.open("a") as out:
        for row in rows:
            out.write(
                f"{row['seed_id']}\t{row.get('seed_source', 'fasta')}\t{row['iteration']}"
                f"\t{row['n_hits']}\t{row['n_included']}\t{row['converged']}\n"
            )


def _write_hits_table(
    path: Path, hits: list[Hit], sources: Optional[dict[str, str]] = None
) -> None:
    path.write_text("contig\tprotein\tseed_id\tseed_source\tbitscore\tevalue\n")
    if not hits:
        return
    sources = sources or {}
    with path.open("a") as out:
        for hit in hits:
            out.write(
                f"{hit.contig}\t{hit.prot_id}\t{hit.model}\t{sources.get(hit.model, 'fasta')}"
                f"\t{hit.bitscore:.4f}\t{hit.evalue:.6g}\n"
            )


def _write_query_record(
    path: Path, cfg: "JackConfig", phrog_seeds: list[Resolved], notes: list[str]
) -> None:
    """Record what the seeds were, as `phu screen` does in .phu/query.json."""
    record: dict[str, object] = {
        "phu_version": __version__,
        "command": "jack",
        "seed_fasta": str(cfg.seed_marker) if cfg.seed_marker else None,
        "phrog_queries": list(cfg.phrog_queries),
        "phrog_seeds": [
            {
                "id": r.phrog.phrog_id,
                "annotation": r.phrog.annotation,
                "category": r.phrog.category,
                "selected_by": r.query,
            }
            for r in phrog_seeds
        ],
        "resolve_notes": notes,
        "iterations": cfg.iterations,
        "inc_evalue": cfg.inc_evalue,
        "max_evalue": cfg.max_evalue,
        "combine_mode": cfg.combine_mode,
    }
    if phrog_seeds:
        record["phu_db"] = bundle_provenance()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n")


def _safe_seed_filename(seed_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", seed_id)


def _jack(cfg: JackConfig) -> None:
    """Run seed-based iterative jackhmmer screening on predicted proteins."""
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")
    if cfg.seed_marker is not None and not cfg.seed_marker.exists():
        raise FileNotFoundError(f"Seed marker file not found: {cfg.seed_marker}")

    # PHROG queries are checked before any slow work, so a typo or an
    # oversized selection fails at once.
    phrog_seeds: list[Resolved] = []
    phrog_notes: list[str] = []
    if cfg.phrog_queries:
        phrog_seeds, phrog_notes = resolve_phrog_seeds(
            cfg.phrog_queries, cfg.max_phrog_seeds
        )
        for note in phrog_notes:
            print(f"Note: {note}")

    cfg.outdir.mkdir(parents=True, exist_ok=True)
    _write_query_record(
        cfg.outdir / ".phu" / "query.json", cfg, phrog_seeds, phrog_notes
    )

    proteins_fa = cfg.outdir / PREDICTED_PROTEINS
    kept_ids = cfg.outdir / MATCHED_CONTIG_IDS
    out_contigs = cfg.outdir / MATCHED_CONTIGS
    iter_tsv = cfg.outdir / "jackhmmer-iterations.tsv"
    hits_tsv = cfg.outdir / "jackhmmer-hits.tsv"
    final_hmm_path = cfg.outdir / "final-iteration.hmm"
    hmm_dir = cfg.outdir / "final-iteration-hmms"

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
        "Predicting proteins with pyrodigal…"
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
        if cache_artifact.temp_dir is not None:
            shutil.rmtree(cache_artifact.temp_dir, ignore_errors=True)
        print("No proteins predicted. Exiting with empty outputs.")
        return

    alphabet = pyhmmer.easel.Alphabet.amino()
    seeds: list[tuple[str, object]] = []
    sources: dict[str, str] = {}
    if cfg.seed_marker is not None:
        for seed_id, query, _ in _read_seed_queries(cfg.seed_marker):
            seeds.append((seed_id, query))
            sources[seed_id] = "fasta"
    if phrog_seeds:
        clash = [r.phrog.phrog_id for r in phrog_seeds if r.phrog.phrog_id in sources]
        if clash:
            raise ValueError(
                f"Seed FASTA records named like PHROG seeds: {', '.join(clash)}. "
                "Rename them so every seed ID is unique."
            )
        print(f"Loading {len(phrog_seeds)} PHROG seed model(s)…")
        for phrog_id, hmm in _load_phrog_hmms([r.phrog.phrog_id for r in phrog_seeds]):
            seeds.append((phrog_id, hmm))
            sources[phrog_id] = "phrog"
    n_seeds = len(seeds)
    if cfg.combine_mode == "threshold" and cfg.min_seed_hits > n_seeds:
        raise ValueError(
            f"min_seed_hits ({cfg.min_seed_hits}) cannot exceed total seeds ({n_seeds})"
        )

    n_phrog = sum(1 for source in sources.values() if source == "phrog")
    print(
        f"Running iterative searches for {n_seeds} seed(s)"
        + (f" ({n_phrog} PHROG)" if n_phrog else "")
        + f", up to {cfg.iterations} iteration(s) each (combine_mode={cfg.combine_mode})…"
    )

    all_hits: list[Hit] = []
    iter_rows: list[dict[str, object]] = []

    targets = _read_targets(proteins_fa_actual, alphabet)
    for seed_id, query in seeds:
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
            targets=targets,
            seed_id=seed_id,
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
            row["seed_source"] = sources[seed_id]
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
    _write_hits_table(hits_tsv, kept_hits, sources)

    print(f"Extracting {len(contig_ids)} contig(s)…")
    _extract_contigs(cfg.input_contigs, contig_ids, out_contigs)

    # Output handling: copy proteins to output folder if requested
    if cfg.keep_proteins:
        shutil.copy(proteins_fa_actual, proteins_fa)
        write_prediction_metadata(
            cfg.outdir / PREDICTION_METADATA,
            cache_hit=cache_artifact.cache_hit,
            cache_key=cache_artifact.cache_key,
            cache_dir=cache_artifact.cache_dir,
        )

    # Clean up temp prediction directory (only set when caching is disabled)
    if cache_artifact.temp_dir is not None:
        shutil.rmtree(cache_artifact.temp_dir, ignore_errors=True)

    print(f"Done. Output FASTA: {out_contigs}")
    files_msg = f"Also wrote: {kept_ids.name}, {hits_tsv.name}, {iter_tsv.name}"
    if cfg.keep_proteins:
        files_msg += f", {proteins_fa.name}"
        if cache_artifact.cache_hit:
            files_msg += " (cached)"
    if cfg.save_hmm and final_hmm_path.exists():
        files_msg += f", {final_hmm_path.name}"
    elif cfg.save_hmm and hmm_dir.exists():
        files_msg += f", {hmm_dir.name}/*.hmm"
    print(f"{files_msg}.")
