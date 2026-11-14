from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

from phu import __version__
from phu.avger.reference import load_avg_reference
from phu.avger.tracks import AvgAnnotationResults, annotate_avg_tracks
from phu.avger.decisions import (
    AvgEvidence,
    FilterEvidence,
    PositiveEvidence,
    apply_amg_weight,
    apply_class_filters,
    calculate_scaffold_averages,
    collect_positive_evidence,
    evaluate_candidates,
    index_positive_evidence,
    resolve_final_class,
)
from phu.core.gene_prediction import (
    CacheArtifact,
    PredictionInputs,
    get_or_predict_proteins,
)
from phu.db.binary import build_binary_models, has_binary_models
from phu.db.bundle import bundle_provenance, resolve_db
from phu.db.phrog import ANNOTATIONS_FILE, read_phrog_annotations
from phu.core.progress import ProgressReporter


@dataclass(frozen=True)
class AvgConfig:
    """Validated configuration for the provisional AVG workflow."""

    input_contigs: Path
    output_folder: Path = Path("phu-avger")
    threads: int = 1
    mode: str = "meta"
    min_gene_len: int = 90
    min_protein_len_aa: int = 30
    translation_table: int | None = None
    min_amg_weight: float = 0.6
    filter_mode: str = "standard"
    keep_hits: bool = False
    scaffold_avl_cutoff: float = 3.0
    gene_vl_cutoff: float = 3.0
    gene_v_cutoff: float = 10.0
    scoring_evalue: float = 1e-5
    phrog: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "input_contigs", Path(self.input_contigs))
        object.__setattr__(self, "output_folder", Path(self.output_folder))
        if not self.input_contigs.exists():
            raise FileNotFoundError(f"Input contigs not found: {self.input_contigs}")
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        if self.min_gene_len < 1:
            raise ValueError("min_gene_len must be >= 1")
        if self.min_protein_len_aa < 1:
            raise ValueError("min_protein_len_aa must be >= 1")
        if self.translation_table is not None and self.translation_table < 1:
            raise ValueError("translation_table must be >= 1")
        if not 0.0 <= self.min_amg_weight <= 1.0:
            raise ValueError("min_amg_weight must be between 0 and 1")
        if self.filter_mode not in {"standard", "strict", "none"}:
            raise ValueError("filter_mode must be 'standard', 'strict', or 'none'")
        if self.scaffold_avl_cutoff < 0:
            raise ValueError("scaffold_avl_cutoff must be >= 0")
        if self.gene_vl_cutoff < 0:
            raise ValueError("gene_vl_cutoff must be >= 0")
        if self.gene_v_cutoff < 0:
            raise ValueError("gene_v_cutoff must be >= 0")
        if self.scoring_evalue < 0:
            raise ValueError("scoring_evalue must be >= 0")


@dataclass(frozen=True)
class AvgRunResult:
    prediction: CacheArtifact
    annotations: AvgAnnotationResults
    outputs: tuple[Path, ...] = ()


def require_databases(phrog: bool = True) -> None:
    """Install the phu-db bundle unless the databases avger reads are in place."""
    for name in ("pfam", "kofam", "avg") + (("phrog",) if phrog else ()):
        resolve_db(name)


def run_avg(
    config: AvgConfig, reporter: ProgressReporter | None = None
) -> AvgRunResult:
    """Run prediction and both annotation tracks for a validated config."""
    reporter = reporter or ProgressReporter()
    reporter.start_phase("AVG curation")
    database_task = reporter.start_task("Checking AVG databases")
    run_started_at = datetime.now(timezone.utc).isoformat()
    try:
        require_databases(config.phrog)
        reporter.succeed_task(database_task)
        # Normally done by `phu dbs prepare`; a first run without it does it
        # here, as a step of its own so the wait is explained.
        searched = ("pfam", "kofam") + (("phrog",) if config.phrog else ())
        if not all(has_binary_models(name) for name in searched):
            binary_task = reporter.start_task("Building fast-loading models (one time)")
            for name in searched:
                build_binary_models(name)
            reporter.succeed_task(binary_task)
        prediction_task = reporter.start_task("Predicting proteins")
        prediction = get_or_predict_proteins(
            PredictionInputs(
                input_contigs=config.input_contigs,
                mode=config.mode,
                min_gene_len=config.min_gene_len,
                min_protein_len_aa=config.min_protein_len_aa,
                translation_table=config.translation_table,
            ),
            use_cache=True,
            threads=config.threads,
        )
        reporter.succeed_task(prediction_task)
        annotation_task = reporter.start_task("Annotating proteins")
        annotations = annotate_avg_tracks(
            prediction.proteins_path,
            threads=config.threads,
            scoring_evalue=config.scoring_evalue,
            keep_all_hits=config.keep_hits,
            relaxed_all_hits_path=(
                config.output_folder / "detailed" / "hmm-hits-relaxed.tsv.gz"
            )
            if config.keep_hits
            else None,
            strict_all_hits_path=(
                config.output_folder / "detailed" / "hmm-hits-strict.tsv.gz"
            )
            if config.keep_hits
            else None,
            phrog=config.phrog,
        )
        reporter.succeed_task(annotation_task)
        output_task = reporter.start_task("Writing AVG results")
        outputs = write_avg_outputs(
            config, prediction, annotations, run_started_at=run_started_at
        )
        reporter.succeed_task(output_task)
        return AvgRunResult(
            prediction=prediction, annotations=annotations, outputs=outputs
        )
    except BaseException as exc:
        reporter.fail_running_tasks(str(exc))
        raise
    finally:
        reporter.finish()


def _float(value: str | None) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def _input_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", delete=False, dir=path.parent
    ) as handle:
        temporary_path = Path(handle.name)
        handle.write(content)
    temporary_path.replace(path)


DATABASE_ORDER = ("pfam", "kofam", "phrog")

# Column -> meaning, written as the data dictionary in output-info.txt.
PREDICTION_COLUMNS = {
    "contig": "Contig the gene sits on.",
    "protein_id": "Identifier of the predicted protein.",
    "gene_number": "Position of the gene on its contig (1 is the first gene).",
    "avg_class": (
        "putative_amg (metabolic), putative_apg (physiology) or "
        "putative_areg (regulatory)."
    ),
    "function": (
        "Function of the reference entry that supports the class, or, when it "
        "names none, of the gene's best strict annotation."
    ),
    "supporting_databases": (
        "Databases whose reference entries support the class, separated by ';'."
    ),
    "supporting_models": "Accessions of those reference entries.",
    "amg_weight": (
        "Highest AMG weight among the supporting metabolic entries; empty for "
        "the other classes. A weight of the reference function, not a "
        "probability that this gene is an AMG."
    ),
    "candidate_databases": (
        "Databases in which the candidate rule held: contig AVL above the "
        "cutoff while the gene's VL and V scores are below theirs."
    ),
    "gene_v_score": "V-score of the gene's scoring hit in the scoring database.",
    "gene_vl_score": "VL-score of that hit.",
    "contig_avl": (
        "Average VL-score of the contig in the scoring database. The scoring "
        "database is the first candidate database that also supports the "
        "class, else the first candidate database (pfam, kofam, phrog order)."
    ),
    "pfam_id": "Best strict Pfam hit of the gene; empty when none.",
    "kofam_id": "Best strict KOfam hit of the gene; empty when none.",
    "phrog_id": "Best strict PHROG hit of the gene; empty when none.",
}

SUMMARY_COLUMNS = (
    "contig",
    "genes",
    "genes_scored",
    "pfam_avl",
    "pfam_scored_genes",
    "kofam_avl",
    "kofam_scored_genes",
    "phrog_avl",
    "phrog_scored_genes",
    "candidates",
    "predictions",
    "putative_amg",
    "putative_apg",
    "putative_areg",
)

FILTERED_COLUMNS = (
    "contig",
    "protein_id",
    "gene_number",
    "candidate_databases",
    "proposed_classes",
    "amg_weight",
    "status",
    "reason",
)

_FILE_NOTES = {
    "avg-predictions.tsv": "One row per predicted auxiliary viral gene. The main result.",
    "contig-summary.tsv": "One row per contig: genes, AVL per database, candidates and predictions.",
    "avg-proteins.faa": "Amino-acid sequences of the predicted genes.",
    "detailed/gene-annotations.tsv": "One row per predicted gene, with every annotation, score and rule value.",
    "detailed/filtered-candidates.tsv": "Candidates that did not become predictions, with the reason.",
    "detailed/avg-evidence.tsv": "Every piece of evidence behind the tables above, one row each.",
    "detailed/hmm-hits-relaxed.tsv.gz": "With --keep-hits: every relaxed-track HMM hit.",
    "detailed/hmm-hits-strict.tsv.gz": "With --keep-hits: every strict-track HMM hit.",
    ".phu/run.json": "Versions, parameters, database provenance and output counts of this run.",
}


def _output_info(
    output_counts: dict[str, int], contig_stats: dict[str, dict[str, int]]
) -> str:
    """The text of output-info.txt: run counts, files and the data dictionary."""
    totals = {
        key: sum(stats[key] for stats in contig_stats.values())
        for key in (
            "genes",
            "candidates",
            "predictions",
            "putative_amg",
            "putative_apg",
            "putative_areg",
        )
    }
    lines = [
        "phu avger results",
        "=================",
        "",
        f"phu version: {__version__}",
        "",
        "This run",
        "--------",
        f"Contigs analyzed:    {len(contig_stats):,}",
        f"Genes analyzed:      {totals['genes']:,}",
        f"Candidate genes:     {totals['candidates']:,}",
        f"Predicted AVGs:      {totals['predictions']:,} "
        f"({totals['putative_amg']} metabolic, {totals['putative_apg']} physiological, "
        f"{totals['putative_areg']} regulatory)",
        "",
        "Files",
        "-----",
    ]
    for name, note in _FILE_NOTES.items():
        if name in output_counts or name.startswith(".phu") or "hmm-hits" in name:
            lines.append(f"{name:36s} {note}")
    lines += [
        "",
        "avg-predictions.tsv columns",
        "---------------------------",
    ]
    width = max(len(name) for name in PREDICTION_COLUMNS)
    for name, meaning in PREDICTION_COLUMNS.items():
        lines.append(f"{name:{width}s}  {meaning}")
    lines += [
        "",
        "A prediction is computational evidence for prioritization, not a",
        "confirmed auxiliary viral gene; see the phu documentation.",
        "",
    ]
    return "\n".join(lines)


def write_avg_outputs(
    config: AvgConfig,
    prediction: CacheArtifact,
    annotations: AvgAnnotationResults,
    *,
    run_started_at: str | None = None,
) -> tuple[Path, ...]:
    """Write the initial deterministic AVG tables and run manifest."""
    config.output_folder.mkdir(parents=True, exist_ok=True)
    reference = load_avg_reference(resolve_db("avg"))
    scores = {(row["database"], row["accession"]): row for row in reference.scores}
    positive_rows = reference.positive
    filter_rows = reference.filters

    relaxed_hits = (
        list(annotations.relaxed.best_pfam_by_protein.values())
        + list(annotations.relaxed.best_kofam_by_protein.values())
        + list(annotations.relaxed.best_phrog_by_protein.values())
    )
    strict_hits = (
        list(annotations.strict.best_pfam_by_protein.values())
        + list(annotations.strict.best_kofam_by_protein.values())
        + list(annotations.strict.best_phrog_by_protein.values())
    )
    # Only the category column needs the annotation table.
    phrog_annotations = (
        read_phrog_annotations(resolve_db("phrog") / ANNOTATIONS_FILE)
        if annotations.strict.best_phrog_by_protein
        else {}
    )
    scoring_evidence: list[AvgEvidence] = []
    for hit in relaxed_hits:
        score = scores.get((hit.database, hit.model_id))
        if (
            score is not None
            and _float(score.get("v_score")) is not None
            and _float(score.get("vl_score")) is not None
        ):
            scoring_evidence.append(
                AvgEvidence(
                    protein_id=hit.protein_id,
                    contig_id=hit.contig_id,
                    database=hit.database,
                    accession=hit.model_id,
                    track="relaxed",
                    v_score=float(score["v_score"]),
                    vl_score=float(score["vl_score"]),
                    model_name=hit.model_name or "",
                    model_description=hit.model_description or "",
                )
            )
    decisions = evaluate_candidates(
        scoring_evidence,
        scaffold_avl_cutoff=config.scaffold_avl_cutoff,
        gene_vl_cutoff=config.gene_vl_cutoff,
        gene_v_cutoff=config.gene_v_cutoff,
    )
    decisions_by_protein: dict[str, list] = {}
    for decision in decisions:
        decisions_by_protein.setdefault(decision.protein_id, []).append(decision)

    strict_evidence = [
        AvgEvidence(
            protein_id=hit.protein_id,
            contig_id=hit.contig_id,
            database=hit.database,
            accession=hit.model_id,
            track="strict",
            v_score=0.0,
            vl_score=0.0,
            model_name=hit.model_name or "",
            model_description=hit.model_description or "",
        )
        for hit in strict_hits
    ]
    positive_evidence = [
        PositiveEvidence(
            database=row["database"],
            accession=row["accession"],
            proposed_class=row["proposed_class"],
            name=row.get("name", ""),
            amg_weight=_float(row.get("amg_weight")),
            amg_level=row.get("amg_level") or None,
        )
        for row in positive_rows
    ]
    filter_evidence = [
        FilterEvidence(
            database=row["database"],
            accession=row["accession"],
            proposed_class=row["proposed_class"],
            category=row["filter_category"],
        )
        for row in filter_rows
    ]

    # Looked up once per gene: index them instead of scanning the full lists,
    # which hold tens of thousands of rows once PHROG is included.
    positive_index = index_positive_evidence(positive_evidence)
    filters_by_key: dict[tuple[str, str], list[tuple[int, FilterEvidence]]] = {}
    for position, item in enumerate(filter_evidence):
        filters_by_key.setdefault((item.database, item.accession), []).append(
            (position, item)
        )

    rows = []
    predictions: list[dict[str, object]] = []
    filtered_rows: list[dict[str, object]] = []
    contig_stats: dict[str, dict[str, int]] = {}
    predicted_headers: dict[str, str] = {}
    audit_rows: list[dict[str, object]] = []
    for record in scoring_evidence:
        audit_rows.append({**record.__dict__, "evidence_type": "relaxed_best"})
    for record in strict_evidence:
        audit_rows.append({**record.__dict__, "evidence_type": "strict_hit"})

    scaffold_averages = calculate_scaffold_averages(scoring_evidence)
    strict_by_protein: dict[str, dict[str, object]] = {}
    for hit in strict_hits:
        strict_by_protein.setdefault(hit.protein_id, {})[hit.database] = hit
    strict_evidence_by_protein: dict[str, list[object]] = {}
    for record in strict_evidence:
        strict_evidence_by_protein.setdefault(record.protein_id, []).append(record)
    scoring_by_protein: dict[str, dict[str, object]] = {}
    for record in scoring_evidence:
        scoring_by_protein.setdefault(record.protein_id, {})[record.database] = record

    for gene in prediction.genes or []:
        gene_decisions = decisions_by_protein.get(gene.gene_id, [])
        database_decisions = {
            decision.database: decision for decision in gene_decisions
        }
        pfam_decision = database_decisions.get("pfam")
        kofam_decision = database_decisions.get("kofam")
        phrog_decision = database_decisions.get("phrog")
        is_candidate = any(decision.candidate for decision in gene_decisions)
        gene_strict_evidence = strict_evidence_by_protein.get(gene.gene_id, [])
        matched_positive = (
            collect_positive_evidence(gene_strict_evidence, positive_index)
            if is_candidate
            else ()
        )
        strict_keys = {
            (record.database, record.accession) for record in gene_strict_evidence
        }
        gene_filter_evidence = [
            item
            for _, item in sorted(
                (entry for key in strict_keys for entry in filters_by_key.get(key, ())),
                key=lambda entry: entry[0],
            )
        ]
        for item in matched_positive:
            audit_rows.append(
                {
                    **item.__dict__,
                    "evidence_type": "positive_reference",
                    "protein_id": gene.gene_id,
                }
            )
        for item in gene_filter_evidence:
            audit_rows.append(
                {
                    **item.__dict__,
                    "evidence_type": "filter_reference",
                    "protein_id": gene.gene_id,
                }
            )
        weighted = apply_amg_weight(
            matched_positive, minimum_weight=config.min_amg_weight
        )
        filtered = apply_class_filters(
            weighted.supported_classes,
            gene_filter_evidence,
            filter_mode=config.filter_mode,
        )
        final = resolve_final_class(
            zhou_candidate=is_candidate,
            weighted_classes=weighted.supported_classes,
            filter_decision=filtered,
            below_weight=weighted.below_weight,
        )
        strict_by_database = strict_by_protein.get(gene.gene_id, {})
        scoring_by_database = scoring_by_protein.get(gene.gene_id, {})
        row = {
            "contig": gene.contig_id,
            "protein_id": gene.gene_id,
            "gene_number": gene.ordinal,
            "pfam_id": ""
            if (hit := strict_by_database.get("pfam")) is None
            else hit.model_id,
            "pfam_name": "" if hit is None else hit.model_name or "",
            "pfam_description": "" if hit is None else hit.model_description or "",
            "pfam_bitscore": "" if hit is None else hit.full_score,
            "pfam_evalue": "" if hit is None else hit.evalue,
            "pfam_threshold_method": "" if hit is None else hit.threshold_source,
            "kofam_id": ""
            if (hit := strict_by_database.get("kofam")) is None
            else hit.model_id,
            "kofam_name": "" if hit is None else hit.model_name or "",
            "kofam_description": "" if hit is None else hit.model_description or "",
            "kofam_bitscore": "" if hit is None else hit.full_score,
            "kofam_evalue": "" if hit is None else hit.evalue,
            "kofam_threshold": "" if hit is None else hit.threshold_value,
            "kofam_score_type": "" if hit is None else hit.score_type,
            "phrog_id": ""
            if (hit := strict_by_database.get("phrog")) is None
            else hit.model_id,
            "phrog_annotation": "" if hit is None else hit.model_description or "",
            "phrog_category": ""
            if hit is None or (phrog := phrog_annotations.get(hit.model_id)) is None
            else phrog.category,
            "phrog_bitscore": "" if hit is None else hit.full_score,
            "phrog_evalue": "" if hit is None else hit.evalue,
            "phrog_model_coverage": ""
            if hit is None or hit.model_coverage is None
            else round(hit.model_coverage, 4),
            "phrog_protein_coverage": ""
            if hit is None or hit.protein_coverage is None
            else round(hit.protein_coverage, 4),
            "pfam_scoring_id": ""
            if (record := scoring_by_database.get("pfam")) is None
            else record.accession,
            "pfam_v_score": "" if record is None else record.v_score,
            "pfam_vl_score": "" if record is None else record.vl_score,
            "kofam_scoring_id": ""
            if (record := scoring_by_database.get("kofam")) is None
            else record.accession,
            "kofam_v_score": "" if record is None else record.v_score,
            "kofam_vl_score": "" if record is None else record.vl_score,
            "phrog_scoring_id": ""
            if (record := scoring_by_database.get("phrog")) is None
            else record.accession,
            "phrog_v_score": "" if record is None else record.v_score,
            "phrog_vl_score": "" if record is None else record.vl_score,
            "pfam_scaffold_avl": ""
            if pfam_decision is None
            else pfam_decision.scaffold_avl,
            "pfam_scored_gene_count": ""
            if pfam_decision is None
            or (average := scaffold_averages.get((gene.contig_id, "pfam"))) is None
            else average.denominator,
            "kofam_scaffold_avl": ""
            if kofam_decision is None
            else kofam_decision.scaffold_avl,
            "kofam_scored_gene_count": ""
            if kofam_decision is None
            or (average := scaffold_averages.get((gene.contig_id, "kofam"))) is None
            else average.denominator,
            "phrog_scaffold_avl": ""
            if phrog_decision is None
            else phrog_decision.scaffold_avl,
            "phrog_scored_gene_count": ""
            if phrog_decision is None
            or (average := scaffold_averages.get((gene.contig_id, "phrog"))) is None
            else average.denominator,
            "rule_pfam_candidate": ""
            if pfam_decision is None
            else str(pfam_decision.candidate).lower(),
            "rule_kofam_candidate": ""
            if kofam_decision is None
            else str(kofam_decision.candidate).lower(),
            "rule_phrog_candidate": ""
            if phrog_decision is None
            else str(phrog_decision.candidate).lower(),
            "rule_avg_candidate": str(is_candidate).lower(),
            "proposed_classes": ";".join(sorted(weighted.supported_classes)),
            "amg_weight": ""
            if weighted.maximum_amg_weight is None
            else weighted.maximum_amg_weight,
            "classification_status": final.status,
            "avg_class": final.avg_class,
        }
        rows.append(row)
        stats = contig_stats.setdefault(
            gene.contig_id,
            dict.fromkeys(
                (
                    "genes",
                    "genes_scored",
                    "candidates",
                    "predictions",
                    "putative_amg",
                    "putative_apg",
                    "putative_areg",
                ),
                0,
            ),
        )
        stats["genes"] += 1
        if scoring_by_database:
            stats["genes_scored"] += 1
        if is_candidate:
            stats["candidates"] += 1
        candidate_databases = [
            name
            for name in DATABASE_ORDER
            if (decision := database_decisions.get(name)) is not None
            and decision.candidate
        ]
        if final.status == "classified" and final.avg_class:
            stats["predictions"] += 1
            stats[final.avg_class] += 1
            supporting = [
                item
                for item in weighted.evidence
                if item.proposed_class == final.avg_class
                and (
                    item.proposed_class != "putative_amg"
                    or (
                        item.amg_weight is not None
                        and item.amg_weight >= config.min_amg_weight
                    )
                )
            ]
            supporting_databases = {item.database for item in supporting}
            # The rule values reported are those of one database: the first
            # candidate database that also supports the class, else the first
            # candidate database (in pfam, kofam, phrog order).
            scoring_database = next(
                (n for n in candidate_databases if n in supporting_databases),
                candidate_databases[0],
            )
            decision = database_decisions[scoring_database]
            best_strict = max(
                strict_by_database.values(),
                key=lambda hit: hit.full_score,
                default=None,
            )
            function = next((item.name for item in supporting if item.name), "") or (
                ""
                if best_strict is None
                else best_strict.model_description or best_strict.model_name or ""
            )
            predictions.append(
                {
                    "contig": gene.contig_id,
                    "protein_id": gene.gene_id,
                    "gene_number": gene.ordinal,
                    "avg_class": final.avg_class,
                    "function": function,
                    "supporting_databases": ";".join(sorted(supporting_databases)),
                    "supporting_models": ";".join(
                        sorted({item.accession for item in supporting})
                    ),
                    "amg_weight": row["amg_weight"],
                    "candidate_databases": ";".join(candidate_databases),
                    "gene_v_score": decision.gene_v_score,
                    "gene_vl_score": decision.gene_vl_score,
                    "contig_avl": decision.scaffold_avl,
                    "pfam_id": row["pfam_id"],
                    "kofam_id": row["kofam_id"],
                    "phrog_id": row["phrog_id"],
                }
            )
            predicted_headers[gene.gene_id] = f"{final.avg_class} {function}".strip()
        elif is_candidate:
            filtered_rows.append(
                {
                    "contig": gene.contig_id,
                    "protein_id": gene.gene_id,
                    "gene_number": gene.ordinal,
                    "candidate_databases": ";".join(candidate_databases),
                    "proposed_classes": row["proposed_classes"],
                    "amg_weight": row["amg_weight"],
                    "status": final.status,
                    "reason": ";".join(final.reasons),
                }
            )

    fields = (
        list(rows[0])
        if rows
        else [
            "contig",
            "protein_id",
            "gene_number",
            "rule_avg_candidate",
            "classification_status",
            "avg_class",
        ]
    )
    detailed = config.output_folder / "detailed"
    paths: list[Path] = []
    output_counts: dict[str, int] = {}

    def write_table(path: Path, fieldnames: list[str], selected: list) -> None:
        buffer = StringIO()
        writer = csv.DictWriter(
            buffer,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(selected)
        _atomic_write_text(path, buffer.getvalue())
        output_counts[path.relative_to(config.output_folder).as_posix()] = len(selected)
        paths.append(path)

    write_table(
        config.output_folder / "avg-predictions.tsv",
        list(PREDICTION_COLUMNS),
        predictions,
    )

    summary_rows = []
    for contig, stats in contig_stats.items():
        summary = {"contig": contig, **stats}
        for name in DATABASE_ORDER:
            average = scaffold_averages.get((contig, name))
            summary[f"{name}_avl"] = "" if average is None else average.score
            summary[f"{name}_scored_genes"] = (
                "" if average is None else average.denominator
            )
        summary_rows.append(summary)
    write_table(
        config.output_folder / "contig-summary.tsv",
        list(SUMMARY_COLUMNS),
        summary_rows,
    )

    fasta_path = config.output_folder / "avg-proteins.faa"
    fasta_lines = []
    for gene in prediction.genes or []:
        header = predicted_headers.get(gene.gene_id)
        if header is not None and gene.amino_acid_sequence:
            fasta_lines.append(
                f">{gene.gene_id} {header}\n{gene.amino_acid_sequence}\n"
            )
    _atomic_write_text(fasta_path, "".join(fasta_lines))
    output_counts["avg-proteins.faa"] = len(fasta_lines)
    paths.append(fasta_path)

    write_table(detailed / "gene-annotations.tsv", fields, rows)
    write_table(
        detailed / "filtered-candidates.tsv", list(FILTERED_COLUMNS), filtered_rows
    )
    evidence_fields = [
        "protein_id",
        "contig_id",
        "database",
        "accession",
        "track",
        "evidence_type",
        "model_name",
        "model_description",
        "v_score",
        "vl_score",
        "proposed_class",
        "category",
        "amg_weight",
        "amg_level",
    ]
    write_table(detailed / "avg-evidence.tsv", evidence_fields, audit_rows)

    info_path = config.output_folder / "output-info.txt"
    _atomic_write_text(info_path, _output_info(output_counts, contig_stats))
    paths.append(info_path)

    manifest_path = config.output_folder / ".phu" / "run.json"
    completed_at = datetime.now(timezone.utc).isoformat()
    run_manifest = {
        "phu_version": __version__,
        "command": "phu avger",
        "command_line": sys.argv,
        "input": {
            "path": str(config.input_contigs),
            "sha256": _input_sha256(config.input_contigs),
        },
        "prediction_cache": {
            "hit": prediction.cache_hit,
            "key": prediction.cache_key,
        },
        "parameters": {
            "mode": config.mode,
            "threads": config.threads,
            "min_gene_len": config.min_gene_len,
            "min_protein_len_aa": config.min_protein_len_aa,
            "translation_table": config.translation_table,
            "min_amg_weight": config.min_amg_weight,
            "filter_mode": config.filter_mode,
            "keep_hits": config.keep_hits,
            "scaffold_avl_cutoff": config.scaffold_avl_cutoff,
            "gene_vl_cutoff": config.gene_vl_cutoff,
            "gene_v_cutoff": config.gene_v_cutoff,
            "scoring_evalue": config.scoring_evalue,
            "phrog": config.phrog,
        },
        "phu_db": bundle_provenance(),
        "output_schema_version": 3,
        "output_counts": output_counts,
        "run_started_at": run_started_at or completed_at,
        "run_completed_at": completed_at,
        "completion_state": "completed",
    }
    _atomic_write_text(manifest_path, json.dumps(run_manifest, indent=2) + "\n")
    return tuple(paths)
