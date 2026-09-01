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

from . import __version__
from .avg_reference_db import get_avg_database_status
from .avg_annotation import AvgAnnotationResults, annotate_avg_tracks
from .avg_decisions import (
    AvgEvidence,
    FilterEvidence,
    PositiveEvidence,
    apply_amg_weight,
    apply_class_filters,
    calculate_scaffold_averages,
    collect_positive_evidence,
    evaluate_candidates,
    resolve_final_class,
)
from .gene_prediction_core import (
    CacheArtifact,
    PredictionInputs,
    get_or_predict_proteins,
)
from .kofam_db import get_kofam_database_status
from .pfam_db import get_pfam_database_status
from ._click import ProgressReporter


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


def database_readiness() -> dict[str, bool]:
    """Return readiness for the three databases required by AVG."""
    statuses = {
        "pfam": get_pfam_database_status(),
        "kofam": get_kofam_database_status(),
        "avg": get_avg_database_status(),
    }
    return {
        name: bool(status.get("downloaded") and status.get("indexed"))
        for name, status in statuses.items()
    }


def missing_databases() -> list[str]:
    """Return required database names that are not ready."""
    return [name for name, ready in database_readiness().items() if not ready]


def require_databases() -> None:
    """Raise an actionable error when AVG databases are not ready."""
    missing = missing_databases()
    if missing:
        names = ", ".join(missing)
        raise FileNotFoundError(
            f"Required databases are not ready: {names}. "
            "Run: phu dbs prepare pfam kofam avg"
        )


def run_avg(
    config: AvgConfig, reporter: ProgressReporter | None = None
) -> AvgRunResult:
    """Run prediction and both annotation tracks for a validated config."""
    reporter = reporter or ProgressReporter()
    reporter.start_phase("AVG curation")
    database_task = reporter.start_task("Checking AVG databases")
    run_started_at = datetime.now(timezone.utc).isoformat()
    try:
        require_databases()
        reporter.succeed_task(database_task)
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
            relaxed_all_hits_path=(config.output_folder / "relaxed_hits.tsv.gz")
            if config.keep_hits
            else None,
            strict_all_hits_path=(config.output_folder / "strict_hits.tsv.gz")
            if config.keep_hits
            else None,
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


def _reference_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def write_avg_outputs(
    config: AvgConfig,
    prediction: CacheArtifact,
    annotations: AvgAnnotationResults,
    *,
    run_started_at: str | None = None,
) -> tuple[Path, ...]:
    """Write the initial deterministic AVG tables and run manifest."""
    config.output_folder.mkdir(parents=True, exist_ok=True)
    root = get_avg_database_status()["root"]
    reference_root = Path(str(root))
    scores = {}
    for row in _reference_rows(reference_root / "v_scores.tsv"):
        scores[(row["database"], row["accession"])] = row
    positive_rows = _reference_rows(reference_root / "avg_positive.tsv")
    filter_rows = _reference_rows(reference_root / "avg_filters.tsv")

    relaxed_hits = list(annotations.relaxed.best_pfam_by_protein.values()) + list(
        annotations.relaxed.best_kofam_by_protein.values()
    )
    strict_hits = list(annotations.strict.best_pfam_by_protein.values()) + list(
        annotations.strict.best_kofam_by_protein.values()
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

    rows = []
    candidates = []
    predictions = []
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
        is_candidate = any(decision.candidate for decision in gene_decisions)
        gene_strict_evidence = strict_evidence_by_protein.get(gene.gene_id, [])
        matched_positive = (
            collect_positive_evidence(gene_strict_evidence, positive_evidence)
            if is_candidate
            else ()
        )
        strict_keys = {
            (record.database, record.accession) for record in gene_strict_evidence
        }
        gene_filter_evidence = [
            item
            for item in filter_evidence
            if (item.database, item.accession) in strict_keys
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
            "zhou_pfam_candidate": ""
            if pfam_decision is None
            else str(pfam_decision.candidate).lower(),
            "zhou_kofam_candidate": ""
            if kofam_decision is None
            else str(kofam_decision.candidate).lower(),
            "zhou_avg_candidate": str(is_candidate).lower(),
            "proposed_classes": ";".join(sorted(weighted.supported_classes)),
            "amg_weight": ""
            if weighted.maximum_amg_weight is None
            else weighted.maximum_amg_weight,
            "classification_status": final.status,
            "avg_class": final.avg_class,
        }
        rows.append(row)
        if is_candidate:
            candidates.append(row)
        if final.status == "classified" and final.avg_class:
            predictions.append(row)

    fields = (
        list(rows[0])
        if rows
        else [
            "contig",
            "protein_id",
            "gene_number",
            "zhou_avg_candidate",
            "classification_status",
            "avg_class",
        ]
    )
    paths = []
    output_counts = {}
    for filename, selected, selected_fields in (
        ("genes.tsv", rows, fields),
        ("avg_candidates.tsv", candidates, fields),
        ("avg_predictions.tsv", predictions, fields),
    ):
        path = config.output_folder / filename
        buffer = StringIO()
        writer = csv.DictWriter(
            buffer,
            fieldnames=selected_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected)
        _atomic_write_text(path, buffer.getvalue())
        output_counts[filename] = len(selected)
        paths.append(path)

    evidence_path = config.output_folder / "evidence.tsv"
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
    buffer = StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=evidence_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(audit_rows)
    _atomic_write_text(evidence_path, buffer.getvalue())
    output_counts["evidence.tsv"] = len(audit_rows)
    paths.append(evidence_path)

    manifest_path = config.output_folder / ".phu" / "run.json"
    completed_at = datetime.now(timezone.utc).isoformat()
    avg_status = get_avg_database_status()
    reference_manifest_path = reference_root / "manifest.json"
    reference_manifest = json.loads(reference_manifest_path.read_text())
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
        },
        "databases": {
            "pfam": get_pfam_database_status(),
            "kofam": get_kofam_database_status(),
            "avg": avg_status,
        },
        "avg_reference": reference_manifest,
        "output_schema_version": 1,
        "output_counts": output_counts,
        "run_started_at": run_started_at or completed_at,
        "run_completed_at": completed_at,
        "completion_state": "completed",
    }
    _atomic_write_text(manifest_path, json.dumps(run_manifest, indent=2) + "\n")
    return tuple(paths)
