from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pyhmmer

from .kofam_db import KOFamMetadata, ensure_kofam_database, get_all_kofam_metadata
from .pfam_db import ensure_pfam_database, normalize_pfam_id
from .vscore_db import VScoreRecord
from .avger_classification import ClassificationRules, classify_protein_annotations


@dataclass(frozen=True)
class AnnotationHit:
    protein_id: str
    contig_id: str
    database: str
    model_id: str
    model_accession: str
    score_type: str
    effective_score: float
    full_score: float
    domain_score: Optional[float]
    evalue: float
    hmm_from: Optional[int]
    hmm_to: Optional[int]
    target_from: Optional[int]
    target_to: Optional[int]
    threshold_source: str
    threshold_value: Optional[float]
    gene_start: Optional[int] = None
    gene_end: Optional[int] = None


@dataclass(frozen=True)
class AnnotationConfig:
    threads: int = 1
    max_evalue: float = 1e-5
    protein_batch_size: int = 10_000
    keep_all_hits: bool = False
    all_hits_path: Optional[Path] = None
    pfam_require_ga: bool = True
    pfam_missing_ga_policy: str = "skip_model"  # skip_model | include_without_ga

    def __post_init__(self) -> None:
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.max_evalue < 0:
            raise ValueError("max_evalue must be >= 0")
        if self.protein_batch_size < 1:
            raise ValueError("protein_batch_size must be >= 1")
        if self.pfam_missing_ga_policy not in {"skip_model", "include_without_ga"}:
            raise ValueError(
                "pfam_missing_ga_policy must be 'skip_model' or 'include_without_ga'"
            )


@dataclass
class AnnotationResults:
    best_pfam_by_protein: dict[str, AnnotationHit]
    best_kofam_by_protein: dict[str, AnnotationHit]
    passing_hit_count: int
    scanned_model_count: int
    skipped_pfam_models_missing_ga: int


def _resolve_vscore_for_row(
    row: AnnotationHit,
    results: AnnotationResults,
    vscore_by_accession: Optional[dict[str, VScoreRecord]],
) -> Optional[VScoreRecord]:
    if not vscore_by_accession:
        return None

    if row.model_id in vscore_by_accession:
        return vscore_by_accession[row.model_id]

    ko_hit = results.best_kofam_by_protein.get(row.protein_id)
    if ko_hit is not None and ko_hit.model_id in vscore_by_accession:
        return vscore_by_accession[ko_hit.model_id]

    return None


def write_best_hits_tsv(
    results: AnnotationResults,
    output_path: Path,
    vscore_by_accession: Optional[dict[str, VScoreRecord]] = None,
    classification_rules: Optional[ClassificationRules] = None,
    candidate_evaluations=None,
) -> int:
    """Write one deterministic best annotation row per protein and database."""
    rows = list(results.best_pfam_by_protein.values()) + list(
        results.best_kofam_by_protein.values()
    )
    rows.sort(key=lambda row: (row.protein_id, row.database, row.model_accession))
    rows_by_protein: dict[str, list[AnnotationHit]] = {}
    for row in rows:
        rows_by_protein.setdefault(row.protein_id, []).append(row)
    evaluations_by_key = {
        (item.protein_id, item.database): item for item in (candidate_evaluations or [])
    }
    avg_candidate_by_protein: dict[str, bool] = {}
    for item in candidate_evaluations or []:
        avg_candidate_by_protein[item.protein_id] = (
            avg_candidate_by_protein.get(item.protein_id, False) or item.candidate
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "protein_id",
                "contig_id",
                "database",
                "model_id",
                "score_type",
                "effective_score",
                "full_score",
                "domain_score",
                "evalue",
                "hmm_from",
                "hmm_to",
                "target_from",
                "target_to",
                "threshold_source",
                "threshold_value",
                "v_score",
                "vl_score",
                "v_score_function",
                "v_score_log10_hit_number",
                "v_score_database_origin",
                "contig_avl_score",
                "database_candidate",
                "avg_candidate",
                "evidence_state",
                "upstream_support",
                "downstream_support",
                "nearest_upstream_distance",
                "nearest_downstream_distance",
                "flank_supported",
                "flank_reason_codes",
                "classification",
                "classification_rule_id",
                "classification_rule_version",
            ]
        )
        for row in rows:
            vscore = _resolve_vscore_for_row(row, results, vscore_by_accession)
            classification, rule_id, rule_version = classify_protein_annotations(
                rows_by_protein[row.protein_id], classification_rules, vscore_by_accession
            )
            evaluation = evaluations_by_key.get((row.protein_id, row.database))
            writer.writerow(
                [
                    row.protein_id,
                    row.contig_id,
                    row.database,
                    row.model_id,
                    row.score_type,
                    f"{row.effective_score:.6f}",
                    f"{row.full_score:.6f}",
                    "" if row.domain_score is None else f"{row.domain_score:.6f}",
                    f"{row.evalue:.6g}",
                    "" if row.hmm_from is None else row.hmm_from,
                    "" if row.hmm_to is None else row.hmm_to,
                    "" if row.target_from is None else row.target_from,
                    "" if row.target_to is None else row.target_to,
                    row.threshold_source,
                    "" if row.threshold_value is None else f"{row.threshold_value:.6f}",
                    "" if vscore is None else f"{vscore.v_score:.6f}",
                    "" if vscore is None else f"{vscore.vl_score:.6f}",
                    "" if vscore is None else vscore.protein_function,
                    "" if vscore is None else f"{vscore.log10_hit_number:.6f}",
                    "" if vscore is None else vscore.database_origin,
                    "" if evaluation is None or evaluation.contig_avl_score is None else f"{evaluation.contig_avl_score:.6f}",
                    "" if evaluation is None else ("true" if evaluation.candidate else "false"),
                    "true" if avg_candidate_by_protein.get(row.protein_id, False) else "false",
                    "" if evaluation is None else evaluation.evidence_state,
                    "" if evaluation is None else ("true" if evaluation.flank.upstream_supported else "false"),
                    "" if evaluation is None else ("true" if evaluation.flank.downstream_supported else "false"),
                    "" if evaluation is None else evaluation.flank.nearest_upstream_distance,
                    "" if evaluation is None else evaluation.flank.nearest_downstream_distance,
                    "" if evaluation is None else ("true" if evaluation.flank.flank_supported else "false"),
                    "" if evaluation is None else ";".join(evaluation.flank.reason_codes),
                    classification,
                    "" if rule_id is None else rule_id,
                    "" if rule_version is None else rule_version,
                ]
            )
    return len(rows)


def parse_contig_id_from_protein_id(protein_id: str) -> str:
    marker = "|gene"
    idx = protein_id.rfind(marker)
    if idx != -1:
        suffix = protein_id[idx + len(marker) :]
        if suffix.isdigit():
            return protein_id[:idx]
    if "|" in protein_id:
        return protein_id.rsplit("|", 1)[0]
    return protein_id


def _domain_metrics(hit) -> tuple[Optional[float], Optional[int], Optional[int], Optional[int], Optional[int]]:
    best_score: Optional[float] = None
    hmm_from: Optional[int] = None
    hmm_to: Optional[int] = None
    target_from: Optional[int] = None
    target_to: Optional[int] = None

    for domain in getattr(hit, "domains", []):
        if not getattr(domain, "included", True):
            continue
        score = float(domain.score)
        if best_score is None or score > best_score:
            best_score = score
            alignment = getattr(domain, "alignment", None)
            if alignment is not None:
                hmm_from = int(alignment.hmm_from)
                hmm_to = int(alignment.hmm_to)
                target_from = int(alignment.target_from)
                target_to = int(alignment.target_to)

    return best_score, hmm_from, hmm_to, target_from, target_to


def _better_hit(candidate: AnnotationHit, current: Optional[AnnotationHit]) -> bool:
    if current is None:
        return True
    if candidate.effective_score != current.effective_score:
        return candidate.effective_score > current.effective_score
    if candidate.evalue != current.evalue:
        return candidate.evalue < current.evalue
    return candidate.model_accession < current.model_accession


def _normalize_pfam_accession(hmm) -> Optional[str]:
    accession = getattr(hmm, "accession", None)
    if isinstance(accession, bytes):
        accession = accession.decode()
    if not accession:
        return None
    try:
        return normalize_pfam_id(str(accession))
    except ValueError:
        return None


def _kofam_model_id(hmm) -> Optional[str]:
    name = getattr(hmm, "name", None)
    accession = getattr(hmm, "accession", None)

    for token in (name, accession):
        if isinstance(token, bytes):
            token = token.decode()
        if not token:
            continue
        token = str(token).strip().upper()
        if len(token) == 6 and token.startswith("K") and token[1:].isdigit():
            return token
    return None


def _pfam_model_ga(hmm) -> Optional[tuple[float, float]]:
    cutoffs = getattr(hmm, "cutoffs", None)
    if cutoffs is None:
        return None
    if not cutoffs.gathering_available():
        return None
    return float(cutoffs.gathering1), float(cutoffs.gathering2)


def _open_all_hits_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        handle = gzip.open(path, "wt", newline="")
    else:
        handle = path.open("w", newline="")

    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "protein_id",
            "contig_id",
            "database",
            "model_id",
            "model_accession",
            "score_type",
            "effective_score",
            "full_score",
            "domain_score",
            "evalue",
            "hmm_from",
            "hmm_to",
            "target_from",
            "target_to",
            "threshold_source",
            "threshold_value",
        ]
    )
    return handle, writer


def _search_and_collect(
    *,
    database: str,
    hmm_path: Path,
    proteins_path: Path,
    cfg: AnnotationConfig,
    kofam_meta_by_model: Optional[dict[str, KOFamMetadata]] = None,
    all_hits_writer=None,
) -> tuple[dict[str, AnnotationHit], int, int, int]:
    best_by_protein: dict[str, AnnotationHit] = {}
    passing = 0
    models = 0
    seen_models: set[str] = set()
    skipped_missing_ga = 0
    kofam_meta_by_model = kofam_meta_by_model or {}

    with pyhmmer.easel.SequenceFile(proteins_path, digital=True) as seq_file:
        while True:
            proteins = seq_file.read_block(sequences=cfg.protein_batch_size)
            if not proteins:
                break

            # HMMFile is consumed by hmmsearch, so reopen it for each bounded
            # protein batch rather than retaining all proteins in memory.
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hit_batches = pyhmmer.hmmsearch(hmm_file, proteins, cpus=cfg.threads)
                for top_hits in hit_batches:
                    hmm = top_hits.query
                    raw_model = getattr(hmm, "accession", None) or getattr(hmm, "name", "")
                    if isinstance(raw_model, bytes):
                        raw_model = raw_model.decode()
                    model_key = f"{database}:{raw_model}"
                    if model_key not in seen_models:
                        seen_models.add(model_key)
                        models += 1

                    if database == "pfam":
                        model_accession = _normalize_pfam_accession(hmm)
                        if model_accession is None:
                            continue
                        model_id = model_accession

                        ga_pair = _pfam_model_ga(hmm)
                        if cfg.pfam_require_ga and ga_pair is None:
                            if cfg.pfam_missing_ga_policy == "skip_model":
                                skipped_missing_ga += 1
                                continue
                        seq_ga = ga_pair[0] if ga_pair is not None else None
                        dom_ga = ga_pair[1] if ga_pair is not None else None
                    else:
                        model_id = _kofam_model_id(hmm)
                        if model_id is None:
                            continue
                        model_accession = model_id
                        ko_meta = kofam_meta_by_model.get(model_id)
                        if ko_meta is None:
                            continue

                    for hit in top_hits:
                        evalue = float(hit.evalue)
                        if evalue > cfg.max_evalue:
                            continue

                        protein_id = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                        contig_id = parse_contig_id_from_protein_id(protein_id)

                        full_score = float(hit.score)
                        domain_score, hmm_from, hmm_to, target_from, target_to = _domain_metrics(hit)

                        if database == "pfam":
                            threshold_source = "pfam_ga"
                            threshold_value = seq_ga
                            if cfg.pfam_require_ga and seq_ga is not None and full_score < seq_ga:
                                continue
                            if cfg.pfam_require_ga and dom_ga is not None and domain_score is not None:
                                if domain_score < dom_ga:
                                    continue
                            score_type = "full"
                            effective_score = full_score
                        else:
                            ko_meta = kofam_meta_by_model[model_id]
                            score_type = ko_meta.score_type
                            threshold_source = "kofam_ko_list"
                            threshold_value = ko_meta.threshold
                            if score_type == "domain":
                                if domain_score is None:
                                    continue
                                effective_score = domain_score
                            else:
                                effective_score = full_score
                            if threshold_value is None:
                                continue
                            if effective_score < float(threshold_value):
                                continue

                        record = AnnotationHit(
                            protein_id=protein_id,
                            contig_id=contig_id,
                            database=database,
                            model_id=model_id,
                            model_accession=model_accession,
                            score_type=score_type,
                            effective_score=effective_score,
                            full_score=full_score,
                            domain_score=domain_score,
                            evalue=evalue,
                            hmm_from=hmm_from,
                            hmm_to=hmm_to,
                            target_from=target_from,
                            target_to=target_to,
                            threshold_source=threshold_source,
                            threshold_value=threshold_value,
                        )

                        passing += 1
                        if all_hits_writer is not None:
                            all_hits_writer.writerow(
                                [
                                    record.protein_id,
                                    record.contig_id,
                                    record.database,
                                    record.model_id,
                                    record.model_accession,
                                    record.score_type,
                                    f"{record.effective_score:.6f}",
                                    f"{record.full_score:.6f}",
                                    "" if record.domain_score is None else f"{record.domain_score:.6f}",
                                    f"{record.evalue:.6g}",
                                    "" if record.hmm_from is None else str(record.hmm_from),
                                    "" if record.hmm_to is None else str(record.hmm_to),
                                    "" if record.target_from is None else str(record.target_from),
                                    "" if record.target_to is None else str(record.target_to),
                                    record.threshold_source,
                                    "" if record.threshold_value is None else f"{record.threshold_value:.6f}",
                                ]
                            )

                        current = best_by_protein.get(record.protein_id)
                        if _better_hit(record, current):
                            best_by_protein[record.protein_id] = record

    return best_by_protein, passing, models, skipped_missing_ga


def annotate_proteins_complete_databases(
    proteins_faa: Path,
    cfg: AnnotationConfig,
) -> AnnotationResults:
    if cfg.keep_all_hits and cfg.all_hits_path is None:
        raise ValueError("all_hits_path is required when keep_all_hits is True")

    pfam_meta = ensure_pfam_database(force_refresh=False)
    kofam_meta = ensure_kofam_database(force_refresh=False)

    pfam_hmm_path = Path(pfam_meta["hmm_path"])
    kofam_hmm_path = Path(kofam_meta["hmm_path"])
    all_kofam_meta = get_all_kofam_metadata()

    all_hits_handle = None
    all_hits_writer = None
    if cfg.keep_all_hits:
        all_hits_handle, all_hits_writer = _open_all_hits_writer(cfg.all_hits_path)

    try:
        best_pfam, pass_pfam, models_pfam, skipped_ga = _search_and_collect(
            database="pfam",
            hmm_path=pfam_hmm_path,
            proteins_path=proteins_faa,
            cfg=cfg,
            all_hits_writer=all_hits_writer,
        )
        best_kofam, pass_kofam, models_kofam, _ = _search_and_collect(
            database="kofam",
            hmm_path=kofam_hmm_path,
            proteins_path=proteins_faa,
            cfg=cfg,
            kofam_meta_by_model=all_kofam_meta,
            all_hits_writer=all_hits_writer,
        )
    finally:
        if all_hits_handle is not None:
            all_hits_handle.close()

    return AnnotationResults(
        best_pfam_by_protein=best_pfam,
        best_kofam_by_protein=best_kofam,
        passing_hit_count=pass_pfam + pass_kofam,
        scanned_model_count=models_pfam + models_kofam,
        skipped_pfam_models_missing_ga=skipped_ga,
    )