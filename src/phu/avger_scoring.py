from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

from .avger_annotation import AnnotationHit
from .vscore_db import VScoreRecord

FLANK_DISTANCE_BP = 10_000


@dataclass(frozen=True)
class DatabaseAVL:
    database: str
    contig_id: str
    score: Optional[float]
    denominator: int
    scored_proteins: int


@dataclass(frozen=True)
class FlankEvidence:
    upstream_supported: bool
    downstream_supported: bool
    nearest_upstream_distance: Optional[int]
    nearest_downstream_distance: Optional[int]
    upstream_edge_incomplete: bool
    downstream_edge_incomplete: bool
    flank_supported: bool
    reason_codes: tuple[str, ...]


@dataclass(frozen=True)
class CandidateEvaluation:
    protein_id: str
    contig_id: str
    database: str
    gene_v_score: Optional[float]
    gene_vl_score: Optional[float]
    contig_avl_score: Optional[float]
    candidate: bool
    evidence_state: str
    flank: FlankEvidence


def calculate_database_avl(
    hits: Iterable[AnnotationHit],
    vscore_by_accession: dict[str, VScoreRecord],
) -> dict[tuple[str, str], DatabaseAVL]:
    """Calculate AVL per contig/database using only scored significant best hits."""
    grouped: dict[tuple[str, str], list[float]] = {}
    for hit in hits:
        score = vscore_by_accession.get(hit.model_id)
        if score is not None:
            grouped.setdefault((hit.contig_id, hit.database), []).append(score.vl_score)

    output: dict[tuple[str, str], DatabaseAVL] = {}
    for key, values in grouped.items():
        contig_id, database = key
        output[key] = DatabaseAVL(
            database=database,
            contig_id=contig_id,
            score=sum(values) / len(values),
            denominator=len(values),
            scored_proteins=len(values),
        )
    return output


def _flank_evidence(
    candidate: AnnotationHit,
    database_hits: list[AnnotationHit],
    vscore_by_accession: dict[str, VScoreRecord],
    contig_length: Optional[int],
    require_flank_support: bool,
) -> FlankEvidence:
    upstream: list[int] = []
    downstream: list[int] = []
    for hit in database_hits:
        if hit.protein_id == candidate.protein_id:
            continue
        score = vscore_by_accession.get(hit.model_id)
        if score is None or score.v_score != 10:
            continue
        if hit.target_to is None or hit.target_from is None:
            continue
        candidate_start = candidate.gene_start if candidate.gene_start is not None else candidate.target_from
        candidate_end = candidate.gene_end if candidate.gene_end is not None else candidate.target_to
        hit_start = hit.gene_start if hit.gene_start is not None else hit.target_from
        hit_end = hit.gene_end if hit.gene_end is not None else hit.target_to
        if candidate_start is None or candidate_end is None or hit_start is None or hit_end is None:
            continue
        if hit_end <= candidate_start:
            distance = candidate_start - hit_end
            if distance <= FLANK_DISTANCE_BP:
                upstream.append(distance)
        elif hit_start >= candidate_end:
            distance = hit_start - candidate_end
            if distance <= FLANK_DISTANCE_BP:
                downstream.append(distance)

    candidate_start = candidate.gene_start if candidate.gene_start is not None else candidate.target_from
    candidate_end = candidate.gene_end if candidate.gene_end is not None else candidate.target_to
    upstream_edge_incomplete = contig_length is not None and candidate_start is not None and candidate_start <= FLANK_DISTANCE_BP
    downstream_edge_incomplete = (
        contig_length is not None
        and candidate_end is not None
        and contig_length - candidate_end <= FLANK_DISTANCE_BP
    )
    upstream_supported = bool(upstream)
    downstream_supported = bool(downstream)
    reasons: list[str] = []
    if upstream_edge_incomplete:
        reasons.append("upstream_edge_incomplete")
    if downstream_edge_incomplete:
        reasons.append("downstream_edge_incomplete")
    if not upstream_supported:
        reasons.append("missing_upstream_support")
    if not downstream_supported:
        reasons.append("missing_downstream_support")
    flank_supported = upstream_supported and downstream_supported
    if not flank_supported and not require_flank_support:
        reasons.append("flank_not_required")

    return FlankEvidence(
        upstream_supported=upstream_supported,
        downstream_supported=downstream_supported,
        nearest_upstream_distance=min(upstream) if upstream else None,
        nearest_downstream_distance=min(downstream) if downstream else None,
        upstream_edge_incomplete=upstream_edge_incomplete,
        downstream_edge_incomplete=downstream_edge_incomplete,
        flank_supported=flank_supported,
        reason_codes=tuple(reasons),
    )


def evaluate_database_candidates(
    hits: Iterable[AnnotationHit],
    vscore_by_accession: dict[str, VScoreRecord],
    contig_lengths: Optional[dict[str, int]] = None,
    require_flank_support: bool = False,
) -> list[CandidateEvaluation]:
    """Evaluate Pfam and KOfam candidates independently, then add flank evidence."""
    hits_by_database: dict[str, list[AnnotationHit]] = {"pfam": [], "kofam": []}
    for hit in hits:
        if hit.database in hits_by_database:
            hits_by_database[hit.database].append(hit)
    avl = calculate_database_avl(hits, vscore_by_accession)
    evaluations: list[CandidateEvaluation] = []

    for database in ("pfam", "kofam"):
        for hit in sorted(hits_by_database[database], key=lambda item: item.protein_id):
            score = vscore_by_accession.get(hit.model_id)
            database_avl = avl.get((hit.contig_id, database))
            flank = _flank_evidence(
                hit,
                hits_by_database[database],
                vscore_by_accession,
                (contig_lengths or {}).get(hit.contig_id),
                require_flank_support,
            )
            if score is None or database_avl is None:
                candidate = False
                state = "not_assessable"
            else:
                candidate = score.vl_score < 3 and score.v_score < 10 and database_avl.score > 3
                if not candidate:
                    state = "not_candidate"
                elif require_flank_support and not flank.flank_supported:
                    state = "avg_candidate"
                elif flank.flank_supported:
                    state = "context_supported_avg_candidate"
                else:
                    state = "avg_candidate"
            evaluations.append(
                CandidateEvaluation(
                    protein_id=hit.protein_id,
                    contig_id=hit.contig_id,
                    database=database,
                    gene_v_score=None if score is None else score.v_score,
                    gene_vl_score=None if score is None else score.vl_score,
                    contig_avl_score=None if database_avl is None else database_avl.score,
                    candidate=candidate,
                    evidence_state=state,
                    flank=flank,
                )
            )
    return evaluations