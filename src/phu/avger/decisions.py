from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Literal

Database = Literal["pfam", "kofam"]
Track = Literal["relaxed", "strict"]
FilterMode = Literal["standard", "strict", "none"]
DatabaseKey = tuple[Database, str]


@dataclass(frozen=True)
class AvgEvidence:
    """One database-scoped evidence record for a predicted protein."""

    protein_id: str
    contig_id: str
    database: Database
    accession: str
    track: Track
    v_score: float
    vl_score: float
    model_name: str = ""
    model_description: str = ""


@dataclass(frozen=True)
class ScaffoldAverage:
    contig_id: str
    database: Database
    score: float
    denominator: int


@dataclass(frozen=True)
class CandidateDecision:
    protein_id: str
    contig_id: str
    database: Database
    accession: str
    scaffold_avl: float | None
    gene_vl_score: float
    gene_v_score: float
    candidate: bool
    reason: str


@dataclass(frozen=True)
class CurationDecision:
    classification: str
    weight: float
    positive: bool
    conflicting: bool
    filtered: bool
    reason: str


@dataclass(frozen=True)
class PositiveEvidence:
    database: Database
    accession: str
    proposed_class: str
    name: str = ""
    amg_weight: float | None = None
    amg_level: str | None = None


@dataclass(frozen=True)
class FilterEvidence:
    database: Database
    accession: str
    proposed_class: str
    category: str


@dataclass(frozen=True)
class WeightDecision:
    evidence: tuple[PositiveEvidence, ...]
    supported_classes: frozenset[str]
    maximum_amg_weight: float | None
    below_weight: bool
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class FilterDecision:
    surviving_classes: frozenset[str]
    blocking_matches: tuple[FilterEvidence, ...]
    warning_matches: tuple[FilterEvidence, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class FinalClassDecision:
    avg_class: str
    status: str
    reasons: tuple[str, ...]


SUPPORTED_CLASSES = frozenset({"putative_amg", "putative_apg", "putative_areg"})
FILTER_CATEGORIES = frozenset(
    {"filter_essential", "filter_glucan", "filter_nucleotide", "filter_methyl", "filter_lipid"}
)


def database_key(database: str, accession: str) -> DatabaseKey:
    """Return a normalized database/accession key for joins."""
    normalized_database = database.strip().lower()
    if normalized_database in {"kegg", "ko", "kofam"}:
        normalized_database = "kofam"
    if normalized_database not in {"pfam", "kofam"}:
        raise ValueError(f"Unsupported AVG database: {database}")
    normalized_accession = accession.strip().upper()
    if not normalized_accession:
        raise ValueError("AVG accession must not be empty")
    return normalized_database, normalized_accession  # type: ignore[return-value]


def calculate_scaffold_averages(
    evidence: Iterable[AvgEvidence],
) -> dict[tuple[str, Database], ScaffoldAverage]:
    """Average VL-scores independently for each contig and database."""
    values: dict[tuple[str, Database], list[float]] = defaultdict(list)
    for record in evidence:
        if record.track != "relaxed":
            continue
        values[(record.contig_id, record.database)].append(record.vl_score)

    return {
        key: ScaffoldAverage(
            contig_id=key[0],
            database=key[1],
            score=sum(scores) / len(scores),
            denominator=len(scores),
        )
        for key, scores in values.items()
    }


def passes_candidate_thresholds(
    scaffold_avl: float,
    gene_vl_score: float,
    gene_v_score: float,
    *,
    scaffold_avl_cutoff: float = 3.0,
    gene_vl_cutoff: float = 3.0,
    gene_v_cutoff: float = 10.0,
) -> bool:
    """Apply the strict AVG candidate predicate with exclusive boundaries."""
    return (
        scaffold_avl > scaffold_avl_cutoff
        and gene_vl_score < gene_vl_cutoff
        and gene_v_score < gene_v_cutoff
    )


def evaluate_candidates(
    evidence: Iterable[AvgEvidence],
    *,
    scaffold_avl_cutoff: float = 3.0,
    gene_vl_cutoff: float = 3.0,
    gene_v_cutoff: float = 10.0,
) -> list[CandidateDecision]:
    """Evaluate candidates using only same-database scaffold averages."""
    records = list(evidence)
    averages = calculate_scaffold_averages(records)
    decisions: list[CandidateDecision] = []
    for record in sorted(records, key=lambda item: (item.protein_id, item.database, item.accession)):
        average = averages.get((record.contig_id, record.database))
        candidate = passes_candidate_thresholds(
            average.score,
            record.vl_score,
            record.v_score,
            scaffold_avl_cutoff=scaffold_avl_cutoff,
            gene_vl_cutoff=gene_vl_cutoff,
            gene_v_cutoff=gene_v_cutoff,
        ) if average is not None else False
        decisions.append(
            CandidateDecision(
                protein_id=record.protein_id,
                contig_id=record.contig_id,
                database=record.database,
                accession=record.accession,
                scaffold_avl=None if average is None else average.score,
                gene_vl_score=record.vl_score,
                gene_v_score=record.v_score,
                candidate=candidate,
                reason="passes_thresholds" if candidate else "fails_thresholds",
            )
        )
    return decisions


def curate_candidate(
    *,
    weight: float,
    positive: bool,
    conflicting: bool = False,
    filter_mode: FilterMode = "standard",
    minimum_weight: float = 0.6,
) -> CurationDecision:
    """Preserve candidate states while applying weight and filter policy."""
    if not 0.0 <= weight <= 1.0:
        raise ValueError("weight must be between 0 and 1")
    if not 0.0 <= minimum_weight <= 1.0:
        raise ValueError("minimum_weight must be between 0 and 1")
    if filter_mode not in {"standard", "strict", "none"}:
        raise ValueError("filter_mode must be 'standard', 'strict', or 'none'")

    below_weight = weight < minimum_weight
    filtered = filter_mode != "none" and below_weight
    if conflicting:
        classification = "conflicting"
        reason = "conflicting_evidence"
    elif filtered:
        classification = "filtered"
        reason = "below_weight" if filter_mode == "standard" else "strict_filter"
    elif not positive:
        classification = "unclassified_avg_candidate"
        reason = "no_positive_evidence"
    else:
        classification = "avg_candidate"
        reason = "positive_evidence"

    return CurationDecision(
        classification=classification,
        weight=weight,
        positive=positive,
        conflicting=conflicting,
        filtered=filtered,
        reason=reason,
    )


def collect_positive_evidence(
    evidence: Iterable[AvgEvidence],
    references: Iterable[PositiveEvidence],
) -> tuple[PositiveEvidence, ...]:
    """Match positive references against strict evidence using composite keys."""
    strict_keys = {
        database_key(record.database, record.accession)
        for record in evidence
        if record.track == "strict"
    }
    matches = [
        reference
        for reference in references
        if reference.proposed_class in SUPPORTED_CLASSES
        and database_key(reference.database, reference.accession) in strict_keys
    ]
    return tuple(
        sorted(
            matches,
            key=lambda item: (item.database, item.accession, item.proposed_class),
        )
    )


def apply_amg_weight(
    evidence: Iterable[PositiveEvidence],
    *,
    minimum_weight: float = 0.6,
) -> WeightDecision:
    """Retain supported classes after applying reference-level AMG weights."""
    if not 0.0 <= minimum_weight <= 1.0:
        raise ValueError("minimum_weight must be between 0 and 1")
    records = tuple(evidence)
    amg_records = tuple(
        item for item in records if item.proposed_class == "putative_amg"
    )
    maximum_weight = max(
        (item.amg_weight for item in amg_records if item.amg_weight is not None),
        default=None,
    )
    supported = {
        item.proposed_class
        for item in records
        if item.proposed_class in {"putative_apg", "putative_areg"}
        or (
            item.proposed_class == "putative_amg"
            and item.amg_weight is not None
            and item.amg_weight >= minimum_weight
        )
    }
    reasons: list[str] = []
    below_weight = bool(amg_records) and not any(
        item.amg_weight is not None and item.amg_weight >= minimum_weight
        for item in amg_records
    )
    if below_weight:
        reasons.append("amg_weight_below_cutoff")
    if not records:
        reasons.append("no_positive_evidence")
    return WeightDecision(
        evidence=records,
        supported_classes=frozenset(supported),
        maximum_amg_weight=maximum_weight,
        below_weight=below_weight,
        reasons=tuple(reasons),
    )


def apply_class_filters(
    proposed_classes: Iterable[str],
    filters: Iterable[FilterEvidence],
    *,
    filter_mode: FilterMode = "standard",
) -> FilterDecision:
    """Apply class-specific filter references while retaining audit evidence."""
    if filter_mode not in {"standard", "strict", "none"}:
        raise ValueError("filter_mode must be 'standard', 'strict', or 'none'")
    classes = frozenset(proposed_classes)
    if not classes.issubset(SUPPORTED_CLASSES):
        raise ValueError("proposed_classes contains an unsupported class")
    matches = tuple(
        sorted(
            (
                item
                for item in filters
                if item.proposed_class in classes
                and item.category in FILTER_CATEGORIES
            ),
            key=lambda item: (item.proposed_class, item.category, item.database, item.accession),
        )
    )
    enforced = {"filter_essential"}
    if filter_mode == "strict":
        enforced = set(FILTER_CATEGORIES)
    blocking = tuple(item for item in matches if item.category in enforced)
    warnings = tuple(item for item in matches if item.category not in enforced)
    blocked_classes = {item.proposed_class for item in blocking} if filter_mode != "none" else set()
    surviving = classes - blocked_classes
    reasons: list[str] = []
    if blocking:
        reasons.append("enforced_filter_match")
    if warnings:
        reasons.append("warning_filter_match")
    return FilterDecision(
        surviving_classes=frozenset(surviving),
        blocking_matches=blocking,
        warning_matches=warnings,
        reasons=tuple(reasons),
    )


def resolve_final_class(
    *,
    zhou_candidate: bool,
    weighted_classes: Iterable[str],
    filter_decision: FilterDecision | None = None,
    below_weight: bool = False,
) -> FinalClassDecision:
    """Resolve one documented status without imposing class priority."""
    if not zhou_candidate:
        return FinalClassDecision("", "not_candidate", ("fails_zhou_rule",))
    classes = frozenset(weighted_classes)
    if filter_decision is not None:
        classes = filter_decision.surviving_classes
    if len(classes) > 1:
        return FinalClassDecision("", "class_conflict", ("multiple_classes_survive",))
    if not classes:
        if filter_decision is not None and filter_decision.blocking_matches:
            return FinalClassDecision("", "filtered", filter_decision.reasons)
        if below_weight:
            return FinalClassDecision("", "below_amg_weight", ("amg_weight_below_cutoff",))
        return FinalClassDecision("", "unclassified_candidate", ("no_supported_class",))
    return FinalClassDecision(next(iter(classes)), "classified", ("one_class_survives",))
