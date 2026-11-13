from __future__ import annotations

import csv
import gzip
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pyhmmer

from phu.core.log import log_call
from phu.db.kofam import KOFamMetadata, ensure_kofam_database, get_all_kofam_metadata
from phu.db.pfam import ensure_pfam_database, normalize_pfam_id
from phu.db.phrog import (
    PhrogAnnotation,
    ensure_phrog_database,
    get_all_phrog_annotations,
    is_phrog_id,
    normalize_phrog_id,
)


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
    model_name: Optional[str] = None
    model_description: Optional[str] = None
    track: str = "strict"
    # Fractions of the model and of the protein covered by the reported
    # domain alignments together (see _coverage).
    model_coverage: Optional[float] = None
    protein_coverage: Optional[float] = None


@dataclass(frozen=True)
class AnnotationConfig:
    threads: int = 1
    max_evalue: float = 1e-5
    protein_batch_size: int = 10_000
    keep_all_hits: bool = False
    all_hits_path: Optional[Path] = None
    pfam_require_ga: bool = True
    pfam_missing_ga_policy: str = "skip_model"  # skip_model | include_without_ga
    kofam_require_thresholds: bool = True
    # PHROG has no curated per-model cutoffs. The strict track uses fixed ones,
    # as CheckAMG does: sequence bit score, and coverage of both model and
    # protein. CheckAMG's coverage default is 0.4; phu asks for half.
    phrog_require_thresholds: bool = True
    phrog_min_score: float = 60.0
    phrog_min_coverage: float = 0.5
    track: str = "strict"

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
        if not 0.0 <= self.phrog_min_coverage <= 1.0:
            raise ValueError("phrog_min_coverage must be between 0 and 1")
        if self.track not in {"relaxed", "strict"}:
            raise ValueError("track must be 'relaxed' or 'strict'")


@dataclass
class AnnotationResults:
    best_pfam_by_protein: dict[str, AnnotationHit]
    best_kofam_by_protein: dict[str, AnnotationHit]
    passing_hit_count: int
    scanned_model_count: int
    skipped_pfam_models_missing_ga: int
    best_phrog_by_protein: dict[str, AnnotationHit] = field(default_factory=dict)

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


def _covered(intervals: list[tuple[int, int]]) -> int:
    """Positions covered by 1-based inclusive intervals, overlaps counted once."""
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        if stop > end:
            total += stop - max(start, end + 1) + 1
            end = stop
    return total


def _coverage(hit, model_length: int) -> tuple[Optional[float], Optional[float]]:
    """Model and protein fractions covered by the reported domains of a hit.

    This is CheckAMG's definition: every reported domain counts, not only the
    best one, so a protein matched in two pieces is credited for both.
    """
    domains = getattr(hit, "domains", [])
    model: list[tuple[int, int]] = []
    protein: list[tuple[int, int]] = []
    for domain in getattr(domains, "reported", domains):
        alignment = getattr(domain, "alignment", None)
        if alignment is None:
            continue
        model.append((int(alignment.hmm_from), int(alignment.hmm_to)))
        protein.append((int(alignment.target_from), int(alignment.target_to)))
    if not model:
        return None, None
    protein_length = getattr(hit, "length", None)
    return (
        _covered(model) / model_length if model_length else None,
        _covered(protein) / protein_length if protein_length else None,
    )


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


def _hmm_text(value) -> Optional[str]:
    if isinstance(value, bytes):
        value = value.decode(errors="replace")
    if not value:
        return None
    return str(value).strip() or None


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
    """Search one database for one track (see _search_tracks)."""
    return _search_tracks(
        database=database,
        hmm_path=hmm_path,
        proteins_path=proteins_path,
        cfgs=[cfg],
        kofam_meta_by_model=kofam_meta_by_model,
        all_hits_writers=[all_hits_writer],
    )[0]


def _score_thresholds(database: str, cfg: AnnotationConfig) -> str:
    if database == "pfam":
        return "pfam_ga" if cfg.pfam_require_ga else "none"
    if database == "phrog":
        return "phrog_score_coverage" if cfg.phrog_require_thresholds else "none"
    return "kofam_ko_list" if cfg.kofam_require_thresholds else "none"


def _phrog_model_id(hmm) -> Optional[str]:
    name = _hmm_text(getattr(hmm, "name", None))
    if name is None or not is_phrog_id(name):
        return None
    return normalize_phrog_id(name)


def _search_tracks(
    *,
    database: str,
    hmm_path: Path,
    proteins_path: Path,
    cfgs: Sequence[AnnotationConfig],
    kofam_meta_by_model: Optional[dict[str, KOFamMetadata]] = None,
    phrog_annotations: Optional[dict[str, PhrogAnnotation]] = None,
    all_hits_writers: Optional[Sequence] = None,
) -> list[tuple[dict[str, AnnotationHit], int, int, int]]:
    """Search one database once and collect the hits of every track in ``cfgs``.

    The tracks differ only in how hits are filtered afterwards (E-value limit,
    Pfam GA, KOfam thresholds, PHROG score and coverage), never in how the search is run, so one search
    serves them all. Returns, per track and in the order of ``cfgs``: best hit
    per protein, passing hit count, models seen, Pfam models skipped for
    lacking GA cutoffs.
    """
    if not cfgs:
        raise ValueError("at least one annotation track is required")
    first = cfgs[0]
    if any(
        cfg.threads != first.threads or cfg.protein_batch_size != first.protein_batch_size
        for cfg in cfgs
    ):
        raise ValueError("tracks searched together must share threads and batch size")
    writers = list(all_hits_writers) if all_hits_writers is not None else [None] * len(cfgs)
    if len(writers) != len(cfgs):
        raise ValueError("all_hits_writers must match the tracks")

    best_by_protein: list[dict[str, AnnotationHit]] = [{} for _ in cfgs]
    passing = [0] * len(cfgs)
    skipped_missing_ga = [0] * len(cfgs)
    models = 0
    seen_models: set[str] = set()
    kofam_meta_by_model = kofam_meta_by_model or {}
    phrog_annotations = phrog_annotations or {}
    log_call(
        "pyhmmer.hmmsearch",
        database=database,
        track="+".join(cfg.track for cfg in cfgs),
        hmm_file=hmm_path,
        proteins=proteins_path,
        protein_batch_size=first.protein_batch_size,
        cpus=first.threads,
        parallel="queries",
        bit_cutoffs="none",
        max_evalue="+".join(str(cfg.max_evalue) for cfg in cfgs),
        score_thresholds="+".join(_score_thresholds(database, cfg) for cfg in cfgs),
    )

    with pyhmmer.easel.SequenceFile(proteins_path, digital=True) as seq_file:
        while True:
            proteins = seq_file.read_block(sequences=first.protein_batch_size)
            if not proteins:
                break

            # HMMFile is consumed by hmmsearch, so reopen it for each bounded
            # protein batch rather than retaining all proteins in memory.
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hit_batches = pyhmmer.hmmsearch(
                    hmm_file,
                    proteins,
                    cpus=first.threads,
                    parallel="queries",
                )
                for top_hits in hit_batches:
                    hmm = top_hits.query
                    raw_model = getattr(hmm, "accession", None) or getattr(hmm, "name", "")
                    if isinstance(raw_model, bytes):
                        raw_model = raw_model.decode()
                    model_key = f"{database}:{raw_model}"
                    if model_key not in seen_models:
                        seen_models.add(model_key)
                        models += 1

                    seq_ga = dom_ga = None
                    ko_meta = None
                    # Tracks that take this model at all.
                    active = list(range(len(cfgs)))
                    if database == "pfam":
                        model_accession = _normalize_pfam_accession(hmm)
                        if model_accession is None:
                            continue
                        model_id = model_accession
                        model_name = _hmm_text(getattr(hmm, "name", None))
                        model_description = _hmm_text(getattr(hmm, "description", None))

                        ga_pair = _pfam_model_ga(hmm)
                        if ga_pair is None:
                            active = []
                            for index, cfg in enumerate(cfgs):
                                if (
                                    cfg.pfam_require_ga
                                    and cfg.pfam_missing_ga_policy == "skip_model"
                                ):
                                    skipped_missing_ga[index] += 1
                                else:
                                    active.append(index)
                            if not active:
                                continue
                        else:
                            seq_ga, dom_ga = ga_pair
                    elif database == "phrog":
                        model_id = _phrog_model_id(hmm)
                        if model_id is None:
                            continue
                        model_accession = model_id
                        model_name = model_id
                        phrog = phrog_annotations.get(model_id)
                        # Most PHROGs have no known function: keep the hit anyway.
                        model_description = phrog.annotation if phrog else ""
                        model_length = int(hmm.M)
                    else:
                        model_id = _kofam_model_id(hmm)
                        if model_id is None:
                            continue
                        model_accession = model_id
                        model_name = model_id
                        ko_meta = kofam_meta_by_model.get(model_id)
                        if ko_meta is None:
                            continue
                        model_description = ko_meta.definition

                    for hit in top_hits:
                        evalue = float(hit.evalue)
                        tracks = [i for i in active if evalue <= cfgs[i].max_evalue]
                        if not tracks:
                            continue

                        protein_id = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                        contig_id = parse_contig_id_from_protein_id(protein_id)

                        full_score = float(hit.score)
                        domain_score, hmm_from, hmm_to, target_from, target_to = _domain_metrics(hit)
                        model_coverage = protein_coverage = None

                        if database == "phrog":
                            threshold_source = "phrog_score_coverage"
                            score_type = "full"
                            effective_score = full_score
                            model_coverage, protein_coverage = _coverage(hit, model_length)
                            threshold_value = None
                            strict_tracks = []
                            for i in tracks:
                                cfg = cfgs[i]
                                if cfg.phrog_require_thresholds and (
                                    full_score < cfg.phrog_min_score
                                    or model_coverage is None
                                    or model_coverage < cfg.phrog_min_coverage
                                    or protein_coverage is None
                                    or protein_coverage < cfg.phrog_min_coverage
                                ):
                                    continue
                                strict_tracks.append(i)
                            tracks = strict_tracks
                        elif database == "pfam":
                            threshold_source = "pfam_ga"
                            threshold_value = seq_ga
                            score_type = "full"
                            effective_score = full_score
                            below_ga = (seq_ga is not None and full_score < seq_ga) or (
                                dom_ga is not None
                                and domain_score is not None
                                and domain_score < dom_ga
                            )
                            tracks = [
                                i for i in tracks if not (cfgs[i].pfam_require_ga and below_ga)
                            ]
                        else:
                            score_type = ko_meta.score_type
                            threshold_source = "kofam_ko_list"
                            threshold_value = ko_meta.threshold
                            if score_type == "domain":
                                if domain_score is None:
                                    continue
                                effective_score = domain_score
                            else:
                                effective_score = full_score
                            below_threshold = threshold_value is None or effective_score < float(
                                threshold_value
                            )
                            tracks = [
                                i
                                for i in tracks
                                if not (cfgs[i].kofam_require_thresholds and below_threshold)
                            ]

                        for index in tracks:
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
                                threshold_value=(
                                    cfgs[index].phrog_min_score
                                    if database == "phrog"
                                    and cfgs[index].phrog_require_thresholds
                                    else threshold_value
                                ),
                                model_name=model_name,
                                model_description=model_description,
                                track=cfgs[index].track,
                                model_coverage=model_coverage,
                                protein_coverage=protein_coverage,
                            )

                            passing[index] += 1
                            if writers[index] is not None:
                                _write_hit(writers[index], record)

                            current = best_by_protein[index].get(record.protein_id)
                            if _better_hit(record, current):
                                best_by_protein[index][record.protein_id] = record

    return [
        (best_by_protein[index], passing[index], models, skipped_missing_ga[index])
        for index in range(len(cfgs))
    ]


def _write_hit(writer, record: AnnotationHit) -> None:
    writer.writerow(
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


def annotate_proteins_complete_databases(
    proteins_faa: Path,
    cfg: AnnotationConfig,
    *,
    phrog: bool = False,
) -> AnnotationResults:
    """Annotate against the complete databases for one track."""
    return annotate_proteins_tracks(proteins_faa, [cfg], phrog=phrog)[0]


def annotate_proteins_tracks(
    proteins_faa: Path,
    cfgs: Sequence[AnnotationConfig],
    *,
    phrog: bool = False,
) -> list[AnnotationResults]:
    """Annotate for several tracks with one search of each database.

    Pfam and KOfam are always searched; PHROG only when ``phrog`` is true.
    """
    for cfg in cfgs:
        if cfg.keep_all_hits and cfg.all_hits_path is None:
            raise ValueError("all_hits_path is required when keep_all_hits is True")

    pfam_meta = ensure_pfam_database()
    kofam_meta = ensure_kofam_database()

    pfam_hmm_path = Path(pfam_meta["hmm_path"])
    kofam_hmm_path = Path(kofam_meta["hmm_path"])
    all_kofam_meta = get_all_kofam_metadata()

    handles = []
    writers = []
    try:
        for cfg in cfgs:
            if cfg.keep_all_hits:
                handle, writer = _open_all_hits_writer(cfg.all_hits_path)
                handles.append(handle)
                writers.append(writer)
            else:
                writers.append(None)

        pfam = _search_tracks(
            database="pfam",
            hmm_path=pfam_hmm_path,
            proteins_path=proteins_faa,
            cfgs=cfgs,
            all_hits_writers=writers,
        )
        kofam = _search_tracks(
            database="kofam",
            hmm_path=kofam_hmm_path,
            proteins_path=proteins_faa,
            cfgs=cfgs,
            kofam_meta_by_model=all_kofam_meta,
            all_hits_writers=writers,
        )
        phrog_results = None
        if phrog:
            phrog_meta = ensure_phrog_database()
            phrog_results = _search_tracks(
                database="phrog",
                hmm_path=Path(phrog_meta["hmm_path"]),
                proteins_path=proteins_faa,
                cfgs=cfgs,
                phrog_annotations=get_all_phrog_annotations(),
                all_hits_writers=writers,
            )
    finally:
        for handle in handles:
            handle.close()

    return [
        AnnotationResults(
            best_pfam_by_protein=pfam[index][0],
            best_kofam_by_protein=kofam[index][0],
            passing_hit_count=pfam[index][1]
            + kofam[index][1]
            + (phrog_results[index][1] if phrog_results else 0),
            scanned_model_count=pfam[index][2]
            + kofam[index][2]
            + (phrog_results[index][2] if phrog_results else 0),
            skipped_pfam_models_missing_ga=pfam[index][3],
            best_phrog_by_protein=phrog_results[index][0] if phrog_results else {},
        )
        for index in range(len(cfgs))
    ]
