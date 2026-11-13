from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from phu.avger.annotation import (
    AnnotationConfig,
    AnnotationResults,
    annotate_proteins_tracks,
)


@dataclass(frozen=True)
class AvgAnnotationResults:
    """Independent relaxed scoring and strict functional annotation tracks."""

    relaxed: AnnotationResults
    strict: AnnotationResults


def annotate_avg_tracks(
    proteins_faa: Path,
    *,
    threads: int = 1,
    scoring_evalue: float = 1e-5,
    keep_all_hits: bool = False,
    relaxed_all_hits_path: Path | None = None,
    strict_all_hits_path: Path | None = None,
    phrog: bool = True,
) -> AvgAnnotationResults:
    """Annotate complete Pfam, KOfam and (unless ``phrog`` is false) PHROG.

    The relaxed track uses only the independent E-value cutoff. The strict
    track retains Pfam GA and KOfam adaptive thresholds, and for PHROG a bit
    score of at least 60 with at least 40% of both model and protein aligned,
    for functional calls.
    """
    if keep_all_hits and (
        relaxed_all_hits_path is None or strict_all_hits_path is None
    ):
        raise ValueError(
            "both all-hit paths are required when keep_all_hits is True"
        )

    # One search of each database serves both tracks: they differ only in the
    # thresholds applied to the hits.
    relaxed, strict = annotate_proteins_tracks(
        proteins_faa,
        [
            AnnotationConfig(
                threads=threads,
                max_evalue=scoring_evalue,
                keep_all_hits=keep_all_hits,
                all_hits_path=relaxed_all_hits_path,
                pfam_require_ga=False,
                kofam_require_thresholds=False,
                phrog_require_thresholds=False,
                track="relaxed",
            ),
            AnnotationConfig(
                threads=threads,
                max_evalue=scoring_evalue,
                keep_all_hits=keep_all_hits,
                all_hits_path=strict_all_hits_path,
                track="strict",
            ),
        ],
        phrog=phrog,
    )
    return AvgAnnotationResults(relaxed=relaxed, strict=strict)
