"""Alignment coverage of a pyhmmer hit, computed the way CheckAMG computes it.

Every reported domain of the hit counts, and overlapping domains are counted
once, on the model side and on the protein side.
"""

from __future__ import annotations

from typing import Optional


def covered_positions(intervals: list[tuple[int, int]]) -> int:
    """Positions covered by 1-based inclusive intervals, overlaps counted once."""
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        if stop > end:
            total += stop - max(start, end + 1) + 1
            end = stop
    return total


def hit_coverage(hit, model_length: int) -> tuple[Optional[float], Optional[float]]:
    """Fractions of the model and of the protein covered by a hit's domains."""
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
        covered_positions(model) / model_length if model_length else None,
        covered_positions(protein) / protein_length if protein_length else None,
    )
