import pytest

from phu.avger_annotation import AnnotationHit
from phu.avger_scoring import evaluate_database_candidates
from phu.vscore_db import VScoreRecord


def hit(protein_id, database, model_id, start, end):
    return AnnotationHit(
        protein_id=protein_id, contig_id="ctg", database=database,
        model_id=model_id, model_accession=model_id, score_type="full",
        effective_score=50, full_score=50, domain_score=None, evalue=1e-10,
        hmm_from=None, hmm_to=None, target_from=start, target_to=end,
        threshold_source="test", threshold_value=1,
        gene_start=start, gene_end=end,
    )


def test_candidate_uses_strict_three_and_ten_boundaries():
    hits = [hit("gene1", "pfam", "PF1", 20_000, 20_100), hit("gene2", "pfam", "PF2", 1, 100)]
    scores = {
        "PF1": VScoreRecord("PF1", "candidate", 9.999, 2.999, 1, "Pfam"),
        "PF2": VScoreRecord("PF2", "flank", 10, 5, 1, "Pfam"),
    }
    result = evaluate_database_candidates(hits, scores, {"ctg": 30_000})
    candidate = result[0]
    assert candidate.candidate is True
    assert candidate.contig_avl_score == pytest.approx(3.9995)

    scores["PF1"] = VScoreRecord("PF1", "boundary", 10, 3, 1, "Pfam")
    assert evaluate_database_candidates(hits, scores, {"ctg": 30_000})[0].candidate is False


def test_flanks_are_database_consistent_and_not_required_by_default():
    hits = [
        hit("candidate", "pfam", "PF1", 20_000, 20_100),
        hit("up", "pfam", "PF2", 15_000, 15_100),
        hit("down", "kofam", "K1", 25_000, 25_100),
    ]
    scores = {
        "PF1": VScoreRecord("PF1", "candidate", 9, 1, 1, "Pfam"),
        "PF2": VScoreRecord("PF2", "viral", 10, 6, 1, "Pfam"),
        "K1": VScoreRecord("K1", "viral", 10, 1, 1, "KEGG"),
    }
    result = evaluate_database_candidates(hits, scores, {"ctg": 30_000})[0]
    assert result.candidate is True
    assert result.flank.upstream_supported is True
    assert result.flank.downstream_supported is False
    assert result.evidence_state == "avg_candidate"
    assert evaluate_database_candidates(hits, scores, {"ctg": 30_000}, True)[0].evidence_state == "avg_candidate"