import pytest

from phu.avg_decisions import (
    AvgEvidence,
    FilterEvidence,
    PositiveEvidence,
    apply_amg_weight,
    apply_class_filters,
    calculate_scaffold_averages,
    collect_positive_evidence,
    curate_candidate,
    database_key,
    evaluate_candidates,
    passes_candidate_thresholds,
    resolve_final_class,
)


def evidence(protein_id, database, accession, vl_score, v_score):
    return AvgEvidence(
        protein_id=protein_id,
        contig_id="contig-1",
        database=database,
        accession=accession,
        track="relaxed",
        v_score=v_score,
        vl_score=vl_score,
    )


def test_database_key_keeps_database_origin_in_join_key():
    assert database_key(" Pfam ", "pf00001") == ("pfam", "PF00001")
    assert database_key("KEGG", "k00001") == ("kofam", "K00001")


def test_database_key_rejects_unknown_or_empty_values():
    with pytest.raises(ValueError, match="Unsupported AVG database"):
        database_key("vscore", "PF00001")
    with pytest.raises(ValueError, match="must not be empty"):
        database_key("pfam", "")


def test_scaffold_averages_are_database_specific():
    records = [
        evidence("gene-1", "pfam", "PF00001", 4.0, 10.0),
        evidence("gene-2", "pfam", "PF00002", 2.0, 10.0),
        evidence("gene-3", "kofam", "K00001", 1.0, 10.0),
    ]

    averages = calculate_scaffold_averages(records)

    assert averages[("contig-1", "pfam")].score == pytest.approx(3.0)
    assert averages[("contig-1", "pfam")].denominator == 2
    assert averages[("contig-1", "kofam")].score == 1.0


def test_candidate_predicate_uses_exclusive_boundaries():
    assert passes_candidate_thresholds(3.01, 2.99, 9.99)
    assert not passes_candidate_thresholds(3.0, 2.99, 9.99)
    assert not passes_candidate_thresholds(3.01, 3.0, 9.99)
    assert not passes_candidate_thresholds(3.01, 2.99, 10.0)


def test_candidates_do_not_mix_databases_and_are_deterministic():
    records = [
        evidence("gene-2", "pfam", "PF00002", 6.0, 10.0),
        evidence("gene-1", "pfam", "PF00001", 1.0, 9.0),
        evidence("gene-3", "kofam", "K00001", 1.0, 9.0),
    ]

    decisions = evaluate_candidates(records)

    assert [decision.protein_id for decision in decisions] == [
        "gene-1",
        "gene-2",
        "gene-3",
    ]
    assert decisions[0].candidate is True
    assert decisions[1].candidate is False
    assert decisions[2].candidate is False


def test_curation_preserves_conflicting_and_filtered_states():
    assert curate_candidate(weight=0.59, positive=True).classification == "filtered"
    assert curate_candidate(
        weight=0.59, positive=True, filter_mode="none"
    ).classification == "avg_candidate"
    assert curate_candidate(
        weight=0.9, positive=False
    ).classification == "unclassified_avg_candidate"
    assert curate_candidate(
        weight=0.9, positive=True, conflicting=True
    ).classification == "conflicting"


def test_curation_validates_modes_and_weights():
    with pytest.raises(ValueError, match="filter_mode"):
        curate_candidate(weight=0.8, positive=True, filter_mode="unknown")
    with pytest.raises(ValueError, match="weight"):
        curate_candidate(weight=1.1, positive=True)


def test_positive_matching_requires_strict_composite_key_evidence():
    records = [evidence("gene-1", "pfam", "PF00001", 2.0, 9.0)]
    records[0] = AvgEvidence(**{**records[0].__dict__, "track": "relaxed"})
    references = [
        PositiveEvidence("pfam", "PF00001", "putative_amg", amg_weight=0.9),
        PositiveEvidence("kofam", "PF00001", "putative_apg"),
    ]

    assert collect_positive_evidence(records, references) == ()


def test_amg_weight_uses_maximum_and_does_not_affect_other_classes():
    decision = apply_amg_weight(
        [
            PositiveEvidence("kofam", "K00001", "putative_amg", amg_weight=0.5),
            PositiveEvidence("kofam", "K00002", "putative_amg", amg_weight=0.6),
            PositiveEvidence("pfam", "PF00001", "putative_apg"),
        ]
    )

    assert decision.maximum_amg_weight == 0.6
    assert decision.supported_classes == {"putative_amg", "putative_apg"}
    assert decision.below_weight is False


def test_class_filters_have_standard_strict_and_none_modes():
    filters = [
        FilterEvidence("pfam", "PF1", "putative_amg", "filter_essential"),
        FilterEvidence("pfam", "PF2", "putative_amg", "filter_glucan"),
    ]

    standard = apply_class_filters({"putative_amg"}, filters)
    strict = apply_class_filters({"putative_amg"}, filters, filter_mode="strict")
    none = apply_class_filters({"putative_amg"}, filters, filter_mode="none")

    assert standard.surviving_classes == set()
    assert len(standard.warning_matches) == 1
    assert strict.surviving_classes == set()
    assert len(strict.blocking_matches) == 2
    assert none.surviving_classes == {"putative_amg"}
    assert len(none.blocking_matches) + len(none.warning_matches) == 2


def test_final_class_resolution_preserves_candidate_statuses():
    empty_filters = apply_class_filters({"putative_amg"}, [])
    assert resolve_final_class(
        zhou_candidate=False, weighted_classes=[]
    ).status == "not_candidate"
    assert resolve_final_class(
        zhou_candidate=True, weighted_classes=[], below_weight=True
    ).status == "below_amg_weight"
    assert resolve_final_class(
        zhou_candidate=True, weighted_classes=["putative_amg"], filter_decision=empty_filters
    ).status == "classified"
    assert resolve_final_class(
        zhou_candidate=True,
        weighted_classes=["putative_amg", "putative_apg"],
    ).status == "class_conflict"
