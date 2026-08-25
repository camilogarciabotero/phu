from pathlib import Path

import pytest

from phu.avger_annotation import AnnotationHit
from phu.avger_classification import (
    UNCLASSIFIED_AVG_CANDIDATE,
    classify_protein_annotations,
    load_default_classification_rules,
    load_classification_rules,
)
from phu.vscore_db import VScoreRecord


def _hit(database: str, model_id: str) -> AnnotationHit:
    return AnnotationHit(
        protein_id="ctg|gene1",
        contig_id="ctg",
        database=database,
        model_id=model_id,
        model_accession=model_id,
        score_type="full",
        effective_score=50.0,
        full_score=50.0,
        domain_score=None,
        evalue=1e-10,
        hmm_from=None,
        hmm_to=None,
        target_from=None,
        target_to=None,
        threshold_source="test",
        threshold_value=1.0,
    )


def test_default_rules_are_bundled_and_conservative():
    rules = load_default_classification_rules()
    assert rules.version == "builtin-2026-08-24"
    assert rules.rules == ()


def test_rules_are_versioned_and_applied_in_declared_order(tmp_path: Path):
    path = tmp_path / "rules.json"
    path.write_text(
        '{"version":"2026-08-24","rules":['
        '{"rule_id":"broad","classification":"broad"},'
        '{"rule_id":"specific","classification":"specific","required_kofam":["K00001"]}'
        ']}'
    )
    rules = load_classification_rules(path)

    classification = classify_protein_annotations([_hit("kofam", "K00001")], rules)
    assert classification == ("broad", "broad", "2026-08-24")


def test_unmatched_annotations_remain_explicitly_unclassified(tmp_path: Path):
    path = tmp_path / "rules.json"
    path.write_text(
        '{"version":"v1","rules":[{"rule_id":"ko","classification":"known","required_kofam":["K00001"]}]}'
    )
    rules = load_classification_rules(path)

    assert classify_protein_annotations([_hit("pfam", "PF00001")], rules) == (
        UNCLASSIFIED_AVG_CANDIDATE,
        None,
        "v1",
    )


def test_vscore_threshold_can_gate_a_rule(tmp_path: Path):
    path = tmp_path / "rules.json"
    path.write_text(
        '{"version":"v1","rules":[{"rule_id":"high-v","classification":"high-v","min_v_score":8}]}'
    )
    rules = load_classification_rules(path)
    scores = {
        "K00001": VScoreRecord("K00001", "function", 9.0, 2.0, 4.0, "KEGG")
    }

    assert classify_protein_annotations([_hit("kofam", "K00001")], rules, scores)[0] == "high-v"
    scores["K00001"] = VScoreRecord("K00001", "function", 7.0, 2.0, 4.0, "KEGG")
    assert classify_protein_annotations([_hit("kofam", "K00001")], rules, scores)[0] == UNCLASSIFIED_AVG_CANDIDATE


def test_duplicate_rule_ids_are_rejected(tmp_path: Path):
    path = tmp_path / "rules.json"
    path.write_text(
        '{"version":"v1","rules":[{"rule_id":"same","classification":"a"},{"rule_id":"same","classification":"b"}]}'
    )
    with pytest.raises(ValueError, match="Duplicate classification"):
        load_classification_rules(path)