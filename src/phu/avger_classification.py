from __future__ import annotations

import json
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from .vscore_db import VScoreRecord

if TYPE_CHECKING:
    from .avger_annotation import AnnotationHit

UNCLASSIFIED_AVG_CANDIDATE = "unclassified_avg_candidate"


@dataclass(frozen=True)
class ClassificationRule:
    rule_id: str
    classification: str
    required_kofam: frozenset[str] = frozenset()
    required_pfam: frozenset[str] = frozenset()
    min_v_score: Optional[float] = None


@dataclass(frozen=True)
class ClassificationRules:
    version: str
    rules: tuple[ClassificationRule, ...]


def load_default_classification_rules() -> ClassificationRules:
    """Load the bundled, versioned rule set used by the avger workflow."""
    resource = files("phu").joinpath("data", "avger_classification_rules.json")
    try:
        payload = json.loads(resource.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("Bundled avger classification rules are unavailable") from exc
    return _parse_classification_rules(payload)


def _parse_classification_rules(payload: object) -> ClassificationRules:
    if not isinstance(payload, dict) or not isinstance(payload.get("version"), str):
        raise TypeError("Classification rules require a string 'version'")
    raw_rules = payload.get("rules")
    if not isinstance(raw_rules, list):
        raise TypeError("Classification rules require a list 'rules'")

    rules: list[ClassificationRule] = []
    seen_ids: set[str] = set()
    for index, raw_rule in enumerate(raw_rules, start=1):
        if not isinstance(raw_rule, dict):
            raise TypeError(f"Classification rule {index} must be an object")
        rule_id = raw_rule.get("rule_id")
        classification = raw_rule.get("classification")
        if not isinstance(rule_id, str) or not rule_id:
            raise ValueError(f"Classification rule {index} requires 'rule_id'")
        if rule_id in seen_ids:
            raise ValueError(f"Duplicate classification rule_id: {rule_id}")
        if not isinstance(classification, str) or not classification:
            raise ValueError(f"Classification rule {rule_id} requires 'classification'")

        def identifiers(key: str) -> frozenset[str]:
            value = raw_rule.get(key, [])
            if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
                raise TypeError(f"Classification rule {rule_id} field '{key}' must be a list of strings")
            return frozenset(item.strip().upper() for item in value if item.strip())

        min_v_score = raw_rule.get("min_v_score")
        if min_v_score is not None:
            if not isinstance(min_v_score, (int, float)) or isinstance(min_v_score, bool):
                raise TypeError(f"Classification rule {rule_id} min_v_score must be numeric")
            min_v_score = float(min_v_score)

        rules.append(
            ClassificationRule(
                rule_id=rule_id,
                classification=classification,
                required_kofam=identifiers("required_kofam"),
                required_pfam=identifiers("required_pfam"),
                min_v_score=min_v_score,
            )
        )
        seen_ids.add(rule_id)
    return ClassificationRules(version=payload["version"], rules=tuple(rules))


def load_classification_rules(path: Path) -> ClassificationRules:
    """Load a versioned JSON rule set with strict structural validation."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read classification rules: {path}") from exc

    return _parse_classification_rules(payload)


def classify_protein_annotations(
    hits: list["AnnotationHit"],
    rules: Optional[ClassificationRules],
    vscore_by_accession: Optional[dict[str, VScoreRecord]] = None,
) -> tuple[str, Optional[str], Optional[str]]:
    """Apply rules in declared order, retaining an explicit unresolved state."""
    if rules is None:
        return UNCLASSIFIED_AVG_CANDIDATE, None, None

    kofam_ids = {hit.model_id.upper() for hit in hits if hit.database == "kofam"}
    pfam_ids = {hit.model_id.upper() for hit in hits if hit.database == "pfam"}
    best_v_score = max(
        (
            vscore_by_accession[hit.model_id].v_score
            for hit in hits
            if vscore_by_accession is not None and hit.model_id in vscore_by_accession
        ),
        default=None,
    )

    for rule in rules.rules:
        if not rule.required_kofam.issubset(kofam_ids):
            continue
        if not rule.required_pfam.issubset(pfam_ids):
            continue
        if rule.min_v_score is not None and (best_v_score is None or best_v_score < rule.min_v_score):
            continue
        return rule.classification, rule.rule_id, rules.version

    return UNCLASSIFIED_AVG_CANDIDATE, None, rules.version