from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from phu.normalize_vcontact_lineage import parse_candidate

FIXTURE = Path(__file__).parent / "fixtures" / "normalize_vcontact_lineage_cases.tsv"


def fixture_rows() -> list[dict[str, str]]:
    with FIXTURE.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def raw_value(row: dict[str, str]) -> str:
    value = row["raw_value"]
    return value


@pytest.mark.parametrize("row", fixture_rows(), ids=lambda row: row["case_id"])
def test_fixture_cases(row: dict[str, str]) -> None:
    value = raw_value(row)
    if "||" in value:
        results = [
            parse_candidate(part, row["column_rank"] or None)
            for part in value.split("||")
        ]
        status = (
            "unparsed"
            if any(item.status == "unparsed" for item in results)
            else "normalized"
        )
        assert status == row["expected_status"]
        assert "||".join(item.value or "" for item in results) == row["expected_value"]
        return
    result = parse_candidate(value, row["column_rank"] or None)
    assert result.status == row["expected_status"]
    assert result.rank_mismatch is (row["expected_rank_mismatch"] == "true")
    if row["expected_status"] == "missing":
        assert result.value is None
    else:
        assert result.value == row["expected_value"]


@pytest.mark.parametrize("value", [None, float("nan"), pd.NA])
def test_python_missing_values_are_table_layer_inputs(value: object) -> None:
    assert pd.isna(value)


@pytest.mark.parametrize("row", fixture_rows(), ids=lambda row: row["case_id"])
def test_normalized_values_are_idempotent(row: dict[str, str]) -> None:
    result = parse_candidate(raw_value(row), row["column_rank"] or None)
    if result.status == "normalized":
        again = parse_candidate(result.value or "", row["column_rank"] or None)
        assert again.status == "unchanged"
        assert again.value == result.value


@pytest.mark.parametrize("row", fixture_rows(), ids=lambda row: row["case_id"])
def test_unparsed_values_are_preserved(row: dict[str, str]) -> None:
    result = parse_candidate(raw_value(row), row["column_rank"] or None)
    if result.status == "unparsed":
        assert result.value == raw_value(row).strip()
