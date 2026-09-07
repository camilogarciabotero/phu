from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from typer.testing import CliRunner

from phu.cli import app
from phu.normalize_vcontact_lineage import QualitySummary, normalize_cell, parse_candidate

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


def test_quality_summary_counts_one_status_per_cell() -> None:
    summary = QualitySummary()

    value = normalize_cell("novel_genus_1_of_Viruses||novel_genus_2_of_Viruses", "genus", summary)

    assert value == "Viruses:NG1||Viruses:NG2"
    assert summary.cells_examined == 1
    assert summary.multi_candidate == 1
    assert summary.normalized == 1
    assert summary.unchanged == 0
    assert summary.unparsed == 0


def test_quiet_mode_hides_qa_summary() -> None:
    runner = CliRunner()
    input_file = FIXTURE.parent / "normalize_vcontact_lineage_quiet.tsv"
    input_file.write_text("genus_prediction\nnovel_genus_1_of_Viruses\n")
    output_file = input_file.with_name("normalize_vcontact_lineage_quiet_out.tsv")

    result = runner.invoke(
        app,
        [
            "normalize-lineage",
            "--input-file",
            str(input_file),
            "--output-file",
            str(output_file),
            "--quiet",
        ],
    )

    assert result.exit_code == 0
    assert "QA Summary" not in result.output
    assert output_file.exists()

    input_file.unlink(missing_ok=True)
    output_file.unlink(missing_ok=True)


@pytest.mark.parametrize("row", fixture_rows(), ids=lambda row: row["case_id"])
def test_unparsed_values_are_preserved(row: dict[str, str]) -> None:
    result = parse_candidate(raw_value(row), row["column_rank"] or None)
    if result.status == "unparsed":
        assert result.value == raw_value(row).strip()
