from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from typer.testing import CliRunner

from phu.cli import app
from phu.commands.normalize_lineage import (
    QualitySummary,
    normalize_cell,
    normalize_dataframe,
    parse_candidate,
)

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


@pytest.mark.parametrize(
    ("raw", "rank", "expected"),
    [
        ("unplaced_phylum_within_unplaced_kingdom_within_default", "phylum", "DEF:UK:UP"),
        (
            "unplaced_class_within_unplaced_phylum_within_unplaced_kingdom_within_default",
            "class",
            "DEF:UK:UP:UC",
        ),
        ("unplaced_order_within_Caudoviricetes", "order", "Caudoviricetes:UO"),
        (
            "novel_genus_1_of_novel_family_2_of_unplaced_order_within_Caudoviricetes",
            "genus",
            "Caudoviricetes:UO:NF2:NG1",
        ),
    ],
)
def test_unplaced_chains_are_normalized(raw: str, rank: str, expected: str) -> None:
    result = parse_candidate(raw, rank)

    assert result.status == "normalized"
    assert result.value == expected
    assert result.deepest_rank == rank
    assert result.rank_mismatch is False

    again = parse_candidate(expected, rank)
    assert again.status == "unchanged"
    assert again.value == expected
    assert again.deepest_rank == rank
    assert again.rank_mismatch is False


@pytest.mark.parametrize(
    "raw",
    [
        "unplaced_kingdom_within_unplaced_phylum_within_default",
        "unplaced_phylum_within_unplaced_kingdom",
        "unplaced_phylum",
    ],
)
def test_invalid_unplaced_chains_are_unparsed(raw: str) -> None:
    result = parse_candidate(raw, "phylum")

    assert result.status == "unparsed"
    assert result.value == raw


def test_unplaced_multi_candidate_cell_counts_as_normalized() -> None:
    summary = QualitySummary()

    value = normalize_cell(
        "unplaced_phylum_within_unplaced_kingdom_within_default"
        "||novel_phylum_3_of_Viruses",
        "phylum",
        summary,
    )

    assert value == "DEF:UK:UP||Viruses:NP3"
    assert summary.cells_examined == 1
    assert summary.multi_candidate == 1
    assert summary.normalized == 1
    assert summary.unchanged == 0
    assert summary.unparsed == 0
    assert summary.rank_mismatch == 0


def test_unplaced_dataframe_second_run_is_a_no_op() -> None:
    dataframe = pd.DataFrame(
        {
            "contig": ["a", "b"],
            "kingdom_prediction": ["unplaced_kingdom_within_default", "Bamfordvirae"],
            "phylum_prediction": [
                "unplaced_phylum_within_unplaced_kingdom_within_default",
                "",
            ],
            "class_prediction": [
                "unplaced_class_within_unplaced_phylum_within_unplaced_kingdom_within_default",
                "novel_class_1_of_unplaced_phylum_within_Bamfordvirae",
            ],
        }
    )

    first, first_summary = normalize_dataframe(dataframe, add_lineage=True)

    assert first["kingdom_prediction"].tolist() == ["DEF:UK", "Bamfordvirae"]
    assert first["phylum_prediction"].iloc[0] == "DEF:UK:UP"
    assert first["class_prediction"].tolist() == ["DEF:UK:UP:UC", "Bamfordvirae:UP:NC1"]
    assert first["compact_lineage"].tolist() == ["DEF:UK:UP:UC", "Bamfordvirae:UP:NC1"]
    assert first_summary.normalized == 4
    assert first_summary.unchanged == 1
    assert first_summary.missing == 1
    assert first_summary.unparsed == 0
    assert first_summary.rank_mismatch == 0
    assert first_summary.lineage_skipped_unparsed == 0

    second, second_summary = normalize_dataframe(
        first.drop(columns="compact_lineage"), add_lineage=True
    )

    pd.testing.assert_frame_equal(second, first)
    assert second_summary.normalized == 0
    assert second_summary.unchanged == 5
    assert second_summary.unparsed == 0
    assert second_summary.rank_mismatch == 0
    assert second_summary.lineage_skipped_unparsed == 0
