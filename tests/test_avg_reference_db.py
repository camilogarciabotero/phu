import json
from pathlib import Path

import pytest

import phu.avg_reference_db as reference


SOURCE_CONTENT = {
    "AMGs.tsv": "id\tdb\tname\tmetabolic_paths\ttotal_paths\tpath_source\tmetabolic_ratio\tVL-score\tAMG_weight\tamg_level\nPF00123.17\tPfam\tenergy enzyme\t1\t2\tPfam\t0.5\t2.0\t0.8\tvery high\nK00001\tKEGG\tko enzyme\t1\t1\tKEGG\t1\t1.0\t0.7\thigh\n",
    "APGs.tsv": "id\tdb\tname\nPF00123\tPfam\tenergy enzyme\n",
    "AReGs.tsv": "id\tdb\tname\nK00002\tKEGG\tregulator\n",
    "AMG_filters.tsv": "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl\tfilter_nucleotide\nPF00123\tPfam\tenergy enzyme\tTrue\tFalse\tFalse\tFalse\tFalse\n",
    "APG_filters.tsv": "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl\tfilter_nucleotide\nK00001\tKEGG\tko enzyme\tFalse\tTrue\tFalse\tFalse\tFalse\n",
    "AReG_filters.tsv": "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl\tfilter_nucleotide\nK00002\tKEGG\tregulator\tFalse\tFalse\tFalse\tFalse\tTrue\n",
    "vscores.tsv": "id\tname\tV-score\tVL-score\tdb\nPF00123.17\tenergy enzyme\t9.0\t2.0\tPfam\nK00001\tko enzyme\t8.0\t4.0\tKEGG\n",
}


def install_fixture(monkeypatch, tmp_path: Path) -> Path:
    root = tmp_path / "avg"
    root.mkdir()
    monkeypatch.setattr(reference, "_avg_root", lambda: root)
    for name, content in SOURCE_CONTENT.items():
        (root / name).write_text(content)
    return root


def test_prepare_compiles_normalized_tables_and_manifest(monkeypatch, tmp_path):
    root = install_fixture(monkeypatch, tmp_path)

    manifest = reference.prepare_avg_database()

    scores = (root / "v_scores.tsv").read_text()
    positive = (root / "avg_positive.tsv").read_text()
    filters = (root / "avg_filters.tsv").read_text()
    assert "pfam\tPF00123\tPF00123.17" in scores
    assert "kofam\tK00001" in scores
    assert "putative_amg" in positive
    assert "filter_essential" in filters
    assert manifest["release"] == "v1.1.1"
    assert manifest["normalized_row_counts"] == {"positive": 4, "filters": 3, "scores": 2}
    assert json.loads((root / "manifest.json").read_text())["files"]["v_scores.tsv"]["sha256"]


def test_prepare_skips_other_database_rows(monkeypatch, tmp_path):
    root = install_fixture(monkeypatch, tmp_path)
    (root / "vscores.tsv").write_text(
        SOURCE_CONTENT["vscores.tsv"] + "X1\tother\t1\t1\tCAMPER\n"
    )

    manifest = reference.prepare_avg_database()

    assert manifest["normalized_row_counts"]["scores"] == 2


def test_conflicting_duplicate_scores_fail(monkeypatch, tmp_path):
    root = install_fixture(monkeypatch, tmp_path)
    (root / "vscores.tsv").write_text(
        "id\tname\tV-score\tVL-score\tdb\n"
        "PF00123\tfirst\t9\t2\tPfam\n"
        "PF00123\tsecond\t9\t2\tPfam\n"
    )

    with pytest.raises(ValueError, match="Conflicting AVG score rows"):
        reference.prepare_avg_database()


def test_status_reports_ready_reference(monkeypatch, tmp_path):
    install_fixture(monkeypatch, tmp_path)
    reference.prepare_avg_database()

    status = reference.get_avg_database_status()

    assert status["ready"] is True
    assert status["release"] == "v1.1.1"
    assert status["normalized_row_counts"]["scores"] == 2


def test_ensure_repairs_hash_invalid_file(monkeypatch, tmp_path):
    root = install_fixture(monkeypatch, tmp_path)
    monkeypatch.setattr(
        reference,
        "_download",
        lambda _url, destination: destination.write_text(SOURCE_CONTENT[destination.name]),
    )
    reference.prepare_avg_database()
    (root / "v_scores.tsv").write_text("corrupt\n")

    repaired = reference.ensure_avg_database()

    assert repaired["files"]["v_scores.tsv"]["size"] == (root / "v_scores.tsv").stat().st_size
    assert "PF00123" in (root / "v_scores.tsv").read_text()
