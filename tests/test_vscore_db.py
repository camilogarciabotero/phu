from pathlib import Path

import pytest

from phu.vscore_db import parse_vscore_csv


def test_parse_vscore_csv_validates_required_columns(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("Accession,V-Score\nK00001,10\n")

    with pytest.raises(ValueError, match="missing required columns"):
        parse_vscore_csv(path)


def test_parse_vscore_csv_normalizes_accession_and_reads_values(tmp_path: Path):
    path = tmp_path / "scores.csv"
    path.write_text(
        "Accession,Protein Function,V-Score,Normalized VL-score,Log10[Hit Number],Database Origin\n"
        "k00001,integrase,10,2.5,4.2,KEGG\n"
    )

    records = parse_vscore_csv(path)
    assert records["K00001"].v_score == 10.0
    assert records["K00001"].vl_score == 2.5
    assert records["K00001"].protein_function == "integrase"


def test_parse_vscore_csv_accepts_legacy_normalized_v_score_column(tmp_path: Path):
    path = tmp_path / "legacy.csv"
    path.write_text(
        "Accession,Protein Function,V-Score,Normalized V-score,Log10[Hit Number],Database Origin\n"
        "k00001,integrase,10,2.5,4.2,KEGG\n"
    )

    records = parse_vscore_csv(path)
    assert records["K00001"].vl_score == 2.5


def test_get_vscore_map_refreshes_stale_csv(monkeypatch, tmp_path: Path):
    csv_path = tmp_path / "VScoreData.csv"
    csv_path.write_text(
        "Accession,Protein Function,V-Score,Log10[Hit Number],Database Origin\n"
        "K00001,integrase,10,4.2,KEGG\n"
    )

    called = {"refresh": False}

    def fake_ensure(force_refresh: bool = False):
        called["refresh"] = force_refresh
        csv_path.write_text(
            "Accession,Protein Function,V-Score,Normalized VL-score,Log10[Hit Number],Database Origin\n"
            "K00001,integrase,10,2.5,4.2,KEGG\n"
        )
        return {"csv_path": str(csv_path)}

    monkeypatch.setattr("phu.vscore_db._vscore_csv_path", lambda: csv_path)
    monkeypatch.setattr("phu.vscore_db.ensure_vscore_database", fake_ensure)

    records = __import__("phu.vscore_db", fromlist=["get_vscore_map"]).get_vscore_map()

    assert called["refresh"] is True
    assert records["K00001"].vl_score == 2.5


def test_parse_vscore_csv_rejects_malformed_numeric_values(tmp_path: Path):
    path = tmp_path / "bad-value.csv"
    path.write_text(
        "Accession,Protein Function,V-Score,Normalized VL-score,Log10[Hit Number],Database Origin\n"
        "K00001,integrase,not-a-number,2.5,4.2,KEGG\n"
    )

    with pytest.raises(ValueError, match="numeric fields are malformed"):
        parse_vscore_csv(path)
