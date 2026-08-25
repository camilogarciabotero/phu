from pathlib import Path

import pytest

from phu.avger_annotation import AnnotationHit, AnnotationResults, write_best_hits_tsv
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


def test_best_hit_writer_adds_vscore_columns_and_sorts_rows(tmp_path: Path):
    rows = [
        AnnotationHit(
            protein_id="ctgB|gene1",
            contig_id="ctgB",
            database="kofam",
            model_id="K00001",
            model_accession="K00001",
            score_type="full",
            effective_score=80,
            full_score=80,
            domain_score=None,
            evalue=1e-20,
            hmm_from=None,
            hmm_to=None,
            target_from=None,
            target_to=None,
            threshold_source="kofam_ko_list",
            threshold_value=50,
        ),
        AnnotationHit(
            protein_id="ctgA|gene1",
            contig_id="ctgA",
            database="pfam",
            model_id="PF00001",
            model_accession="PF00001",
            score_type="full",
            effective_score=40,
            full_score=40,
            domain_score=35,
            evalue=1e-10,
            hmm_from=1,
            hmm_to=20,
            target_from=2,
            target_to=21,
            threshold_source="pfam_ga",
            threshold_value=30,
        ),
    ]
    result = AnnotationResults(
        best_pfam_by_protein={rows[1].protein_id: rows[1]},
        best_kofam_by_protein={rows[0].protein_id: rows[0]},
        passing_hit_count=2,
        scanned_model_count=2,
        skipped_pfam_models_missing_ga=0,
    )
    output = tmp_path / "best.tsv"

    assert write_best_hits_tsv(
        result,
        output,
        {"K00001": type("Score", (), {
            "v_score": 10.0,
            "vl_score": 2.0,
            "protein_function": "integrase",
            "log10_hit_number": 4.2,
            "database_origin": "KEGG",
        })()},
    ) == 2
    lines = output.read_text().splitlines()
    assert lines[0].startswith("protein_id\tcontig_id")
    assert lines[1].startswith("ctgA|gene1\tctgA\tpfam")
    assert "10.000000\t2.000000\tintegrase\t4.200000\tKEGG" in lines[2]
    assert "unclassified_avg_candidate" in lines[2]