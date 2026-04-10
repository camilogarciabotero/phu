from pathlib import Path

import pytest

import phu.pfam_db as pfam_db
from phu.pfam_db import extract_pfam_models, is_pfam_id, normalize_pfam_id


def test_is_pfam_id_accepts_basic_and_versioned_ids():
    assert is_pfam_id("PF00001")
    assert is_pfam_id("pf12345")
    assert is_pfam_id("PF54321.12")


def test_is_pfam_id_rejects_invalid_ids():
    assert not is_pfam_id("PF0001")
    assert not is_pfam_id("PFABCDE")
    assert not is_pfam_id("PFX0001")
    assert not is_pfam_id("PF12345.")


def test_normalize_pfam_id_removes_version():
    assert normalize_pfam_id("pf00001.22") == "PF00001"
    assert normalize_pfam_id("PF54321") == "PF54321"


def test_normalize_pfam_id_raises_for_invalid_values():
    with pytest.raises(ValueError, match="Invalid PFAM accession"):
        normalize_pfam_id("PF12")


def test_extract_pfam_models_extracts_requested_models(tmp_path: Path):
    hmm_db = tmp_path / "Pfam-A.hmm"
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  ModelA\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
        "NAME  ModelB\n"
        "ACC   PF00002.3\n"
        "LENG  33\n"
        "//\n"
    )

    out_dir = tmp_path / "resolved"
    extracted, missing = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001", "PF00002.99"],
        output_dir=out_dir,
    )

    assert missing == []
    assert len(extracted) == 2
    assert (out_dir / "PF00001.hmm").exists()
    assert (out_dir / "PF00002.hmm").exists()


def test_extract_pfam_models_reports_missing_ids(tmp_path: Path):
    hmm_db = tmp_path / "Pfam-A.hmm"
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  ModelA\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
    )

    out_dir = tmp_path / "resolved"
    extracted, missing = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001", "PF99999"],
        output_dir=out_dir,
    )

    assert len(extracted) == 1
    assert missing == ["PF99999"]


def test_extract_pfam_models_reuses_split_cache_without_rescanning(tmp_path: Path, monkeypatch):
    hmm_db = tmp_path / "Pfam-A.hmm"
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  ModelA\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
        "NAME  ModelB\n"
        "ACC   PF00002.3\n"
        "LENG  33\n"
        "//\n"
    )

    first_out = tmp_path / "resolved_first"
    extracted_1, missing_1 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001"],
        output_dir=first_out,
    )
    assert missing_1 == []
    assert len(extracted_1) == 1

    index_path = tmp_path / "offsets.json"
    assert index_path.exists()

    def _fail_rebuild(*args, **kwargs):
        raise AssertionError("unexpected offsets-index rebuild")

    monkeypatch.setattr(pfam_db, "_build_offsets_index", _fail_rebuild)

    second_out = tmp_path / "resolved_second"
    extracted_2, missing_2 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00002"],
        output_dir=second_out,
    )

    assert missing_2 == []
    assert len(extracted_2) == 1
    assert (second_out / "PF00002.hmm").exists()


def test_extract_pfam_models_rebuilds_split_cache_when_source_changes(tmp_path: Path):
    hmm_db = tmp_path / "Pfam-A.hmm"
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  ModelA\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
    )

    first_out = tmp_path / "resolved_first"
    extracted_1, missing_1 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001"],
        output_dir=first_out,
    )
    assert missing_1 == []
    assert len(extracted_1) == 1

    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  ModelA\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
        "NAME  ModelZ\n"
        "ACC   PF99999.4\n"
        "LENG  17\n"
        "//\n"
    )

    second_out = tmp_path / "resolved_second"
    extracted_2, missing_2 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF99999"],
        output_dir=second_out,
    )

    assert missing_2 == []
    assert len(extracted_2) == 1
    assert (second_out / "PF99999.hmm").exists()


def test_extract_pfam_models_invalidates_stale_sparse_cache_when_source_changes(tmp_path: Path):
    hmm_db = tmp_path / "Pfam-A.hmm"
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  OldModel\n"
        "ACC   PF00001.1\n"
        "LENG  42\n"
        "//\n"
    )

    out_first = tmp_path / "resolved_first"
    extracted_1, missing_1 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001"],
        output_dir=out_first,
    )
    assert missing_1 == []
    assert len(extracted_1) == 1

    # Change source DB for same accession; stale sparse model cache must not be reused.
    hmm_db.write_text(
        "HMMER3/f\n"
        "NAME  NewModel\n"
        "ACC   PF00001.2\n"
        "LENG  55\n"
        "DESC  updated\n"
        "//\n"
    )

    out_second = tmp_path / "resolved_second"
    extracted_2, missing_2 = extract_pfam_models(
        hmm_db_path=hmm_db,
        requested_ids=["PF00001"],
        output_dir=out_second,
    )

    assert missing_2 == []
    assert len(extracted_2) == 1
    contents = (out_second / "PF00001.hmm").read_text()
    assert "NAME  NewModel" in contents
    assert "DESC  updated" in contents


def test_prepare_pfam_database_missing_hmm_message_points_to_dbs_command(tmp_path: Path, monkeypatch):
    missing_hmm = tmp_path / "pfam" / "Pfam-A.hmm"

    monkeypatch.setattr(pfam_db, "_pfam_hmm_path", lambda: missing_hmm)

    with pytest.raises(FileNotFoundError, match="phu dbs prepare pfam"):
        pfam_db.prepare_pfam_database(download=False, index=False, force_refresh=False)
