from pathlib import Path

import pytest

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
