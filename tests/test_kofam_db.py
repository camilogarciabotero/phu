from pathlib import Path

import pytest

from phu.kofam_db import (
    extract_kofam_models,
    is_kofam_id,
    normalize_kofam_id,
)


def test_is_kofam_id_accepts_valid_ids():
    assert is_kofam_id("K00001")
    assert is_kofam_id("k12345")


def test_is_kofam_id_rejects_invalid_ids():
    assert not is_kofam_id("K0001")
    assert not is_kofam_id("KO0001")
    assert not is_kofam_id("K12A45")


def test_normalize_kofam_id_normalizes_case():
    assert normalize_kofam_id("k00077") == "K00077"


def test_normalize_kofam_id_raises_on_invalid_id():
    with pytest.raises(ValueError, match="Invalid KO identifier"):
        normalize_kofam_id("K123")


def test_extract_kofam_models_extracts_from_profiles(tmp_path: Path, monkeypatch):
    hmm_db_path = tmp_path / "kofam.hmm"
    hmm_db_path.write_text("HMMER3/f\nNAME  K00001\nLENG  100\n//\n")

    models_dir = tmp_path / "kofam" / "models"
    out_dir = tmp_path / "resolved"

    import phu.kofam_db as kofam_db

    monkeypatch.setattr(kofam_db, "_kofam_hmm_path", lambda: hmm_db_path)
    monkeypatch.setattr(kofam_db, "_kofam_models_dir", lambda: models_dir)
    monkeypatch.setattr(
        kofam_db, "_kofam_offsets_index_path", lambda: tmp_path / "offsets.json"
    )

    extracted, missing = extract_kofam_models(["K00001"], out_dir)

    assert missing == []
    assert len(extracted) == 1
    assert (out_dir / "K00001.hmm").exists()
    assert (models_dir / "K00001.hmm").exists()
