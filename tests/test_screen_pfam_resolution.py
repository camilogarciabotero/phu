from pathlib import Path

import pytest

import phu.screen as screen
from phu.kofam_db import KOFamMetadata


def test_resolve_hmm_inputs_keeps_existing_local_files(tmp_path: Path):
    local_hmm = tmp_path / "custom.hmm"
    local_hmm.write_text("HMMER3/f\n//\n")

    resolved, ko_meta = screen._resolve_hmm_inputs([local_hmm], tmp_path / "out")
    assert resolved == [local_hmm]
    assert ko_meta == {}


def test_resolve_hmm_inputs_resolves_pfam_ids(monkeypatch, tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    pfam_db_hmm = tmp_path / "Pfam-A.hmm"
    pfam_db_hmm.write_text("dummy")

    def _fake_ensure(force_refresh: bool = False):
        return {"hmm_path": str(pfam_db_hmm)}

    def _fake_extract(hmm_db_path: Path, requested_ids, output_dir: Path):
        assert hmm_db_path == pfam_db_hmm
        assert requested_ids == ["PF00001"]
        extracted = output_dir / "PF00001.hmm"
        output_dir.mkdir(parents=True, exist_ok=True)
        extracted.write_text("HMMER3/f\nACC   PF00001.1\n//\n")
        return [extracted], []

    monkeypatch.setattr(screen, "ensure_pfam_database", _fake_ensure)
    monkeypatch.setattr(screen, "extract_pfam_models", _fake_extract)

    resolved, ko_meta = screen._resolve_hmm_inputs([Path("PF00001")], outdir)

    assert len(resolved) == 1
    assert resolved[0].name == "PF00001.hmm"
    assert ko_meta == {}


def test_resolve_hmm_inputs_resolves_kofam_ids(monkeypatch, tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    def _fake_ensure(force_refresh: bool = False):
        return {"profiles_path": str(tmp_path / "profiles")}

    def _fake_extract(requested_ids, output_dir: Path):
        assert requested_ids == ["K00001"]
        extracted = output_dir / "K00001.hmm"
        output_dir.mkdir(parents=True, exist_ok=True)
        extracted.write_text("HMMER3/f\nNAME  K00001\n//\n")
        return [extracted], []

    def _fake_metadata(ko_ids):
        assert ko_ids == ["K00001"]
        return {
            "K00001": KOFamMetadata(
                ko_id="K00001",
                threshold=55.0,
                score_type="full",
                profile_type="trim",
                definition="hexokinase",
            )
        }

    monkeypatch.setattr(screen, "ensure_kofam_database", _fake_ensure)
    monkeypatch.setattr(screen, "extract_kofam_models", _fake_extract)
    monkeypatch.setattr(screen, "get_kofam_metadata_map", _fake_metadata)

    resolved, ko_meta = screen._resolve_hmm_inputs([Path("K00001")], outdir)

    assert len(resolved) == 1
    assert resolved[0].name == "K00001.hmm"
    assert ko_meta["K00001"].definition == "hexokinase"


def test_resolve_hmm_inputs_raises_for_missing_pfam(monkeypatch, tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    pfam_db_hmm = tmp_path / "Pfam-A.hmm"
    pfam_db_hmm.write_text("dummy")

    def _fake_ensure(force_refresh: bool = False):
        return {"hmm_path": str(pfam_db_hmm)}

    def _fake_extract(hmm_db_path: Path, requested_ids, output_dir: Path):
        return [], ["PF99999"]

    monkeypatch.setattr(screen, "ensure_pfam_database", _fake_ensure)
    monkeypatch.setattr(screen, "extract_pfam_models", _fake_extract)

    with pytest.raises(FileNotFoundError, match=r"PFAM accession\(s\) not found"):
        screen._resolve_hmm_inputs([Path("PF99999")], outdir)


def test_resolve_hmm_inputs_rejects_unknown_missing_tokens(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="HMM file not found"):
        screen._resolve_hmm_inputs([Path("not_a_file_or_pfam")], tmp_path / "out")
