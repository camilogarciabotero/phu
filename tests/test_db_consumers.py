"""The database modules, driven through a real install of the fixture bundle."""

from pathlib import Path

import pyhmmer
import pytest

from phu import avg, db_bundle, dbcan_db, kofam_db, pfam_db
from phu.avg_reference_db import load_avg_reference


@pytest.fixture
def installed(monkeypatch, tmp_path, fixture_bundle, fixture_pins, fixture_manifest_bytes):
    """phu pinned to the fixture bundle, which installs from PHU_DB_BUNDLE on first use."""
    db_folder = tmp_path / "db"
    monkeypatch.setenv("PHU_DB_FOLDER", str(db_folder))
    monkeypatch.setenv("PHU_DB_BUNDLE", str(fixture_bundle))
    monkeypatch.delenv("PHU_DB_VERIFY", raising=False)
    monkeypatch.setattr(
        db_bundle, "_production_pins", lambda: (fixture_pins, fixture_manifest_bytes)
    )
    monkeypatch.setattr(db_bundle, "_log_stderr", lambda message: None)

    def fetch(*args, **kwargs):
        raise AssertionError("a local bundle needs no download")

    monkeypatch.setattr(db_bundle, "download_verified", fetch)
    kofam_db._kofam_root.cache_clear()
    dbcan_db._dbcan_root.cache_clear()
    yield db_folder / "phu-db" / fixture_pins.bundle_version
    kofam_db._kofam_root.cache_clear()
    dbcan_db._dbcan_root.cache_clear()


def _model_names(path):
    with pyhmmer.plan7.HMMFile(path) as handle:
        return [hmm.name for hmm in handle]


def test_pfam_models_come_from_the_installed_bundle(installed, tmp_path):
    meta = pfam_db.ensure_pfam_database()

    assert meta["hmm_path"] == str(installed / "pfam" / "profiles.hmm")
    paths, missing = pfam_db.extract_pfam_models(
        installed / "pfam" / "profiles.hmm", ["PF00002.7", "PF09999"], tmp_path / "out"
    )
    assert [path.name for path in paths] == ["PF00002.hmm"]
    assert _model_names(paths[0]) == ["FixtureB"]
    assert missing == ["PF09999"]


def test_kofam_models_and_thresholds_come_from_the_installed_bundle(installed, tmp_path):
    meta = kofam_db.ensure_kofam_database()

    assert meta["hmm_path"] == str(installed / "kofam" / "profiles.hmm")
    assert len(_model_names(meta["hmm_path"])) == 3
    paths, missing = kofam_db.extract_kofam_models(["k00003", "K09999"], tmp_path / "out")
    assert _model_names(paths[0]) == ["K00003"]
    assert missing == ["K09999"]
    metadata = kofam_db.get_all_kofam_metadata()
    assert set(metadata) == {"K00001", "K00002", "K00003"}
    assert metadata["K00001"].threshold == 120.5
    assert metadata["K00001"].score_type == "domain"
    assert metadata["K00002"].threshold is None


def test_dbcan_models_and_pul_rules_come_from_the_installed_bundle(
    installed, tmp_path, fixture_manifest
):
    meta = dbcan_db.ensure_dbcan_database()

    dbcan = fixture_manifest.databases["dbcan"]
    hashes = {item.path: item.sha256 for item in dbcan.prepared_files}
    assert meta["hmm_path"] == str(installed / "dbcan" / "profiles.hmm")
    assert meta["release"] == dbcan.upstream_version
    assert meta["hmm_sha256"] == hashes["profiles.hmm"]
    assert dbcan_db.get_dbcan_model_inventory() == (("GH5", "CBM48"), ())
    paths, missing = dbcan_db.extract_dbcan_models(["CBM48", "GH9"], tmp_path / "out")
    assert _model_names(paths[0]) == ["CBM48.hmm"]
    assert missing == ["GH9"]
    rules = dbcan_db.get_dbcan_pul_rules()
    assert rules["PUL0001"].families == ("GH5", "CBM48")
    assert rules["PUL0001"].resolved
    assert rules["PUL0003"].unresolved_tokens == ("GH999",)
    assert dbcan_db.get_dbcan_pul("PUL0002").substrate == "starch glycogen"


def test_avg_reference_comes_from_the_installed_bundle(installed):
    avg.require_databases()

    reference = load_avg_reference(db_bundle.resolve_db("avg"))
    assert [(row["database"], row["accession"]) for row in reference.scores][:2] == [
        ("kofam", "K00001"),
        ("kofam", "K00002"),
    ]
    assert db_bundle.bundle_provenance()["source"] == "local"


def test_index_files_beside_the_databases_do_not_invalidate_the_install(installed, tmp_path):
    hmm_path = pfam_db.ensure_pfam_database()["hmm_path"]
    pfam_db.extract_pfam_models(Path(hmm_path), ["PF00001"], tmp_path / "pfam-out")
    kofam_db.extract_kofam_models(["K00001"], tmp_path / "kofam-out")
    kofam_db.get_all_kofam_metadata()
    dbcan_db.extract_dbcan_models(["GH5"], tmp_path / "dbcan-out")

    assert (installed / "pfam" / "offsets.json").is_file()
    assert (installed / "kofam" / "offsets.json").is_file()
    assert (installed / "dbcan" / "family_offsets.json").is_file()
    pins, manifest, _ = db_bundle._pinned(None, None)
    assert db_bundle.installed_problems(installed.parent.parent, pins, manifest, deep=True) == []
