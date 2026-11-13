"""The database modules, driven through a real install of the fixture bundle."""

import os
from pathlib import Path

import pyhmmer
import pytest
from typer.testing import CliRunner

from phu.commands import avger as avg
from phu.db import bundle as db_bundle
from phu.db import dbcan as dbcan_db
from phu.db import kofam as kofam_db
from phu.db import pfam as pfam_db
from phu.commands import screen
from phu.avger.reference import load_avg_reference
from phu.cli import app


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
    monkeypatch.setattr(db_bundle, "_RESOLVED", {})
    monkeypatch.setattr(db_bundle, "_DEEP_VERIFIED", set())
    return db_folder / "phu-db" / fixture_pins.bundle_version


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


def test_avger_searches_phrog_from_the_installed_bundle(installed, tmp_path):
    from phu.avger.annotation import AnnotationConfig, annotate_proteins_tracks
    from phu.db import phrog as phrog_db

    proteins = tmp_path / "proteins.faa"
    # The consensus of fixture PHROG 2, and a protein that matches nothing.
    proteins.write_text(">ctg|gene1\nDMPQIFEWGYSHTNR\n>ctg|gene2\nAAAAAAAAAAAAAAA\n")
    # The fixture models are 15 residues long, far too short to reach 60 bits.
    cfg = AnnotationConfig(max_evalue=10.0, phrog_min_score=5.0)

    with_phrog = annotate_proteins_tracks(proteins, [cfg], phrog=True)[0]
    without = annotate_proteins_tracks(proteins, [cfg])[0]

    hit = with_phrog.best_phrog_by_protein["ctg|gene1"]
    assert hit.model_id == "phrog_2"
    assert hit.model_description == "fixture anti-restriction protein"
    assert hit.model_coverage == pytest.approx(1.0)
    assert hit.protein_coverage == pytest.approx(1.0)
    assert "ctg|gene2" not in with_phrog.best_phrog_by_protein
    assert without.best_phrog_by_protein == {}
    annotations = phrog_db.get_all_phrog_annotations()
    assert annotations["phrog_4"].annotation == ""
    assert annotations["phrog_4"].category == "unknown function"
    assert phrog_db.normalize_phrog_id("PHROG_007") == "phrog_7"


def test_screen_resolves_phrog_identifiers_from_the_installed_bundle(installed, tmp_path):
    from phu.db import phrog as phrog_db

    paths, missing = phrog_db.extract_phrog_models(["PHROG_2", "phrog_1", "phrog_4"], tmp_path / "x")
    assert [path.name for path in paths] == ["phrog_2.hmm", "phrog_1.hmm"]
    assert missing == ["phrog_4"]  # annotated upstream, but has no alignment
    assert _model_names(paths[0]) == ["phrog_2"]
    assert (installed / "phrog" / "offsets.json").is_file()

    hmms, ko_metadata, query = screen._resolve_hmm_inputs(
        [Path("phrog_2"), Path("PF00001")], tmp_path / "out", return_query=True
    )
    assert [path.name for path in hmms] == ["PF00001.hmm", "phrog_2.hmm"]
    assert query.normalized_models == ("PF00001", "phrog_2")
    assert query.databases == ("pfam", "phrog")
    manifest = screen._query_manifest(query, "any")
    assert manifest["phrog_annotations"] == {
        "phrog_2": {
            "annotation": "fixture anti-restriction protein",
            "category": "moron, auxiliary metabolic gene and host takeover",
        }
    }
    assert any(name.startswith("phrog/") for name in manifest["database_hashes"])

    with pytest.raises(FileNotFoundError, match="phrog_4"):
        screen._resolve_hmm_inputs([Path("phrog_4")], tmp_path / "out2")
    pins, bundle_manifest, _ = db_bundle._pinned(None, None)
    assert db_bundle.installed_problems(installed.parent.parent, pins, bundle_manifest, deep=True) == []


def test_avger_evalues_do_not_depend_on_the_protein_batch_size(installed, tmp_path):
    from phu.avger.annotation import AnnotationConfig, annotate_proteins_tracks

    proteins = tmp_path / "proteins.faa"
    proteins.write_text(
        ">ctg|gene1\nMKVLAAGIVGLLLAQ\n>ctg|gene2\nGSWTKHEPLRCDNYF\n>ctg|gene3\nDMPQIFEWGYSHTNR\n"
    )

    def evalues(batch_size):
        cfg = AnnotationConfig(
            max_evalue=10.0,
            protein_batch_size=batch_size,
            pfam_require_ga=False,
            kofam_require_thresholds=False,
            phrog_require_thresholds=False,
        )
        result = annotate_proteins_tracks(proteins, [cfg], phrog=True)[0]
        return {
            (protein, hit.model_id): hit.evalue
            for best in (result.best_pfam_by_protein, result.best_phrog_by_protein)
            for protein, hit in best.items()
        }

    whole = evalues(10_000)
    assert len(whole) >= 3
    assert evalues(1) == whole
    # The searches read binary copies built beside the databases.
    for name in ("pfam", "kofam", "phrog"):
        assert (installed / name / "profiles.h3m").is_file()
    pins, manifest, _ = db_bundle._pinned(None, None)
    assert db_bundle.installed_problems(installed.parent.parent, pins, manifest, deep=True) == []


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


def test_dbs_prepare_installs_the_whole_bundle(installed, fixture_pins):
    result = CliRunner().invoke(app, ["dbs", "prepare"])

    assert result.exit_code == 0, result.output
    assert f"{fixture_pins.bundle_version} is ready: {installed}" in result.stdout
    for name in ("pfam", "kofam", "dbcan", "avg"):
        assert (installed / name).is_dir()
        assert f"  {name}: " in result.stdout
    # The model files avger searches get their binary copies; dbCAN does not.
    for name in ("pfam", "kofam", "phrog"):
        assert f"Building fast-loading models for {name}" in result.output
        assert (installed / name / "profiles.h3m").is_file()
    assert not (installed / "dbcan" / "profiles.h3m").exists()

    again = CliRunner().invoke(app, ["dbs", "prepare"])
    assert again.exit_code == 0, again.output
    assert "Building fast-loading models" not in again.output


def test_dbs_prepare_takes_a_bundle_path_and_old_arguments(
    installed, monkeypatch, fixture_bundle
):
    monkeypatch.delenv("PHU_DB_BUNDLE")

    result = CliRunner().invoke(
        app, ["dbs", "prepare", "pfam", "--force-refresh", "--bundle", str(fixture_bundle)]
    )

    assert result.exit_code == 0, result.output
    assert (installed / "pfam" / "profiles.hmm").is_file()
    assert "all databases are installed together" in result.stderr
    assert "--force-refresh is ignored" in result.stderr


def test_dbs_prepare_verify_reinstalls_a_corrupt_file(installed):
    runner = CliRunner()
    assert runner.invoke(app, ["dbs", "prepare"]).exit_code == 0
    target = installed / "avg" / "v_scores.tsv"
    good = target.read_bytes()
    target.write_bytes(good[:-1] + (b"0" if good[-1:] != b"0" else b"1"))

    assert runner.invoke(app, ["dbs", "prepare"]).exit_code == 0
    assert target.read_bytes() != good
    assert runner.invoke(app, ["dbs", "prepare", "--verify"]).exit_code == 0
    assert target.read_bytes() == good


def test_dbs_prepare_reports_a_wrong_bundle(installed, monkeypatch, tmp_path):
    wrong = tmp_path / "wrong.tar.gz"
    wrong.write_bytes(b"not the bundle")

    result = CliRunner().invoke(app, ["dbs", "prepare", "--bundle", str(wrong)])

    assert result.exit_code == 1
    assert "different or corrupt bundle" in result.stderr


@pytest.fixture
def index_locks(monkeypatch):
    """Record each directory whose index lock is taken."""
    taken = []
    real = db_bundle.derived_files_lock

    def recording(directory, *args, **kwargs):
        taken.append(Path(directory).name)
        return real(directory, *args, **kwargs)

    for module in (pfam_db, kofam_db, dbcan_db):
        monkeypatch.setattr(module, "derived_files_lock", recording)
    return taken


def test_indexes_are_built_under_a_lock_and_only_once(installed, tmp_path, index_locks):
    pfam_hmm = Path(pfam_db.ensure_pfam_database()["hmm_path"])
    for _ in range(2):
        pfam_db.extract_pfam_models(pfam_hmm, ["PF00001"], tmp_path / "pfam")
        kofam_db.extract_kofam_models(["K00003"], tmp_path / "kofam")
        dbcan_db.extract_dbcan_models(["CBM48"], tmp_path / "dbcan")

    # The second round finds each index valid and takes no lock.
    assert index_locks == ["pfam", "kofam", "dbcan"]
    for name, model in (("pfam", "PF00001"), ("kofam", "K00003"), ("dbcan", "CBM48")):
        models = installed / name / "models"
        assert [path.name for path in models.iterdir()] == [f"{model}.hmm"]


def test_a_job_that_waited_for_the_index_does_not_clear_the_cache(
    installed, tmp_path, monkeypatch
):
    dbcan_db.extract_dbcan_models(["CBM48"], tmp_path / "first")
    cached = installed / "dbcan" / "models" / "CBM48.hmm"
    assert cached.is_file()
    # This job saw no index; by the time it holds the lock another job has built it.
    real = dbcan_db._read_offsets
    answers = iter([None])
    monkeypatch.setattr(
        dbcan_db, "_read_offsets", lambda source: next(answers, None) or real(source)
    )
    monkeypatch.setattr(
        dbcan_db,
        "_clear_sparse_models_cache",
        lambda: pytest.fail("the cache another job is using was cleared"),
    )

    paths, missing = dbcan_db.extract_dbcan_models(["GH5"], tmp_path / "second")

    assert missing == [] and _model_names(paths[0]) == ["GH5.hmm"]
    assert cached.is_file()


def test_query_manifest_records_the_bundle_for_every_database_query(installed, tmp_path):
    _, _, query = screen._resolve_hmm_inputs(
        [Path("PF00001"), Path("K00003"), Path("CBM48")], tmp_path / "out", return_query=True
    )

    manifest = screen._query_manifest(query, "any")

    assert query.databases == ("pfam", "kofam", "dbcan")
    assert query.normalized_models == ("PF00001", "K00003", "CBM48")
    assert manifest["phu_db"] == db_bundle.bundle_provenance()
    assert manifest["database_hashes"]["pfam/profiles.hmm"] == (
        db_bundle.database_hashes("pfam")["profiles.hmm"]
    )
    assert {key.split("/")[0] for key in manifest["database_hashes"]} == {
        "pfam",
        "kofam",
        "dbcan",
    }


def test_query_manifest_for_local_hmms_needs_no_bundle(monkeypatch):
    monkeypatch.setattr(
        screen, "bundle_provenance", lambda: pytest.fail("local HMMs do not use the bundle")
    )
    query = screen.ScreenQuery(original_inputs=("fam.hmm",), normalized_models=())

    manifest = screen._query_manifest(query, "any")

    assert manifest["database_hashes"] == {}
    assert "phu_db" not in manifest


def test_a_read_only_install_keeps_its_indexes_in_the_user_cache(
    installed, tmp_path, monkeypatch
):
    if os.geteuid() == 0:
        pytest.skip("root can write to a read-only directory")
    cache = tmp_path / "cache"
    monkeypatch.setenv("PHU_CACHE_DIR", str(cache))
    pfam_hmm = Path(pfam_db.ensure_pfam_database()["hmm_path"])
    kofam_db.ensure_kofam_database()
    dbcan_db.ensure_dbcan_database()
    locked = [installed, *(installed / name for name in ("pfam", "kofam", "dbcan"))]
    for path in locked:
        path.chmod(0o555)
    try:
        pfam, _ = pfam_db.extract_pfam_models(pfam_hmm, ["PF00001"], tmp_path / "pfam")
        kofam, _ = kofam_db.extract_kofam_models(["K00003"], tmp_path / "kofam")
        dbcan, _ = dbcan_db.extract_dbcan_models(["CBM48"], tmp_path / "dbcan")
        metadata = kofam_db.get_all_kofam_metadata()
        before = {str(path.relative_to(installed)) for path in installed.rglob("*")}
    finally:
        for path in locked:
            path.chmod(0o755)

    assert [len(paths) for paths in (pfam, kofam, dbcan)] == [1, 1, 1]
    assert set(metadata) == {"K00001", "K00002", "K00003"}
    # Nothing was added to the install; the indexes and models are in the cache.
    assert not {name for name in before if "models" in name or name.endswith(".json")} - {
        "installed.json",
        "phu-db.json",
    }
    derived = sorted(path.name.rsplit("-", 1)[0] for path in (cache / "phu-db-derived").iterdir())
    assert derived == ["dbcan", "kofam", "pfam"]
    assert len(list((cache / "phu-db-derived").rglob("*.hmm"))) == 3


def test_dbs_prepare_connections_option_sets_the_download_setting(installed, monkeypatch):
    monkeypatch.delenv("PHU_DB_CONNECTIONS", raising=False)
    monkeypatch.setattr(os, "environ", dict(os.environ))

    result = CliRunner().invoke(app, ["dbs", "prepare", "--connections", "2"])

    assert result.exit_code == 0, result.output
    assert os.environ["PHU_DB_CONNECTIONS"] == "2"
    assert CliRunner().invoke(app, ["dbs", "prepare", "--connections", "0"]).exit_code == 2
