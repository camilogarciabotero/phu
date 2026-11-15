"""`phu clean`: list, then remove the cache and/or the installed databases."""

from __future__ import annotations

from typer.testing import CliRunner

from phu.cli import app

runner = CliRunner()


def _setup(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    (cache / "proteins").mkdir(parents=True)
    (cache / "proteins" / "p.faa").write_text(">a\nM\n")
    db = tmp_path / "db"
    (db / "phu-db" / "phu-db-0.2.0").mkdir(parents=True)
    (db / "phu-db" / "phu-db-0.2.0" / "profiles.hmm").write_bytes(b"x" * 2048)
    (db / ".downloads").mkdir()
    (db / "mine.txt").write_text("not phu's")  # never touched
    monkeypatch.setenv("PHU_CACHE_DIR", str(cache))
    monkeypatch.setenv("PHU_DB_FOLDER", str(db))
    return cache, db


def test_without_options_lists_and_removes_nothing(tmp_path, monkeypatch):
    cache, db = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["clean"])
    assert result.exit_code == 0
    assert "cache" in result.stdout and "2.0 KB" in result.stdout
    assert "Nothing removed" in result.stdout
    assert cache.exists() and (db / "phu-db").exists()


def test_cache_only(tmp_path, monkeypatch):
    cache, db = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["clean", "--cache"])
    assert result.exit_code == 0
    assert not cache.exists()
    assert (db / "phu-db").exists()


def test_dbs_asks_and_keeps_on_no(tmp_path, monkeypatch):
    _, db = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["clean", "--dbs"], input="n\n")
    assert result.exit_code == 0
    assert "Databases kept" in result.stdout
    assert (db / "phu-db").exists()


def test_all_with_yes_removes_only_phus_files(tmp_path, monkeypatch):
    cache, db = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["clean", "--all", "--yes"])
    assert result.exit_code == 0, result.output
    assert not cache.exists()
    assert not (db / "phu-db").exists() and not (db / ".downloads").exists()
    assert (db / "mine.txt").read_text() == "not phu's"


def test_dry_run(tmp_path, monkeypatch):
    cache, db = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["clean", "--all", "--dry-run"])
    assert result.exit_code == 0 and "Dry run" in result.stdout
    assert cache.exists() and (db / "phu-db").exists()


def test_nothing_to_remove(tmp_path, monkeypatch):
    monkeypatch.setenv("PHU_CACHE_DIR", str(tmp_path / "none"))
    monkeypatch.setenv("PHU_DB_FOLDER", str(tmp_path / "nodb"))
    result = runner.invoke(app, ["clean", "--all", "--yes"])
    assert result.exit_code == 0 and "Nothing to remove" in result.stdout


def test_old_flag_still_works_with_a_notice(tmp_path, monkeypatch):
    cache, _ = _setup(tmp_path, monkeypatch)
    result = runner.invoke(app, ["--clean-cache"])
    assert result.exit_code == 0
    assert "phu clean --cache" in result.stderr
    assert not cache.exists()
