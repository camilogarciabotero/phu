from typer.testing import CliRunner
from phu.cli import app
import phu.cli as cli_module

runner = CliRunner()


def test_root_help_runs():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "cluster" in result.stdout
    assert "screen" in result.stdout
    assert "jack" in result.stdout
    assert "dbs" in result.stdout
    assert "simplify-taxa" in result.stdout
    assert "--clean-cache" in result.stdout


def test_clean_cache_removes_prediction_cache(tmp_path, monkeypatch):
    cache_dir = tmp_path / "phu-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "dummy.txt").write_text("cached")
    monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

    result = runner.invoke(app, ["--clean-cache"])

    assert result.exit_code == 0
    assert "Removed cache directory" in result.stdout
    assert not cache_dir.exists()


def test_dbs_help_runs():
    result = runner.invoke(app, ["dbs", "--help"])
    assert result.exit_code == 0
    assert "list" in result.stdout
    assert "status" in result.stdout
    assert "prepare" in result.stdout
    assert "refresh" in result.stdout
    assert "remove" in result.stdout


def test_dbs_prepare_calls_pfam_prepare(monkeypatch, tmp_path):
    hmm_path = tmp_path / "Pfam-A.hmm"
    offsets_path = tmp_path / "offsets.json"

    def fake_prepare_pfam_database(*, download, index, force_refresh):
        assert download is True
        assert index is True
        assert force_refresh is False
        hmm_path.write_text("dummy")
        return {"hmm_path": str(hmm_path), "offsets_path": str(offsets_path), "downloaded": True, "indexed": True}

    monkeypatch.setattr(cli_module, "prepare_pfam_database", fake_prepare_pfam_database)

    result = runner.invoke(app, ["dbs", "prepare", "pfam"])

    assert result.exit_code == 0
    assert "Prepared pfam" in result.stdout
    assert "Index ready" in result.stdout


def test_dbs_refresh_calls_pfam_refresh(monkeypatch, tmp_path):
    hmm_path = tmp_path / "Pfam-A.hmm"
    hmm_path.write_text("dummy")
    offsets_path = tmp_path / "offsets.json"

    def fake_refresh_pfam_database():
        offsets_path.write_text("{}")
        return {"hmm_path": str(hmm_path), "offsets_path": str(offsets_path), "downloaded": True, "indexed": True}

    monkeypatch.setattr(cli_module, "refresh_pfam_database", fake_refresh_pfam_database)

    result = runner.invoke(app, ["dbs", "refresh", "pfam"])

    assert result.exit_code == 0
    assert "Refreshed pfam" in result.stdout
    assert "Index ready" in result.stdout


def test_dbs_status_shows_pfam_state(monkeypatch):
    def fake_status():
        return {
            "name": "pfam",
            "downloaded": True,
            "indexed": True,
            "manifest_exists": True,
            "model_count": 100,
            "sparse_cached_count": 5,
            "root": "/tmp/pfam",
        }

    monkeypatch.setattr(cli_module, "get_pfam_database_status", fake_status)

    result = runner.invoke(app, ["dbs", "status", "pfam"])

    assert result.exit_code == 0
    assert "pfam:" in result.stdout
    assert "downloaded: True" in result.stdout
    assert "indexed: True" in result.stdout


def test_dbs_remove_requires_yes():
    result = runner.invoke(app, ["dbs", "remove", "pfam"])
    assert result.exit_code == 1


def test_dbs_remove_calls_pfam_remove(monkeypatch):
    monkeypatch.setattr(cli_module, "remove_pfam_database", lambda: True)

    result = runner.invoke(app, ["dbs", "remove", "pfam", "--yes"])

    assert result.exit_code == 0
    assert "Removed pfam database" in result.stdout


def test_cluster_short_options_present_in_help():
    result = runner.invoke(app, ["cluster", "--help"])
    assert result.exit_code == 0
    assert "--input-contigs" in result.stdout
    assert "-i" in result.stdout
    assert "--output-folder" in result.stdout
    assert "-o" in result.stdout
    assert "--threads" in result.stdout
    assert "-t" in result.stdout


def test_screen_short_options_present_in_help():
    result = runner.invoke(app, ["screen", "--help"], terminal_width=200)
    assert result.exit_code == 0
    assert "--input-contigs" in result.stdout
    assert "-i" in result.stdout
    assert "--mode" in result.stdout
    assert "-m" in result.stdout
    assert "-g" in result.stdout
    assert "--combine-mode" in result.stdout
    assert "-c" in result.stdout
    assert "--cut-ga" in result.stdout


def test_jack_short_options_present_in_help():
    result = runner.invoke(app, ["jack", "--help"], terminal_width=200)
    assert result.exit_code == 0
    assert "--input-contigs" in result.stdout
    assert "-i" in result.stdout
    assert "--iterations" in result.stdout
    assert "-I" in result.stdout
    assert "--max-evalue" in result.stdout
    assert "-e" in result.stdout
    assert "--combine-mode" in result.stdout
    assert "-c" in result.stdout
    assert "--min-seed-hits" in result.stdout
    assert "-k" in result.stdout


def test_simplify_taxa_short_options_present_in_help():
    result = runner.invoke(app, ["simplify-taxa", "--help"])
    assert result.exit_code == 0
    assert "--input-file" in result.stdout
    assert "-i" in result.stdout
    assert "--output-file" in result.stdout
    assert "-o" in result.stdout
    assert "--sep" in result.stdout
    assert "-s" in result.stdout
