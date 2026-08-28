import re

from typer.testing import CliRunner

import phu.cli as cli_module
from phu.cli import app

runner = CliRunner()


def plain_output(result) -> str:
    return re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", result.stdout)


def test_root_help_runs():
    result = runner.invoke(app, ["--help"])
    output = plain_output(result)
    assert result.exit_code == 0
    assert "Workflow" in output
    assert "Database Management" in output
    assert "cluster" in output
    assert "screen" in output
    assert "avger" in output
    assert "jack" in output
    assert "dbs" in output
    assert "simplify-taxa" in output
    assert "--clean-cache" in output


def test_clean_cache_removes_prediction_cache(tmp_path, monkeypatch):
    cache_dir = tmp_path / "phu-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "dummy.txt").write_text("cached")
    monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

    result = runner.invoke(app, ["--clean-cache"])

    assert result.exit_code == 0
    assert "Removed cache directory" in result.stdout
    assert not cache_dir.exists()


def test_avg_help_shows_core_options():
    result = runner.invoke(app, ["avger", "--help"], terminal_width=400)
    output = re.sub(
        r"\x1b\[[0-?]*[ -/]*[@-~]", "", result.stdout + result.stderr
    )
    assert result.exit_code == 0
    assert "Usage: root avger" in output
    assert "Predict and curate putative auxiliary viral genes" in output


def test_avg_fails_before_prediction_when_databases_are_missing(tmp_path, monkeypatch):
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c1\nATG\n")
    monkeypatch.setattr(cli_module, "run_avg", lambda config: (_ for _ in ()).throw(
        FileNotFoundError(
            "Required databases are not ready: pfam, avg. "
            "Run: phu dbs prepare pfam kofam avg"
        )
    ))

    result = runner.invoke(app, ["avger", "-i", str(contigs)])

    assert result.exit_code == 1
    assert "Required databases are not ready: pfam, avg" in result.stderr
    assert "phu dbs prepare pfam kofam avg" in result.stderr


def test_dbs_help_runs():
    result = runner.invoke(app, ["dbs", "--help"])
    output = plain_output(result)
    assert result.exit_code == 0
    assert "list" in output
    assert "status" in output
    assert "prepare" in output
    assert "refresh" in output
    assert "remove" in output


def test_dbs_list_includes_avg(monkeypatch):
    monkeypatch.setattr(
        cli_module,
        "get_avg_database_status",
        lambda: {"downloaded": True, "indexed": True},
    )

    result = runner.invoke(app, ["dbs", "list"])

    assert result.exit_code == 0
    assert "avg\tready" in result.stdout


def test_dbs_prepare_calls_avg_prepare(monkeypatch):
    called = {}

    def fake_prepare(*, force_refresh):
        called["force_refresh"] = force_refresh
        return {"root": "/tmp/avg", "release": "v1.1.1"}

    monkeypatch.setattr(cli_module, "ensure_avg_database", fake_prepare)

    result = runner.invoke(app, ["dbs", "prepare", "avg"])

    assert result.exit_code == 0
    assert called["force_refresh"] is False
    assert "Prepared avg" in result.stdout


def test_dbs_remove_calls_avg_remove(monkeypatch):
    monkeypatch.setattr(cli_module, "remove_avg_database", lambda: True)

    result = runner.invoke(app, ["dbs", "remove", "avg", "--yes"])

    assert result.exit_code == 0
    assert "Removed avg database" in result.stdout


def test_dbs_prepare_calls_pfam_prepare(monkeypatch, tmp_path):
    hmm_path = tmp_path / "Pfam-A.hmm"
    offsets_path = tmp_path / "offsets.json"

    def fake_prepare_pfam_database(*, download, index, force_refresh):
        assert download is True
        assert index is True
        assert force_refresh is False
        hmm_path.write_text("dummy")
        return {
            "hmm_path": str(hmm_path),
            "offsets_path": str(offsets_path),
            "downloaded": True,
            "indexed": True,
        }

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
        return {
            "hmm_path": str(hmm_path),
            "offsets_path": str(offsets_path),
            "downloaded": True,
            "indexed": True,
        }

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
    output = plain_output(result)
    assert result.exit_code == 0
    assert "--input-contigs" in output
    assert "-i" in output
    assert "--output-folder" in output
    assert "-o" in output
    assert "--threads" in output
    assert "-t" in output


def test_screen_short_options_present_in_help():
    result = runner.invoke(app, ["screen", "--help"], terminal_width=200)
    output = plain_output(result)
    assert result.exit_code == 0
    assert "--input-contigs" in output
    assert "-i" in output
    assert "--mode" in output
    assert "-m" in output
    assert "-g" in output
    assert "--combine-mode" in output
    assert "-c" in output
    assert "--cut-ga" in output


def test_screen_cut_ga_enabled_by_default(tmp_path, monkeypatch):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATGATGATG\n")
    hmm = tmp_path / "model.hmm"
    hmm.write_text("HMMER3/f\n")

    captured = {}

    def fake_screen(cfg):
        captured["cut_ga"] = cfg.cut_ga

    monkeypatch.setattr(cli_module, "_screen", fake_screen)

    result = runner.invoke(app, ["screen", "-i", str(contigs), str(hmm)])

    assert result.exit_code == 0
    assert captured["cut_ga"] is True


def test_screen_value_error_exits_cleanly(tmp_path, monkeypatch):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATGATGATG\n")
    hmm = tmp_path / "model.hmm"
    hmm.write_text("HMMER3/f\n")

    def fake_screen(cfg):
        raise ValueError("No HMM models found in the supplied files")

    monkeypatch.setattr(cli_module, "_screen", fake_screen)

    result = runner.invoke(app, ["screen", "-i", str(contigs), str(hmm)])

    assert result.exit_code == 1
    assert "No HMM models found in the supplied files" in result.stderr


def test_jack_short_options_present_in_help():
    result = runner.invoke(app, ["jack", "--help"], terminal_width=200)
    output = plain_output(result)
    assert result.exit_code == 0
    assert "--input-contigs" in output
    assert "-i" in output
    assert "--iterations" in output
    assert "-I" in output
    assert "--max-evalue" in output
    assert "-e" in output
    assert "--combine-mode" in output
    assert "-c" in output
    assert "--min-seed-hits" in output
    assert "-k" in output


def test_simplify_taxa_short_options_present_in_help():
    result = runner.invoke(app, ["simplify-taxa", "--help"])
    output = plain_output(result)
    assert result.exit_code == 0
    assert "--input-file" in output
    assert "-i" in output
    assert "--output-file" in output
    assert "-o" in output
    assert "--sep" in output
    assert "-s" in output
