import re

import pytest
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
    assert "cluster" in output
    assert "screen" in output
    assert "avger" in output
    assert "jack" in output
    assert "dbs" in output
    assert "simplify-taxa" not in output  # deprecated: hidden from the listing
    assert "--clean-cache" in output


def test_deprecated_simplify_taxa_help_still_runs():
    result = runner.invoke(app, ["simplify-taxa", "--help"])
    output = plain_output(result)
    assert result.exit_code == 0
    assert "Usage:" in output
    assert "simplify-taxa" in output
    assert "Deprecated" in output


def test_simplify_taxa_prints_deprecation_notice_on_stderr(tmp_path):
    table = tmp_path / "in.tsv"
    table.write_text("genome\tfamily\ng1\tnovel_family_2_of_Caudoviricetes\n")
    out = tmp_path / "out.tsv"

    result = runner.invoke(app, ["simplify-taxa", "-i", str(table), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert out.exists()
    assert "deprecated" in result.stderr
    assert "deprecated" not in result.stdout


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
            "Required databases are not ready: pfam, avg."
        )
    ))

    result = runner.invoke(app, ["avger", "-i", str(contigs)])

    assert result.exit_code == 1
    assert "Required databases are not ready: pfam, avg" in result.stderr


@pytest.mark.parametrize(
    "args",
    [["list"], ["status", "pfam"], ["refresh", "--all"], ["remove", "--all", "--yes"]],
)
def test_removed_dbs_commands_point_to_prepare(args):
    result = runner.invoke(app, ["dbs", *args])

    assert result.exit_code == 0
    assert result.stdout == ""
    assert "phu dbs prepare" in result.stderr


def test_dbs_help_lists_only_prepare():
    result = runner.invoke(app, ["dbs", "--help"])
    output = plain_output(result)

    assert result.exit_code == 0
    assert "prepare" in output
    assert "refresh" not in output and "remove" not in output


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


@pytest.mark.parametrize("command", [["avger"], ["screen", "PF00001"]])
def test_database_errors_are_reported_without_a_traceback(tmp_path, monkeypatch, command):
    from phu.db_fetch import DatabaseError

    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c1\nATG\n")

    def fail(*args, **kwargs):
        raise DatabaseError("phu-db could not be installed: no space")

    monkeypatch.setattr(cli_module, "run_avg", fail)
    monkeypatch.setattr(cli_module, "_screen", fail)

    result = runner.invoke(app, [*command, "-i", str(contigs), "-o", str(tmp_path / "out")])

    assert result.exit_code == 1
    assert "phu-db could not be installed: no space" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "command",
    ["screen", "jack", "cluster", "avger", "normalize-lineage", "simplify-taxa"],
)
def test_commands_have_quiet_and_verbose(command: str) -> None:
    result = runner.invoke(app, [command, "--help"], terminal_width=400)
    output = plain_output(result)
    assert result.exit_code == 0
    assert "--quiet" in output
    assert "--verbose" in output


def test_screen_reports_an_os_error_from_the_database_install(tmp_path, monkeypatch):
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c1\nATG\n")

    def full_disk(config):
        raise OSError(122, "Disk quota exceeded", str(tmp_path / "db"))

    monkeypatch.setattr(cli_module, "_screen", full_disk)

    result = runner.invoke(app, ["screen", "-i", str(contigs), "PF00001"])

    assert result.exit_code == 1
    assert "Disk quota exceeded" in result.stderr
    assert "Traceback" not in result.output
