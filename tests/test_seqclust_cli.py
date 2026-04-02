from typer.testing import CliRunner
from phu.cli import app

runner = CliRunner()


def test_root_help_runs():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "cluster" in result.stdout
    assert "screen" in result.stdout
    assert "jack" in result.stdout
    assert "simplify-taxa" in result.stdout


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
    result = runner.invoke(app, ["screen", "--help"])
    assert result.exit_code == 0
    assert "--input-contigs" in result.stdout
    assert "-i" in result.stdout
    assert "--mode" in result.stdout
    assert "-m" in result.stdout
    assert "--min-protein-len-aa" in result.stdout
    assert "-g" in result.stdout
    assert "--combine-mode" in result.stdout
    assert "-c" in result.stdout


def test_jack_short_options_present_in_help():
    result = runner.invoke(app, ["jack", "--help"])
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
