from typer.testing import CliRunner
from phu.cli import app

runner = CliRunner()


def test_help_runs():
    result = runner.invoke(app, ["seqclust", "--help"])
    assert result.exit_code == 0
    assert "seqclust" in result.stdout
