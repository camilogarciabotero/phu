from __future__ import annotations

import re
from datetime import datetime
from io import StringIO

import pytest
from rich.console import Console
from typer.testing import CliRunner

from phu.core import log as _log
from phu.core import progress as _click
from phu.core.progress import ProgressReporter, download_progress, run_click_task
from phu.core.exec import run
from phu.cli import app

ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
runner = CliRunner()


@pytest.fixture(autouse=True)
def default_verbosity():
    _log.set_verbosity()
    yield
    _log.set_verbosity()


class FakeClock:
    def __init__(self, start: float = 100.0) -> None:
        self.value = start

    def __call__(self) -> float:
        return self.value


def test_format_duration() -> None:
    assert _log.format_duration(0) == "0.0s"
    assert _log.format_duration(42.04) == "42.0s"
    assert _log.format_duration(134) == "2m14s"
    assert _log.format_duration(3723) == "1h02m03s"


def test_elapsed_time_is_real_with_injected_clock() -> None:
    clock = FakeClock()
    stream = StringIO()
    reporter = ProgressReporter(
        console=Console(file=stream, force_terminal=False),
        clock=clock,
        now=lambda: datetime(2026, 1, 2, 14, 2, 11),
    )
    task = reporter.start_task("Predicting proteins")
    clock.value += 60
    assert reporter.elapsed(task) == 60  # still running
    clock.value += 74
    reporter.succeed_task(task)
    clock.value += 1000  # time after the end does not count
    assert reporter.elapsed(task) == 134
    reporter.finish()

    lines = stream.getvalue().splitlines()
    assert lines[0] == "[14:02:11] Predicting proteins ..."
    assert lines[1] == "[14:02:11] Predicting proteins done in 2m14s"


def test_long_plain_lines_are_not_wrapped() -> None:
    stream = StringIO()
    reporter = ProgressReporter(
        console=Console(file=stream, force_terminal=False, width=40)
    )
    reporter.fail_task(reporter.start_task("Step"), "x" * 200)
    assert len(stream.getvalue().splitlines()) == 2


def test_run_click_task_follows_process_verbosity(capsys) -> None:
    assert run_click_task("Step", lambda: 7) == 7
    out = capsys.readouterr().out
    assert re.search(r"^\[\d\d:\d\d:\d\d\] Step \.\.\.$", out, re.M)
    assert re.search(r"^\[\d\d:\d\d:\d\d\] Step done in \d+\.\ds$", out, re.M)

    _log.set_verbosity(quiet=True)
    assert run_click_task("Step", lambda: 7) == 7
    assert capsys.readouterr().out == ""


def test_run_click_task_reports_failure_and_reraises(capsys) -> None:
    def boom() -> None:
        raise ValueError("bad input")

    with pytest.raises(ValueError):
        run_click_task("Step", boom)
    assert re.search(r"Step failed after \d+\.\ds: bad input", capsys.readouterr().out)


def test_verbose_prints_exact_argv_for_external_command(capfd) -> None:
    argv = ["true", "two words", "--flag=1", "it's"]

    run(argv)
    assert capfd.readouterr().out == ""

    _log.set_verbosity(verbose=True)
    run(argv)
    assert capfd.readouterr().out == "[cmd] true 'two words' --flag=1 'it'\"'\"'s'\n"

    _log.set_verbosity(quiet=True, verbose=True)  # quiet wins
    run(argv)
    assert capfd.readouterr().out == ""


def test_verbose_command_line_shows_redirect_and_in_process_calls(capsys) -> None:
    _log.log_command(
        ["seqkit", "grep", "-f", "ids.txt", "in.fa"], stdout="out dir/x.fa"
    )
    _log.log_call("pyhmmer.hmmsearch", models=3, cpus=4)
    assert capsys.readouterr().out == ""

    _log.set_verbosity(verbose=True)
    _log.log_command(
        ["seqkit", "grep", "-f", "ids.txt", "in.fa"], stdout="out dir/x.fa"
    )
    _log.log_call("pyhmmer.hmmsearch", models=3, cpus=4, bit_cutoffs="gathering")
    assert capsys.readouterr().out.splitlines() == [
        "[cmd] seqkit grep -f ids.txt in.fa > 'out dir/x.fa'",
        "[cmd] pyhmmer.hmmsearch models=3 cpus=4 bit_cutoffs=gathering",
    ]


def test_gene_prediction_states_its_parameters_when_verbose(tmp_path, capsys) -> None:
    from phu.core.gene_prediction import PredictionInputs, predict_genes_pyrodigal

    contigs = tmp_path / "contigs.fna"
    contigs.write_text(
        ">c1\n" + "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACC" * 20 + "\n"
    )
    inputs = PredictionInputs(input_contigs=contigs, mode="meta")

    predict_genes_pyrodigal(inputs, threads=1)
    assert "[cmd]" not in capsys.readouterr().out

    _log.set_verbosity(verbose=True)
    predict_genes_pyrodigal(inputs, threads=1)
    line = capsys.readouterr().out.strip()
    assert line.startswith("[cmd] pyrodigal_gv.ViralGeneFinder contigs=1 mode=meta ")
    assert "translation_table=auto" in line
    assert "threads=1" in line


def test_cli_logs_step_and_quiet_silences_it(tmp_path) -> None:
    table = tmp_path / "in.tsv"
    table.write_text("genome\tfamily\ng1\tnovel_family_2_of_Caudoviricetes\n")
    args = ["simplify-taxa", "-i", str(table), "-o", str(tmp_path / "out.tsv")]

    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert re.search(
        r"^\[\d\d:\d\d:\d\d\] Simplifying taxonomy done in \d+\.\ds$",
        result.stdout,
        re.M,
    )
    assert "\x1b" not in result.stdout

    quiet = runner.invoke(app, [*args, "--quiet"])
    assert quiet.exit_code == 0, quiet.output
    assert "Simplifying taxonomy" not in quiet.stdout
    assert not _log.is_quiet()  # the flag does not leak into later invocations


def test_download_progress_logs_every_tenth_when_not_a_terminal() -> None:
    lines: list[str] = []
    console = Console(file=StringIO(), force_terminal=False)

    with download_progress("downloading", 2_000_000_000, log=lines.append, console=console) as update:
        for done in range(0, 2_000_000_001, 50_000_000):
            update(done)

    assert len(lines) == 11
    assert lines[0] == "downloading: 0% (0.0 MB of 2.00 GB)"
    assert lines[5] == "downloading: 50% (1.00 GB of 2.00 GB)"
    assert lines[-1] == "downloading: 100% (2.00 GB of 2.00 GB)"


def test_download_progress_draws_a_spinner_and_bar_on_a_terminal() -> None:
    out = StringIO()
    console = Console(file=out, force_terminal=True, width=100)
    lines: list[str] = []

    with download_progress("downloading phu-db", 1000, log=lines.append, console=console) as update:
        update(250)
        update(1000)

    assert lines == []
    assert "downloading phu-db" in out.getvalue()
    assert "1.0/1.0 kB" in out.getvalue()
    assert "MB/s" in out.getvalue()


def test_download_progress_writes_into_the_step_already_on_screen(monkeypatch) -> None:
    reporter = ProgressReporter(console=Console(file=StringIO(), force_terminal=False))
    task_id = reporter.start_task("Preparing PFAM database")
    monkeypatch.setattr(_click, "_LIVE_REPORTERS", [reporter])
    lines: list[str] = []

    with download_progress("downloading", 2_000_000_000, log=lines.append) as update:
        update(500_000_000)
        assert reporter._tasks[task_id].detail == "downloading 500.0 MB of 2.00 GB (25%)"

    assert reporter._tasks[task_id].detail == ""
    assert lines == []


def test_download_progress_is_silent_when_quiet() -> None:
    _log.set_verbosity(quiet=True)
    lines: list[str] = []
    out = StringIO()

    with download_progress(
        "downloading", 100, log=lines.append, console=Console(file=out, force_terminal=True)
    ) as update:
        update(100)

    assert lines == [] and out.getvalue() == ""
