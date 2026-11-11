from __future__ import annotations

import re
import time
from io import StringIO

from rich.console import Console

from phu._click import ProgressReporter

ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def test_reporter_plain_summary_and_states() -> None:
    stream = StringIO()
    reporter = ProgressReporter(console=Console(file=stream, force_terminal=False))
    reporter.start_phase("Pipeline")
    completed = reporter.start_task("Measured", total=10)
    running = reporter.start_task("Waiting")
    reporter.update_task(completed, completed=10)
    reporter.succeed_task(completed)
    reporter.skip_task(running, "not requested")
    reporter.finish()

    output = stream.getvalue()
    assert "Pipeline" in output
    assert "Completed: 1" in output
    assert "Skipped: 1" in output
    assert "Failed: 0" in output
    assert "\\x1b[" not in output


def test_reporter_failure_and_quiet_mode() -> None:
    stream = StringIO()
    reporter = ProgressReporter(
        console=Console(file=stream, force_terminal=False), quiet=True
    )
    task_id = reporter.start_task("Command")
    reporter.fail_task(task_id, "exit code 2")
    reporter.finish()

    assert stream.getvalue() == ""


def test_terminal_final_line_shows_real_elapsed_time() -> None:
    # Uses only the constructor arguments that exist on main, where the final
    # line always read 0:00:00.
    stream = StringIO()
    reporter = ProgressReporter(
        console=Console(file=stream, force_terminal=True, width=80)
    )
    reporter.start_phase("Pipeline")
    task = reporter.start_task("Measured")
    time.sleep(0.25)
    reporter.succeed_task(task)
    reporter.finish()

    final = [
        line
        for line in ANSI.sub("", stream.getvalue()).replace("\r", "\n").splitlines()
        if "✓ Measured" in line
    ][-1]
    match = re.search(r"(\d+\.\d)s", final)
    assert match, final
    assert 0.2 <= float(match.group(1)) < 5


def test_non_terminal_prints_start_and_end_lines_without_ansi() -> None:
    stream = StringIO()
    reporter = ProgressReporter(console=Console(file=stream, force_terminal=False))
    ok = reporter.start_task("Predicting proteins")
    reporter.succeed_task(ok)
    bad = reporter.start_task("Annotating proteins")
    reporter.fail_task(bad, "exit code 2\nsecond line")
    reporter.finish()

    output = stream.getvalue()
    lines = output.splitlines()
    stamp = r"\[\d\d:\d\d:\d\d\] "
    assert re.fullmatch(stamp + r"Predicting proteins \.\.\.", lines[0])
    assert re.fullmatch(stamp + r"Predicting proteins done in \d+\.\ds", lines[1])
    assert re.fullmatch(stamp + r"Annotating proteins \.\.\.", lines[2])
    assert re.fullmatch(
        stamp + r"Annotating proteins failed after \d+\.\ds: exit code 2 second line",
        lines[3],
    )
    assert "\x1b" not in output
