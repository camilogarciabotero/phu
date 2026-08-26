from __future__ import annotations

from io import StringIO

from rich.console import Console

from phu._click import ProgressReporter


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
