from __future__ import annotations

import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from threading import Lock
from typing import Any

from rich.console import Console
from rich.live import Live
from rich.markup import escape
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)
from rich.spinner import Spinner
from rich.table import Table

from ._log import format_duration, is_quiet, is_verbose

_SYMBOLS = {"success": "✓", "failed": "✗", "skipped": "○", "cancelled": "!"}
_PLAIN_ENDINGS = {
    "success": "done in",
    "failed": "failed after",
    "skipped": "skipped after",
    "cancelled": "cancelled after",
}


# Reporters whose live display is on screen now; a nested step writes into the
# innermost one instead of starting a second live display on the same terminal.
_LIVE_REPORTERS: list["ProgressReporter"] = []


@dataclass
class _Task:
    label: str
    started: float
    status: str = "running"
    detail: str = ""
    completed: int = 0
    total: int | None = None
    ended: float | None = None


class ProgressReporter:
    """Small, thread-safe reporter for blocking CLI work.

    On a terminal it shows a live line per step with its elapsed time. When
    the output is not a terminal (a SLURM log, a pipe) it prints one plain
    line when a step starts and one when it ends::

        [14:02:11] Predicting proteins ...
        [14:04:25] Predicting proteins done in 2m14s

    ``quiet`` and ``verbose`` default to the process-wide setting made by the
    CLI (see ``phu._log``). ``clock`` (monotonic seconds) and ``now``
    (wall-clock ``datetime``) can be replaced in tests.
    """

    def __init__(
        self,
        *,
        console: Console | None = None,
        quiet: bool | None = None,
        verbose: bool | None = None,
        clock: Callable[[], float] = time.monotonic,
        now: Callable[[], datetime] = datetime.now,
    ) -> None:
        self.console = console or Console()
        self.quiet = is_quiet() if quiet is None else quiet
        self.verbose = is_verbose() if verbose is None else verbose
        self._clock = clock
        self._now = now
        self._tasks: dict[TaskID, _Task] = {}
        self._spinner = Spinner("dots")
        self._lock = Lock()
        self._next_id = 0
        self._live: Live | None = None

    def elapsed(self, task_id: TaskID) -> float:
        """Seconds the task has run so far, or its total once it has ended."""
        task = self._tasks[task_id]
        end = task.ended if task.ended is not None else self._clock()
        return max(0.0, end - task.started)

    def _plain(self, message: str) -> None:
        stamp = self._now().strftime("%H:%M:%S")
        self.console.print(
            f"[{stamp}] {message}", markup=False, highlight=False, soft_wrap=True
        )

    def _render(self) -> Table:
        # Called by the Live refresh thread, so elapsed times keep ticking.
        # It must not take self._lock: Live holds its own lock while rendering.
        grid = Table.grid(padding=(0, 1))
        for task_id, task in list(self._tasks.items()):
            mark: Any = _SYMBOLS.get(task.status, self._spinner)
            text = escape(task.label)
            if task.total and task.status == "running":
                text += f"  {task.completed}/{task.total}"
            if task.detail:
                text += f"  {escape(task.detail)}"
            grid.add_row(mark, text, format_duration(self.elapsed(task_id)))
        return grid

    def start_phase(self, label: str, *, announce: bool = True) -> None:
        if self.quiet:
            return
        if self.console.is_terminal:
            self._live = Live(
                get_renderable=self._render,
                console=self.console,
                refresh_per_second=8,
            )
            self._live.start()
            _LIVE_REPORTERS.append(self)
        elif announce:
            self._plain(label)

    def start_task(self, label: str, *, total: int | None = None) -> TaskID:
        with self._lock:
            self._next_id += 1
            task_id = TaskID(self._next_id)
            self._tasks[task_id] = _Task(label, started=self._clock(), total=total)
            if not self.quiet and not self.console.is_terminal:
                self._plain(f"{label} ...")
            return task_id

    def update_task(
        self, task_id: TaskID, *, completed: int | None = None, total: int | None = None
    ) -> None:
        with self._lock:
            task = self._tasks[task_id]
            if completed is not None:
                task.completed = completed
            if total is not None:
                task.total = total

    def set_running_detail(self, detail: str) -> None:
        """Show ``detail`` beside every step that is still running."""
        with self._lock:
            for task in self._tasks.values():
                if task.status == "running":
                    task.detail = detail

    def _finish(self, task_id: TaskID, status: str, detail: str = "") -> None:
        with self._lock:
            task = self._tasks[task_id]
            if task.ended is None:
                task.ended = self._clock()
            task.status, task.detail = status, detail
            if not self.quiet and not self.console.is_terminal:
                line = (
                    f"{task.label} {_PLAIN_ENDINGS[status]} "
                    f"{format_duration(self.elapsed(task_id))}"
                )
                if detail:
                    line += f": {' '.join(detail.split())}"
                self._plain(line)

    def succeed_task(self, task_id: TaskID) -> None:
        self._finish(task_id, "success")

    def fail_task(self, task_id: TaskID, detail: str = "") -> None:
        self._finish(task_id, "failed", detail)

    def skip_task(self, task_id: TaskID, detail: str = "skipped") -> None:
        self._finish(task_id, "skipped", detail)

    def cancel_task(self, task_id: TaskID, detail: str = "cancelled") -> None:
        self._finish(task_id, "cancelled", detail)

    def fail_running_tasks(self, detail: str) -> None:
        for task_id, task in list(self._tasks.items()):
            if task.status == "running":
                self.fail_task(task_id, detail)

    def finish(self, *, summary: bool = True) -> None:
        with self._lock:
            live, self._live = self._live, None
        if live:
            if self in _LIVE_REPORTERS:
                _LIVE_REPORTERS.remove(self)
            live.stop()  # renders the final state, with the real durations
        if self.quiet or not summary:
            return
        counts = {
            status: sum(task.status == status for task in self._tasks.values())
            for status in ("success", "skipped", "failed", "cancelled")
        }
        self.console.print(f"Completed: {counts['success']}")
        self.console.print(f"Skipped: {counts['skipped']}")
        self.console.print(f"Failed: {counts['failed'] + counts['cancelled']}")


def run_click_task(
    label: str,
    func: Callable[..., Any],
    *args: Any,
    quiet: bool | None = None,
    verbose: bool | None = None,
    **kwargs: Any,
) -> Any:
    """Run a blocking task as one timed step.

    ``quiet`` and ``verbose`` default to the process-wide setting, so nested
    steps follow the command's ``--quiet`` / ``--verbose`` flags. The step's
    own final line carries its status and duration, so no summary is printed.
    """
    reporter = ProgressReporter(quiet=quiet, verbose=verbose)
    reporter.start_phase(label, announce=False)
    task_id = reporter.start_task(label)
    try:
        result = func(*args, **kwargs)
    except BaseException as exc:
        reporter.fail_task(task_id, str(exc))
        reporter.finish(summary=False)
        raise
    reporter.succeed_task(task_id)
    reporter.finish(summary=False)
    return result


def _megabytes(size: int) -> str:
    return f"{size / 1e9:.2f} GB" if size >= 1e9 else f"{size / 1e6:.1f} MB"


@contextmanager
def download_progress(
    label: str,
    total: int,
    *,
    log: Callable[[str], None],
    console: Console | None = None,
    step_percent: int = 10,
) -> Iterator[Callable[[int], None]]:
    """Show a download of ``total`` bytes; yields ``update(bytes_so_far)``.

    On a terminal: a rich spinner with a bar, size, speed and time left, on
    stderr, removed when the download ends. Inside a step that already shows a
    spinner, the percentage is added to that step's line instead. When stderr
    is not a terminal (a SLURM log, a pipe), one ``log`` line per
    ``step_percent`` of the file. Nothing is shown with ``--quiet``.
    """
    if is_quiet():
        yield lambda done: None
        return

    if _LIVE_REPORTERS:
        reporter = _LIVE_REPORTERS[-1]

        def update(done: int) -> None:
            percent = 100 * done // total if total else 100
            reporter.set_running_detail(
                f"{label} {_megabytes(done)} of {_megabytes(total)} ({percent}%)"
            )

        try:
            yield update
        finally:
            reporter.set_running_detail("")
        return

    console = console or Console(stderr=True)
    if console.is_terminal:
        progress = Progress(
            SpinnerColumn(),
            TextColumn("{task.description}"),
            BarColumn(),
            DownloadColumn(),
            TransferSpeedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=True,
        )
        with progress:
            task_id = progress.add_task(escape(label), total=total)
            yield lambda done: progress.update(task_id, completed=done)
        return

    reported = -1

    def update(done: int) -> None:
        nonlocal reported
        step = (100 * done // total if total else 100) // step_percent
        if step != reported:
            reported = step
            log(
                f"{label}: {step * step_percent}% "
                f"({_megabytes(done)} of {_megabytes(total)})"
            )

    yield update
