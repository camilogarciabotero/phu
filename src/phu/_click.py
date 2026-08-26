from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from threading import Lock
from typing import Any

from rich.console import Console
from rich.live import Live
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.progress import TaskID


@dataclass
class _Task:
    label: str
    status: str = "running"
    detail: str = ""
    completed: int = 0
    total: int | None = None


class ProgressReporter:
    """Small, thread-safe Rich reporter for blocking CLI work."""

    def __init__(self, *, console: Console | None = None, quiet: bool = False, verbose: bool = False) -> None:
        self.console = console or Console()
        self.quiet = quiet
        self.verbose = verbose
        self._tasks: dict[TaskID, _Task] = {}
        self._lock = Lock()
        self._next_id = 0
        self._live: Live | None = None

    def _render(self) -> Progress:
        progress = Progress(
            SpinnerColumn(), TextColumn("{task.description}"), TimeElapsedColumn(),
            console=self.console,
        )
        for task in self._tasks.values():
            symbol = {"success": "✓", "failed": "✗", "skipped": "○", "cancelled": "!"}.get(task.status, "⠋")
            progress.add_task(f"{symbol} {task.label}" + (f"  {task.detail}" if task.detail else ""), total=task.total, completed=task.completed)
        return progress

    def start_phase(self, label: str) -> None:
        if self.quiet:
            return
        if self.console.is_terminal:
            self._live = Live(self._render(), console=self.console, refresh_per_second=8)
            self._live.start()
        else:
            self.console.print(label)

    def start_task(self, label: str, *, total: int | None = None) -> TaskID:
        with self._lock:
            self._next_id += 1
            task_id = TaskID(self._next_id)
            self._tasks[task_id] = _Task(label, total=total)
            if self._live:
                self._live.update(self._render())
            return task_id

    def update_task(self, task_id: TaskID, *, completed: int | None = None, total: int | None = None) -> None:
        with self._lock:
            task = self._tasks[task_id]
            if completed is not None:
                task.completed = completed
            if total is not None:
                task.total = total
            if self._live:
                self._live.update(self._render())

    def _finish(self, task_id: TaskID, status: str, detail: str = "") -> None:
        with self._lock:
            task = self._tasks[task_id]
            task.status, task.detail = status, detail
            if self._live:
                self._live.update(self._render())

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

    def finish(self) -> None:
        with self._lock:
            if self._live:
                self._live.update(self._render())
                self._live.stop()
                self._live = None
            if self.quiet:
                return
            counts = {status: sum(task.status == status for task in self._tasks.values()) for status in ("success", "skipped", "failed", "cancelled")}
            self.console.print(f"Completed: {counts['success']}")
            self.console.print(f"Skipped: {counts['skipped']}")
            self.console.print(f"Failed: {counts['failed'] + counts['cancelled']}")


def run_click_task(
    label: str,
    func: Callable[..., Any],
    *args: Any,
    quiet: bool = False,
    verbose: bool = False,
    **kwargs: Any,
) -> Any:
    """Run a blocking task with a Rich progress indicator."""
    reporter = ProgressReporter(quiet=quiet, verbose=verbose)
    reporter.start_phase(label)
    task_id = reporter.start_task(label)
    try:
        result = func(*args, **kwargs)
    except BaseException as exc:
        reporter.fail_task(task_id, str(exc))
        reporter.finish()
        raise
    reporter.succeed_task(task_id)
    reporter.finish()
    return result
