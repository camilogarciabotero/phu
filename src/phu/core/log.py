"""Process-wide verbosity for console output.

The CLI calls :func:`set_verbosity` once, before any work starts. Code that
runs a tool then calls :func:`log_command` (external programs) or
:func:`log_call` (in-process libraries); both print nothing unless
``--verbose`` was given. Every such line starts with ``[cmd] `` so a run log
can be filtered with ``grep '^\\[cmd\\]'``.

The setting is written once and only read afterwards, so worker threads can
call the logging helpers freely; a lock keeps each line whole.
"""

from __future__ import annotations

import shlex
import sys
from collections.abc import Sequence
from pathlib import Path
from threading import Lock

CMD_PREFIX = "[cmd] "

_quiet = False
_verbose = False
_write_lock = Lock()


def set_verbosity(*, quiet: bool = False, verbose: bool = False) -> None:
    """Set the verbosity for this process. ``quiet`` wins over ``verbose``."""
    global _quiet, _verbose
    _quiet = bool(quiet)
    _verbose = bool(verbose) and not quiet


def is_quiet() -> bool:
    return _quiet


def is_verbose() -> bool:
    return _verbose


def format_duration(seconds: float) -> str:
    """Compact duration: ``0.4s``, ``42.0s``, ``2m14s``, ``1h02m03s``."""
    seconds = max(0.0, float(seconds))
    if seconds < 59.95:
        return f"{seconds:.1f}s"
    whole = int(round(seconds))
    hours, rest = divmod(whole, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{secs:02d}s"
    return f"{minutes}m{secs:02d}s"


def _emit(text: str) -> None:
    # Flush so the line lands before the output of a child process that
    # writes to the same file descriptor.
    with _write_lock:
        sys.stdout.write(CMD_PREFIX + text + "\n")
        sys.stdout.flush()


def log_command(
    argv: Sequence[object],
    *,
    stdout: str | Path | None = None,
    cwd: str | Path | None = None,
) -> None:
    """With ``--verbose``, print the exact argv of an external command."""
    if not _verbose:
        return
    line = shlex.join(str(arg) for arg in argv)
    if stdout is not None:
        line += f" > {shlex.quote(str(stdout))}"
    if cwd is not None:
        line = f"(cd {shlex.quote(str(cwd))} && {line})"
    _emit(line)


def log_call(name: str, **params: object) -> None:
    """With ``--verbose``, print the effective parameters of an in-process step."""
    if not _verbose:
        return
    _emit(" ".join([name, *(f"{key}={value}" for key, value in params.items())]))
