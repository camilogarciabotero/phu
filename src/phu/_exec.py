from __future__ import annotations

import subprocess
from collections.abc import Sequence
from pathlib import Path
from shutil import which


class CmdNotFound(RuntimeError):
    pass


def _executable(candidates: Sequence[str]) -> str:
    """
    Find an executable in the PATH.
    """
    for name in candidates:
        path = which(name)
        if path:
            return path
    raise CmdNotFound(f"None of the executables found in PATH: {', '.join(candidates)}")


def run(
    argv: Sequence[str],
    cwd: str | Path | None = None,
    check: bool = True,
    timeout: float | None = None,
) -> subprocess.CompletedProcess[str]:
    """
    Run a command in a subprocess.
    """
    proc = subprocess.run(
        list(map(str, argv)),
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=False,  # stream to console
        check=check,
        timeout=timeout,
    )
    return proc
