from __future__ import annotations

import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from shutil import which

from phu.core.log import log_command


class CmdNotFound(RuntimeError):
    pass


def _executable(candidates: Sequence[str]) -> str:
    """
    Find an executable beside the running Python, then on PATH.

    Tools installed as phu's dependencies (vclust) land in the environment's
    bin/ directory, which `uv tool install` and pipx do not put on PATH.
    """
    env_bin = str(Path(sys.executable).parent)
    for name in candidates:
        path = which(name, path=env_bin) or which(name)
        if path:
            return path
    raise CmdNotFound(
        f"None of the executables found beside Python or in PATH: {', '.join(candidates)}"
    )


def run(
    argv: Sequence[str],
    cwd: str | Path | None = None,
    check: bool = True,
    timeout: float | None = None,
) -> subprocess.CompletedProcess[str]:
    """
    Run a command in a subprocess.

    With ``--verbose`` the exact argv is printed first (see ``phu.core.log``).
    """
    command = list(map(str, argv))
    log_command(command, cwd=cwd)
    proc = subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=False,  # stream to console
        check=check,
        timeout=timeout,
    )
    return proc
