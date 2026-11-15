"""How much memory a step may use before it should read from disk instead.

pyhmmer searches fastest when its targets are loaded into memory, but loading
a file larger than the machine (or the batch job) allows ends the run. These
helpers make that choice with the rule from pyhmmer's performance notes: load
a file only when it takes a small share of the memory still available.

Under a batch scheduler the job's own limit (its cgroup) is usually much
smaller than the node's free memory, so that limit is honoured when present.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import psutil

# Share of the available memory one loaded file may take.
LOAD_FRACTION = 0.1
# A digital sequence block takes up to about twice its FASTA's size (many
# short sequences carry more per-sequence overhead).
SEQUENCE_OVERHEAD = 2.0

_CGROUP_V2 = ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current")
_CGROUP_V1 = (
    "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    "/sys/fs/cgroup/memory/memory.usage_in_bytes",
)


def _read_int(path: str) -> Optional[int]:
    try:
        text = Path(path).read_text().strip()
    except OSError:
        return None
    return int(text) if text.isdigit() else None  # "max" means no limit


def _cgroup_room() -> Optional[int]:
    """Bytes left under this process's cgroup memory limit, if it has one."""
    for limit_path, usage_path in (_CGROUP_V2, _CGROUP_V1):
        limit = _read_int(limit_path)
        usage = _read_int(usage_path)
        # cgroup v1 reports "no limit" as a huge number.
        if limit is not None and usage is not None and limit < 1 << 60:
            return max(limit - usage, 0)
    return None


def available_memory() -> Optional[int]:
    """Bytes this process can still allocate, or None when it cannot be told."""
    try:
        room = psutil.virtual_memory().available
    except Exception:  # psutil can fail in unusual containers
        room = None
    cgroup = _cgroup_room()
    known = [value for value in (room, cgroup) if value is not None]
    return min(known) if known else None


def fits_in_memory(
    path: Path | str,
    *,
    overhead: float = SEQUENCE_OVERHEAD,
    fraction: float = LOAD_FRACTION,
) -> bool:
    """Whether loading ``path`` stays within ``fraction`` of available memory.

    True when the available memory cannot be measured, so behaviour is then
    what it was before this check existed.
    """
    available = available_memory()
    if available is None:
        return True
    size = os.stat(path).st_size
    return size * overhead <= available * fraction


def protein_batch_size(path: Path | str, total: int, minimum: int) -> int:
    """Proteins to search per batch: all of them when they fit, else ``minimum``."""
    if total > minimum and fits_in_memory(path):
        return total
    return minimum
