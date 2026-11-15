"""`phu clean`: remove what phu keeps between runs.

Two kinds of data, removed separately:

- cache: predicted proteins reused across runs, and, for read-only shared
  database installs, the models and indexes built from them. Rebuilt on the
  next run that needs them.
- dbs: the installed phu-db bundle and any partial download. Getting it back
  means downloading and unpacking about 2 GB again (`phu dbs prepare`).
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

from phu.core.gene_prediction import get_cache_dir
from phu.db.bundle import DOWNLOADS_DIR, default_db_folder

BUNDLE_DIR = "phu-db"  # install_path(): <db folder>/phu-db/<bundle version>


@dataclass(frozen=True)
class Target:
    kind: str  # "cache" or "dbs"
    path: Path
    size: int  # bytes


def _size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).lstat().st_size
            except OSError:
                pass
    return total


def find_targets(kinds: tuple[str, ...] = ("cache", "dbs")) -> list[Target]:
    """What `phu clean` would remove for ``kinds``, existing paths only."""
    paths: list[tuple[str, Path]] = []
    if "cache" in kinds:
        paths.append(("cache", get_cache_dir()))
    if "dbs" in kinds:
        folder = default_db_folder()
        paths += [
            ("dbs", folder / BUNDLE_DIR),
            ("dbs", folder / DOWNLOADS_DIR),
        ]
    return [Target(kind, path, _size(path)) for kind, path in paths if path.exists()]


def remove(target: Target) -> None:
    if target.path.is_dir() and not target.path.is_symlink():
        shutil.rmtree(target.path)
    else:
        target.path.unlink()


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GB"
