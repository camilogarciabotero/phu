"""Binary copies of the bundle's HMM files, for fast loading.

The bundle ships HMMs as text, which is portable but slow to parse: reading
Pfam, KOfam and PHROG as text takes minutes, as binary a few seconds. The
binary format depends on the machine and on the HMMER version, so the copies
are built here on first use, kept with the other derived files of the
database, and rebuilt whenever the text file or pyhmmer changes.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pyhmmer

from phu.core.log import log_call
from phu.db.bundle import derived_dir, derived_files_lock, resolve_db, write_atomically

# Databases whose whole model file is searched (by avger), so worth a copy.
SEARCHED_DATABASES = ("pfam", "kofam", "phrog")
TEXT_NAME = "profiles.hmm"

BINARY_NAME = "profiles.h3m"
STAMP_NAME = "profiles.h3m.json"
_SCHEMA_VERSION = 1


def _stamp(text_path: Path) -> dict[str, object]:
    stat = text_path.stat()
    return {
        "schema_version": _SCHEMA_VERSION,
        "source_path": str(text_path),
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
        "pyhmmer": pyhmmer.__version__,
    }


def _is_current(binary: Path, stamp_path: Path, expected: dict[str, object]) -> bool:
    if not binary.is_file() or not stamp_path.is_file():
        return False
    try:
        return json.loads(stamp_path.read_text()) == expected
    except (OSError, ValueError):
        return False


def has_binary_models(name: str) -> bool:
    """Whether database ``name`` of the installed bundle has a current binary copy."""
    text_path = resolve_db(name) / TEXT_NAME
    place = derived_dir(text_path.parent)
    try:
        return _is_current(place / BINARY_NAME, place / STAMP_NAME, _stamp(text_path))
    except OSError:
        return False


def build_binary_models(name: str) -> Path:
    """Build (if needed) the binary copy of database ``name``; returns the file to read."""
    return binary_hmm_path(resolve_db(name) / TEXT_NAME)


def _build(text_path: Path, binary: Path) -> int:
    """Write every model of ``text_path`` to ``binary`` through a rename."""
    binary.parent.mkdir(parents=True, exist_ok=True)
    temporary = binary.with_name(f".{binary.name}.{os.urandom(6).hex()}.tmp")
    models = 0
    try:
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
        with os.fdopen(fd, "wb") as output, pyhmmer.plan7.HMMFile(text_path) as hmm_file:
            for hmm in hmm_file:
                hmm.write(output, binary=True)
                models += 1
        os.replace(temporary, binary)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return models


def binary_hmm_path(text_path: Path) -> Path:
    """The binary copy of an installed ``profiles.hmm``, built if needed.

    Returns ``text_path`` itself when the copy cannot be written (no space, no
    writable place): the search is then slower, never wrong.
    """
    text_path = Path(text_path)
    place = derived_dir(text_path.parent)
    binary = place / BINARY_NAME
    stamp_path = place / STAMP_NAME
    try:
        expected = _stamp(text_path)
        if _is_current(binary, stamp_path, expected):
            return binary
        place.mkdir(parents=True, exist_ok=True)
        with derived_files_lock(place):
            # Another job may have built it while this one waited.
            if not _is_current(binary, stamp_path, expected):
                stamp_path.unlink(missing_ok=True)
                models = _build(text_path, binary)
                write_atomically(stamp_path, json.dumps(expected, indent=2).encode())
                log_call("phu.db.binary", source=text_path, binary=binary, models=models)
    except OSError:
        return text_path
    return binary
