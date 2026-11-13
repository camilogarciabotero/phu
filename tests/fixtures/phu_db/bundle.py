"""Deterministic phu-db bundles built from the committed fixture tree.

The ``fixture_bundle`` fixture in tests/conftest.py packs prepared/ under
``phu-db/`` plus phu-db.json, the way phu-db-prep's pack does: regular files
only, sorted, mtime 0, owner root. Negative tests start from bundle_files(),
drop or replace entries and add hand-made members through write_bundle.
"""

from __future__ import annotations

import gzip
import io
import tarfile
from collections.abc import Iterable, Mapping
from pathlib import Path

FIXTURE_DIR = Path(__file__).parent
BUNDLE_ROOT = "phu-db"
URL = "https://phu-db-fixture.invalid/phu-db.tar.gz"


def bundle_files() -> dict[str, bytes]:
    """Member name -> bytes of the full fixture bundle (16 members)."""
    files = {f"{BUNDLE_ROOT}/phu-db.json": (FIXTURE_DIR / "phu-db.json").read_bytes()}
    prepared = FIXTURE_DIR / "prepared"
    for path in sorted(prepared.rglob("*")):
        if path.is_file():
            files[f"{BUNDLE_ROOT}/{path.relative_to(prepared).as_posix()}"] = path.read_bytes()
    return dict(sorted(files.items()))


def member(name: str, type: bytes = tarfile.REGTYPE, *, linkname: str = "") -> tarfile.TarInfo:
    """A header with the bundle's fixed metadata; size is set by write_bundle."""
    info = tarfile.TarInfo(name)
    info.type = type
    info.linkname = linkname
    info.mtime = 0
    info.uid = info.gid = 0
    info.uname = info.gname = "root"
    info.mode = 0o755 if type == tarfile.DIRTYPE else 0o644
    return info


def write_bundle(
    path: Path,
    files: Mapping[str, bytes],
    extra: Iterable[tuple[tarfile.TarInfo, bytes | None]] = (),
) -> Path:
    """Write ``files`` as regular members plus ``extra`` members, sorted by name."""
    entries = [(member(name), data) for name, data in files.items()] + list(extra)
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.GNU_FORMAT) as archive:
        for info, data in sorted(entries, key=lambda entry: entry[0].name):
            if data is None:
                archive.addfile(info)
            else:
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
    path = Path(path)
    path.write_bytes(gzip.compress(buffer.getvalue(), mtime=0))
    return path
