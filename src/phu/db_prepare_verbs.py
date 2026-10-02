"""Prepare verbs for phu-db upstream sources, and the tar member checks they share.

Every tar archive phu unpacks, the phu-db bundle and upstream archives alike,
passes check_tar_members before a single byte is written: no links, devices or
fifos, and no member name that could land outside the extraction directory.
"""

from __future__ import annotations

import os
import shutil
import tarfile
from collections import Counter
from pathlib import Path

from phu.db_fetch import DatabaseError

COPY_CHUNK = 1024 * 1024


class ArchiveError(DatabaseError):
    """A tar archive could not be read or extracted."""


class UnsafeArchive(ArchiveError):
    """A tar archive holds a member that must not be extracted."""


def _check_name(name: str) -> None:
    """Reject a member name that is absolute, climbs out, or is ambiguous."""
    if not name or "\\" in name or "\0" in name or name.startswith("/"):
        raise UnsafeArchive(f"unsafe archive member name: {name!r}")
    if any(part in ("", "..") for part in name.split("/")):
        raise UnsafeArchive(f"unsafe archive member name: {name!r}")


def _kind(member: tarfile.TarInfo) -> str:
    if member.issym():
        return "symbolic link"
    if member.islnk():
        return "hard link"
    if member.ischr() or member.isblk():
        return "device"
    if member.isfifo():
        return "fifo"
    if member.isdir():
        return "directory"
    return f"member of tar type {member.type!r}"


def check_member(member: tarfile.TarInfo, *, regular_only: bool) -> None:
    """Raise UnsafeArchive unless ``member`` is a safe regular file (or directory)."""
    _check_name(member.name)
    if member.isreg() or (member.isdir() and not regular_only):
        return
    allowed = "regular files" if regular_only else "regular files and directories"
    raise UnsafeArchive(
        f"archive member {member.name!r} is a {_kind(member)}; only {allowed} may be extracted"
    )


def check_tar_members(
    tar: tarfile.TarFile, expected: frozenset[str] | None = None
) -> list[tarfile.TarInfo]:
    """Check every member of ``tar`` and return them, in archive order.

    With ``expected`` (a bundle), members must be regular files, each named in
    ``expected`` and none twice; which subsets of ``expected`` are acceptable is
    the caller's rule. Without it (an upstream archive), directories are also
    allowed, as in the contract's ``extract`` verb.
    """
    members = tar.getmembers()
    for member in members:
        check_member(member, regular_only=expected is not None)
    if expected is not None:
        counts = Counter(member.name for member in members)
        duplicates = sorted(name for name, count in counts.items() if count > 1)
        if duplicates:
            raise UnsafeArchive(f"archive members appear more than once: {duplicates}")
        extra = sorted(set(counts) - expected)
        if extra:
            raise UnsafeArchive(f"unexpected archive members: {extra}")
    return members


def extract_member(
    tar: tarfile.TarFile, member: tarfile.TarInfo, dest: Path, *, regular_only: bool = False
) -> None:
    """Check one member with check_member, then write it below ``dest``.

    The check runs here too, so both write paths refuse the same members whatever
    the caller checked before. Uses tarfile's ``data`` filter where this Python
    has it; older 3.10 and 3.11 patch releases get an equivalent manual write.
    Either way files get the current time and default permissions, not the
    archive's: an epoch mtime is the first thing scratch-filesystem purges delete.
    """
    check_member(member, regular_only=regular_only)
    # Looked up per call so tests can take the filter away.
    if hasattr(tarfile, "data_filter"):
        tar.extract(member, dest, set_attrs=False, filter="data")
        return
    target = Path(dest).joinpath(*member.name.split("/"))
    if member.isdir():
        target.mkdir(parents=True, exist_ok=True)
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    # Every directory above target was made here from checked names, so only the
    # final component could be a link; O_NOFOLLOW refuses to write through one.
    flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC | getattr(os, "O_NOFOLLOW", 0)
    source = tar.extractfile(member)
    with source, os.fdopen(os.open(target, flags, 0o666), "wb") as handle:
        shutil.copyfileobj(source, handle, COPY_CHUNK)
