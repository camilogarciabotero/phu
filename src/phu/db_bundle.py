"""Install and locate the pinned phu-db bundle.

Pins (tarball digest and size, bundle version, manifest digest and the full
member set) are passed in as a BundlePins value; nothing here trusts what a
tarball says about itself. extract_bundle reads the whole archive once to check
it against the pins and the vendored manifest, and only then writes anything.
"""

from __future__ import annotations

import hashlib
import shutil
import tarfile
import tempfile
import zlib
from dataclasses import dataclass
from pathlib import Path

from phu.db_fetch import DatabaseError
from phu.db_manifest import BUNDLE_ROOT, MANIFEST_MEMBER, Manifest
from phu.db_prepare_verbs import ArchiveError, UnsafeArchive, check_tar_members, extract_member

# Extraction directories inside the db folder start with this.
TEMP_PREFIX = ".phu-db-tmp-"

# What a truncated, corrupt or non-gzip tarball raises while it is read.
_READ_ERRORS = (tarfile.TarError, OSError, EOFError, zlib.error)


@dataclass(frozen=True)
class BundlePins:
    url: str
    sha256: str
    size: int
    bundle_version: str
    manifest_sha256: str
    expected_members: frozenset[str]  # paths like "phu-db/pfam/profiles.hmm"


class BundleMismatch(DatabaseError):
    """The bundle, the vendored manifest and the pins do not agree."""


@dataclass(frozen=True)
class ExtractedBundle:
    """A checked bundle unpacked into ``tempdir``, which the caller must remove."""

    tempdir: Path
    databases: frozenset[str]

    @property
    def root(self) -> Path:
        """The unpacked ``phu-db/`` directory: phu-db.json plus one directory per database."""
        return self.tempdir / BUNDLE_ROOT


def _check_pins(pins: BundlePins, manifest: Manifest, vendored_manifest_bytes: bytes) -> None:
    """The vendored manifest must be the one the pins were cut from."""
    if hashlib.sha256(vendored_manifest_bytes).hexdigest() != pins.manifest_sha256:
        raise BundleMismatch("the vendored phu-db.json does not match its pinned sha256")
    if manifest.bundle_version != pins.bundle_version:
        raise BundleMismatch(
            f"the vendored manifest is {manifest.bundle_version!r}, "
            f"the pins are {pins.bundle_version!r}"
        )
    if manifest.expected_members() != pins.expected_members:
        raise BundleMismatch("the vendored manifest and the pins list different bundle members")


def _scan(
    tarball: Path, expected: frozenset[str], manifest_size: int
) -> tuple[list[tarfile.TarInfo], bytes | None]:
    """Read every member header, and the bundled phu-db.json if its size is right."""
    bundled_manifest = None
    try:
        with tarfile.open(tarball, "r|gz") as tar:
            for member in tar:
                if (
                    member.name == MANIFEST_MEMBER
                    and member.isreg()
                    and member.size == manifest_size
                    and bundled_manifest is None
                ):
                    bundled_manifest = tar.extractfile(member).read()
            members = check_tar_members(tar, expected)
    except _READ_ERRORS as error:
        raise ArchiveError(f"cannot read {tarball}: {error}") from error
    return members, bundled_manifest


def _databases_present(names: set[str], manifest: Manifest) -> frozenset[str]:
    """Contract section 1: each database is entirely present or entirely absent."""
    present = set()
    for name, entry in manifest.databases.items():
        members = entry.member_paths()
        found = members & names
        if found == members:
            present.add(name)
        elif found:
            raise UnsafeArchive(
                f"bundle holds only part of the {name} database; "
                f"missing {sorted(members - found)}"
            )
    return frozenset(present)


def _check_sizes(
    members: dict[str, tarfile.TarInfo], manifest: Manifest, databases: frozenset[str]
) -> None:
    for name in sorted(databases):
        for item in manifest.databases[name].prepared_files:
            member = members[manifest.member_path(name, item.path)]
            if member.size != item.bytes:
                raise BundleMismatch(
                    f"bundle member {member.name} is {member.size} bytes, "
                    f"the manifest says {item.bytes}"
                )


def _identity(member: tarfile.TarInfo) -> tuple[str, bytes, int]:
    return member.name, member.type, member.size


def _unpack(tarball: Path, members: list[tarfile.TarInfo], dest: Path) -> None:
    """Extract ``members`` in one streaming pass, refusing anything that changed."""
    expected = [_identity(member) for member in members]
    seen = 0
    try:
        with tarfile.open(tarball, "r|gz") as tar:
            for member in tar:
                # The file could be swapped between the scan and this pass.
                if seen >= len(expected) or _identity(member) != expected[seen]:
                    raise UnsafeArchive(f"{tarball} changed while it was being extracted")
                extract_member(tar, member, dest, regular_only=True)
                seen += 1
    except _READ_ERRORS as error:
        raise ArchiveError(f"cannot extract {tarball}: {error}") from error
    if seen != len(expected):
        raise UnsafeArchive(f"{tarball} changed while it was being extracted")


def extract_bundle(
    tarball: Path,
    db_folder: Path,
    pins: BundlePins,
    manifest: Manifest,
    vendored_manifest_bytes: bytes,
) -> ExtractedBundle:
    """Check ``tarball`` against the pins and the vendored manifest, then unpack it.

    Nothing is written until every member has been checked: regular files only,
    each in ``pins.expected_members``, ``phu-db/phu-db.json`` byte-identical to
    the vendored copy, and every database either complete (with the manifest's
    byte sizes) or absent. The bundle is unpacked into a fresh ``.phu-db-tmp-*``
    directory inside ``db_folder``, which is removed again on any failure.
    Databases absent from the bundle are left to the caller.

    ``tarball`` is not hashed here: the caller must already have checked it
    against ``pins.size`` and ``pins.sha256`` (download_verified, or the same
    check on PHU_DB_BUNDLE), and must deep-validate the extracted tree.
    """
    tarball = Path(tarball)
    db_folder = Path(db_folder)
    _check_pins(pins, manifest, vendored_manifest_bytes)

    members, bundled_manifest = _scan(tarball, pins.expected_members, len(vendored_manifest_bytes))
    by_name = {member.name: member for member in members}
    if MANIFEST_MEMBER not in by_name:
        raise UnsafeArchive(f"bundle has no {MANIFEST_MEMBER}")
    if bundled_manifest != vendored_manifest_bytes:
        raise BundleMismatch(f"{MANIFEST_MEMBER} in the bundle differs from the vendored copy")
    databases = _databases_present(set(by_name), manifest)
    _check_sizes(by_name, manifest, databases)

    db_folder.mkdir(parents=True, exist_ok=True)
    tempdir = Path(tempfile.mkdtemp(dir=db_folder, prefix=TEMP_PREFIX))
    try:
        _unpack(tarball, members, tempdir)
        if (tempdir / MANIFEST_MEMBER).read_bytes() != vendored_manifest_bytes:
            raise BundleMismatch(f"extracted {MANIFEST_MEMBER} differs from the vendored copy")
    except BaseException:
        shutil.rmtree(tempdir, ignore_errors=True)
        raise
    return ExtractedBundle(tempdir=tempdir, databases=databases)
