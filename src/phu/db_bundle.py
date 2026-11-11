"""Install and locate the pinned phu-db bundle.

Pins (tarball digest and size, bundle version, manifest digest and the full
member set) are passed in as a BundlePins value; nothing here trusts what a
tarball says about itself. extract_bundle reads the whole archive once to check
it against the pins and the vendored manifest, and only then writes anything.

install() runs under an exclusive lock on the db folder: it verifies the tarball
(downloaded, or the user's PHU_DB_BUNDLE) like a download, extracts it, deep
validates the result against the manifest and moves the finished tree into
``phu-db/<bundle_version>`` in one rename, so no reader ever sees half a tree.
A bundle that leaves a database out is refused.

resolve_db() is what commands call: it checks the installed tree (shallowly, or
deeply with PHU_DB_VERIFY=1), installs when that fails, and returns the path of
one database. bundle_provenance() reports what is installed, for run.json.
"""

from __future__ import annotations

import errno
import hashlib
import importlib
import importlib.resources
import json
import os
import shutil
import socket
import sys
import tarfile
import tempfile
import time
import zlib
from collections.abc import Callable, Iterable, Iterator
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from phu.db_fetch import (
    CHUNK_SIZE,
    DatabaseError,
    IntegrityError,
    download_verified,
    ensure_free_space,
    required_install_bytes,
)
from phu.db_manifest import (
    BUNDLE_ROOT,
    DATABASE_NAMES,
    MANIFEST_MEMBER,
    MANIFEST_NAME,
    Manifest,
    PreparedFile,
    loads_manifest,
)
from phu.db_archive import (
    ArchiveError,
    UnsafeArchive,
    check_tar_members,
    extract_member,
)

try:
    import fcntl
except ImportError:  # not POSIX: only the lock directory is available
    fcntl = None

# Extraction directories inside the db folder start with this.
TEMP_PREFIX = ".phu-db-tmp-"
# Inside the db folder: installed trees, downloads, and the install lock.
INSTALL_DIR = BUNDLE_ROOT
DOWNLOADS_DIR = ".downloads"
TARBALL_NAME = "phu-db.tar.gz"
INSTALLED_NAME = "installed.json"
LOCK_NAME = ".phu-db.lock"
LOCK_DIR_NAME = ".phu-db.lock.d"
LOCK_OWNER_NAME = "owner.json"
# A lock directory this old was left by a killed install and may be broken.
STALE_LOCK_SECONDS = 6 * 60 * 60
LOCK_POLL_SECONDS = 2.0
BUNDLE_ENV = "PHU_DB_BUNDLE"
FOLDER_ENV = "PHU_DB_FOLDER"
VERIFY_ENV = "PHU_DB_VERIFY"

# flock errors meaning "this filesystem cannot lock", e.g. NFS without lockd.
_NO_FLOCK = frozenset({errno.ENOLCK, errno.EOPNOTSUPP, errno.ENOTSUP, errno.EINVAL})

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


class TarballMismatch(IntegrityError):
    """The tarball (downloaded or PHU_DB_BUNDLE) is not the pinned bundle."""


class TreeInvalid(DatabaseError):
    """A phu-db tree does not match its manifest."""

    def __init__(self, problems: list[str], root: Path | None = None) -> None:
        self.problems = list(problems)
        self.root = root
        shown = "; ".join(self.problems[:5])
        more = f"; and {len(self.problems) - 5} more" if len(self.problems) > 5 else ""
        where = f" {root}" if root is not None else ""
        super().__init__(f"phu-db tree{where} is invalid: {shown}{more}")


class LockTimeout(DatabaseError):
    """Another install held the lock for longer than the caller would wait."""


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


# Validation.


class _Counter:
    """Contract section 5 record count, fed in chunks.

    Lines end at b"\n" only. A ``tsv`` line counts when its strip() is non-empty
    (header included); an ``hmm`` line counts when it starts with b"NAME ".
    """

    def __init__(self, file_format: str) -> None:
        self.hmm = file_format == "hmm"
        self.count = 0
        # hmm: the last bytes seen, starting with a virtual b"\n" before line one.
        self._tail = b"\n"
        # tsv: whether the unfinished line so far is all whitespace.
        self._blank = True

    def feed(self, chunk: bytes) -> None:
        if self.hmm:
            # Every NAME line is one b"\nNAME "; the carried five bytes let a
            # match split across chunks be seen, and only once.
            data = self._tail + chunk
            self.count += data.count(b"\nNAME ")
            self._tail = data[-5:]
            return
        lines = chunk.split(b"\n")
        for line in lines[:-1]:
            if not self._blank or line.strip():
                self.count += 1
            self._blank = True
        if lines[-1].strip():
            self._blank = False

    def total(self) -> int:
        if not self.hmm and not self._blank:
            return self.count + 1
        return self.count


def _scan_file(path: Path, *, digest: bool, file_format: str | None) -> tuple[str | None, int | None]:
    """Stream ``path`` once; return its sha256 (if ``digest``) and record count."""
    sha256 = hashlib.sha256() if digest else None
    counter = _Counter(file_format) if file_format else None
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            if sha256 is not None:
                sha256.update(chunk)
            if counter is not None:
                counter.feed(chunk)
    return (
        sha256.hexdigest() if sha256 is not None else None,
        counter.total() if counter is not None else None,
    )


def _file_problems(path: Path, label: str, item: PreparedFile, deep: bool) -> list[str]:
    try:
        size = path.stat().st_size if path.is_file() else None
    except FileNotFoundError:
        # Moved or removed by an install while this unlocked check ran.
        size = None
    except OSError as error:
        return [f"{label} cannot be read: {error}"]
    if size is None:
        return [f"{label} is missing"]
    if size != item.bytes:
        return [f"{label} is {size} bytes, expected {item.bytes}"]
    # Shallow checks count TSV records only; HMMs are large, so they wait for deep.
    counted = item.format == "tsv" or deep
    try:
        sha256, count = _scan_file(path, digest=deep, file_format=item.format if counted else None)
    except OSError as error:
        return [f"{label} cannot be read: {error}"]
    problems = []
    if sha256 is not None and sha256 != item.sha256:
        problems.append(f"{label} has sha256 {sha256}, expected {item.sha256}")
    if count is not None and count != item.count:
        problems.append(f"{label} has {count} records, expected {item.count}")
    return problems


def tree_problems(
    root: Path, manifest: Manifest, *, deep: bool, databases: Iterable[str] | None = None
) -> list[str]:
    """What is wrong with the prepared files under ``root``; empty when valid.

    ``root`` holds ``<layout_dir>/<path>`` for every prepared file of
    ``databases`` (default: all). Shallow checks existence, exact byte size and
    TSV record counts; ``deep`` adds every sha256 and HMM model count.
    """
    root = Path(root)
    names = list(manifest.databases) if databases is None else list(databases)
    unknown = sorted(set(names) - set(manifest.databases))
    if unknown:
        raise ValueError(f"unknown phu-db databases: {unknown}")
    problems = []
    for name in names:
        entry = manifest.databases[name]
        for item in entry.prepared_files:
            label = f"{entry.layout_dir}/{item.path}"
            path = root.joinpath(entry.layout_dir, *item.path.split("/"))
            problems.extend(_file_problems(path, label, item, deep))
    return problems


def validate_tree(
    root: Path, manifest: Manifest, *, deep: bool, databases: Iterable[str] | None = None
) -> None:
    """Raise TreeInvalid unless every prepared file under ``root`` matches the manifest."""
    problems = tree_problems(root, manifest, deep=deep, databases=databases)
    if problems:
        raise TreeInvalid(problems, root=Path(root))


def installed_problems(
    db_folder: Path,
    pins: BundlePins,
    manifest: Manifest,
    *,
    deep: bool,
    databases: Iterable[str] | None = None,
) -> list[str]:
    """What is wrong with the installed ``phu-db/<bundle_version>``; empty when valid.

    Checks installed.json and the prepared files of ``databases`` (default: all).
    """
    target = install_path(db_folder, pins)
    if not target.is_dir():
        return [f"{target} is not installed"]
    problems = []
    try:
        record = json.loads((target / INSTALLED_NAME).read_bytes())
    except (OSError, ValueError) as error:
        problems.append(f"{INSTALLED_NAME} cannot be read: {error}")
    else:
        if not isinstance(record, dict) or record.get("bundle_version") != pins.bundle_version:
            problems.append(f"{INSTALLED_NAME} does not record {pins.bundle_version}")
    return problems + tree_problems(target, manifest, deep=deep, databases=databases)


def install_path(db_folder: Path, pins: BundlePins) -> Path:
    return Path(db_folder) / INSTALL_DIR / pins.bundle_version


# Locking.


def _log_stderr(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def _flock(
    fd: int,
    path: Path,
    *,
    log: Callable[[str], None],
    timeout: float | None,
    sleep: Callable[[float], None],
) -> bool:
    """Take an exclusive flock on ``fd``; False if this filesystem cannot lock."""
    if fcntl is None:
        return False
    deadline = None if timeout is None else time.monotonic() + timeout
    waited = False
    try:
        while True:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return True
            except BlockingIOError:
                pass
            if not waited:
                log(f"phu-db: waiting for another install to finish (lock {path})")
                waited = True
            if deadline is None:
                # The kernel releases a flock when its holder exits, even if killed.
                fcntl.flock(fd, fcntl.LOCK_EX)
                return True
            if time.monotonic() >= deadline:
                raise LockTimeout(f"{path} is still locked after {timeout} seconds")
            sleep(0.05)
    except OSError as error:
        if error.errno in _NO_FLOCK:
            return False
        raise


def _read_owner(lock_dir: Path) -> dict:
    try:
        owner = json.loads((lock_dir / LOCK_OWNER_NAME).read_bytes())
    except (OSError, ValueError):
        return {}
    return owner if isinstance(owner, dict) else {}


def _lock_age(lock_dir: Path, clock: Callable[[], float]) -> float | None:
    """Seconds since ``lock_dir`` was taken, or None if it is gone."""
    created = _read_owner(lock_dir).get("created")
    if not isinstance(created, (int, float)) or isinstance(created, bool):
        # Killed between mkdir and writing the owner file: use the directory's age.
        try:
            created = lock_dir.stat().st_mtime
        except FileNotFoundError:
            return None
    return clock() - created


def _publish_lock_dir(lock_dir: Path, owner: dict) -> bool:
    """Create ``lock_dir`` holding ``owner``; False if someone else holds it.

    The directory is built under a private name and renamed into place, so a
    lock directory never exists without its owner file: rename(2) refuses to
    replace a non-empty directory, which makes the rename the atomic step that
    os.mkdir would otherwise be.
    """
    staging = lock_dir.with_name(f"{lock_dir.name}.new-{owner['token']}")
    os.mkdir(staging)
    try:
        (staging / LOCK_OWNER_NAME).write_text(json.dumps(owner))
        os.rename(staging, lock_dir)
        return True
    except OSError as error:
        if error.errno in (errno.EEXIST, errno.ENOTEMPTY):
            return False
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _break_stale_lock(
    lock_dir: Path, clock: Callable[[], float], log: Callable[[str], None]
) -> None:
    """Remove ``lock_dir`` if it is still a stale lock once it has been moved aside.

    Several waiters may judge the same lock stale. Only one can rename it away;
    another may instead rename the fresh lock the first one has just taken, so
    the directory actually moved is judged again and a live one is put back.
    """
    aside = lock_dir.with_name(f"{lock_dir.name}.stale-{os.urandom(6).hex()}")
    try:
        os.rename(lock_dir, aside)
    except FileNotFoundError:
        return
    age = _lock_age(aside, clock)
    if age is not None and age <= STALE_LOCK_SECONDS:
        try:
            os.rename(aside, lock_dir)
            return
        except OSError:
            # A third install took the free lock in the moment it was aside.
            log(
                f"phu-db: install lock {lock_dir} was taken twice while a stale lock "
                f"was broken; one of the installs may fail"
            )
    else:
        owner = _read_owner(aside)
        hours = f"{age / 3600:.1f}" if age is not None else "?"
        log(
            f"phu-db: breaking stale install lock {lock_dir} "
            f"(pid {owner.get('pid', '?')} on {owner.get('hostname', '?')}, "
            f"{hours} hours old)"
        )
    shutil.rmtree(aside, ignore_errors=True)


@contextmanager
def _dir_lock(
    lock_dir: Path,
    *,
    log: Callable[[str], None],
    timeout: float | None,
    poll: float,
    clock: Callable[[], float],
    sleep: Callable[[float], None],
) -> Iterator[None]:
    """Hold the lock directory ``lock_dir`` for the block."""
    deadline = None if timeout is None else time.monotonic() + timeout
    waited = False
    token = os.urandom(8).hex()
    while True:
        owner = {
            "pid": os.getpid(),
            "hostname": socket.gethostname(),
            "created": clock(),
            "token": token,
        }
        if _publish_lock_dir(lock_dir, owner):
            break
        age = _lock_age(lock_dir, clock)
        if age is None:
            continue
        if age > STALE_LOCK_SECONDS:
            _break_stale_lock(lock_dir, clock, log)
            continue
        if deadline is not None and time.monotonic() >= deadline:
            raise LockTimeout(f"{lock_dir} is still held after {timeout} seconds")
        if not waited:
            log(f"phu-db: waiting for another install to finish (lock {lock_dir})")
            waited = True
        sleep(poll)
    try:
        yield
    finally:
        # A holder slower than STALE_LOCK_SECONDS may have lost the lock to a
        # waiter; never remove a lock directory someone else now owns.
        if _read_owner(lock_dir).get("token") == token:
            shutil.rmtree(lock_dir, ignore_errors=True)


@contextmanager
def install_lock(
    db_folder: Path,
    *,
    log: Callable[[str], None] = _log_stderr,
    timeout: float | None = None,
    poll: float = LOCK_POLL_SECONDS,
    clock: Callable[[], float] = time.time,
    sleep: Callable[[float], None] = time.sleep,
) -> Iterator[None]:
    """Hold the exclusive install lock on ``db_folder`` for the block.

    An flock on ``db_folder/.phu-db.lock``; where the filesystem cannot flock
    (for example NFS without lockd), a lock directory ``.phu-db.lock.d``,
    created atomically with an owner file holding the owner's pid, host and
    start time, which a waiter breaks once it is more than six hours old. Waiting is logged once.
    ``timeout`` (seconds, default forever) raises LockTimeout.
    """
    db_folder = Path(db_folder)
    db_folder.mkdir(parents=True, exist_ok=True)
    path = db_folder / LOCK_NAME
    with ExitStack() as stack:
        # The lock file is never deleted: a waiter could lock an unlinked inode.
        try:
            fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o666)
        except PermissionError:
            # Another user's lock file in a shared folder; flock needs only read.
            fd = os.open(path, os.O_RDONLY)
        stack.callback(os.close, fd)
        if not _flock(fd, path, log=log, timeout=timeout, sleep=sleep):
            log(f"phu-db: {db_folder} does not support flock; using a lock directory")
            stack.enter_context(
                _dir_lock(
                    db_folder / LOCK_DIR_NAME,
                    log=log,
                    timeout=timeout,
                    poll=poll,
                    clock=clock,
                    sleep=sleep,
                )
            )
        yield


# Install.


def _human(size: int) -> str:
    if size < 1000:
        return f"{size} bytes"
    value = float(size)
    for unit in ("KB", "MB", "GB"):
        value /= 1000
        if value < 1000 or unit == "GB":
            break
    return f"{value:.1f} {unit}"


def _expected(pins: BundlePins) -> str:
    return (
        f"phu needs {pins.bundle_version} ({TARBALL_NAME}, {pins.size} bytes, "
        f"sha256 {pins.sha256})"
    )


def _verify_tarball(path: Path, pins: BundlePins) -> None:
    """Check size, then sha256, of ``path`` against the pins, reading it only."""
    size = path.stat().st_size
    if size != pins.size:
        raise TarballMismatch(
            f"{path} is a different or corrupt bundle: {_expected(pins)}, "
            f"but the file is {size} bytes"
        )
    sha256, _ = _scan_file(path, digest=True, file_format=None)
    if sha256 != pins.sha256:
        raise TarballMismatch(
            f"{path} is a different or corrupt bundle: {_expected(pins)}, "
            f"but the file has sha256 {sha256}"
        )


def _local_bundle(bundle_path: Path, pins: BundlePins, log: Callable[[str], None]) -> Path:
    """The user's PHU_DB_BUNDLE, verified like a download and never modified."""
    if not bundle_path.is_file():
        raise DatabaseError(f"{BUNDLE_ENV}={bundle_path} is not a file")
    log(f"phu-db: verifying {pins.bundle_version} bundle {bundle_path} ({BUNDLE_ENV})")
    _verify_tarball(bundle_path, pins)
    return bundle_path


def _download_bundle(
    db_folder: Path, pins: BundlePins, fetch: Callable[..., Path], log: Callable[[str], None]
) -> Path:
    dest = db_folder / DOWNLOADS_DIR / TARBALL_NAME
    if dest.is_file():
        # Left by an install that failed after its download; reuse it if it is ours.
        try:
            _verify_tarball(dest, pins)
            log(f"phu-db: reusing downloaded {pins.bundle_version} bundle {dest}")
            return dest
        except TarballMismatch:
            dest.unlink()
    log(
        f"phu-db: downloading {pins.bundle_version} ({_human(pins.size)}) "
        f"from {pins.url} to {dest}"
    )
    try:
        fetch(pins.url, dest, pins.sha256, pins.size)
    except IntegrityError as error:
        raise TarballMismatch(
            f"the download from {pins.url} is a different or corrupt bundle: "
            f"{_expected(pins)} ({error})"
        ) from error
    # A fetch is injectable; never extract a file it did not leave verified.
    try:
        _verify_tarball(dest, pins)
    except TarballMismatch:
        dest.unlink(missing_ok=True)
        raise
    return dest


def _write_json(path: Path, data: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.urandom(6).hex()}.tmp")
    try:
        temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _remove_leftovers(db_folder: Path, log: Callable[[str], None]) -> None:
    """Remove temp directories of installs that died; the caller holds the lock."""
    for path in sorted(db_folder.glob(f"{TEMP_PREFIX}*")):
        if path.is_dir() and not path.is_symlink():
            log(f"phu-db: removing leftover {path}")
            shutil.rmtree(path, ignore_errors=True)


def _move_into_place(staged: Path, target: Path, db_folder: Path) -> None:
    """os.replace ``staged`` onto ``target``, first moving a stale ``target`` aside."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if not (target.exists() or target.is_symlink()):
        os.replace(staged, target)
        return
    aside = Path(tempfile.mkdtemp(dir=db_folder, prefix=TEMP_PREFIX))
    try:
        os.replace(target, aside / "stale")
        os.replace(staged, target)
    finally:
        shutil.rmtree(aside, ignore_errors=True)


def _remove_other_versions(target: Path, log: Callable[[str], None]) -> None:
    for path in sorted(target.parent.iterdir()):
        if path.name != target.name and path.is_dir() and not path.is_symlink():
            log(f"phu-db: removing previous install {path}")
            shutil.rmtree(path, ignore_errors=True)


def _same_file(path: Path, other: Path | None) -> bool:
    try:
        return other is not None and os.path.samefile(path, other)
    except OSError:
        return False


def _clean_up(
    db_folder: Path, target: Path, bundle_path: Path | None, log: Callable[[str], None]
) -> None:
    """Remove what installs leave behind once ``target`` is valid; under the lock.

    Run on every install() call, so whatever a run killed after its replace
    left (older versions, the downloaded tarball, temp directories) is removed
    by the next one. The user's own tarball is never touched.
    """
    _remove_leftovers(db_folder, log)
    _remove_other_versions(target, log)
    downloaded = db_folder / DOWNLOADS_DIR / TARBALL_NAME
    for path in (downloaded, downloaded.with_name(downloaded.name + ".part")):
        if path.is_file() and not _same_file(path, bundle_path):
            path.unlink(missing_ok=True)


def install(
    *,
    db_folder: Path,
    pins: BundlePins,
    manifest: Manifest,
    vendored_manifest_bytes: bytes,
    bundle_path: Path | None = None,
    fetch: Callable[..., Path] = download_verified,
    log: Callable[[str], None] = _log_stderr,
    deep: bool = False,
) -> Path:
    """Install the pinned bundle under ``db_folder``; returns ``phu-db/<bundle_version>``.

    Under the install lock: if another process has installed a valid tree in the
    meantime (checked shallowly, or with ``deep`` as the caller's own check was),
    return it. Otherwise check free space, take the tarball (``bundle_path``,
    else $PHU_DB_BUNDLE, else a download to ``db_folder/.downloads``) only after
    its size and sha256 match the pins, extract it, deep validate every file and
    replace the installed tree in one rename. Older versions and the downloaded
    tarball are removed after, and again by any later call, so an interrupted
    clean-up is finished. The user's own tarball is never modified or deleted.
    Progress goes to ``log`` (stderr by default).
    """
    db_folder = Path(db_folder)
    if bundle_path is None and os.environ.get(BUNDLE_ENV):
        bundle_path = Path(os.environ[BUNDLE_ENV])
    target = install_path(db_folder, pins)

    with install_lock(db_folder, log=log):
        if not installed_problems(db_folder, pins, manifest, deep=deep):
            log(f"phu-db: {pins.bundle_version} is already installed at {target}")
            _clean_up(db_folder, target, bundle_path, log)
            return target
        _remove_leftovers(db_folder, log)

        prepared_bytes = sum(
            item.bytes for entry in manifest.databases.values() for item in entry.prepared_files
        )
        # A local tarball is read in place; only a download needs room for it.
        tarball_bytes = 0 if bundle_path is not None else pins.size
        ensure_free_space(db_folder, required_install_bytes(tarball_bytes, prepared_bytes))

        if bundle_path is not None:
            tarball = _local_bundle(Path(bundle_path), pins, log)
        else:
            tarball = _download_bundle(db_folder, pins, fetch, log)

        log(f"phu-db: extracting {pins.bundle_version} into {db_folder}")
        extracted = extract_bundle(tarball, db_folder, pins, manifest, vendored_manifest_bytes)
        try:
            missing = sorted(set(manifest.databases) - extracted.databases)
            if missing:
                raise BundleMismatch(
                    f"the {pins.bundle_version} bundle does not include {', '.join(missing)}"
                )
            validate_tree(extracted.root, manifest, deep=True)
            _write_json(
                extracted.root / INSTALLED_NAME,
                {
                    "bundle_version": pins.bundle_version,
                    "source": "local" if bundle_path is not None else "zenodo",
                    "tarball_sha256": pins.sha256,
                    "installed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                },
            )
            _move_into_place(extracted.root, target, db_folder)
        finally:
            shutil.rmtree(extracted.tempdir, ignore_errors=True)

        _clean_up(db_folder, target, bundle_path, log)
        log(f"phu-db: installed {pins.bundle_version} at {target}")
    return target


# Locating databases.


def default_db_folder() -> Path:
    """$PHU_DB_FOLDER, else $XDG_DATA_HOME/phu/db, else ~/.local/share/phu/db."""
    if db_env := os.environ.get(FOLDER_ENV):
        return Path(db_env)

    if xdg_data_home := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg_data_home) / "phu" / "db"
    return Path.home() / ".local" / "share" / "phu" / "db"


_NOT_PINNED = (
    "phu-db bundle is not pinned yet: this phu has no phu._bundle pins or vendored "
    "phu/data/phu-db.json, so it cannot install or locate its databases"
)


def _production_pins() -> tuple[BundlePins, bytes]:
    """The pins and vendored manifest shipped with phu."""
    try:
        bundle = importlib.import_module("phu._bundle")
    except ModuleNotFoundError as error:
        # Only a missing pins module means "not pinned"; a broken one is a bug.
        if error.name != "phu._bundle":
            raise
        raise DatabaseError(_NOT_PINNED) from None
    # One joinpath step at a time: Traversable takes a single child before 3.11.
    vendored = importlib.resources.files("phu").joinpath("data").joinpath(MANIFEST_NAME)
    if not vendored.is_file():
        raise DatabaseError(_NOT_PINNED)
    return bundle.PRODUCTION_PINS, vendored.read_bytes()


def _pinned(
    pins: BundlePins | None, vendored_manifest_bytes: bytes | None
) -> tuple[BundlePins, Manifest, bytes]:
    if pins is None and vendored_manifest_bytes is None:
        pins, vendored_manifest_bytes = _production_pins()
    elif pins is None or vendored_manifest_bytes is None:
        raise ValueError("pins and vendored_manifest_bytes must be given together")
    manifest = loads_manifest(vendored_manifest_bytes)
    # Every entry point, not only an install, uses a manifest the pins vouch for.
    _check_pins(pins, manifest, vendored_manifest_bytes)
    return pins, manifest, vendored_manifest_bytes


def resolve_db(
    name: str,
    *,
    pins: BundlePins | None = None,
    db_folder: Path | None = None,
    vendored_manifest_bytes: bytes | None = None,
    fetch: Callable[..., Path] = download_verified,
    log: Callable[[str], None] = _log_stderr,
) -> Path:
    """Return the directory of database ``name``, installing the bundle first if needed.

    ``name`` is one of pfam, kofam, dbcan, avg. The installed tree under
    ``db_folder`` (default: default_db_folder()) is checked for that database,
    shallowly, or by sha256 with PHU_DB_VERIFY=1, and install() runs when the
    check fails.
    ``pins`` and ``vendored_manifest_bytes`` default to the ones shipped with
    phu; tests pass both.
    """
    if name not in DATABASE_NAMES:
        raise ValueError(
            f"unknown phu-db database {name!r}; expected one of {sorted(DATABASE_NAMES)}"
        )
    pins, manifest, vendored_manifest_bytes = _pinned(pins, vendored_manifest_bytes)
    db_folder = default_db_folder() if db_folder is None else Path(db_folder)
    deep = os.environ.get(VERIFY_ENV) == "1"
    target = install_path(db_folder, pins)
    # Only the requested database: a run resolving several hashes each one once.
    if installed_problems(db_folder, pins, manifest, deep=deep, databases=[name]):
        target = install(
            db_folder=db_folder,
            pins=pins,
            manifest=manifest,
            vendored_manifest_bytes=vendored_manifest_bytes,
            fetch=fetch,
            log=log,
            deep=deep,
        )
    return target / manifest.databases[name].layout_dir


def prepare_bundle(
    *,
    pins: BundlePins | None = None,
    db_folder: Path | None = None,
    vendored_manifest_bytes: bytes | None = None,
    bundle_path: Path | None = None,
    deep: bool = False,
    fetch: Callable[..., Path] = download_verified,
    log: Callable[[str], None] = _log_stderr,
) -> Path:
    """Install the whole pinned bundle ahead of time; returns ``phu-db/<bundle_version>``.

    What ``phu dbs prepare`` runs. Every database is checked, shallowly, or by
    sha256 with ``deep`` or PHU_DB_VERIFY=1, and install() runs when the check
    fails. ``bundle_path`` (else $PHU_DB_BUNDLE) is a pre-downloaded tarball.
    """
    pins, manifest, vendored_manifest_bytes = _pinned(pins, vendored_manifest_bytes)
    db_folder = default_db_folder() if db_folder is None else Path(db_folder)
    deep = deep or os.environ.get(VERIFY_ENV) == "1"
    if not installed_problems(db_folder, pins, manifest, deep=deep):
        return install_path(db_folder, pins)
    return install(
        db_folder=db_folder,
        pins=pins,
        manifest=manifest,
        vendored_manifest_bytes=vendored_manifest_bytes,
        bundle_path=bundle_path,
        fetch=fetch,
        log=log,
        deep=deep,
    )


def database_hashes(
    name: str,
    *,
    pins: BundlePins | None = None,
    vendored_manifest_bytes: bytes | None = None,
) -> dict[str, str]:
    """Prepared file -> pinned sha256 for database ``name``; never installs."""
    _, manifest, _ = _pinned(pins, vendored_manifest_bytes)
    if name not in manifest.databases:
        raise ValueError(f"unknown phu-db database {name!r}")
    return {item.path: item.sha256 for item in manifest.databases[name].prepared_files}


def bundle_provenance(
    *,
    pins: BundlePins | None = None,
    db_folder: Path | None = None,
    vendored_manifest_bytes: bytes | None = None,
) -> dict[str, object]:
    """Where the installed databases came from, for .phu/run.json; never installs.

    Returns ``source`` ("zenodo" or "local"), ``bundle_version``, ``sha256``
    of the tarball and ``upstream_versions`` (database -> its manifest
    upstream_version).
    """
    pins, manifest, _ = _pinned(pins, vendored_manifest_bytes)
    db_folder = default_db_folder() if db_folder is None else Path(db_folder)
    path = install_path(db_folder, pins) / INSTALLED_NAME
    try:
        record = json.loads(path.read_bytes())
    except (OSError, ValueError) as error:
        raise DatabaseError(f"phu-db {pins.bundle_version} is not installed ({error})") from error
    if not isinstance(record, dict) or record.get("bundle_version") != pins.bundle_version:
        raise DatabaseError(f"{path} does not record {pins.bundle_version}")
    return {
        "source": record.get("source"),
        "bundle_version": pins.bundle_version,
        "sha256": record.get("tarball_sha256"),
        "upstream_versions": {
            name: entry.upstream_version for name, entry in manifest.databases.items()
        },
    }
