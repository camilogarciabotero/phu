"""Prepare verbs for phu-db upstream sources, and the tar member checks they share.

Every tar archive phu unpacks, the phu-db bundle and upstream archives alike,
passes check_tar_members before a single byte is written: no links, devices or
fifos, and no member name that could land outside the extraction directory.

The verbs turn a database's verified downloads into its prepared files with the
exact semantics of the phu-db 0.1.1 contract (section 3), which phu-db-prep's
prepare step shares, so a database fetched upstream hashes like the bundled one.
Every verb writes its output atomically: on failure no partial target is left.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import tarfile
import tempfile
import zipfile
import zlib
from collections import Counter
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO
from xml.etree import ElementTree

from phu.db_fetch import DatabaseError
from phu.db_manifest import ManifestError, Step, parse_step

COPY_CHUNK = 1024 * 1024


class ArchiveError(DatabaseError):
    """A tar archive could not be read or extracted."""


class UnsafeArchive(ArchiveError):
    """A tar archive holds a member that must not be extracted."""


class PrepareError(DatabaseError):
    """A prepare step could not turn its source into its target."""


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


# Prepare verbs. Each takes and returns plain paths and raises PrepareError.


@contextmanager
def _atomic_output(dst: Path) -> Iterator[BinaryIO]:
    """Yield a binary handle whose bytes replace ``dst`` only if the block succeeds.

    The temporary file sits beside ``dst`` (same filesystem, so os.replace is
    atomic) and is created with default permissions, unlike mkstemp's 0600.
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary = dst.with_name(f".{dst.name}.{os.urandom(6).hex()}.tmp")
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
    try:
        with os.fdopen(fd, "wb") as handle:
            yield handle
        os.replace(temporary, dst)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def gunzip(src: Path, dst: Path) -> Path:
    """Decompress the gzip file ``src`` to ``dst``."""
    try:
        with gzip.open(src, "rb") as stream, _atomic_output(dst) as sink:
            shutil.copyfileobj(stream, sink, COPY_CHUNK)
    # BadGzipFile is an OSError; a truncated stream raises EOFError.
    except (OSError, EOFError, zlib.error) as error:
        raise PrepareError(f"cannot gunzip {Path(src).name}: {error}") from error
    return Path(dst)


def _merge(staging: Path, workdir: Path) -> None:
    """Move everything below ``staging`` to the same place below ``workdir``."""
    for directory, subdirectories, files in os.walk(staging):
        relative = Path(directory).relative_to(staging)
        for name in subdirectories:
            (workdir / relative / name).mkdir(exist_ok=True)
        for name in files:
            os.replace(Path(directory, name), workdir / relative / name)


def extract(archive: Path, workdir: Path, target_dir: Path) -> Path:
    """Unpack the tar ``archive`` into the workdir root; ``target_dir`` must result.

    Archives carry their own top-level directory (kofam's ``profiles/``), so
    the archive is not unpacked into ``target_dir`` but must produce it. Members
    pass check_tar_members first and are unpacked into a staging directory, so
    a bad member or a missing target leaves the workdir untouched.
    """
    workdir = Path(workdir)
    target = workdir / target_dir
    try:
        relative = target.relative_to(workdir)
    except ValueError:
        relative = None
    if relative is None or not relative.parts or ".." in relative.parts:
        raise PrepareError(f"extract target {target_dir} is not inside {workdir}")
    try:
        staging = Path(tempfile.mkdtemp(dir=workdir, prefix=".phu-extract-"))
    except OSError as error:
        raise PrepareError(f"cannot extract {Path(archive).name}: {error}") from error
    try:
        try:
            with tarfile.open(archive, "r:*") as tar:
                for member in check_tar_members(tar, None):
                    extract_member(tar, member, staging)
        except (OSError, EOFError, zlib.error, tarfile.TarError, DatabaseError) as error:
            raise PrepareError(f"cannot extract {Path(archive).name}: {error}") from error
        if not ((staging / relative).is_dir() or target.is_dir()):
            raise PrepareError(
                f"extracting {Path(archive).name} did not produce {relative.as_posix()}/"
            )
        try:
            _merge(staging, workdir)
        except OSError as error:
            raise PrepareError(f"cannot extract {Path(archive).name}: {error}") from error
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return target


def _glob(base: Path, pattern: str, dst: Path) -> list[Path]:
    try:
        found = [path for path in base.glob(pattern) if path.absolute() != dst.absolute()]
    # pathlib refuses empty and absolute patterns.
    except (ValueError, NotImplementedError) as error:
        raise PrepareError(f"bad concatenate pattern {pattern!r}: {error}") from error
    # By relative path string, as phu-db-prep sorts; Path ordering compares
    # segment by segment and puts "a/x" before "a-b/x".
    return sorted(found, key=lambda path: path.relative_to(base).as_posix())


def concatenate(pattern: str, workdir: Path, downloads: Path, dst: Path) -> Path:
    """Join every file matching the glob ``pattern`` into ``dst``.

    Matches come from the workdir (earlier step targets) and only if there are
    none there, from the downloads. A ``.tsv`` target keeps the first file's
    header and drops every later file's first line, inserting a newline where a
    file did not end with one; anything else is a byte-exact append.
    """
    dst = Path(dst)
    if ".." in pattern.split("/"):
        raise PrepareError(f"bad concatenate pattern {pattern!r}")
    matches = _glob(Path(workdir), pattern, dst) or _glob(Path(downloads), pattern, dst)
    if not matches:
        raise PrepareError(f"concatenate matched no files: {pattern}")
    skip_header = dst.name.endswith(".tsv")
    try:
        with _atomic_output(dst) as output:
            ends_with_newline = True
            for index, path in enumerate(matches):
                with path.open("rb") as stream:
                    if skip_header and index:
                        stream.readline()
                    if skip_header and not ends_with_newline:
                        output.write(b"\n")
                        ends_with_newline = True
                    for chunk in iter(lambda: stream.read(COPY_CHUNK), b""):
                        output.write(chunk)
                        ends_with_newline = chunk.endswith(b"\n")
    except OSError as error:
        raise PrepareError(f"cannot concatenate {pattern} into {dst.name}: {error}") from error
    return dst


def rename(src: Path, dst: Path) -> Path:
    """Copy ``src`` to ``dst`` with its bytes unchanged."""
    try:
        with Path(src).open("rb") as stream, _atomic_output(dst) as sink:
            shutil.copyfileobj(stream, sink, COPY_CHUNK)
    except OSError as error:
        raise PrepareError(f"cannot copy {Path(src).name} to {Path(dst).name}: {error}") from error
    return Path(dst)


_SHEET = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_ID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
_PUL_ID_RE = re.compile(r"PUL\d+")
_PUL_COLUMNS = ("id", "pul", "pul_id", "pul id")
_RULE_COLUMN = "cazymes_predicted_dbcan"
PUL_TABLE_HEADER = "pul_id\tsubstrate\tcazymes_predicted_dbcan"
# What a damaged workbook raises: a missing part or relationship is a KeyError,
# a bad shared-string index an IndexError or ValueError; zipfile raises
# RuntimeError for an encrypted member, NotImplementedError for an unknown
# compression method.
_XLSX_ERRORS = (
    OSError,
    EOFError,
    zlib.error,
    zipfile.BadZipFile,
    ElementTree.ParseError,
    KeyError,
    IndexError,
    ValueError,
    RuntimeError,
    NotImplementedError,
)


def _xlsx_rows(path: Path) -> list[list[str]]:
    """Cell text of the first worksheet, decoded like phu main's dbcan_db._xlsx_values.

    Shared strings, inlineStr and plain <v> text; the column comes from the
    letters of the cell ref; rows without cells are skipped and missing cells
    are "". As in phu-db-prep, an absolute rels target ("/xl/...") resolves
    where phu main's lookup would fail; every workbook main reads decodes the same.
    """
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in root.findall(f"{_SHEET}si")]
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {item.attrib["Id"]: item.attrib["Target"] for item in rels}
        sheet = workbook.find(f".//{_SHEET}sheet")
        if sheet is None:
            raise PrepareError(f"{Path(path).name} has no worksheets")
        target = targets[sheet.attrib[_REL_ID]].removeprefix("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ElementTree.fromstring(archive.read(target))

    rows: list[list[str]] = []
    for row in root.findall(f".//{_SHEET}row"):
        values: dict[int, str] = {}
        for cell in row.findall(f"{_SHEET}c"):
            letters = re.search(r"[A-Z]+", cell.attrib.get("r", ""))
            if not letters:
                continue
            column = 0
            for character in letters.group():
                column = column * 26 + ord(character) - 64
            value = cell.find(f"{_SHEET}v")
            text = "" if value is None else value.text or ""
            if cell.attrib.get("t") == "s" and text:
                text = shared[int(text)]
            elif cell.attrib.get("t") == "inlineStr":
                inline = cell.find(f"{_SHEET}is")
                text = "" if inline is None else "".join(inline.itertext())
            values[column - 1] = text
        if values:
            rows.append([values.get(index, "") for index in range(max(values) + 1)])
    return rows


def _clean(value: str) -> str:
    # One record per line: a tab or line break inside a cell would split it.
    return re.sub(r"[\t\r\n]", " ", value.strip())


def pul_table(xlsx: Path, dst: Path) -> Path:
    """Write the PUL rows of the dbCAN workbook as pul_id/substrate/rule TSV.

    No family resolution happens here; phu tokenizes the rule column against
    the families in profiles.hmm when it reads the table. Duplicate PUL ids are
    all kept, in workbook order.
    """
    name = Path(xlsx).name
    try:
        rows = _xlsx_rows(Path(xlsx))
    # PrepareError is a RuntimeError; let _xlsx_rows' own message through.
    except PrepareError:
        raise
    except _XLSX_ERRORS as error:
        raise PrepareError(f"{name} is not a readable xlsx workbook: {error!r}") from error
    if not rows:
        raise PrepareError(f"{name} is empty")
    # A name -> index dict, as phu main's parse_pul_workbook and phu-db-prep
    # build it, so a repeated header resolves the same way in all three: its
    # last column wins. Distinct names still go in first-column order.
    headers = {value.strip().lower(): index for index, value in enumerate(rows[0])}
    pul = next((headers[column] for column in _PUL_COLUMNS if column in headers), None)
    if pul is None:
        raise PrepareError(f"{name} is missing columns: pul")
    rule = headers.get(_RULE_COLUMN)
    if rule is None:
        raise PrepareError(f"{name} is missing columns: {_RULE_COLUMN}")
    substrate = next((index for header, index in headers.items() if "substrate" in header), None)

    def cell(row: list[str], index: int | None) -> str:
        return row[index] if index is not None and index < len(row) else ""

    lines = [PUL_TABLE_HEADER]
    for row in rows[1:]:
        pul_id = cell(row, pul).strip().upper()
        if _PUL_ID_RE.fullmatch(pul_id):
            lines.append(f"{pul_id}\t{_clean(cell(row, substrate))}\t{_clean(cell(row, rule))}")
    try:
        with _atomic_output(dst) as output:
            output.write(("\n".join(lines) + "\n").encode("utf-8"))
    except OSError as error:
        raise PrepareError(f"cannot write {Path(dst).name}: {error}") from error
    return Path(dst)


_FILE_VERBS = {"gunzip": gunzip, "rename": rename, "pul-table": pul_table}


def _resolve(name: str, workdir: Path, downloads: Path) -> Path:
    """An earlier step's output in the workdir wins over a download of that name."""
    for base in (workdir, downloads):
        if (base / name).exists():
            return base / name
    raise PrepareError(f"prepare input {name} is neither an earlier target nor a download")


def _step_text(step: Step) -> str:
    return f"{step.verb}:{step.source} -> {step.target}{'/' if step.is_directory else ''}"


def run_steps(
    steps: Iterable[str | Step], *, workdir: Path, downloads: Path
) -> list[Path]:
    """Run prepare steps in order inside ``workdir``; returns each step's target path.

    ``downloads`` holds the database's verified sources under their manifest
    names. Steps are strings in the contract grammar (or Steps).
    """
    workdir = Path(workdir)
    downloads = Path(downloads)
    workdir.mkdir(parents=True, exist_ok=True)
    produced = []
    for text in steps:
        # A Step can be built by hand; it goes through the same verb and path
        # checks as a manifest string.
        try:
            step = parse_step(_step_text(text) if isinstance(text, Step) else text)
        except ManifestError as error:
            raise PrepareError(str(error)) from error
        target = workdir / step.target
        if step.verb == "concatenate":
            concatenate(step.source, workdir, downloads, target)
        elif step.verb == "extract":
            extract(_resolve(step.source, workdir, downloads), workdir, Path(step.target))
        else:
            _FILE_VERBS[step.verb](_resolve(step.source, workdir, downloads), target)
        produced.append(target)
    return produced
