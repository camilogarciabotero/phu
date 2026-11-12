"""Model of phu-db.json, the manifest shipped inside the phu-db bundle.

The manifest is trusted only after load_manifest has checked it: every path it
names is a safe relative POSIX path, every digest is well formed, and every
prepared file is produced by a step that the phu-db 0.1.1 contract allows.
"""

from __future__ import annotations

import fnmatch
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any

LAYOUT_VERSION = 1
DATABASE_NAMES = frozenset({"pfam", "kofam", "dbcan", "avg"})
BUNDLE_ROOT = "phu-db"
MANIFEST_NAME = "phu-db.json"
MANIFEST_MEMBER = f"{BUNDLE_ROOT}/{MANIFEST_NAME}"
PREPARED_FORMATS = frozenset({"hmm", "tsv"})

STEP_RE = re.compile(r"(gunzip|extract|concatenate|rename|pul-table):(\S+) -> (\S+)")
_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_GLOB_CHARS = frozenset("*?[")


class ManifestError(ValueError):
    """The manifest is malformed or breaks the bundle contract."""


@dataclass(frozen=True)
class Source:
    name: str
    url: str
    sha256: str
    bytes: int


@dataclass(frozen=True)
class PreparedFile:
    path: str
    format: str
    bytes: int
    sha256: str
    count: int


@dataclass(frozen=True)
class Step:
    """One parsed prepare step. ``target`` carries no trailing slash."""

    verb: str
    source: str
    target: str

    @property
    def is_directory(self) -> bool:
        return self.verb == "extract"

    @property
    def is_glob(self) -> bool:
        return self.verb == "concatenate" and bool(_GLOB_CHARS & set(self.source))


@dataclass(frozen=True)
class DatabaseEntry:
    name: str
    provenance: str
    redistribution: str
    upstream_version: str | None
    license: str | None
    citation: str | None
    layout_dir: str
    sources: tuple[Source, ...]
    steps: tuple[str, ...]
    prepared_files: tuple[PreparedFile, ...]

    @property
    def parsed_steps(self) -> tuple[Step, ...]:
        return tuple(parse_step(step) for step in self.steps)

    def member_paths(self) -> frozenset[str]:
        """Tarball member names of this database's prepared files."""
        return frozenset(
            f"{BUNDLE_ROOT}/{self.layout_dir}/{item.path}" for item in self.prepared_files
        )


@dataclass(frozen=True)
class Manifest:
    bundle_version: str
    layout_version: int
    # A mapping proxy is unhashable; hash on the version fields alone.
    databases: Mapping[str, DatabaseEntry] = field(hash=False)

    def member_path(self, database: str, path: str) -> str:
        return f"{BUNDLE_ROOT}/{self.databases[database].layout_dir}/{path}"

    def expected_members(self) -> frozenset[str]:
        """Every member of a full bundle, phu-db/phu-db.json included."""
        members = {MANIFEST_MEMBER}
        for entry in self.databases.values():
            members |= entry.member_paths()
        return frozenset(members)


def _check_path(value: str, what: str) -> None:
    """Reject anything but a plain relative POSIX path."""
    if not value or "\\" in value or "\0" in value or value.startswith("/"):
        raise ManifestError(f"unsafe {what}: {value!r}")
    if any(part in ("", ".", "..") for part in value.split("/")):
        raise ManifestError(f"unsafe {what}: {value!r}")


def parse_step(text: str) -> Step:
    """Parse ``verb:source -> target`` with the contract grammar and path rules."""
    if not isinstance(text, str):
        raise ManifestError(f"prepare step must be a string: {text!r}")
    match = STEP_RE.fullmatch(text)
    if not match:
        raise ManifestError(f"malformed prepare step: {text!r}")
    verb, source, target = match.groups()
    if target.endswith("/"):
        if verb != "extract":
            raise ManifestError(f"only extract may target a directory: {text!r}")
        target = target[:-1]
    _check_path(source, f"step source in {text!r}")
    _check_path(target, f"step target in {text!r}")
    return Step(verb, source, target)


def _require(mapping: Mapping[str, Any], key: str, kind: type | tuple, where: str) -> Any:
    if key not in mapping:
        raise ManifestError(f"{where} is missing {key!r}")
    value = mapping[key]
    # bool is an int subclass; a JSON true is never a valid size or count.
    if not isinstance(value, kind) or (isinstance(value, bool) and kind is not bool):
        raise ManifestError(f"{where}.{key} has the wrong type: {value!r}")
    return value


def _size(mapping: Mapping[str, Any], key: str, where: str) -> int:
    value = _require(mapping, key, int, where)
    if value < 0:
        raise ManifestError(f"{where}.{key} is negative: {value}")
    return value


def _sha256(mapping: Mapping[str, Any], where: str) -> str:
    value = _require(mapping, "sha256", str, where)
    if not _SHA256_RE.fullmatch(value):
        raise ManifestError(f"{where}.sha256 is not 64 lowercase hex digits: {value!r}")
    return value


def _optional_text(mapping: Mapping[str, Any], key: str, where: str) -> str | None:
    value = mapping.get(key)
    if value is not None and not isinstance(value, str):
        raise ManifestError(f"{where}.{key} has the wrong type: {value!r}")
    return value


def _check_steps(name: str, steps: tuple[Step, ...], source_names: set[str]) -> set[str]:
    """Check step references in order; returns the file targets the steps produce."""
    files: set[str] = set()
    directories: set[str] = set()
    for step in steps:
        known = step.source in source_names or step.source in files
        if not known and step.is_glob:
            # A glob may only range over earlier targets: files it could match
            # directly, or files an earlier extract unpacked under a directory.
            depth = len(PurePosixPath(step.source).parts)
            known = any(
                len(PurePosixPath(target).parts) == depth
                and fnmatch.fnmatchcase(target, step.source)
                for target in files
            ) or any(step.source.startswith(f"{directory}/") for directory in directories)
        if not known:
            raise ManifestError(
                f"databases.{name} step {step.verb}:{step.source} refers to neither a "
                f"source nor an earlier target"
            )
        _check_target_kind(name, step, files, directories)
        (directories if step.is_directory else files).add(step.target)
    return files


def _check_target_kind(name: str, step: Step, files: set[str], directories: set[str]) -> None:
    """Reject a target that would need one path to be both a file and a directory."""
    target = step.target
    clash = target in (files if step.is_directory else directories)
    # Nothing can sit below a file, and a new file cannot sit above an earlier target.
    clash = clash or any(target.startswith(f"{earlier}/") for earlier in files)
    if not step.is_directory:
        clash = clash or any(
            earlier.startswith(f"{target}/") for earlier in files | directories
        )
    if clash:
        raise ManifestError(
            f"databases.{name} step target {target!r} conflicts with an earlier "
            f"file or directory target"
        )


def _database(name: str, raw: Any) -> DatabaseEntry:
    where = f"databases.{name}"
    if not isinstance(raw, Mapping):
        raise ManifestError(f"{where} must be an object")
    layout_dir = _require(raw, "layout_dir", str, where)
    _check_path(layout_dir, f"{where}.layout_dir")
    if layout_dir != name:
        raise ManifestError(f"{where}.layout_dir must be {name!r}, got {layout_dir!r}")

    sources = []
    for index, item in enumerate(_require(raw, "sources", list, where)):
        here = f"{where}.sources[{index}]"
        if not isinstance(item, Mapping):
            raise ManifestError(f"{here} must be an object")
        source_name = _require(item, "name", str, here)
        # Names become download file names, so they must be a single segment.
        _check_path(source_name, f"{here}.name")
        if "/" in source_name:
            raise ManifestError(f"unsafe {here}.name: {source_name!r}")
        sources.append(
            Source(
                name=source_name,
                url=_require(item, "url", str, here),
                sha256=_sha256(item, here),
                bytes=_size(item, "bytes", here),
            )
        )
    source_names = {source.name for source in sources}
    if len(source_names) != len(sources):
        raise ManifestError(f"{where}.sources has duplicate names")

    step_texts = _require(raw, "steps", list, where)
    steps = tuple(parse_step(step) for step in step_texts)
    produced = _check_steps(name, steps, source_names)

    prepared_files = []
    for index, item in enumerate(_require(raw, "prepared_files", list, where)):
        here = f"{where}.prepared_files[{index}]"
        if not isinstance(item, Mapping):
            raise ManifestError(f"{here} must be an object")
        path = _require(item, "path", str, here)
        _check_path(path, f"{here}.path")
        if path not in produced:
            raise ManifestError(f"{here}.path {path!r} is not the target of any step")
        file_format = _require(item, "format", str, here)
        if file_format not in PREPARED_FORMATS:
            raise ManifestError(f"{here}.format must be one of hmm, tsv: {file_format!r}")
        prepared_files.append(
            PreparedFile(
                path=path,
                format=file_format,
                bytes=_size(item, "bytes", here),
                sha256=_sha256(item, here),
                count=_size(item, "count", here),
            )
        )
    if not prepared_files:
        raise ManifestError(f"{where}.prepared_files is empty")
    if len({item.path for item in prepared_files}) != len(prepared_files):
        raise ManifestError(f"{where}.prepared_files has duplicate paths")

    return DatabaseEntry(
        name=name,
        provenance=_require(raw, "provenance", str, where),
        redistribution=_require(raw, "redistribution", str, where),
        upstream_version=_optional_text(raw, "upstream_version", where),
        license=_optional_text(raw, "license", where),
        citation=_optional_text(raw, "citation", where),
        layout_dir=layout_dir,
        sources=tuple(sources),
        steps=tuple(step_texts),
        prepared_files=tuple(prepared_files),
    )


def _no_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ManifestError(f"duplicate key in manifest: {key!r}")
        result[key] = value
    return result


def loads_manifest(data: bytes | str) -> Manifest:
    """Parse and validate manifest JSON. Unknown keys are ignored."""
    try:
        raw = json.loads(data, object_pairs_hook=_no_duplicate_keys)
    # RecursionError: absurdly nested arrays or objects overflow the parser.
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise ManifestError(f"manifest is not valid JSON: {error}") from error
    if not isinstance(raw, Mapping):
        raise ManifestError("manifest must be a JSON object")

    layout_version = _require(raw, "layout_version", int, "manifest")
    if layout_version != LAYOUT_VERSION:
        raise ManifestError(
            f"unsupported layout_version {layout_version}; this phu reads {LAYOUT_VERSION}"
        )
    bundle_version = _require(raw, "bundle_version", str, "manifest")
    # The version names the install directory, so it must be one path segment.
    _check_path(bundle_version, "bundle_version")
    if "/" in bundle_version:
        raise ManifestError(f"unsafe bundle_version: {bundle_version!r}")

    databases = _require(raw, "databases", Mapping, "manifest")
    if set(databases) != DATABASE_NAMES:
        raise ManifestError(
            f"manifest databases must be {sorted(DATABASE_NAMES)}, got {sorted(databases)}"
        )
    return Manifest(
        bundle_version=bundle_version,
        layout_version=layout_version,
        databases=MappingProxyType(
            {name: _database(name, entry) for name, entry in databases.items()}
        ),
    )


def load_manifest(path: Path) -> Manifest:
    return loads_manifest(Path(path).read_bytes())
