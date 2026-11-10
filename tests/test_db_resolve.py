import hashlib
import importlib.util
import json
import socket
import sys
import types
from pathlib import Path

import pytest

from phu import db_bundle
from phu.db_bundle import (
    bundle_provenance,
    default_db_folder,
    resolve_db,
    tree_problems,
)
from phu.db_fetch import DatabaseError, download_verified
from phu.db_manifest import DATABASE_NAMES

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"
# Every URL in the fixture manifest and pins lives under this host.
FIXTURE_HOST = "https://phu-db-fixture.invalid/"


@pytest.fixture(autouse=True)
def loopback_only(monkeypatch):
    """Fail any test that tries to open a connection off this machine."""
    real = socket.create_connection

    def guarded(address, *args, **kwargs):
        if address[0] not in ("127.0.0.1", "::1", "localhost"):
            raise AssertionError(f"test tried to reach {address[0]}")
        return real(address, *args, **kwargs)

    monkeypatch.setattr(socket, "create_connection", guarded)


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for name in ("PHU_DB_BUNDLE", "PHU_DB_VERIFY", "PHU_DB_FOLDER", "XDG_DATA_HOME"):
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def db_folder(tmp_path):
    return tmp_path / "db"


class LoopbackFetch:
    """download_verified, with phu-db-fixture.invalid URLs served by the loopback server."""

    def __init__(self, server) -> None:
        self.server = server
        self.urls = []

    def __call__(self, url, dest, sha256, size):
        assert url.startswith(FIXTURE_HOST), url
        self.urls.append(url)
        return download_verified(
            self.server.url(url[len(FIXTURE_HOST) :]),
            dest,
            sha256,
            size,
            timeout=5,
            retries=0,
            sleep=lambda seconds: None,
        )


@pytest.fixture
def server(loopback_server, fixture_bundle):
    """Serves the full fixture bundle."""
    loopback_server.files["phu-db.tar.gz"] = fixture_bundle.read_bytes()
    return loopback_server


@pytest.fixture
def fetch(server):
    return LoopbackFetch(server)


@pytest.fixture
def resolve(db_folder, fetch, fixture_pins, fixture_manifest_bytes):
    def run(name, **kwargs):
        kwargs.setdefault("pins", fixture_pins)
        kwargs.setdefault("db_folder", db_folder)
        kwargs.setdefault("vendored_manifest_bytes", fixture_manifest_bytes)
        kwargs.setdefault("fetch", fetch)
        kwargs.setdefault("log", lambda message: None)
        return resolve_db(name, **kwargs)

    return run


def _target(db_folder, pins):
    return db_folder / "phu-db" / pins.bundle_version


def _leftovers(db_folder):
    return sorted(db_folder.glob(".phu-db-tmp-*")) if db_folder.exists() else []


def _tree(root):
    """Relative path -> sha256 (None for directories) of everything under ``root``."""
    return {
        path.relative_to(root).as_posix(): (
            hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        )
        for path in sorted(root.rglob("*"))
    }


def _flip_byte(path: Path, *, avoid=b"\n \t\r") -> None:
    data = bytearray(path.read_bytes())
    index = next(i for i in range(len(data) // 2, len(data)) if data[i] not in avoid)
    data[index] = ord("0") if data[index] != ord("0") else ord("1")
    path.write_bytes(bytes(data))


def _installed(db_folder, pins):
    return json.loads((_target(db_folder, pins) / "installed.json").read_text())


# resolve_db.


@pytest.mark.parametrize("name", sorted(DATABASE_NAMES))
def test_resolve_db_returns_the_database_directory(resolve, db_folder, fixture_pins, name):
    path = resolve(name)

    assert path == db_folder / "phu-db" / fixture_pins.bundle_version / name
    assert path.is_dir()


def test_resolve_db_installs_only_once(resolve, db_folder, fetch, fixture_pins):
    logs = []
    first = resolve("pfam", log=logs.append)
    assert fetch.urls == [fixture_pins.url]
    record = (_target(db_folder, fixture_pins) / "installed.json").read_bytes()
    fetch.urls.clear()
    logs.clear()

    assert resolve("pfam", log=logs.append) == first
    assert resolve("kofam", log=logs.append) == first.parent / "kofam"

    assert fetch.urls == []
    assert logs == []
    assert (_target(db_folder, fixture_pins) / "installed.json").read_bytes() == record


def test_phu_db_verify_catches_what_shallow_validation_misses(
    resolve, monkeypatch, db_folder, fetch, fixture_pins, fixture_manifest
):
    path = resolve("pfam")
    hmm = path / "profiles.hmm"
    original = hmm.read_bytes()
    _flip_byte(hmm)
    assert len(hmm.read_bytes()) == len(original)
    target = _target(db_folder, fixture_pins)
    assert tree_problems(target, fixture_manifest, deep=False) == []
    fetch.urls.clear()

    # Shallow: same size, so the altered tree is accepted.
    resolve("pfam")
    assert fetch.urls == []
    assert hmm.read_bytes() != original

    monkeypatch.setenv("PHU_DB_VERIFY", "1")
    assert resolve("pfam") == path

    assert fetch.urls == [fixture_pins.url]
    assert hmm.read_bytes() == original
    assert tree_problems(target, fixture_manifest, deep=True) == []


@pytest.mark.parametrize("name", ["cazy", "", "PFAM", "../pfam"])
def test_unknown_database_name_raises(resolve, db_folder, fetch, name):
    with pytest.raises(ValueError, match="unknown phu-db database"):
        resolve(name)

    assert fetch.urls == []
    assert not db_folder.exists()


def test_pins_without_a_manifest_are_refused(db_folder, fixture_pins, fixture_manifest_bytes):
    with pytest.raises(ValueError, match="together"):
        resolve_db("pfam", pins=fixture_pins, db_folder=db_folder)
    with pytest.raises(ValueError, match="together"):
        resolve_db("pfam", db_folder=db_folder, vendored_manifest_bytes=fixture_manifest_bytes)


def test_unpinned_phu_raises_a_database_error(monkeypatch, db_folder):
    # No phu._bundle (until 0.1.1 is published): None in sys.modules makes it unimportable.
    monkeypatch.setitem(sys.modules, "phu._bundle", None)

    with pytest.raises(DatabaseError, match="not pinned yet") as raised:
        resolve_db("pfam", db_folder=db_folder)
    assert not isinstance(raised.value, ImportError)
    assert raised.value.__cause__ is None
    with pytest.raises(DatabaseError, match="not pinned yet"):
        bundle_provenance(db_folder=db_folder)
    assert not db_folder.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("phu._bundle") is not None,
    reason="phu._bundle exists: phu is pinned",
)
def test_this_phu_is_not_pinned_yet(db_folder):
    # Nothing patched: the real state until 0.1.1 is published.
    with pytest.raises(DatabaseError, match="not pinned yet") as raised:
        resolve_db("pfam", db_folder=db_folder)
    assert not isinstance(raised.value, ImportError)
    assert not db_folder.exists()


def _vendored_data(monkeypatch, root: Path) -> None:
    """Make importlib.resources.files("phu") point at ``root``."""
    monkeypatch.setattr(db_bundle.importlib.resources, "files", lambda package: root)


def test_pins_without_a_vendored_manifest_raise_a_database_error(
    monkeypatch, tmp_path, db_folder, fixture_pins
):
    monkeypatch.setitem(
        sys.modules, "phu._bundle", types.SimpleNamespace(PRODUCTION_PINS=fixture_pins)
    )
    _vendored_data(monkeypatch, tmp_path / "package")

    with pytest.raises(DatabaseError, match="not pinned yet"):
        resolve_db("pfam", db_folder=db_folder)
    assert not db_folder.exists()


def test_production_pins_and_vendored_manifest_are_used_by_default(
    monkeypatch, tmp_path, db_folder, fixture_bundle, fixture_pins, fixture_manifest_bytes
):
    monkeypatch.setitem(
        sys.modules, "phu._bundle", types.SimpleNamespace(PRODUCTION_PINS=fixture_pins)
    )
    package = tmp_path / "package"
    (package / "data").mkdir(parents=True)
    (package / "data" / "phu-db.json").write_bytes(fixture_manifest_bytes)
    _vendored_data(monkeypatch, package)
    monkeypatch.setenv("PHU_DB_BUNDLE", str(fixture_bundle))

    path = resolve_db("dbcan", db_folder=db_folder, log=lambda message: None)

    assert path == _target(db_folder, fixture_pins) / "dbcan"
    assert bundle_provenance(db_folder=db_folder)["source"] == "local"


def test_default_db_folder(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    assert default_db_folder() == tmp_path / "home" / ".local" / "share" / "phu" / "db"

    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    assert default_db_folder() == tmp_path / "xdg" / "phu" / "db"

    monkeypatch.setenv("PHU_DB_FOLDER", str(tmp_path / "mine"))
    assert default_db_folder() == tmp_path / "mine"


@pytest.mark.parametrize("variable", ["PHU_DB_FOLDER", "XDG_DATA_HOME"])
def test_resolve_db_installs_into_the_default_db_folder(
    resolve, monkeypatch, tmp_path, fixture_pins, variable
):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv(variable, str(tmp_path / "env"))
    expected = tmp_path / "env" if variable == "PHU_DB_FOLDER" else tmp_path / "env" / "phu" / "db"

    path = resolve("avg", db_folder=None)

    assert path == _target(expected, fixture_pins) / "avg"
    assert path.is_dir()
    assert not (tmp_path / "home").exists()


def test_resolve_db_checks_only_the_requested_database(
    resolve, monkeypatch, db_folder, fetch, fixture_pins
):
    target = resolve("avg").parent
    _flip_byte(target / "pfam" / "profiles.hmm")
    fetch.urls.clear()
    monkeypatch.setenv("PHU_DB_VERIFY", "1")

    resolve("avg")
    assert fetch.urls == []

    resolve("pfam")
    assert fetch.urls == [fixture_pins.url]


def test_vendored_manifest_must_match_its_pin(
    resolve, db_folder, fetch, fixture_pins, fixture_manifest_bytes
):
    resolve("pfam")
    fetch.urls.clear()
    drifted = fixture_manifest_bytes.replace(b"\n}", b"\n}\n", 1)
    assert drifted != fixture_manifest_bytes

    # A valid tree is installed, but the manifest is not the pinned one.
    with pytest.raises(db_bundle.BundleMismatch, match="pinned sha256"):
        resolve("pfam", vendored_manifest_bytes=drifted)
    with pytest.raises(db_bundle.BundleMismatch, match="pinned sha256"):
        bundle_provenance(pins=fixture_pins, db_folder=db_folder, vendored_manifest_bytes=drifted)
    assert fetch.urls == []


# bundle_provenance.


def test_bundle_provenance_reports_the_install(
    resolve, db_folder, fetch, fixture_pins, fixture_manifest, fixture_manifest_bytes
):
    pins = fixture_pins
    resolve("pfam")
    fetch.urls.clear()

    provenance = bundle_provenance(
        pins=pins, db_folder=db_folder, vendored_manifest_bytes=fixture_manifest_bytes
    )

    assert provenance == {
        "source": "zenodo",
        "bundle_version": pins.bundle_version,
        "sha256": pins.sha256,
        "upstream_versions": {
            name: entry.upstream_version for name, entry in fixture_manifest.databases.items()
        },
    }
    assert set(provenance["upstream_versions"]) == DATABASE_NAMES
    assert all(provenance["upstream_versions"].values())
    assert fetch.urls == []


def test_bundle_provenance_never_installs(
    db_folder, fixture_pins, fixture_manifest_bytes
):
    with pytest.raises(DatabaseError, match="not installed"):
        bundle_provenance(
            pins=fixture_pins, db_folder=db_folder, vendored_manifest_bytes=fixture_manifest_bytes
        )

    assert not db_folder.exists()
