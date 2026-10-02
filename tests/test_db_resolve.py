import hashlib
import importlib.util
import json
import socket
import sys
import types
from dataclasses import replace
from pathlib import Path

import pytest

from phu import db_bundle
from phu.db_bundle import (
    UpstreamMismatch,
    bundle_provenance,
    default_db_folder,
    resolve_db,
    tree_problems,
)
from phu.db_fetch import DatabaseError, IntegrityError, download_verified
from phu.db_manifest import DATABASE_NAMES

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"
UPSTREAM = FIXTURE_DIR / "upstream"
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

    @property
    def upstream_urls(self):
        return [url for url in self.urls if not url.endswith("/phu-db.tar.gz")]


@pytest.fixture
def server(loopback_server, fixture_bundle):
    """Serves the full fixture bundle and every fixture upstream source."""
    loopback_server.serve_dir(UPSTREAM)
    loopback_server.files["phu-db.tar.gz"] = fixture_bundle.read_bytes()
    return loopback_server


@pytest.fixture
def fetch(server):
    return LoopbackFetch(server)


@pytest.fixture
def partial_bundle(make_bundle, server, fixture_manifest, fixture_pins):
    """Serve a bundle without the named databases; returns its pins."""

    def build(*missing, files=None):
        drop = set()
        for name in missing:
            drop |= fixture_manifest.databases[name].member_paths()
        tarball = make_bundle(drop=drop, replace=files)
        data = tarball.read_bytes()
        server.files["phu-db.tar.gz"] = data
        return replace(fixture_pins, sha256=hashlib.sha256(data).hexdigest(), size=len(data))

    return build


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


# Upstream fallback.


@pytest.fixture
def full_install(
    tmp_path_factory, fixture_bundle, fixture_pins, fixture_manifest, fixture_manifest_bytes
):
    """The tree a full-bundle install produces, as _tree() of it."""
    folder = tmp_path_factory.mktemp("full") / "db"

    def no_fetch(url, dest, sha256, size):
        raise AssertionError("a full local bundle needs no download")

    db_bundle.install(
        db_folder=folder,
        pins=fixture_pins,
        manifest=fixture_manifest,
        vendored_manifest_bytes=fixture_manifest_bytes,
        bundle_path=fixture_bundle,
        fetch=no_fetch,
        log=lambda message: None,
    )
    return _tree(_target(folder, fixture_pins))


@pytest.mark.parametrize(
    "missing",
    [("dbcan",), ("kofam", "dbcan"), ("pfam", "kofam", "dbcan", "avg")],
    ids=["dbcan", "kofam+dbcan", "all"],
)
def test_databases_missing_from_the_bundle_come_from_upstream(
    resolve, partial_bundle, db_folder, fetch, fixture_manifest, full_install, missing
):
    pins = partial_bundle(*missing)

    path = resolve("dbcan", pins=pins)

    target = _target(db_folder, pins)
    assert path == target / "dbcan"
    tree = _tree(target)
    # installed.json differs by install time and source; everything else is identical.
    assert set(tree) == set(full_install)
    del tree["installed.json"]
    assert tree == {key: value for key, value in full_install.items() if key != "installed.json"}
    expected_urls = [
        source.url
        for name in sorted(missing)
        for source in fixture_manifest.databases[name].sources
    ]
    assert fetch.upstream_urls == expected_urls
    record = _installed(db_folder, pins)
    assert record["upstream_databases"] == sorted(missing)
    assert record["tarball_sha256"] == pins.sha256
    assert _leftovers(db_folder) == []
    assert list((db_folder / ".downloads").iterdir()) == []


def test_upstream_source_with_a_wrong_byte_installs_nothing(
    resolve, partial_bundle, server, db_folder
):
    pins = partial_bundle("dbcan")
    data = bytearray(server.files["dbcan/dbCAN-PUL.xlsx"])
    data[len(data) // 2] ^= 0xFF
    server.files["dbcan/dbCAN-PUL.xlsx"] = bytes(data)

    with pytest.raises(UpstreamMismatch) as raised:
        resolve("dbcan", pins=pins)

    assert isinstance(raised.value, IntegrityError)
    assert "dbCAN-PUL.xlsx" in str(raised.value)
    assert not (db_folder / "phu-db").exists()
    assert _leftovers(db_folder) == []
    # Only the verified bundle download is kept, for the next attempt to reuse.
    assert sorted(path.name for path in (db_folder / ".downloads").iterdir()) == ["phu-db.tar.gz"]


def test_upstream_source_is_verified_after_an_injected_fetch(
    resolve, partial_bundle, server, db_folder, fetch
):
    pins = partial_bundle("avg")

    def careless(url, dest, sha256, size):
        if url.endswith("/phu-db.tar.gz"):
            return fetch(url, dest, sha256, size)
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        Path(dest).write_bytes(b"x" * size)
        return Path(dest)

    with pytest.raises(UpstreamMismatch, match="has sha256"):
        resolve("avg", pins=pins, fetch=careless)

    assert not (db_folder / "phu-db").exists()
    assert _leftovers(db_folder) == []


def test_prepared_upstream_file_must_match_the_manifest(
    resolve, partial_bundle, db_folder, fixture_manifest_bytes, fixture_pins
):
    # Sources match their pins, but a prepared file's digest does not.
    raw = json.loads(fixture_manifest_bytes)
    pul_rules = next(
        item
        for item in raw["databases"]["dbcan"]["prepared_files"]
        if item["path"] == "pul_rules.tsv"
    )
    pul_rules["sha256"] = "0" * 64
    manifest_bytes = (json.dumps(raw, indent=2) + "\n").encode()
    pins = partial_bundle("dbcan", files={"phu-db/phu-db.json": manifest_bytes})
    pins = replace(pins, manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest())

    with pytest.raises(UpstreamMismatch) as raised:
        resolve("dbcan", pins=pins, vendored_manifest_bytes=manifest_bytes)

    assert "pul_rules.tsv has sha256" in str(raised.value)
    assert not (db_folder / "phu-db").exists()
    assert _leftovers(db_folder) == []


def _bytes_under(root: Path) -> int:
    return sum(path.stat().st_size for path in root.rglob("*") if path.is_file())


@pytest.mark.parametrize(
    "missing", [("kofam",), ("pfam", "kofam", "dbcan", "avg")], ids=["kofam", "all"]
)
def test_upstream_space_check_covers_the_peak(
    resolve, partial_bundle, monkeypatch, db_folder, missing
):
    # kofam's extract leaves every per-KO HMM beside the concatenated profiles.hmm.
    pins = partial_bundle(*missing)
    checks = []
    peaks = []
    real_ensure = db_bundle.ensure_free_space
    real_run_steps = db_bundle.run_steps

    def ensure(path, required):
        checks.append((required, _bytes_under(db_folder)))
        return real_ensure(path, required)

    def run_steps(steps, *, workdir, downloads):
        produced = real_run_steps(steps, workdir=workdir, downloads=downloads)
        # Nothing a database leaves on disk is removed before its steps finish.
        peaks.append(_bytes_under(db_folder))
        return produced

    monkeypatch.setattr(db_bundle, "ensure_free_space", ensure)
    monkeypatch.setattr(db_bundle, "run_steps", run_steps)

    resolve("kofam", pins=pins)

    # The install preflight, then the upstream check before any source is fetched.
    assert len(checks) == 2
    required, used_before = checks[1]
    assert len(peaks) == len(missing)
    assert max(peaks) - used_before <= required


def test_install_log_names_upstream_databases_once(
    capsys, partial_bundle, db_folder, fetch, fixture_manifest_bytes, fixture_manifest
):
    pins = partial_bundle("kofam", "dbcan")
    kwargs = dict(
        pins=pins,
        db_folder=db_folder,
        vendored_manifest_bytes=fixture_manifest_bytes,
        fetch=fetch,
    )

    resolve_db("pfam", **kwargs)

    captured = capsys.readouterr()
    assert captured.out == ""
    lines = [line for line in captured.err.splitlines() if "upstream" in line]
    assert len(lines) == 1
    assert "dbcan, kofam" in lines[0]
    # Names where the databases end up, not only the staging directory.
    assert f"into {_target(db_folder, pins)} " in lines[0]
    sources = [
        source for name in ("dbcan", "kofam") for source in fixture_manifest.databases[name].sources
    ]
    assert f"{len(sources)} upstream files" in lines[0]
    assert db_bundle._human(sum(source.bytes for source in sources)) in lines[0]
    for source in sources:
        assert source.url in captured.err

    # Already installed: nothing is downloaded, and nothing logged.
    resolve_db("dbcan", **kwargs)
    assert capsys.readouterr().err == ""


def test_full_bundle_install_logs_no_upstream(
    capsys, db_folder, fetch, fixture_pins, fixture_manifest_bytes
):
    resolve_db(
        "avg",
        pins=fixture_pins,
        db_folder=db_folder,
        vendored_manifest_bytes=fixture_manifest_bytes,
        fetch=fetch,
    )

    assert "upstream" not in capsys.readouterr().err
    assert fetch.upstream_urls == []
    assert _installed(db_folder, fixture_pins)["upstream_databases"] == []


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
    resolve, partial_bundle, db_folder, fetch, fixture_manifest, fixture_manifest_bytes
):
    pins = partial_bundle("dbcan")
    resolve("pfam", pins=pins)
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
        "upstream_databases": ["dbcan"],
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
