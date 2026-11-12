"""pytest configuration and fixtures."""

import hashlib
import importlib.util
import itertools
import sys
from pathlib import Path

import pytest

# Ensure src/ is in sys.path for imports
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

PHU_DB_FIXTURES = Path(__file__).parent / "fixtures" / "phu_db"


def _load_fixture_module(name: str):
    """Import a helper from tests/fixtures/phu_db without making it a package."""
    spec = importlib.util.spec_from_file_location(
        f"phu_db_fixture_{name}", PHU_DB_FIXTURES / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_loopback = _load_fixture_module("loopback")
_bundle = _load_fixture_module("bundle")


@pytest.fixture
def loopback_server():
    """A fresh HTTP server on 127.0.0.1; see tests/fixtures/phu_db/loopback.py."""
    server = _loopback.LoopbackServer().start()
    try:
        yield server
    finally:
        server.stop()


@pytest.fixture(scope="session")
def fixture_manifest_bytes() -> bytes:
    """The fixture phu-db.json, standing in for the manifest vendored with phu."""
    return (PHU_DB_FIXTURES / "phu-db.json").read_bytes()


@pytest.fixture(scope="session")
def fixture_manifest(fixture_manifest_bytes):
    from phu.db_manifest import loads_manifest

    return loads_manifest(fixture_manifest_bytes)


@pytest.fixture(scope="session")
def fixture_bundle(tmp_path_factory) -> Path:
    """The full fixture bundle, phu-db.tar.gz, built once per session; never committed."""
    path = tmp_path_factory.mktemp("phu-db-bundle") / "phu-db.tar.gz"
    return _bundle.write_bundle(path, _bundle.bundle_files())


@pytest.fixture(scope="session")
def fixture_pins(fixture_bundle, fixture_manifest, fixture_manifest_bytes):
    """BundlePins for the fixture bundle; the only place tests get pins from."""
    from phu.db_bundle import BundlePins

    data = fixture_bundle.read_bytes()
    return BundlePins(
        url=_bundle.URL,
        sha256=hashlib.sha256(data).hexdigest(),
        size=len(data),
        bundle_version=fixture_manifest.bundle_version,
        manifest_sha256=hashlib.sha256(fixture_manifest_bytes).hexdigest(),
        expected_members=frozenset(_bundle.bundle_files()),
    )


@pytest.fixture
def make_bundle(tmp_path):
    """Build a variant of the fixture bundle for negative tests.

    ``make_bundle(drop=names, replace={name: bytes}, extra=[(TarInfo, bytes | None)])``
    returns the path of a new tarball; ``make_bundle.member`` builds a TarInfo
    with the bundle's fixed metadata.
    """
    counter = itertools.count()

    def build(*, drop=(), replace=None, extra=()) -> Path:
        files = _bundle.bundle_files()
        for name in drop:
            del files[name]
        files.update(replace or {})
        path = tmp_path / f"variant-{next(counter)}.tar.gz"
        return _bundle.write_bundle(path, files, extra)

    build.member = _bundle.member
    build.files = _bundle.bundle_files
    return build
