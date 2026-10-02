"""pytest configuration and fixtures."""

import importlib.util
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


@pytest.fixture
def loopback_server():
    """A fresh HTTP server on 127.0.0.1; see tests/fixtures/phu_db/loopback.py."""
    server = _loopback.LoopbackServer().start()
    try:
        yield server
    finally:
        server.stop()
