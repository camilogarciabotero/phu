"""Binary copies of the bundle's HMM files."""

import shutil
from pathlib import Path

import pyhmmer
import pytest

from phu.db import binary

FIXTURE_PFAM = Path(__file__).parent / "fixtures" / "phu_db" / "prepared" / "pfam" / "profiles.hmm"


@pytest.fixture
def text_hmm(tmp_path, monkeypatch) -> Path:
    monkeypatch.setenv("PHU_CACHE_DIR", str(tmp_path / "cache"))
    database = tmp_path / "pfam"
    database.mkdir()
    return Path(shutil.copy(FIXTURE_PFAM, database / "profiles.hmm"))


def _models(path):
    with pyhmmer.plan7.HMMFile(path) as handle:
        return [
            (hmm.name, hmm.accession, hmm.description, hmm.M, hmm.cutoffs.gathering, hmm.checksum)
            for hmm in handle
        ]


def test_binary_copy_holds_the_same_models(text_hmm):
    path = binary.binary_hmm_path(text_hmm)

    assert path == text_hmm.parent / "profiles.h3m"
    assert _models(path) == _models(text_hmm)
    assert len(_models(path)) == 2


def test_binary_copy_is_built_once_and_rebuilt_when_the_source_changes(text_hmm, monkeypatch):
    builds = []
    build = binary._build
    monkeypatch.setattr(binary, "_build", lambda *args: builds.append(args) or build(*args))

    binary.binary_hmm_path(text_hmm)
    binary.binary_hmm_path(text_hmm)
    assert len(builds) == 1

    with pyhmmer.plan7.HMMFile(text_hmm) as handle:
        first = next(handle)
    with text_hmm.open("wb") as output:
        first.write(output)
    path = binary.binary_hmm_path(text_hmm)

    assert len(builds) == 2
    assert len(_models(path)) == 1


def test_a_stale_stamp_from_another_pyhmmer_forces_a_rebuild(text_hmm, monkeypatch):
    path = binary.binary_hmm_path(text_hmm)
    monkeypatch.setattr(binary.pyhmmer, "__version__", "0.0.0")
    path.write_bytes(b"not a model file")

    assert _models(binary.binary_hmm_path(text_hmm)) == _models(text_hmm)


def test_text_file_is_used_when_the_copy_cannot_be_written(text_hmm, monkeypatch):
    def fail(*args):
        raise OSError("No space left on device")

    monkeypatch.setattr(binary, "_build", fail)

    assert binary.binary_hmm_path(text_hmm) == text_hmm


def test_read_only_install_keeps_the_copy_in_the_user_cache(text_hmm, tmp_path):
    text_hmm.parent.chmod(0o555)
    try:
        path = binary.binary_hmm_path(text_hmm)
    finally:
        text_hmm.parent.chmod(0o755)

    assert tmp_path / "cache" in path.parents
    assert _models(path) == _models(text_hmm)
