import errno
import fcntl
import hashlib
import json
import multiprocessing
import os
import shutil
import socket
import time
from collections import namedtuple
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from phu import db_bundle
from phu.db_bundle import (
    LockTimeout,
    BundleMismatch,
    TarballMismatch,
    TreeInvalid,
    install,
    install_lock,
    installed_problems,
    tree_problems,
    validate_tree,
)
from phu.db_fetch import (
    InsufficientSpace,
    download_verified,
    required_install_bytes,
)
from phu.db_manifest import DATABASE_NAMES, loads_manifest

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"
PREPARED = FIXTURE_DIR / "prepared"

DiskUsage = namedtuple("DiskUsage", "total used free")


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
def no_bundle_env(monkeypatch):
    monkeypatch.delenv("PHU_DB_BUNDLE", raising=False)


@pytest.fixture
def db_folder(tmp_path):
    return tmp_path / "db"


@pytest.fixture
def logs():
    return []


class FetchSpy:
    """Stands in for download_verified: copies a local tarball to ``dest``."""

    def __init__(self, source: Path | None) -> None:
        self.source = source
        self.calls = []
        self.during = None

    def __call__(self, url, dest, sha256, size):
        self.calls.append((url, Path(dest), sha256, size))
        if self.during is not None:
            self.during(Path(dest))
        if self.source is None:
            raise AssertionError("install must not download")
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.source, dest)
        return Path(dest)


@pytest.fixture
def fetch(fixture_bundle):
    return FetchSpy(fixture_bundle)


@pytest.fixture
def run(db_folder, logs, fetch, fixture_pins, fixture_manifest, fixture_manifest_bytes):
    def run_install(**kwargs):
        kwargs.setdefault("pins", fixture_pins)
        kwargs.setdefault("fetch", fetch)
        kwargs.setdefault("log", logs.append)
        return install(
            db_folder=db_folder,
            manifest=fixture_manifest,
            vendored_manifest_bytes=fixture_manifest_bytes,
            **kwargs,
        )

    return run_install


@pytest.fixture
def user_bundle(tmp_path, fixture_bundle):
    """A copy of the fixture bundle in the user's own directory."""
    path = tmp_path / "mine" / "phu-db.tar.gz"
    path.parent.mkdir()
    shutil.copyfile(fixture_bundle, path)
    return path


def _target(db_folder, pins):
    return db_folder / "phu-db" / pins.bundle_version


def _leftovers(db_folder):
    return sorted(db_folder.glob(".phu-db-tmp-*")) if db_folder.exists() else []


def _assert_installed(target):
    assert {path.name for path in target.iterdir() if path.is_dir()} == DATABASE_NAMES
    for path in PREPARED.rglob("*"):
        if path.is_file():
            relative = path.relative_to(PREPARED)
            assert (target / relative).read_bytes() == path.read_bytes(), relative


def _flip_byte(path: Path, *, avoid=b"\n \t\r") -> None:
    data = bytearray(path.read_bytes())
    index = next(i for i in range(len(data) // 2, len(data)) if data[i] not in avoid)
    data[index] = ord("0") if data[index] != ord("0") else ord("1")
    path.write_bytes(bytes(data))


# Installing.


def test_install_from_phu_db_bundle_env(
    run, monkeypatch, db_folder, fetch, user_bundle, fixture_pins, fixture_manifest, logs
):
    monkeypatch.setenv("PHU_DB_BUNDLE", str(user_bundle))
    before = user_bundle.stat()
    fetch.source = None

    target = run()

    assert target == _target(db_folder, fixture_pins)
    _assert_installed(target)
    record = json.loads((target / "installed.json").read_text())
    assert set(record) == {
        "bundle_version",
        "source",
        "tarball_sha256",
        "installed_at",
    }
    assert record["bundle_version"] == fixture_pins.bundle_version
    assert record["source"] == "local"
    assert record["tarball_sha256"] == fixture_pins.sha256
    installed_at = datetime.fromisoformat(record["installed_at"])
    assert installed_at.utcoffset() == timedelta(0)
    assert abs(datetime.now(timezone.utc) - installed_at) < timedelta(minutes=5)
    # The user's tarball is read, never moved, modified or deleted.
    after = user_bundle.stat()
    assert (after.st_ino, after.st_size, after.st_mtime_ns) == (
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    assert hashlib.sha256(user_bundle.read_bytes()).hexdigest() == fixture_pins.sha256
    assert fetch.calls == []
    assert _leftovers(db_folder) == []
    assert not (db_folder / ".downloads").exists()
    assert installed_problems(db_folder, fixture_pins, fixture_manifest, deep=True) == []
    assert any(str(user_bundle) in message for message in logs)


def test_install_downloads_over_loopback(
    run, loopback_server, db_folder, fixture_bundle, fixture_pins, logs
):
    loopback_server.files["records/phu-db.tar.gz"] = fixture_bundle.read_bytes()

    def fetch(url, dest, sha256, size):
        assert url == fixture_pins.url
        return download_verified(
            loopback_server.url("records/phu-db.tar.gz"),
            dest,
            sha256,
            size,
            timeout=5,
            sleep=lambda seconds: None,
        )

    target = run(fetch=fetch)

    _assert_installed(target)
    assert json.loads((target / "installed.json").read_text())["source"] == "zenodo"
    assert [request.status for request in loopback_server.requests] == [200]
    # The downloaded tarball is deleted once installed.
    assert list((db_folder / ".downloads").iterdir()) == []
    assert _leftovers(db_folder) == []
    download_log = [message for message in logs if "downloading" in message]
    assert len(download_log) == 1
    assert fixture_pins.url in download_log[0]
    assert str(db_folder / ".downloads" / "phu-db.tar.gz") in download_log[0]


def test_install_progress_goes_to_stderr_by_default(
    monkeypatch,
    capsys,
    user_bundle,
    db_folder,
    fixture_pins,
    fixture_manifest,
    fixture_manifest_bytes,
):
    monkeypatch.setenv("PHU_DB_BUNDLE", str(user_bundle))

    install(
        db_folder=db_folder,
        pins=fixture_pins,
        manifest=fixture_manifest,
        vendored_manifest_bytes=fixture_manifest_bytes,
        fetch=FetchSpy(None),
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert f"installed {fixture_pins.bundle_version}" in captured.err


def test_valid_install_is_not_repeated(run, db_folder, fetch, fixture_pins, logs):
    target = run()
    record = (target / "installed.json").read_bytes()
    fetch.calls.clear()
    fetch.source = None

    assert run() == target

    assert fetch.calls == []
    assert (target / "installed.json").read_bytes() == record
    assert "already installed" in logs[-1]


def test_failed_deep_check_reinstalls_when_asked(run, db_folder, fetch, fixture_pins):
    target = run()
    _flip_byte(target / "pfam" / "profiles.hmm")
    fetch.calls.clear()

    # A shallow re-check under the lock would accept the tree; a deep one does not.
    run(deep=True)

    assert len(fetch.calls) == 1
    _assert_installed(target)


def test_other_versions_are_kept(run, db_folder, fixture_pins, logs):
    # Another phu release sharing the folder may be running against its own version.
    old = db_folder / "phu-db" / "phu-db-0.1.0"
    (old / "pfam").mkdir(parents=True)
    (old / "pfam" / "profiles.hmm").write_text("old\n")

    target = run()
    run()

    assert sorted(path.name for path in (db_folder / "phu-db").iterdir()) == [
        "phu-db-0.1.0",
        target.name,
    ]
    assert (old / "pfam" / "profiles.hmm").read_text() == "old\n"
    assert not any(str(old) in message for message in logs)


def test_stale_invalid_target_is_replaced(run, db_folder, fixture_pins):
    target = _target(db_folder, fixture_pins)
    (target / "pfam").mkdir(parents=True)
    (target / "pfam" / "profiles.hmm").write_text("truncated")
    (target / "junk.txt").write_text("left by a crash")
    # A temp directory a killed install left behind.
    (db_folder / ".phu-db-tmp-dead" / "phu-db").mkdir(parents=True)

    assert run() == target

    _assert_installed(target)
    assert not (target / "junk.txt").exists()
    assert _leftovers(db_folder) == []


def test_bundle_missing_a_database_installs_nothing(
    run, make_bundle, db_folder, fixture_pins, fixture_manifest
):
    avg = fixture_manifest.databases["avg"].member_paths()
    tarball = make_bundle(drop=avg)
    data = tarball.read_bytes()
    pins = replace(fixture_pins, sha256=hashlib.sha256(data).hexdigest(), size=len(data))

    def no_fetch(url, dest, sha256, size):
        raise AssertionError("nothing may be downloaded")

    with pytest.raises(BundleMismatch, match="does not include avg"):
        run(pins=pins, bundle_path=tarball, fetch=no_fetch)

    assert _leftovers(db_folder) == []
    assert not (db_folder / "phu-db").exists()
    assert tarball.read_bytes() == data


def test_bundle_with_altered_content_is_not_installed(
    run, make_bundle, db_folder, fixture_pins, fixture_manifest
):
    # Same size, so only install()'s deep validation of the extracted tree sees it.
    original = (PREPARED / "pfam" / "profiles.hmm").read_bytes()
    altered = bytearray(original)
    index = original.index(b"\nHMM ") + 1
    altered[index] = ord("h")
    tarball = make_bundle(replace={"phu-db/pfam/profiles.hmm": bytes(altered)})
    data = tarball.read_bytes()
    pins = replace(fixture_pins, sha256=hashlib.sha256(data).hexdigest(), size=len(data))

    with pytest.raises(TreeInvalid) as raised:
        run(pins=pins, bundle_path=tarball)

    assert [problem.split(" ")[:3] for problem in raised.value.problems] == [
        ["pfam/profiles.hmm", "has", "sha256"]
    ]
    assert not (db_folder / "phu-db").exists()
    assert _leftovers(db_folder) == []
    assert tarball.read_bytes() == data


def test_interrupted_clean_up_is_finished_by_the_next_install(
    run, monkeypatch, db_folder, fetch, fixture_pins, logs
):
    real = db_bundle._clean_up

    def killed(*args):
        raise KeyboardInterrupt("job killed while cleaning up")

    monkeypatch.setattr(db_bundle, "_clean_up", killed)
    with pytest.raises(KeyboardInterrupt):
        run()
    downloads = db_folder / ".downloads"
    assert (downloads / "phu-db.tar.gz").is_file()
    # And a later failed download, and a killed extraction.
    (downloads / "phu-db.tar.gz.part").write_bytes(b"partial")
    (db_folder / ".phu-db-tmp-dead").mkdir()
    monkeypatch.setattr(db_bundle, "_clean_up", real)
    fetch.calls.clear()
    fetch.source = None

    target = run()

    assert fetch.calls == []
    assert any("already installed" in message for message in logs)
    assert sorted(path.name for path in (db_folder / "phu-db").iterdir()) == [target.name]
    assert list(downloads.iterdir()) == []
    assert _leftovers(db_folder) == []
    _assert_installed(target)


def test_clean_up_never_deletes_a_phu_db_bundle_in_downloads(
    run, monkeypatch, db_folder, fetch, fixture_bundle
):
    # The user pointed PHU_DB_BUNDLE at the file phu itself would download to.
    mine = db_folder / ".downloads" / "phu-db.tar.gz"
    mine.parent.mkdir(parents=True)
    shutil.copyfile(fixture_bundle, mine)
    monkeypatch.setenv("PHU_DB_BUNDLE", str(mine))
    fetch.source = None

    run()
    run()

    assert mine.read_bytes() == fixture_bundle.read_bytes()


# Refusing tarballs.


def _assert_mismatch_message(message, pins):
    assert "different or corrupt bundle" in message
    assert pins.bundle_version in message
    assert pins.sha256 in message
    assert str(pins.size) in message


@pytest.mark.parametrize("change", ["flip", "append"])
def test_phu_db_bundle_with_wrong_digest_is_refused(
    run, monkeypatch, user_bundle, db_folder, fetch, fixture_pins, change
):
    data = bytearray(user_bundle.read_bytes())
    if change == "flip":
        data[len(data) // 2] ^= 0xFF
    else:
        data += b"\0"
    user_bundle.write_bytes(bytes(data))
    monkeypatch.setenv("PHU_DB_BUNDLE", str(user_bundle))
    fetch.source = None

    with pytest.raises(TarballMismatch) as raised:
        run()

    _assert_mismatch_message(str(raised.value), fixture_pins)
    assert str(user_bundle) in str(raised.value)
    assert user_bundle.read_bytes() == bytes(data)
    assert fetch.calls == []
    assert not (db_folder / "phu-db").exists()
    assert _leftovers(db_folder) == []


def test_corrupt_download_is_refused(run, loopback_server, db_folder, fixture_bundle, fixture_pins):
    loopback_server.files["phu-db.tar.gz"] = fixture_bundle.read_bytes()
    loopback_server.switches.corrupt_offset = 100

    def fetch(url, dest, sha256, size):
        return download_verified(
            loopback_server.url("phu-db.tar.gz"), dest, sha256, size, timeout=5
        )

    with pytest.raises(TarballMismatch) as raised:
        run(fetch=fetch)

    _assert_mismatch_message(str(raised.value), fixture_pins)
    assert list((db_folder / ".downloads").iterdir()) == []
    assert not (db_folder / "phu-db").exists()


def test_fetch_result_is_verified_again(run, db_folder, fetch, tmp_path, fixture_bundle):
    # An injected fetch that does not verify cannot slip a bad tarball through.
    bad = tmp_path / "bad.tar.gz"
    data = bytearray(fixture_bundle.read_bytes())
    data[50] ^= 0xFF
    bad.write_bytes(bytes(data))
    fetch.source = bad

    with pytest.raises(TarballMismatch):
        run()

    assert list((db_folder / ".downloads").iterdir()) == []
    assert not (db_folder / "phu-db").exists()


def test_insufficient_space_raises_before_download(
    run, monkeypatch, db_folder, fetch, fixture_pins, fixture_manifest
):
    monkeypatch.setattr(shutil, "disk_usage", lambda path: DiskUsage(10**6, 10**6 - 1000, 1000))
    prepared = sum(
        item.bytes for entry in fixture_manifest.databases.values() for item in entry.prepared_files
    )
    required = required_install_bytes(fixture_pins.size, prepared)

    with pytest.raises(InsufficientSpace) as raised:
        run()

    assert (raised.value.required, raised.value.available) == (required, 1000)
    assert f"{required} bytes required" in str(raised.value)
    assert "1000 bytes available" in str(raised.value)
    assert fetch.calls == []
    assert not (db_folder / ".downloads").exists()


def test_insufficient_space_with_phu_db_bundle_counts_no_tarball(
    run, monkeypatch, db_folder, user_bundle, fixture_manifest
):
    # The user's tarball is read where it is; only the prepared tree needs room.
    monkeypatch.setenv("PHU_DB_BUNDLE", str(user_bundle))
    monkeypatch.setattr(shutil, "disk_usage", lambda path: DiskUsage(10**6, 10**6 - 1000, 1000))
    prepared = sum(
        item.bytes for entry in fixture_manifest.databases.values() for item in entry.prepared_files
    )

    with pytest.raises(InsufficientSpace) as raised:
        run()

    assert raised.value.required == required_install_bytes(0, prepared)
    assert not (db_folder / "phu-db").exists()


# Validation.


@pytest.fixture
def installed(run):
    return run()


def test_installed_tree_is_valid(installed, fixture_manifest):
    assert tree_problems(installed, fixture_manifest, deep=False) == []
    assert tree_problems(installed, fixture_manifest, deep=True) == []
    validate_tree(installed, fixture_manifest, deep=True)


def test_deep_validation_catches_one_altered_byte(installed, fixture_manifest):
    _flip_byte(installed / "pfam" / "profiles.hmm")

    validate_tree(installed, fixture_manifest, deep=False)
    with pytest.raises(TreeInvalid) as raised:
        validate_tree(installed, fixture_manifest, deep=True)

    assert len(raised.value.problems) == 1
    assert raised.value.problems[0].startswith("pfam/profiles.hmm has sha256 ")


def test_shallow_validation_catches_a_truncated_file(installed, fixture_manifest):
    path = installed / "avg" / "v_scores.tsv"
    data = path.read_bytes()
    path.write_bytes(data[:-1])

    with pytest.raises(TreeInvalid) as raised:
        validate_tree(installed, fixture_manifest, deep=False)

    assert raised.value.problems == [
        f"avg/v_scores.tsv is {len(data) - 1} bytes, expected {len(data)}"
    ]


def test_shallow_validation_catches_a_tsv_count_mismatch(installed, fixture_manifest):
    path = installed / "kofam" / "ko_list.tsv"
    data = path.read_bytes()
    lines = data.split(b"\n")
    # Blank out the last row, keeping the byte size.
    lines[-2] = b" " * len(lines[-2])
    path.write_bytes(b"\n".join(lines))

    with pytest.raises(TreeInvalid) as raised:
        validate_tree(installed, fixture_manifest, deep=False)

    assert raised.value.problems == ["kofam/ko_list.tsv has 3 records, expected 4"]


def test_hmm_counts_are_checked_deep_only(installed, fixture_manifest):
    path = installed / "dbcan" / "profiles.hmm"
    data = path.read_bytes()
    path.write_bytes(data.replace(b"\nNAME ", b"\nname ", 1))

    assert tree_problems(installed, fixture_manifest, deep=False) == []
    problems = tree_problems(installed, fixture_manifest, deep=True)
    assert problems[-1] == "dbcan/profiles.hmm has 1 records, expected 2"


def test_missing_files_and_database_subsets(installed, fixture_manifest):
    shutil.rmtree(installed / "dbcan")

    assert tree_problems(installed, fixture_manifest, deep=False, databases=["pfam", "avg"]) == []
    assert tree_problems(installed, fixture_manifest, deep=False, databases=["dbcan"]) == [
        "dbcan/profiles.hmm is missing",
        "dbcan/family_substrate.tsv is missing",
        "dbcan/pul_rules.tsv is missing",
    ]
    with pytest.raises(ValueError):
        tree_problems(installed, fixture_manifest, deep=False, databases=["vscore"])


def test_file_removed_during_the_check_is_reported_missing(
    installed, monkeypatch, fixture_manifest
):
    gone = installed / "avg" / "v_scores.tsv"
    real = Path.is_file

    def is_file(path):
        # Moved away by an install between is_file() and stat().
        if path == gone:
            gone.unlink(missing_ok=True)
            return True
        return real(path)

    monkeypatch.setattr(Path, "is_file", is_file)

    assert tree_problems(installed, fixture_manifest, deep=False) == ["avg/v_scores.tsv is missing"]


def test_tree_invalid_takes_problems_first():
    error = TreeInvalid(["pfam/profiles.hmm is missing"])

    assert error.problems == ["pfam/profiles.hmm is missing"]
    assert error.root is None
    assert "pfam/profiles.hmm is missing" in str(error)


def test_installed_check_needs_installed_json(installed, db_folder, fixture_pins, fixture_manifest):
    (installed / "installed.json").unlink()

    assert installed_problems(db_folder, fixture_pins, fixture_manifest, deep=False)


def _reference_count(data: bytes, file_format: str) -> int:
    lines = data.split(b"\n")
    if file_format == "hmm":
        return sum(1 for line in lines if line.startswith(b"NAME "))
    return sum(1 for line in lines if line.strip())


COUNT_SAMPLES = [
    b"",
    b"\n",
    b"a",
    b"a\n\n  \n\t\nb",
    b"header\r\nrow\r\n\r\n",
    b"one\rtwo\nthree",  # a bare \r does not end a line
    b"   \n\x0b\x0c\n x \n",
    b"NAME a\nNAME b\n NAME c\nNAMEd\nxNAME e\nNAME ",
    b"HMMER3/f\nNAME  x\n//\nNAME y",
]


@pytest.mark.parametrize("chunk_size", [1, 2, 3, 5, 6, 7, 1024 * 1024])
@pytest.mark.parametrize("file_format", ["tsv", "hmm"])
def test_counts_follow_contract_across_chunk_boundaries(
    tmp_path, monkeypatch, chunk_size, file_format
):
    monkeypatch.setattr(db_bundle, "CHUNK_SIZE", chunk_size)
    samples = COUNT_SAMPLES + [path.read_bytes() for path in sorted(PREPARED.rglob("*.*"))]

    for index, data in enumerate(samples):
        path = tmp_path / f"{index}.{file_format}"
        path.write_bytes(data)
        sha256, count = db_bundle._scan_file(path, digest=True, file_format=file_format)
        assert count == _reference_count(data, file_format), data[:40]
        assert sha256 == hashlib.sha256(data).hexdigest()


# Locking.


def _no_flock(monkeypatch, code=errno.ENOLCK):
    def flock(fd, operation):
        raise OSError(code, os.strerror(code))

    monkeypatch.setattr(fcntl, "flock", flock)


def _write_lock_dir(db_folder, *, created, owner=True, token="x"):
    lock_dir = db_folder / ".phu-db.lock.d"
    lock_dir.mkdir(parents=True)
    if owner:
        (lock_dir / "owner.json").write_text(
            json.dumps({"pid": 1, "hostname": "elsewhere", "created": created, "token": token})
        )
    return lock_dir


def test_mkdir_lock_is_used_when_flock_is_unsupported(run, monkeypatch, db_folder, fetch, logs):
    _no_flock(monkeypatch)
    seen = {}

    def during(dest):
        owner = json.loads((db_folder / ".phu-db.lock.d" / "owner.json").read_text())
        seen.update(owner)

    fetch.during = during

    target = run()

    _assert_installed(target)
    assert seen["pid"] == os.getpid()
    assert seen["hostname"] == socket.gethostname()
    assert abs(seen["created"] - time.time()) < 600
    assert not (db_folder / ".phu-db.lock.d").exists()
    assert any("does not support flock" in message for message in logs)


@pytest.mark.parametrize("code", [errno.EOPNOTSUPP, errno.EINVAL])
def test_other_unsupported_flock_errors_fall_back(monkeypatch, db_folder, logs, code):
    _no_flock(monkeypatch, code)

    with install_lock(db_folder, log=logs.append):
        assert (db_folder / ".phu-db.lock.d").is_dir()

    assert not (db_folder / ".phu-db.lock.d").exists()


def test_other_flock_errors_are_raised(monkeypatch, db_folder):
    _no_flock(monkeypatch, errno.EBADF)

    with pytest.raises(OSError):
        with install_lock(db_folder):
            pass


def test_fresh_lock_dir_blocks(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    lock_dir = _write_lock_dir(db_folder, created=time.time() - 60)

    with pytest.raises(LockTimeout):
        with install_lock(db_folder, log=logs.append, timeout=0.2, poll=0.01):
            pass

    assert lock_dir.is_dir()
    assert json.loads((lock_dir / "owner.json").read_text())["token"] == "x"
    assert sum("waiting" in message for message in logs) == 1


def test_stale_lock_dir_is_broken(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    _write_lock_dir(db_folder, created=time.time() - 7 * 3600)

    with install_lock(db_folder, log=logs.append, timeout=1, poll=0.01):
        owner = json.loads((db_folder / ".phu-db.lock.d" / "owner.json").read_text())
        assert owner["pid"] == os.getpid()

    assert any("breaking stale install lock" in message for message in logs)
    assert not any("waiting" in message for message in logs)
    assert sorted(path.name for path in db_folder.iterdir()) == [".phu-db.lock"]


def test_lock_dir_age_uses_the_injected_clock(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    _write_lock_dir(db_folder, created=time.time())

    with install_lock(
        db_folder, log=logs.append, timeout=1, clock=lambda: time.time() + 6 * 3600 + 1
    ):
        pass

    assert any("breaking stale install lock" in message for message in logs)


def test_lock_dir_with_unreadable_owner_ages_by_mtime(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    lock_dir = _write_lock_dir(db_folder, created=None, owner=False)
    (lock_dir / "owner.json").write_text("{not json")
    old = time.time() - 7 * 3600
    os.utime(lock_dir, (old, old))

    with install_lock(db_folder, log=logs.append, timeout=1, poll=0.01):
        pass

    assert any("breaking stale install lock" in message for message in logs)


def test_fresh_lock_dir_with_unreadable_owner_blocks(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    lock_dir = _write_lock_dir(db_folder, created=None, owner=False)
    (lock_dir / "owner.json").write_text("{not json")

    with pytest.raises(LockTimeout):
        with install_lock(db_folder, log=logs.append, timeout=0.2, poll=0.01):
            pass

    assert (lock_dir / "owner.json").read_text() == "{not json"


def test_lock_dir_never_appears_without_its_owner(monkeypatch, db_folder):
    _no_flock(monkeypatch)
    seen = []
    real_rename = os.rename

    def rename(source, dest):
        # The rename that publishes the lock already carries the owner file.
        if Path(dest).name == ".phu-db.lock.d":
            seen.append(json.loads((Path(source) / "owner.json").read_text())["pid"])
        real_rename(source, dest)

    monkeypatch.setattr(os, "rename", rename)

    with install_lock(db_folder):
        assert not (db_folder / ".phu-db.lock.d").is_symlink()

    assert seen == [os.getpid()]
    assert sorted(path.name for path in db_folder.iterdir()) == [".phu-db.lock"]


def _racing_lock_age(monkeypatch, lock_dir, effects):
    """Run ``effects[i]`` right after the i-th staleness check, as another waiter would."""
    real = db_bundle._lock_age
    calls = []

    def lock_age(path, clock):
        age = real(path, clock)
        if len(calls) < len(effects):
            effects[len(calls)](Path(path))
        calls.append(Path(path))
        return age

    monkeypatch.setattr(db_bundle, "_lock_age", lock_age)
    return calls


def test_racing_breakers_leave_a_fresh_lock_alone(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    lock_dir = _write_lock_dir(db_folder, created=time.time() - 7 * 3600)

    def other_waiter_breaks_and_takes(path):
        # Between this waiter judging the lock stale and breaking it, another
        # waiter breaks the same stale lock and takes a fresh one.
        shutil.rmtree(lock_dir)
        _write_lock_dir(db_folder, created=time.time(), token="other")

    _racing_lock_age(monkeypatch, lock_dir, [other_waiter_breaks_and_takes])

    with pytest.raises(LockTimeout):
        with install_lock(db_folder, log=logs.append, timeout=0.2, poll=0.01):
            pass

    assert json.loads((lock_dir / "owner.json").read_text())["token"] == "other"
    assert not any("breaking" in message for message in logs)
    assert sorted(path.name for path in db_folder.iterdir()) == [".phu-db.lock", ".phu-db.lock.d"]


def test_lock_taken_while_a_fresh_lock_was_aside_is_reported(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    lock_dir = _write_lock_dir(db_folder, created=time.time() - 7 * 3600)

    def other_waiter_breaks_and_takes(path):
        shutil.rmtree(lock_dir)
        _write_lock_dir(db_folder, created=time.time(), token="other")

    def third_waiter_takes(path):
        # The fresh lock is aside for a moment, and a third waiter takes the free one.
        assert path != lock_dir
        _write_lock_dir(db_folder, created=time.time(), token="third")

    _racing_lock_age(monkeypatch, lock_dir, [other_waiter_breaks_and_takes, third_waiter_takes])

    with pytest.raises(LockTimeout):
        with install_lock(db_folder, log=logs.append, timeout=0.2, poll=0.01):
            pass

    assert json.loads((lock_dir / "owner.json").read_text())["token"] == "third"
    assert any("was taken twice" in message for message in logs)
    assert sorted(path.name for path in db_folder.iterdir()) == [".phu-db.lock", ".phu-db.lock.d"]


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root ignores modes")
def test_read_only_lock_file_is_still_locked(db_folder):
    db_folder.mkdir()
    (db_folder / ".phu-db.lock").touch(0o444)
    (db_folder / ".phu-db.lock").chmod(0o444)

    with install_lock(db_folder):
        with pytest.raises(LockTimeout):
            with install_lock(db_folder, timeout=0.2):
                pass


def test_lock_dir_taken_over_by_another_owner_is_left_alone(monkeypatch, db_folder):
    _no_flock(monkeypatch)
    lock_dir = db_folder / ".phu-db.lock.d"

    with install_lock(db_folder):
        # A waiter broke this lock as stale and took it.
        (lock_dir / "owner.json").write_text(json.dumps({"token": "theirs", "created": 0}))

    assert lock_dir.is_dir()


def test_flock_wait_is_logged_once_and_can_time_out(db_folder, logs):
    with install_lock(db_folder):
        with pytest.raises(LockTimeout):
            with install_lock(db_folder, log=logs.append, timeout=0.2):
                pass

    assert sum("waiting" in message for message in logs) == 1
    with install_lock(db_folder, timeout=0.2):
        pass


def _wait_for(predicate, timeout=60.0):
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError("timed out waiting for the other process")
        time.sleep(0.01)


def _install_worker(name, db_folder, bundle, pins, manifest_bytes, signals):
    """One of two processes installing into the same db folder (see below)."""
    signals = Path(signals)

    def log(message):
        with (signals / f"{name}.log").open("a") as handle:
            handle.write(message + "\n")

    def fetch(url, dest, sha256, size):
        with (signals / "fetches").open("a") as handle:
            handle.write(f"{name}\n")
        (signals / "fetching").touch()
        _wait_for((signals / "go").exists)
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(bundle, dest)
        return Path(dest)

    target = install(
        db_folder=Path(db_folder),
        pins=pins,
        manifest=loads_manifest(manifest_bytes),
        vendored_manifest_bytes=manifest_bytes,
        fetch=fetch,
        log=log,
    )
    stat = (target / "installed.json").stat()
    (signals / f"{name}.done").write_text(json.dumps([stat.st_ino, stat.st_mtime_ns]))


def test_concurrent_installs_wait_and_install_once(
    tmp_path, db_folder, fixture_bundle, fixture_pins, fixture_manifest_bytes
):
    signals = tmp_path / "signals"
    signals.mkdir()
    context = multiprocessing.get_context("spawn")
    args = (str(db_folder), str(fixture_bundle), fixture_pins, fixture_manifest_bytes, str(signals))
    first = context.Process(target=_install_worker, args=("first", *args))
    second = context.Process(target=_install_worker, args=("second", *args))
    try:
        first.start()
        _wait_for((signals / "fetching").exists)
        second.start()
        # The second process is blocked on the lock while the first downloads.
        _wait_for(
            lambda: (signals / "second.log").exists()
            and "waiting" in (signals / "second.log").read_text()
        )
        (signals / "go").touch()
        first.join(60)
        second.join(60)
        assert (first.exitcode, second.exitcode) == (0, 0)
    finally:
        for process in (first, second):
            if process.is_alive():
                process.terminate()
                process.join()

    assert (signals / "fetches").read_text() == "first\n"
    second_log = (signals / "second.log").read_text()
    assert "already installed" in second_log
    assert "extracting" not in second_log
    # installed.json is still the one the first process wrote.
    stat = (_target(db_folder, fixture_pins) / "installed.json").stat()
    assert json.loads((signals / "first.done").read_text()) == [stat.st_ino, stat.st_mtime_ns]
    assert json.loads((signals / "second.done").read_text()) == [stat.st_ino, stat.st_mtime_ns]
    _assert_installed(_target(db_folder, fixture_pins))


def test_shallow_validation_never_opens_an_hmm_file(installed, fixture_manifest, monkeypatch):
    real = db_bundle._scan_file
    scanned = []

    def recording(path, **kwargs):
        scanned.append(path.name)
        return real(path, **kwargs)

    monkeypatch.setattr(db_bundle, "_scan_file", recording)

    assert tree_problems(installed, fixture_manifest, deep=False) == []
    assert scanned and "profiles.hmm" not in scanned
    scanned.clear()
    assert tree_problems(installed, fixture_manifest, deep=True) == []
    assert scanned.count("profiles.hmm") == 3


def test_write_atomically_leaves_one_readable_file(tmp_path):
    path = tmp_path / "models" / "GH5.hmm"

    db_bundle.write_atomically(path, b"first")
    db_bundle.write_atomically(path, b"second")

    assert path.read_bytes() == b"second"
    assert [item.name for item in path.parent.iterdir()] == ["GH5.hmm"]
    umask = os.umask(0)
    os.umask(umask)
    assert path.stat().st_mode & 0o777 == 0o666 & ~umask


def test_a_live_holder_keeps_its_lock_dir_fresh(monkeypatch, db_folder, logs):
    _no_flock(monkeypatch)
    now = [time.time()]

    with install_lock(db_folder, clock=lambda: now[0], refresh=0.01):
        owner_path = db_folder / ".phu-db.lock.d" / "owner.json"
        # Seven hours pass; the holder is slow but alive and re-stamps its lock.
        now[0] += 7 * 3600
        deadline = time.monotonic() + 5
        while "refreshed" not in json.loads(owner_path.read_text()):
            assert time.monotonic() < deadline, "the lock was never re-stamped"
            time.sleep(0.01)

        with pytest.raises(LockTimeout):
            with install_lock(
                db_folder, log=logs.append, timeout=0.2, poll=0.01, clock=lambda: now[0]
            ):
                pass

        assert json.loads(owner_path.read_text())["pid"] == os.getpid()

    assert not any("breaking stale install lock" in message for message in logs)
    assert sorted(path.name for path in db_folder.iterdir()) == [".phu-db.lock"]
