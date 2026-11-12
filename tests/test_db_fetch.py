import hashlib
import shutil
import socket
from collections import namedtuple
from pathlib import Path

import pytest

from phu import db_fetch
from phu.db_fetch import (
    DatabaseError,
    DownloadError,
    InsecureRedirect,
    InsecureURL,
    InsufficientSpace,
    IntegrityError,
    download_verified,
    ensure_free_space,
    redirect_allowed,
    required_install_bytes,
)
from phu.db_manifest import load_manifest

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"

# Deterministic, incompressible-looking payload of a few hundred KiB.
DATA = b"".join(hashlib.sha256(str(i).encode()).digest() for i in range(10_000))
DATA_SHA256 = hashlib.sha256(DATA).hexdigest()

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


@pytest.fixture
def sleeps():
    """An injected sleep that records the delays instead of waiting."""
    return []


@pytest.fixture
def server(loopback_server):
    loopback_server.files["bundle.tar.gz"] = DATA
    return loopback_server


def _fetch(server, dest, sleeps, *, sha256=DATA_SHA256, size=len(DATA), **kwargs):
    kwargs.setdefault("timeout", 5)
    return download_verified(
        server.url("bundle.tar.gz"), dest, sha256, size, sleep=sleeps.append, **kwargs
    )


def _part(dest: Path) -> Path:
    return dest.with_name(dest.name + ".part")


def test_full_download_of_a_fixture_source(loopback_server, tmp_path, sleeps):
    loopback_server.serve_dir(FIXTURE_DIR / "upstream")
    source = load_manifest(FIXTURE_DIR / "phu-db.json").databases["pfam"].sources[0]
    dest = tmp_path / "downloads" / source.name

    result = download_verified(
        loopback_server.url(f"pfam/{source.name}"),
        dest,
        source.sha256,
        source.bytes,
        sleep=sleeps.append,
    )

    assert result == dest
    assert dest.read_bytes() == (FIXTURE_DIR / "upstream" / "pfam" / source.name).read_bytes()
    assert not _part(dest).exists()
    assert [(r.range, r.status) for r in loopback_server.requests] == [(None, 200)]
    assert sleeps == []


def test_resume_after_a_cut_sends_range_and_appends(server, tmp_path, sleeps):
    server.switches.cut_after = 70_000
    server.switches.cut_times = 1
    dest = tmp_path / "phu-db.tar.gz"

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [
        (None, 200),
        ("bytes=70000-", 206),
    ]
    assert sleeps == [1.0]


def test_existing_part_is_rehashed_and_resumed(server, tmp_path, sleeps):
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:5000])

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [("bytes=5000-", 206)]
    assert not _part(dest).exists()


def test_existing_part_with_wrong_bytes_is_fetched_again_from_zero(server, tmp_path, sleeps):
    # The resumed tail is right; the re-hashed prefix (say, of an older bundle) is not.
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(b"\0" * 300_000)

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [
        ("bytes=300000-", 206),
        (None, 200),
    ]
    assert not _part(dest).exists()
    assert sleeps == []


def test_complete_part_with_wrong_bytes_is_fetched_again_from_zero(server, tmp_path, sleeps):
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(b"\0" * len(DATA))

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [(None, 200)]


def test_fresh_pass_after_a_bad_part_is_made_only_once(server, tmp_path, sleeps):
    server.switches.corrupt_offset = 12_345
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:5000])

    with pytest.raises(IntegrityError, match="sha256"):
        _fetch(server, dest, sleeps)

    assert [(r.range, r.status) for r in server.requests] == [
        ("bytes=5000-", 206),
        (None, 200),
    ]
    assert not _part(dest).exists()
    assert not dest.exists()


def test_complete_part_is_verified_without_a_request(server, tmp_path, sleeps):
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA)

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert server.requests == []


def test_part_longer_than_the_pin_is_discarded(server, tmp_path, sleeps):
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA + b"extra")

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [(None, 200)]


def test_200_reply_of_another_size_keeps_the_resumed_part(server, tmp_path, sleeps):
    # A proxy or portal page answering a Range request must not wipe a good prefix.
    server.switches.ignore_range = True
    server.files["bundle.tar.gz"] = b"<html>log in first</html>"
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:200_000])

    with pytest.raises(IntegrityError, match="25 bytes"):
        _fetch(server, dest, sleeps)

    assert _part(dest).read_bytes() == DATA[:200_000]
    assert not dest.exists()


def test_200_reply_to_a_range_request_restarts_from_zero(server, tmp_path, sleeps):
    server.switches.ignore_range = True
    dest = tmp_path / "phu-db.tar.gz"
    # Wrong bytes: only a true restart can still verify.
    _part(dest).write_bytes(b"\xff" * 5000)

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [("bytes=5000-", 200)]


def test_corrupted_byte_raises_and_leaves_nothing(server, tmp_path, sleeps):
    server.switches.corrupt_offset = 12_345
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(IntegrityError, match="sha256"):
        _fetch(server, dest, sleeps)

    assert not _part(dest).exists()
    assert not dest.exists()
    assert sleeps == []


@pytest.mark.parametrize("size", [len(DATA) - 1, len(DATA) + 1])
def test_wrong_size_raises_and_leaves_nothing(server, tmp_path, sleeps, size):
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(IntegrityError):
        _fetch(server, dest, sleeps, size=size)

    assert not _part(dest).exists()
    assert not dest.exists()
    assert sleeps == []


def test_short_file_without_content_length_raises_integrity_error(server, tmp_path, sleeps):
    # Like a PHP download endpoint: no length, so the short end looks like a cut
    # until the resume's 416 reports the server's total.
    server.files["bundle.tar.gz"] = DATA[:1000]
    server.switches.omit_length = True
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(IntegrityError, match="1000 bytes"):
        _fetch(server, dest, sleeps)

    assert [(r.range, r.status) for r in server.requests] == [
        (None, 200),
        ("bytes=1000-", 416),
    ]
    assert not _part(dest).exists()
    assert not dest.exists()


def test_download_without_content_length_verifies(server, tmp_path, sleeps):
    server.switches.omit_length = True
    dest = tmp_path / "phu-db.tar.gz"

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [(r.range, r.status) for r in server.requests] == [(None, 200)]


def test_redirect_to_plain_http_off_loopback_is_refused(server, tmp_path, sleeps):
    # The target never resolves; reaching it would also trip loopback_only.
    server.switches.redirects["bundle.tar.gz"] = "http://phu-db-fixture.invalid/x.tar.gz"
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(InsecureRedirect, match="phu-db-fixture.invalid"):
        _fetch(server, dest, sleeps)

    assert [(r.path, r.status) for r in server.requests] == [("bundle.tar.gz", 302)]
    assert sleeps == []
    assert not dest.exists()


def test_redirect_between_loopback_urls_is_followed(server, tmp_path, sleeps):
    server.files["moved.tar.gz"] = DATA
    server.switches.redirects["bundle.tar.gz"] = server.url("moved.tar.gz")
    dest = tmp_path / "phu-db.tar.gz"

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [r.path for r in server.requests] == ["bundle.tar.gz", "moved.tar.gz"]


def test_redirect_rule():
    assert redirect_allowed("https://zenodo.org/a", "https://cdn.example.org/b")
    assert redirect_allowed("http://127.0.0.1:8000/a", "https://zenodo.org/b")
    assert redirect_allowed("http://127.0.0.1:8000/a", "http://localhost:9000/b")
    assert redirect_allowed("http://[::1]:8000/a", "http://127.0.0.1:9000/b")
    assert not redirect_allowed("https://zenodo.org/a", "http://zenodo.org/b")
    assert not redirect_allowed("https://zenodo.org/a", "http://127.0.0.1/b")
    assert not redirect_allowed("http://127.0.0.1/a", "http://zenodo.org/b")
    assert not redirect_allowed("http://127.0.0.1/a", "ftp://127.0.0.1/b")
    # Only the listed loopback hosts, not all of 127.0.0.0/8.
    assert not redirect_allowed("http://127.0.0.2/a", "http://127.0.0.1/b")
    assert not redirect_allowed("http://[::ffff:127.0.0.1]/a", "http://127.0.0.1/b")


def test_malformed_url_is_not_retried(tmp_path, sleeps):
    with pytest.raises(DownloadError, match="port"):
        download_verified(
            "http://127.0.0.1:abc/b", tmp_path / "b", DATA_SHA256, len(DATA), sleep=sleeps.append
        )

    assert sleeps == []


@pytest.mark.parametrize(
    "url",
    [
        "http://zenodo.org/records/1/files/phu-db.tar.gz",
        "ftp://ftp.ebi.ac.uk/pub/Pfam-A.hmm.gz",
        "file:///etc/passwd",
        "zenodo.org/phu-db.tar.gz",
    ],
)
def test_insecure_initial_url_is_refused(tmp_path, sleeps, url):
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(InsecureRedirect):
        download_verified(url, dest, DATA_SHA256, len(DATA), sleep=sleeps.append)

    assert sleeps == []
    assert not _part(dest).exists()


def test_stall_past_the_timeout_is_retried_then_gives_up(server, tmp_path, sleeps):
    server.switches.stall = 30
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(DownloadError, match="3 attempts"):
        _fetch(server, dest, sleeps, timeout=0.2, retries=2)

    assert len(server.requests) == 3
    assert sleeps == [1.0, 2.0]
    assert not dest.exists()


def test_retries_exhausted_after_five_retries(server, tmp_path, sleeps):
    server.switches.cut_after = 0
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(DownloadError, match="6 attempts"):
        _fetch(server, dest, sleeps, backoff=0.5)

    assert len(server.requests) == 6
    assert sleeps == [0.5, 1.0, 2.0, 4.0, 8.0]
    assert not dest.exists()


def test_transient_failure_then_success(server, tmp_path, sleeps):
    server.switches.fail_status = 503
    dest = tmp_path / "phu-db.tar.gz"

    _fetch(server, dest, sleeps)

    assert dest.read_bytes() == DATA
    assert [r.status for r in server.requests] == [503, 200]
    assert sleeps == [1.0]


def test_client_error_is_not_retried(server, tmp_path, sleeps):
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(DownloadError, match="404"):
        download_verified(
            server.url("missing.tar.gz"), dest, DATA_SHA256, len(DATA), sleep=sleeps.append
        )

    assert len(server.requests) == 1
    assert sleeps == []


def test_preflight_refuses_before_any_request(server, tmp_path, sleeps, monkeypatch):
    monkeypatch.setattr(shutil, "disk_usage", lambda path: DiskUsage(10**9, 10**9, 100))
    dest = tmp_path / "phu-db.tar.gz"

    with pytest.raises(InsufficientSpace) as caught:
        _fetch(server, dest, sleeps)

    assert caught.value.required == len(DATA)
    assert caught.value.available == 100
    assert str(len(DATA)) in str(caught.value)
    assert server.requests == []
    assert not _part(dest).exists()


def test_preflight_counts_only_the_bytes_still_to_fetch(server, tmp_path, sleeps, monkeypatch):
    monkeypatch.setattr(shutil, "disk_usage", lambda path: DiskUsage(10**9, 10**9, 100))
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:5000])

    with pytest.raises(InsufficientSpace) as caught:
        _fetch(server, dest, sleeps)

    assert caught.value.required == len(DATA) - 5000
    assert server.requests == []


def test_ensure_free_space_measures_the_nearest_existing_parent(tmp_path, monkeypatch):
    seen = []

    def disk_usage(path):
        seen.append(Path(path))
        return DiskUsage(10**9, 0, 1000)

    monkeypatch.setattr(shutil, "disk_usage", disk_usage)
    ensure_free_space(tmp_path / "not" / "yet", 1000)
    with pytest.raises(InsufficientSpace) as caught:
        ensure_free_space(tmp_path / "not" / "yet", 1001)

    assert seen == [tmp_path, tmp_path]
    assert (caught.value.required, caught.value.available) == (1001, 1000)


def test_preflight_requires_the_install_total(tmp_path, monkeypatch):
    required = required_install_bytes(1_850_000_000, 4_321_000_017)
    monkeypatch.setattr(shutil, "disk_usage", lambda path: DiskUsage(10**10, 0, required - 1))

    with pytest.raises(InsufficientSpace) as caught:
        ensure_free_space(tmp_path, required)

    assert caught.value.required == 6_479_550_018


def test_required_install_bytes():
    assert required_install_bytes(100, 0) == 105
    assert required_install_bytes(60, 40) == 105
    # Rounded up, never down.
    assert required_install_bytes(1, 0) == 2
    assert required_install_bytes(0, 0) == 0
    # 6_171_000_017 * 1.05 = 6_479_550_017.85
    assert required_install_bytes(1_850_000_000, 4_321_000_017) == 6_479_550_018


def test_errors_share_a_base_class():
    for error in (DownloadError, IntegrityError, InsecureURL, InsecureRedirect, InsufficientSpace):
        assert issubclass(error, DatabaseError)
    assert db_fetch.DEFAULT_TIMEOUT == 60
    assert db_fetch.DEFAULT_RETRIES == 5


def test_416_reporting_another_total_raises_at_once(server, tmp_path, sleeps):
    # The server's file is shorter than the pin and than the .part.
    server.files["bundle.tar.gz"] = DATA[:1000]
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:2000])

    with pytest.raises(IntegrityError, match="1000"):
        _fetch(server, dest, sleeps)

    assert [(r.range, r.status) for r in server.requests] == [("bytes=2000-", 416)]
    assert not _part(dest).exists()
    assert sleeps == []


def test_resume_refused_when_server_reports_another_total(server, tmp_path, sleeps):
    server.files["bundle.tar.gz"] = DATA + b"x"
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(DATA[:5000])

    with pytest.raises(IntegrityError, match=str(len(DATA) + 1)):
        _fetch(server, dest, sleeps)

    assert [r.status for r in server.requests] == [206]
    assert not _part(dest).exists()


def test_untrusted_certificate_fails_at_once_and_names_the_fix(tmp_path, monkeypatch):
    import ssl
    import urllib.error

    from phu import db_fetch

    calls = []

    class Untrusted:
        def open(self, request, timeout=None):
            calls.append(request)
            raise urllib.error.URLError(ssl.SSLCertVerificationError(1, "certificate verify failed"))

    monkeypatch.setattr(db_fetch, "_opener", lambda url: Untrusted())

    with pytest.raises(db_fetch.DownloadError, match="SSL_CERT_FILE"):
        db_fetch.download_verified(
            "https://example.org/a", tmp_path / "a", "0" * 64, 10, sleep=lambda seconds: None
        )
    assert len(calls) == 1


def test_progress_reports_the_bytes_on_disk(server, tmp_path, sleeps):
    seen = []

    _fetch(server, tmp_path / "bundle.tar.gz", sleeps, progress=seen.append)

    assert seen[0] == 0
    assert seen[-1] == len(DATA)
    assert seen == sorted(seen)


def test_progress_starts_from_a_resumed_part(server, tmp_path, sleeps):
    dest = tmp_path / "bundle.tar.gz"
    _part(dest).write_bytes(DATA[:1000])
    seen = []

    _fetch(server, dest, sleeps, progress=seen.append)

    assert seen[0] == 1000
    assert seen[-1] == len(DATA)


def test_progress_returns_to_zero_when_a_download_starts_over(server, tmp_path, sleeps):
    # The resumed tail completes the file, the checksum fails, and it is fetched whole.
    dest = tmp_path / "phu-db.tar.gz"
    _part(dest).write_bytes(b"\0" * 300_000)
    seen = []

    _fetch(server, dest, sleeps, progress=seen.append)

    assert seen[0] == 300_000
    full = seen.index(len(DATA))
    # Zero is reported right after the first pass reaches 100%, before any new chunk.
    assert seen[full + 1] == 0
    assert seen[-1] == len(DATA)
    assert seen[full + 1 :] == sorted(seen[full + 1 :])
