"""Verified, resumable downloads for the phu-db bundle and its upstream sources.

download_verified never trusts what it receives: every byte is hashed as it is
written to ``<dest>.part``, and the file reaches ``dest`` only after its size
and sha256 match the pins the caller passes in. Digests always come from the
caller (the vendored manifest or the bundle pins), never from the server.
"""

from __future__ import annotations

import hashlib
import http.client
import os
import re
import shutil
import ssl
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from pathlib import Path
from urllib.parse import urlsplit

CHUNK_SIZE = 1024 * 1024
DEFAULT_TIMEOUT = 60
DEFAULT_RETRIES = 5
DEFAULT_BACKOFF = 1.0
# Headroom on top of the tarball and the prepared tree it unpacks to.
SPACE_MARGIN_PERCENT = 5
# The only hosts plain HTTP may reach (urlsplit drops the brackets of [::1]).
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "::1", "localhost"})

# Status codes worth another attempt; any other HTTP error is final.
_RETRY_STATUSES = frozenset({408, 429, 500, 502, 503, 504})
_CONTENT_RANGE_RE = re.compile(r"bytes (\d+)-(\d+)/(\d+|\*)")
# Content-Range of a 416 reply: the server's full size.
_UNSATISFIED_RANGE_RE = re.compile(r"bytes \*/(\d+)")


class DatabaseError(RuntimeError):
    """Base class for every phu-db install failure."""


class DownloadError(DatabaseError):
    """The download could not be completed."""


class IntegrityError(DatabaseError):
    """The downloaded bytes do not match the pinned size or sha256."""


class InsecureRedirect(DatabaseError):
    """A URL or redirect target is not HTTPS, and not plain HTTP on loopback."""


class InsecureURL(InsecureRedirect):
    """The URL a download starts from is itself insecure."""


class InsufficientSpace(DatabaseError):
    """Not enough free disk space for the install."""

    def __init__(self, required: int, available: int, path: Path | None = None) -> None:
        self.required = required
        self.available = available
        self.path = path
        where = f" at {path}" if path is not None else ""
        super().__init__(
            f"not enough free disk space{where}: {required} bytes required, "
            f"{available} bytes available"
        )


class _Transient(Exception):
    """One attempt failed in a way that a later attempt may not."""


def _is_loopback(url: str) -> bool:
    return (urlsplit(url).hostname or "").lower() in LOOPBACK_HOSTS


def _check_url(url: str) -> None:
    """Allow HTTPS anywhere, and plain HTTP only to this machine."""
    scheme = urlsplit(url).scheme.lower()
    if scheme == "https" or (scheme == "http" and _is_loopback(url)):
        return
    raise InsecureURL(f"refusing to download over {scheme or 'no scheme'}: {url}")


def redirect_allowed(source: str, target: str) -> bool:
    """A redirect must land on HTTPS, unless it stays on the loopback interface."""
    if urlsplit(target).scheme.lower() == "https":
        return True
    return (
        urlsplit(target).scheme.lower() == "http"
        and _is_loopback(source)
        and _is_loopback(target)
    )


class _GuardedRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Refuse an insecure redirect before any connection to its target."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not redirect_allowed(req.full_url, newurl):
            fp.close()
            raise InsecureRedirect(f"refusing redirect from {req.full_url} to {newurl}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _opener(url: str) -> urllib.request.OpenerDirector:
    # Loopback traffic never goes through a proxy configured for the outside.
    proxies = {} if _is_loopback(url) else None
    return urllib.request.build_opener(
        urllib.request.ProxyHandler(proxies), _GuardedRedirectHandler()
    )


def _hash_file(path: Path) -> hashlib._Hash:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            digest.update(chunk)
    return digest


def required_install_bytes(tarball_bytes: int, prepared_bytes: int) -> int:
    """Disk space an install needs: tarball plus prepared tree, plus 5%."""
    # Integer ceiling division: floats lose bytes at bundle sizes.
    total = tarball_bytes + prepared_bytes
    return -(-total * (100 + SPACE_MARGIN_PERCENT) // 100)


def ensure_free_space(path: Path, required_bytes: int) -> None:
    """Raise InsufficientSpace unless ``path``'s filesystem has ``required_bytes`` free."""
    probe = Path(path)
    # The install directory may not exist yet; measure the nearest ancestor.
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    available = shutil.disk_usage(probe).free
    if available < required_bytes:
        raise InsufficientSpace(required_bytes, available, Path(path))


class _Download:
    """State of one download_verified call, shared across its attempts."""

    def __init__(
        self,
        url: str,
        part: Path,
        sha256: str,
        size: int,
        timeout: float,
        progress: Callable[[int], None] | None = None,
    ) -> None:
        self.progress = progress
        self.url = url
        self.part = part
        self.sha256 = sha256
        self.size = size
        self.timeout = timeout
        self.opener = _opener(url)
        self.offset = 0
        self.digest = hashlib.sha256()
        # Bytes taken over from a .part left by an earlier call.
        self.resumed = 0

    def restart(self) -> None:
        self.part.write_bytes(b"")
        self.offset = 0
        self.digest = hashlib.sha256()
        self.resumed = 0

    def resume_existing(self) -> None:
        """Pick up a .part left by an earlier call, re-hashing what it holds.

        A .part that already holds ``size`` bytes needs no request at all; it is
        verified like a finished download.
        """
        if not self.part.is_file():
            return
        if self.part.stat().st_size > self.size:
            # Longer than the pinned size: it cannot be a prefix of the file.
            self.part.unlink()
            return
        self.digest = _hash_file(self.part)
        self.offset = self.resumed = self.part.stat().st_size

    def fail_integrity(self, message: str) -> IntegrityError:
        self.part.unlink(missing_ok=True)
        return IntegrityError(f"{self.url}: {message}")

    def attempt(self) -> None:
        """Fetch the remaining bytes once; raises _Transient when a retry may help."""
        if self.offset == self.size:
            return
        request = urllib.request.Request(
            self.url, headers={"User-Agent": "phu", "Accept-Encoding": "identity"}
        )
        if self.offset:
            request.add_header("Range", f"bytes={self.offset}-")
        try:
            response = self.opener.open(request, timeout=self.timeout)
        except urllib.error.HTTPError as error:
            error.close()
            if error.code == 416 and self.offset:
                total = _UNSATISFIED_RANGE_RE.fullmatch(error.headers.get("Content-Range", ""))
                if total and int(total.group(1)) != self.size:
                    raise self.fail_integrity(
                        f"server reports {total.group(1)} bytes, expected {self.size}"
                    ) from error
                # The server will not resume from here; start over.
                offset = self.offset
                self.restart()
                raise _Transient(f"HTTP 416 resuming at byte {offset}") from error
            if error.code in _RETRY_STATUSES:
                raise _Transient(f"HTTP {error.code}") from error
            raise DownloadError(f"{self.url}: HTTP {error.code} {error.reason}") from error
        except urllib.error.URLError as error:
            if isinstance(error.reason, ssl.SSLCertVerificationError):
                # A Python without a CA bundle (uv-managed ones on some systems)
                # fails here for every host; retrying cannot help.
                raise DownloadError(
                    f"{self.url}: {error.reason}. If this Python has no CA bundle, "
                    "set SSL_CERT_FILE to the system one, e.g. "
                    "/etc/pki/tls/certs/ca-bundle.crt or /etc/ssl/certs/ca-certificates.crt"
                ) from error
            raise _Transient(str(error.reason)) from error
        except http.client.InvalidURL as error:
            # Malformed (for example a non-numeric port): no retry can fix it.
            raise DownloadError(f"{self.url}: {error}") from error
        except (OSError, http.client.HTTPException) as error:
            raise _Transient(str(error) or type(error).__name__) from error
        with response:
            self._check_response(response)
            self._stream(response)

    def _check_response(self, response: http.client.HTTPResponse) -> None:
        if response.status == 206 and self.offset:
            match = _CONTENT_RANGE_RE.fullmatch(response.headers.get("Content-Range", ""))
            if not match or int(match.group(1)) != self.offset:
                self.restart()
                raise _Transient("server sent a range other than the one requested")
            if match.group(3) != "*" and int(match.group(3)) != self.size:
                raise self.fail_integrity(
                    f"server reports {match.group(3)} bytes, expected {self.size}"
                )
            return
        if response.status != 200:
            raise DownloadError(f"{self.url}: unexpected HTTP {response.status}")
        length = response.headers.get("Content-Length")
        if length is not None and length.isdigit() and int(length) != self.size:
            message = f"server sends {length} bytes, expected {self.size}"
            if self.offset:
                # Possibly a proxy or portal page: keep the prefix for a later run.
                raise IntegrityError(f"{self.url}: {message}; kept {self.part}")
            raise self.fail_integrity(message)
        if self.offset:
            # A 200 to a Range request is the whole file again.
            self.restart()

    def _stream(self, response: http.client.HTTPResponse) -> None:
        with self.part.open("ab") as handle:
            while True:
                try:
                    chunk = response.read(CHUNK_SIZE)
                except (OSError, http.client.HTTPException) as error:
                    raise _Transient(str(error) or type(error).__name__) from error
                if not chunk:
                    break
                if self.offset + len(chunk) > self.size:
                    raise self.fail_integrity(f"server sent more than {self.size} bytes")
                handle.write(chunk)
                self.digest.update(chunk)
                self.offset += len(chunk)
                if self.progress is not None:
                    self.progress(self.offset)
            handle.flush()
            os.fsync(handle.fileno())
        if self.offset < self.size:
            raise _Transient(f"connection closed after {self.offset} of {self.size} bytes")


def download_verified(
    url: str,
    dest: Path,
    sha256: str,
    size: int,
    *,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    backoff: float = DEFAULT_BACKOFF,
    sleep: Callable[[float], None] = time.sleep,
    progress: Callable[[int], None] | None = None,
) -> Path:
    """Download ``url`` to ``dest`` and return ``dest`` once size and sha256 match.

    ``progress``, if given, is called with the bytes on disk so far: once after
    an existing ``.part`` is taken over, then after every chunk. The count goes
    back to zero when a download starts over.

    Bytes stream into ``<dest>.part``; an existing ``.part`` is re-hashed and
    resumed with a Range request. If a resumed file fails its sha256, the
    ``.part`` may be a prefix of some other file, so it is fetched once more
    from zero. Network failures are retried ``retries`` times in all, waiting
    ``backoff * 2**n`` seconds before retry ``n``. A mismatch deletes the
    ``.part`` and raises IntegrityError.
    """
    _check_url(url)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    download = _Download(url, part, sha256, size, timeout, progress)
    download.resume_existing()
    if progress is not None:
        progress(download.offset)
    ensure_free_space(dest.parent, size - download.offset)

    failures = 0
    while True:
        try:
            download.attempt()
        except _Transient as error:
            if failures >= retries:
                raise DownloadError(
                    f"{url}: giving up after {failures + 1} attempts: {error}"
                ) from error
            sleep(backoff * 2**failures)
            failures += 1
            continue
        if download.offset != size:
            raise download.fail_integrity(f"got {download.offset} bytes, expected {size}")
        if download.digest.hexdigest() == sha256:
            break
        if download.resumed:
            # The old .part was not a prefix of this file (another bundle
            # version, or torn by a crash); fetch it whole, once.
            download.restart()
            continue
        raise download.fail_integrity(
            f"sha256 {download.digest.hexdigest()} does not match the pinned {sha256}"
        )
    os.replace(part, dest)
    return dest
