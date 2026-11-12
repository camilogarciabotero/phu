"""Verified, resumable downloads for the phu-db bundle and its upstream sources.

download_verified never trusts what it receives: every byte is hashed before
the file leaves ``<dest>.part``, and the file reaches ``dest`` only after its
size and sha256 match the pins the caller passes in. Digests always come from
the caller (the vendored manifest or the bundle pins), never from the server.

A large file is fetched as several byte ranges at once, because hosts such as
Zenodo cap the speed of each connection; a small one, or one whose server
ignores Range requests, is fetched over a single connection.
"""

from __future__ import annotations

import hashlib
import http.client
import json
import os
import re
import shutil
import ssl
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from pathlib import Path
from urllib.parse import urlsplit

CHUNK_SIZE = 1024 * 1024
# Bytes asked of the network at a time: small enough that a slow link still
# reports its progress several times a second.
NETWORK_CHUNK = 128 * 1024
DEFAULT_TIMEOUT = 60
# Ranges of one file fetched at once, and the smallest file worth splitting.
DEFAULT_CONNECTIONS = 8
MAX_CONNECTIONS = 16
PARALLEL_MIN_BYTES = 64 * 1024 * 1024
# One range request covers this much; the pieces are shared among the connections.
SEGMENT_BYTES = 64 * 1024 * 1024
# Beside <dest>.part during a ranged download: how much of each piece is on disk.
SEGMENTS_SUFFIX = ".segments"
SEGMENTS_SAVE_SECONDS = 5.0
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


class _NoRanges(Exception):
    """The server answered a Range request with the whole file."""


def _open(
    opener: urllib.request.OpenerDirector,
    url: str,
    byte_range: str | None,
    timeout: float,
    on_http_error: Callable[[urllib.error.HTTPError], None] | None = None,
) -> http.client.HTTPResponse:
    """Open ``url`` (``byte_range`` like "bytes=10-19"), sorting failures.

    Raises _Transient when another attempt may succeed and DownloadError when
    none can. ``on_http_error`` sees an HTTP error first and may raise its own.
    """
    request = urllib.request.Request(
        url, headers={"User-Agent": "phu", "Accept-Encoding": "identity"}
    )
    if byte_range:
        request.add_header("Range", byte_range)
    try:
        return opener.open(request, timeout=timeout)
    except urllib.error.HTTPError as error:
        error.close()
        if on_http_error is not None:
            on_http_error(error)
        if error.code in _RETRY_STATUSES:
            raise _Transient(f"HTTP {error.code}") from error
        raise DownloadError(f"{url}: HTTP {error.code} {error.reason}") from error
    except urllib.error.URLError as error:
        if isinstance(error.reason, ssl.SSLCertVerificationError):
            # A Python without a CA bundle (uv-managed ones on some systems)
            # fails here for every host; retrying cannot help.
            raise DownloadError(
                f"{url}: {error.reason}. If this Python has no CA bundle, "
                "set SSL_CERT_FILE to the system one, e.g. "
                "/etc/pki/tls/certs/ca-bundle.crt or /etc/ssl/certs/ca-certificates.crt"
            ) from error
        raise _Transient(str(error.reason)) from error
    except http.client.InvalidURL as error:
        # Malformed (for example a non-numeric port): no retry can fix it.
        raise DownloadError(f"{url}: {error}") from error
    except (OSError, http.client.HTTPException) as error:
        raise _Transient(str(error) or type(error).__name__) from error


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
        # The file is empty again: say so now, not when the next chunk arrives.
        if self.progress is not None:
            self.progress(0)

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

        def unsatisfied(error: urllib.error.HTTPError) -> None:
            if error.code != 416 or not self.offset:
                return
            total = _UNSATISFIED_RANGE_RE.fullmatch(error.headers.get("Content-Range", ""))
            if total and int(total.group(1)) != self.size:
                raise self.fail_integrity(
                    f"server reports {total.group(1)} bytes, expected {self.size}"
                ) from error
            # The server will not resume from here; start over.
            offset = self.offset
            self.restart()
            raise _Transient(f"HTTP 416 resuming at byte {offset}") from error

        byte_range = f"bytes={self.offset}-" if self.offset else None
        response = _open(self.opener, self.url, byte_range, self.timeout, unsatisfied)
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
                    chunk = response.read(NETWORK_CHUNK)
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


class _RangedDownload:
    """One download_verified call that fetches ``part`` as several ranges at once.

    ``part`` is created at its full size and filled in place. ``segments``
    holds ``[start, end, done]`` for each piece, ``done`` being the bytes of it
    on disk; the list is saved beside ``part`` so a later call continues every
    piece where it stopped.
    """

    def __init__(
        self,
        url: str,
        part: Path,
        sha256: str,
        size: int,
        *,
        timeout: float,
        retries: int,
        backoff: float,
        sleep: Callable[[float], None],
        progress: Callable[[int], None] | None,
        connections: int,
        segment_bytes: int,
    ) -> None:
        self.url = url
        self.part = part
        self.state_path = part.with_name(part.name + SEGMENTS_SUFFIX)
        self.sha256 = sha256
        self.size = size
        self.timeout = timeout
        self.retries = retries
        self.backoff = backoff
        self.sleep = sleep
        self.progress = progress
        self.connections = connections
        self.segment_bytes = max(1, segment_bytes)
        self.segments: list[list[int]] = []
        # Whether bytes from an earlier call are part of this download.
        self.resumed = False
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._saved = time.monotonic()

    def done(self) -> int:
        return sum(segment[2] for segment in self.segments)

    def discard(self) -> None:
        self.part.unlink(missing_ok=True)
        self.state_path.unlink(missing_ok=True)

    def fail_integrity(self, message: str) -> IntegrityError:
        self.discard()
        return IntegrityError(f"{self.url}: {message}")

    def _saved_segments(self) -> list[list[int]] | None:
        """The segment list of an earlier call, if it describes this file and ``part``."""
        try:
            state = json.loads(self.state_path.read_bytes())
            segments = state["segments"]
            if state["size"] != self.size or state["sha256"] != self.sha256:
                return None
            position = 0
            for start, end, done in segments:
                ints = all(type(value) is int for value in (start, end, done))
                if not ints or start != position or end <= start or not 0 <= done <= end - start:
                    return None
                position = end
            if position != self.size or self.part.stat().st_size != self.size:
                return None
        except (OSError, ValueError, KeyError, TypeError):
            return None
        return [list(segment) for segment in segments]

    def _save(self) -> None:
        state = {"size": self.size, "sha256": self.sha256, "segments": self.segments}
        temporary = self.state_path.with_name(self.state_path.name + ".tmp")
        temporary.write_text(json.dumps(state))
        os.replace(temporary, self.state_path)
        self._saved = time.monotonic()

    def prepare(self) -> None:
        """Continue an earlier call's ``part`` or lay out a new one."""
        saved = self._saved_segments()
        if saved is not None:
            self.segments = saved
            self.resumed = self.done() > 0
            return
        # A .part with no segment list was written front to back by a single
        # connection: keep it as a finished first piece.
        prefix = 0
        if self.part.is_file() and not self.state_path.exists():
            prefix = self.part.stat().st_size
        if prefix > self.size or self.state_path.exists():
            prefix = 0
        if not prefix:
            self.discard()
        self.segments = [[0, prefix, prefix]] if prefix else []
        for start in range(prefix, self.size, self.segment_bytes):
            self.segments.append([start, min(start + self.segment_bytes, self.size), 0])
        self.resumed = prefix > 0
        ensure_free_space(self.part.parent, self.size - prefix)
        self._save()
        with self.part.open("r+b" if prefix else "wb") as handle:
            handle.truncate(self.size)

    def run(self) -> None:
        """Fetch every unfinished piece; raises what stopped the first failing one."""
        pending = [i for i, (start, end, done) in enumerate(self.segments) if done < end - start]
        errors: list[BaseException] = []
        self._stop.clear()

        def work() -> None:
            while not self._stop.is_set():
                with self._lock:
                    if not pending:
                        return
                    index = pending.pop(0)
                try:
                    self._fetch(index)
                except BaseException as error:
                    with self._lock:
                        errors.append(error)
                    self._stop.set()
                    return

        workers = [
            threading.Thread(target=work, daemon=True)
            for _ in range(min(self.connections, len(pending)))
        ]
        try:
            for worker in workers:
                worker.start()
            while any(worker.is_alive() for worker in workers):
                for worker in workers:
                    worker.join(0.1)
        finally:
            # Also reached on Ctrl-C: stop the workers and keep what they fetched.
            self._stop.set()
            for worker in workers:
                worker.join(self.timeout + 1)
            if self.part.is_file():
                with self._lock:
                    self._save()
        if errors:
            # A server without ranges, then a wrong file, outrank a network failure.
            order = (_NoRanges, IntegrityError)
            errors.sort(key=lambda e: next((i for i, t in enumerate(order) if isinstance(e, t)), 9))
            raise errors[0]
        if self.done() != self.size:
            raise DownloadError(
                f"{self.url}: stopped with {self.done()} of {self.size} bytes fetched"
            )
        with self.part.open("rb") as handle:
            os.fsync(handle.fileno())

    def _fetch(self, index: int) -> None:
        failures = 0
        while not self._stop.is_set():
            try:
                self._attempt(index)
                return
            except _Transient as error:
                if failures >= self.retries:
                    raise DownloadError(
                        f"{self.url}: giving up after {failures + 1} attempts: {error}"
                    ) from error
                self.sleep(self.backoff * 2**failures)
                failures += 1

    def _attempt(self, index: int) -> None:
        start, end, done = self.segments[index]
        position = start + done
        if position == end:
            return

        def unsatisfied(error: urllib.error.HTTPError) -> None:
            if error.code != 416:
                return
            total = _UNSATISFIED_RANGE_RE.fullmatch(error.headers.get("Content-Range", ""))
            if total and int(total.group(1)) != self.size:
                raise self.fail_integrity(
                    f"server reports {total.group(1)} bytes, expected {self.size}"
                ) from error

        response = _open(
            _opener(self.url), self.url, f"bytes={position}-{end - 1}", self.timeout, unsatisfied
        )
        with response:
            if response.status == 200:
                raise _NoRanges
            if response.status != 206:
                raise DownloadError(f"{self.url}: unexpected HTTP {response.status}")
            match = _CONTENT_RANGE_RE.fullmatch(response.headers.get("Content-Range", ""))
            if not match or int(match.group(1)) != position:
                raise _Transient("server sent a range other than the one requested")
            if match.group(3) != "*" and int(match.group(3)) != self.size:
                raise self.fail_integrity(
                    f"server reports {match.group(3)} bytes, expected {self.size}"
                )
            with self.part.open("r+b") as handle:
                handle.seek(position)
                while position < end and not self._stop.is_set():
                    try:
                        chunk = response.read(min(NETWORK_CHUNK, end - position))
                    except (OSError, http.client.HTTPException) as error:
                        raise _Transient(str(error) or type(error).__name__) from error
                    if not chunk:
                        raise _Transient(
                            f"connection closed at byte {position} of piece {start}-{end - 1}"
                        )
                    handle.write(chunk)
                    handle.flush()
                    position += len(chunk)
                    with self._lock:
                        self.segments[index][2] = position - start
                        if time.monotonic() - self._saved >= SEGMENTS_SAVE_SECONDS:
                            self._save()
                        if self.progress is not None:
                            self.progress(self.done())


def _download_ranged(download: _RangedDownload, dest: Path) -> Path:
    """Run ``download`` to a verified ``dest``; raises _NoRanges to ask for one connection."""
    download.prepare()
    if download.progress is not None:
        download.progress(download.done())
    while True:
        download.run()
        digest = _hash_file(download.part).hexdigest()
        if digest == download.sha256:
            break
        if download.resumed:
            # What an earlier call left was not part of this file (another
            # bundle version, or torn by a crash); fetch it whole, once.
            download.discard()
            download.prepare()
            if download.progress is not None:
                download.progress(0)
            continue
        raise download.fail_integrity(
            f"sha256 {digest} does not match the pinned {download.sha256}"
        )
    download.state_path.unlink(missing_ok=True)
    os.replace(download.part, dest)
    return dest


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
    connections: int = DEFAULT_CONNECTIONS,
    parallel_min_bytes: int = PARALLEL_MIN_BYTES,
    segment_bytes: int = SEGMENT_BYTES,
) -> Path:
    """Download ``url`` to ``dest`` and return ``dest`` once size and sha256 match.

    A file of at least ``parallel_min_bytes`` is fetched as pieces of
    ``segment_bytes``, up to ``connections`` of them at once, each with its own
    retries; ``<dest>.part.segments`` records how far each piece got, so a
    later call continues all of them. If the server ignores Range requests,
    the file is fetched from zero over one connection, as smaller files are.

    ``progress``, if given, is called with the bytes on disk so far: once after
    an existing ``.part`` is taken over, then after every chunk. The count is
    reported as zero as soon as a download starts over.

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
    segments = part.with_name(part.name + SEGMENTS_SUFFIX)
    connections = max(1, min(int(connections), MAX_CONNECTIONS))
    if connections > 1 and size >= parallel_min_bytes:
        ranged = _RangedDownload(
            url,
            part,
            sha256,
            size,
            timeout=timeout,
            retries=retries,
            backoff=backoff,
            sleep=sleep,
            progress=progress,
            connections=connections,
            segment_bytes=segment_bytes,
        )
        try:
            return _download_ranged(ranged, dest)
        except _NoRanges:
            ranged.discard()
    elif segments.exists():
        # Left by a ranged download: the .part has holes a single connection cannot fill.
        part.unlink(missing_ok=True)
        segments.unlink()
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
