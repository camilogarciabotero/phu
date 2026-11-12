"""A misbehaving-on-request HTTP server on 127.0.0.1 for download tests.

Each test gets a fresh LoopbackServer (see the ``loopback_server`` fixture in
tests/conftest.py). Files are served from ``server.files``; the switches below
make it fail the way real servers and networks do. Every request is recorded
in ``server.requests`` so tests can assert what the client sent.
"""

from __future__ import annotations

import re
import threading
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

_RANGE_RE = re.compile(r"bytes=(\d+)-(\d*)")


@dataclass
class Request:
    path: str
    range: str | None
    status: int


@dataclass
class Switches:
    # Send only this many body bytes, then drop the connection.
    cut_after: int | None = None
    # Apply cut_after to this many requests (None: every request).
    cut_times: int | None = None
    # Answer every request with the whole file and a 200, Range or not.
    ignore_range: bool = False
    # Flip the bits of the byte at this offset of every file served.
    corrupt_offset: int | None = None
    # Wait this long before answering; tests pass a smaller client timeout.
    stall: float | None = None
    stall_times: int | None = None
    # Answer this status (with no body) to the first fail_times requests.
    fail_status: int | None = None
    fail_times: int = 1
    # Send 200/206 bodies without Content-Length, delimited by closing the connection.
    omit_length: bool = False
    # path -> Location header of a 302 reply.
    redirects: dict[str, str] = field(default_factory=dict)


class LoopbackServer:
    def __init__(self) -> None:
        self.files: dict[str, bytes] = {}
        self.switches = Switches()
        self.requests: list[Request] = []
        self._lock = threading.Lock()
        self._released = threading.Event()
        self._httpd = ThreadingHTTPServer(("127.0.0.1", 0), _handler_for(self))
        self._httpd.daemon_threads = True
        # A short poll keeps shutdown() quick; every test stops its own server.
        self._thread = threading.Thread(
            target=self._httpd.serve_forever, kwargs={"poll_interval": 0.02}, daemon=True
        )

    @property
    def base_url(self) -> str:
        host, port = self._httpd.server_address[:2]
        return f"http://{host}:{port}"

    def url(self, path: str) -> str:
        return f"{self.base_url}/{path.lstrip('/')}"

    def serve_dir(self, root: Path, prefix: str = "") -> None:
        """Serve every file under ``root`` at ``<prefix>/<relative path>``."""
        root = Path(root)
        for path in sorted(root.rglob("*")):
            if path.is_file():
                name = path.relative_to(root).as_posix()
                self.files[f"{prefix.strip('/')}/{name}".lstrip("/")] = path.read_bytes()

    def start(self) -> LoopbackServer:
        self._thread.start()
        return self

    def stop(self) -> None:
        self._released.set()  # wake any stalled handler
        self._httpd.shutdown()
        self._httpd.server_close()
        self._thread.join(timeout=5)

    def _record(self, request: Request) -> int:
        """Record a request; returns how many came before it."""
        with self._lock:
            self.requests.append(request)
            return len(self.requests) - 1


def _handler_for(server: LoopbackServer) -> type[BaseHTTPRequestHandler]:
    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.0"

        def log_message(self, format, *args):
            pass

        def do_GET(self) -> None:
            switches = server.switches
            path = self.path.split("?", 1)[0].lstrip("/")
            range_header = self.headers.get("Range")
            index = server._record(Request(path, range_header, 0))
            record = server.requests[index]

            if switches.stall is not None and (
                switches.stall_times is None or index < switches.stall_times
            ):
                if server._released.wait(switches.stall):
                    return

            if switches.fail_status is not None and index < switches.fail_times:
                self._reply_empty(record, switches.fail_status)
                return
            if path in switches.redirects:
                record.status = 302
                self.send_response(302)
                self.send_header("Location", switches.redirects[path])
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            if path not in server.files:
                self._reply_empty(record, 404)
                return

            data = bytearray(server.files[path])
            if switches.corrupt_offset is not None and switches.corrupt_offset < len(data):
                data[switches.corrupt_offset] ^= 0xFF
            start = 0
            end = len(data)
            match = _RANGE_RE.fullmatch(range_header or "")
            if match and not switches.ignore_range:
                start = int(match.group(1))
                if start >= len(data):
                    self._reply_empty(record, 416, {"Content-Range": f"bytes */{len(data)}"})
                    return
                if match.group(2):
                    end = min(int(match.group(2)) + 1, len(data))
                record.status = 206
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{end - 1}/{len(data)}")
            else:
                record.status = 200
                self.send_response(200)
            body = bytes(data[start:end])
            if not switches.omit_length:
                self.send_header("Content-Length", str(len(body)))
            self.end_headers()

            cut = switches.cut_after is not None and (
                switches.cut_times is None or index < switches.cut_times
            )
            if cut:
                body = body[: switches.cut_after]
            try:
                self.wfile.write(body)
                self.wfile.flush()
            except OSError:
                return
            if cut:
                # Drop the connection with the advertised length unmet.
                self.close_connection = True

        def _reply_empty(self, record: Request, status: int, headers=None) -> None:
            record.status = status
            self.send_response(status)
            for key, value in (headers or {}).items():
                self.send_header(key, value)
            self.send_header("Content-Length", "0")
            self.end_headers()

    return Handler
