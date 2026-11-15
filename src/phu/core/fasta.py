"""Pull records out of a FASTA file by identifier.

This replaces ``seqkit grep -f ids.txt input.fa``: a record matches when the
text after ``>`` up to the first whitespace is one of the requested ids.
Matching records are copied verbatim, in input order. The records, headers
and sequences are the same as seqkit's; only line wrapping can differ, since
seqkit rewraps sequences at 60 characters while this keeps the input's.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import IO

from phu.core.log import log_call


def _open_text(path: Path) -> IO[str]:
    """Open a plain or gzip-compressed FASTA file for reading."""
    path = Path(path)
    with path.open("rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    if gzipped:
        return gzip.open(path, "rt")
    return path.open("r")


def _record_id(header: str) -> str:
    """Return the record id: header text after '>' up to the first whitespace."""
    fields = header[1:].split(maxsplit=1)
    return fields[0] if fields else ""


def extract_records_multi(
    input_fa: Path, targets: Mapping[Path, Iterable[str]]
) -> dict[Path, int]:
    """Read ``input_fa`` once and write each record to every output asking for it.

    ``targets`` maps an output path to the ids it should receive. Every output
    file is created, even when none of its ids is found. Returns the number of
    records written to each output.
    """
    wanted: dict[str, list[Path]] = {}
    for out, ids in targets.items():
        for record_id in dict.fromkeys(ids):
            wanted.setdefault(record_id, []).append(Path(out))
    log_call(
        "phu.extract_records",
        input=input_fa,
        ids=len(wanted),
        outputs=len(targets),
    )

    counts = {Path(out): 0 for out in targets}
    handles = {out: out.open("w") for out in counts}
    try:
        current: list[IO[str]] = []
        with _open_text(input_fa) as src:
            for line in src:
                if line.startswith(">"):
                    outs = wanted.get(_record_id(line.rstrip("\r\n")), [])
                    current = [handles[out] for out in outs]
                    for out in outs:
                        counts[out] += 1
                for handle in current:
                    handle.write(line)
    finally:
        for handle in handles.values():
            handle.close()
    return counts


def extract_records(input_fa: Path, ids: Iterable[str], output_fa: Path) -> int:
    """Write the records of ``input_fa`` whose id is in ``ids`` to ``output_fa``.

    Returns the number of records written.
    """
    output_fa = Path(output_fa)
    return extract_records_multi(input_fa, {output_fa: ids})[output_fa]
