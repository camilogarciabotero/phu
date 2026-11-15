from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from phu.core.fasta import extract_records, extract_records_multi

FASTA = (
    ">c1 first contig\n"
    "ACGTACGTAC\n"
    "GTACGT\n"
    ">c2\n"
    "TTTT\n"
    ">c3\tthird\n"
    "GGGGCCCC\n"
    "GG\n"
    ">c1|gene1 x\n"
    "MPEP\n"
    ">c2\n"
    "AAAA\n"
)


def _records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line, ""))
        else:
            header, seq = records[-1]
            records[-1] = (header, seq + line)
    return records


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "in.fa"
    path.write_text(FASTA)
    return path


def test_matches_id_up_to_whitespace_and_keeps_records_verbatim(
    fasta: Path, tmp_path: Path
) -> None:
    out = tmp_path / "out.fa"
    n = extract_records(fasta, ["c3", "c1"], out)
    # Input order, not request order; wrapping and descriptions unchanged.
    assert out.read_text() == (
        ">c1 first contig\nACGTACGTAC\nGTACGT\n>c3\tthird\nGGGGCCCC\nGG\n"
    )
    assert n == 2


def test_writes_every_record_with_a_repeated_id(fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.fa"
    assert extract_records(fasta, ["c2", "c2"], out) == 2
    assert out.read_text() == ">c2\nTTTT\n>c2\nAAAA\n"


def test_ids_with_pipes_match_whole_id_only(fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.fa"
    extract_records(fasta, ["c1|gene1"], out)
    assert out.read_text() == ">c1|gene1 x\nMPEP\n"


def test_no_ids_or_no_match_gives_empty_file(fasta: Path, tmp_path: Path) -> None:
    for ids in ([], ["missing"]):
        out = tmp_path / "out.fa"
        assert extract_records(fasta, ids, out) == 0
        assert out.exists() and out.read_text() == ""


def test_reads_gzipped_input(fasta: Path, tmp_path: Path) -> None:
    gz = tmp_path / "in.fa.gz"
    gz.write_bytes(gzip.compress(fasta.read_bytes()))
    out = tmp_path / "out.fa"
    assert extract_records(gz, ["c2"], out) == 2


def test_multi_reads_once_and_fills_every_output(fasta: Path, tmp_path: Path) -> None:
    a, b, empty = tmp_path / "a.fa", tmp_path / "b.fa", tmp_path / "e.fa"
    counts = extract_records_multi(fasta, {a: ["c1", "c2"], b: ["c2"], empty: []})
    assert counts == {a: 3, b: 2, empty: 0}
    assert b.read_text() == ">c2\nTTTT\n>c2\nAAAA\n"
    assert empty.read_text() == ""


@pytest.mark.skipif(shutil.which("seqkit") is None, reason="seqkit not installed")
def test_records_match_seqkit_grep(fasta: Path, tmp_path: Path) -> None:
    ids = ["c3", "c1", "c2", "c1|gene1", "missing"]
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("\n".join(ids) + "\n")
    seqkit_out = tmp_path / "seqkit.fa"
    subprocess.run(
        ["seqkit", "grep", "-f", str(ids_file), str(fasta), "-o", str(seqkit_out)],
        capture_output=True,
        check=True,
    )
    out = tmp_path / "out.fa"
    extract_records(fasta, ids, out)
    # Same records, headers and sequences in the same order; seqkit rewraps
    # sequence lines at 60 characters, phu keeps the input's wrapping.
    assert _records(out) == _records(seqkit_out)
