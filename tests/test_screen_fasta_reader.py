import gzip

import phu.screen as screen
from phu.screen import _read_fasta


def test_read_fasta_plain(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">contig1 some description\nATGC\nATGC\n>contig2\nTTAA\n")

    records = list(_read_fasta(fasta))
    assert records == [("contig1", "ATGCATGC"), ("contig2", "TTAA")]


def test_read_fasta_gz(tmp_path):
    fasta_gz = tmp_path / "contigs.fa.gz"
    with gzip.open(fasta_gz, "wt") as out:
        out.write(">contigA\nAAAA\n>contigB comment\nCCCC\nGGGG\n")

    records = list(_read_fasta(fasta_gz))
    assert records == [("contigA", "AAAA"), ("contigB", "CCCCGGGG")]


def test_read_fasta_fallback_when_easel_fails(tmp_path, monkeypatch):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">contig1\nATGC\n")

    def _boom(_):
        raise RuntimeError("simulated easel failure")

    monkeypatch.setattr(screen, "_read_fasta_easel", _boom)

    records = list(_read_fasta(fasta))
    assert records == [("contig1", "ATGC")]


def test_read_fasta_fallback_after_partial_python_stream(tmp_path, monkeypatch):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">contig1\nATGC\n")

    def _python(_):
        yield ("contig1", "ATGC")
        raise ValueError("simulated parser failure")

    def _easel(_):
        yield ("contig2", "TTAA")

    monkeypatch.setattr(screen, "_read_fasta_python", _python)
    monkeypatch.setattr(screen, "_read_fasta_easel", _easel)

    records = list(_read_fasta(fasta))
    assert records == [("contig1", "ATGC"), ("contig2", "TTAA")]
