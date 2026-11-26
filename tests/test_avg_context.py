import json
import random
import shutil
from pathlib import Path

import pytest

from phu.avger import context as ctx
from phu.core.gene_prediction import PredictedGene

CODE = {
    a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
        16 * i + 4 * j + k
    ]
    for i, a in enumerate("TCAG")
    for j, b in enumerate("TCAG")
    for k, c in enumerate("TCAG")
}
SENSE = [codon for codon, aa in CODE.items() if aa != "*"]


def _gene(rng, codons=150):
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(codons)) + "TAA"


def _contigs(tmp_path: Path):
    """Two contigs sharing one 5-gene block, a third with an unrelated block.

    The anchor is the middle gene of each block; contig 'a' starts at its
    block so a flank can be short on one side.
    """
    rng = random.Random(3)
    shared = [_gene(rng) for _ in range(5)]
    other = [_gene(rng) for _ in range(5)]
    spacer = lambda: "".join(rng.choice("ACGT") for _ in range(60))  # noqa: E731
    layout = {"a": shared, "b": [_gene(rng)] + shared, "c": other}
    genes, seqs = [], {}
    for contig, block in layout.items():
        seq = ""
        for i, nt in enumerate(block, start=1):
            seq += spacer()
            start = len(seq) + 1
            seq += nt
            aa = "".join(CODE[nt[p : p + 3]] for p in range(0, len(nt), 3))
            genes.append(
                PredictedGene(contig, f"{contig}|{i}", start, len(seq), 1, i, nt, aa)
            )
        seqs[contig] = seq + spacer()
    fasta = tmp_path / "contigs.fna"
    fasta.write_text("".join(f">{k} desc\n{v}\n" for k, v in seqs.items()))
    return fasta, genes


def test_shared_neighbourhoods_share_a_family_and_flanks_are_summarised(tmp_path):
    fasta, genes = _contigs(tmp_path)
    anchors = {"a|3": "putative_amg", "b|4": "filtered", "c|3": "putative_areg"}
    scores = {"a|2": {"pfam": (2.0, 1.0)}, "a|1": {"pfam": (4.0, 3.0)}}
    work = tmp_path / "work"
    work.mkdir()
    result, records = ctx.build_context(
        contigs_path=fasta,
        genes=genes,
        anchors=anchors,
        scores=scores,
        flank=2,
        threads=1,
        workdir=work,
    )
    rows = {r["protein_id"]: r for r in result.neighborhoods}
    assert rows["a|3"]["gcf_id"] == rows["b|4"]["gcf_id"] != rows["c|3"]["gcf_id"]
    assert rows["a|3"]["gcf_windows"] == 2 and rows["a|3"]["gcf_contigs"] == 2
    assert rows["a|3"]["gcf_predictions"] == 1
    assert (rows["a|3"]["genes_left"], rows["a|3"]["genes_right"]) == (2, 2)
    assert len(records) == 3 and result.gcfs == 2
    assert result.versions.keys() == {"igua", "pymmseqs", "gb-io"}

    flank = {(r["protein_id"], r["side"], r["database"]): r for r in result.flanks}
    left = flank[("a|3", "left", "pfam")]
    assert (left["genes"], left["scored_genes"], left["max_v"], left["mean_vl"]) == (
        2,
        2,
        4.0,
        2.0,
    )
    right = flank[("a|3", "right", "kofam")]
    assert (right["genes"], right["scored_genes"], right["max_v"]) == (2, 0, "")


def test_a_window_at_a_contig_end_is_short_on_that_side(tmp_path):
    fasta, genes = _contigs(tmp_path)
    work = tmp_path / "work"
    work.mkdir()
    result, records = ctx.build_context(
        contigs_path=fasta,
        genes=genes,
        anchors={"a|1": "filtered"},
        scores={},
        flank=3,
        threads=1,
        workdir=work,
    )
    (row,) = result.neighborhoods
    assert (row["genes_left"], row["genes_right"]) == (0, 3)
    assert row["window_start"] == genes[0].start
    assert len(records[0].features) == 4
    assert row["gcf_id"]


def test_context_failure_warns_and_keeps_the_avg_outputs(monkeypatch, tmp_path, capsys):
    from test_avg import _phrog_outputs
    from phu.commands import avger as avg

    def broken(**kwargs):
        raise RuntimeError("no mmseqs")

    monkeypatch.setattr(avg, "build_context", broken)
    _, out = _phrog_outputs(monkeypatch, tmp_path, "phrog_2")
    assert "neighbourhood context skipped" in capsys.readouterr().err
    assert len((out / "avg-predictions.tsv").read_text().splitlines()) == 2
    assert not (out / "detailed" / "context-neighborhoods.tsv").exists()
    run = json.loads((out / ".phu" / "run.json").read_text())
    assert run["context"]["status"] == "failed"
    assert "no mmseqs" in run["context"]["error"]


def test_context_flank_must_be_positive(tmp_path):
    from phu.commands.avger import AvgConfig

    contigs = tmp_path / "c.fna"
    contigs.write_text(">x\nACGT\n")
    with pytest.raises(ValueError, match="context_flank"):
        AvgConfig(input_contigs=contigs, context_flank=0)


def test_a_family_shares_genes_not_necessarily_their_order(tmp_path):
    """IGUA groups windows by shared protein content; gene order is not checked.

    Two windows with the same five genes, the second in reversed order, land in
    one family. Documented so that nobody reads a family as conserved synteny.
    """
    fasta, genes = _contigs(tmp_path)
    a = [g for g in genes if g.contig_id == "a"]
    shuffled, seq, seqs = [], "", {}
    rng = random.Random(5)
    for i, g in enumerate(reversed(a), start=1):
        seq += "".join(rng.choice("ACGT") for _ in range(60))
        start = len(seq) + 1
        seq += g.nucleotide_sequence
        shuffled.append(
            PredictedGene(
                "r",
                f"r|{i}",
                start,
                len(seq),
                1,
                i,
                g.nucleotide_sequence,
                g.amino_acid_sequence,
            )
        )
    fasta.write_text(fasta.read_text() + f">r\n{seq}\n")
    work = tmp_path / "work"
    work.mkdir()
    result, _ = ctx.build_context(
        contigs_path=fasta,
        genes=genes + shuffled,
        anchors={"a|3": "filtered", "r|3": "filtered"},
        scores={},
        flank=2,
        threads=1,
        workdir=work,
    )
    rows = {r["protein_id"]: r["gcf_id"] for r in result.neighborhoods}
    assert rows["a|3"] == rows["r|3"]


CONTEXT_FILES = (
    "context-neighborhoods.tsv",
    "context-flanks.tsv",
    "context-windows.gbk",
)


def test_no_context_removes_context_files_left_by_an_earlier_run(monkeypatch, tmp_path):
    from test_avg import _phrog_outputs

    _, out = _phrog_outputs(monkeypatch, tmp_path, "phrog_2")
    assert all((out / "detailed" / name).exists() for name in CONTEXT_FILES)
    shutil.rmtree(tmp_path / "phrog")  # the helper recreates its fixture folder
    _, out = _phrog_outputs(monkeypatch, tmp_path, "phrog_2", context=False)
    assert not any((out / "detailed" / name).exists() for name in CONTEXT_FILES)
    run = json.loads((out / ".phu" / "run.json").read_text())
    assert run["context"] == {"enabled": False}
    assert not any("context" in key for key in run["output_counts"])
    header = (out / "avg-predictions.tsv").read_text().splitlines()
    assert header[1].split("\t")[header[0].split("\t").index("context_gcf")] == ""


def test_a_failed_context_export_warns_and_still_writes_the_metadata(
    monkeypatch, tmp_path, capsys
):
    from test_avg import _phrog_outputs
    from phu.commands import avger as avg

    def broken(records, path):
        path.write_text("partial")
        raise OSError("disk full")

    monkeypatch.setattr(avg, "write_genbank", broken)
    _, out = _phrog_outputs(monkeypatch, tmp_path, "phrog_2")
    assert "neighbourhood context not written" in capsys.readouterr().err
    assert not any((out / "detailed" / name).exists() for name in CONTEXT_FILES)
    run = json.loads((out / ".phu" / "run.json").read_text())
    assert (
        run["context"]["status"] == "failed" and "disk full" in run["context"]["error"]
    )
    assert (out / "output-info.txt").exists()
    rows = (out / "avg-predictions.tsv").read_text().splitlines()
    assert rows[1].split("\t")[rows[0].split("\t").index("context_gcf")] == ""
