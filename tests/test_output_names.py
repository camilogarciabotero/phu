"""Output file and folder names written by screen, jack and cluster."""

from __future__ import annotations

import json
from pathlib import Path

import pyhmmer
import pytest

from phu import cluster, jack, screen
from phu.gene_prediction_core import CacheArtifact

SEQ_A = "MKKLLAAVAGWAAAPAAACDEFGHIKLMNPQRSTVWY"
SEQ_B = "WYVTSRQPNMLKIHGFEDCWWHHYYFFMMKKRRDDEE"


def _write_hmm(path: Path, name: str, sequence: str) -> Path:
    alphabet = pyhmmer.easel.Alphabet.amino()
    msa = pyhmmer.easel.TextMSA(
        name=name.encode(),
        sequences=[
            pyhmmer.easel.TextSequence(name=f"s{i}".encode(), sequence=sequence)
            for i in range(3)
        ],
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        hmm.write(handle)
    return path


def _fasta_records(path: Path) -> list[tuple[str, str]]:
    return list(screen._read_fasta(path))


class _Completed:
    returncode = 0
    stderr = ""


def _fake_run(cmd, stdout=None, **kwargs):
    """Stand in for `seqkit grep -f ids fasta` and `mafft ... fasta`."""
    if cmd[1] == "grep":
        wanted = set(Path(cmd[3]).read_text().split())
        records = [r for r in _fasta_records(Path(cmd[4])) if r[0] in wanted]
    else:  # aligner: the fixture sequences already have equal lengths
        records = _fasta_records(Path(cmd[-1]))
    stdout.write("".join(f">{name}\n{seq}\n" for name, seq in records))
    return _Completed()


@pytest.fixture
def workspace(tmp_path: Path, monkeypatch) -> Path:
    """Contigs, predicted proteins and stubs for the external binaries."""
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c1\nACGT\n>c2\nACGT\n>c3\nACGT\n")
    proteins = tmp_path / "cache" / "proteins.faa"
    proteins.parent.mkdir()
    proteins.write_text(
        f">c1|gene1\n{SEQ_A}\n>c2|gene1\n{SEQ_A}\n>c1|gene2\n{SEQ_B}\n"
        ">c3|gene1\nMPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n"
    )

    def fake_prediction(inputs, use_cache=True, threads=1):
        return CacheArtifact(proteins, 4, False, "key", cache_dir=proteins.parent)

    for module in (screen, jack):
        monkeypatch.setattr(module, "get_or_predict_proteins", fake_prediction)
        monkeypatch.setattr(module, "_executable", lambda names: names[0])
    monkeypatch.setattr(screen.subprocess, "run", _fake_run)
    return tmp_path


def _listing(folder: Path) -> list[str]:
    return sorted(
        path.relative_to(folder).as_posix()
        for path in folder.rglob("*")
        if path.is_file()
    )


def test_screen_default_output_names(workspace: Path) -> None:
    out = workspace / "out"
    screen._screen(
        screen.ScreenConfig(
            input_contigs=workspace / "contigs.fna",
            hmms=[
                _write_hmm(workspace / "famA.hmm", "model_a", SEQ_A),
                _write_hmm(workspace / "famB.hmm", "model_b", SEQ_B),
            ],
            outdir=out,
            cut_ga=False,
        )
    )

    assert _listing(out) == [
        ".phu/query.json",
        "famA-hmmsearch-domains.tbl",
        "famB-hmmsearch-domains.tbl",
        "matched-contig-ids.txt",
        "matched-contigs.fna",
    ]
    assert (out / "matched-contig-ids.txt").read_text().split() == ["c1", "c2"]
    assert [name for name, _ in _fasta_records(out / "matched-contigs.fna")] == [
        "c1",
        "c2",
    ]


def test_screen_optional_output_names(workspace: Path) -> None:
    out = workspace / "out"
    screen._screen(
        screen.ScreenConfig(
            input_contigs=workspace / "contigs.fna",
            hmms=[
                _write_hmm(workspace / "famA.hmm", "model_a", SEQ_A),
                _write_hmm(workspace / "famB.hmm", "model_b", SEQ_B),
            ],
            outdir=out,
            cut_ga=False,
            keep_domtbl=False,
            keep_proteins=True,
            save_target_proteins=True,
            save_target_hmms=True,
        )
    )

    assert _listing(out) == [
        ".phu/prediction.json",
        ".phu/query.json",
        "matched-contig-ids.txt",
        "matched-contigs.fna",
        "predicted-proteins.faa",
        "target-hmms/famA.hmm",
        "target-hmms/famB.hmm",
        "target-proteins/famA.faa",
        "target-proteins/famB.faa",
    ]
    assert json.loads((out / ".phu" / "prediction.json").read_text())["cache_key"]
    # famA has two proteins, so it was aligned; the alignment is not kept.
    assert len(_fasta_records(out / "target-proteins" / "famA.faa")) == 2


def test_screen_bulk_generated_files_use_new_names() -> None:
    assert screen.ALL_CAZYMES_FILES == [
        "matched-contigs.fna",
        "matched-contig-ids.txt",
        "cazyme-matches.tsv.gz",
        "cazyme-family-summary.tsv",
        "cazyme-class-summary.tsv",
        "evidence/cazyme-hits.tsv.gz",
        ".phu/run.json",
    ]
    assert screen.ALL_PULS_FILES == [
        "matched-contigs.fna",
        "matched-contig-ids.txt",
        "pul-matches.tsv.gz",
        "pul-family-support.tsv.gz",
        "pul-summary.tsv",
        "pul-substrate-summary.tsv",
        "evidence/cazyme-hits.tsv.gz",
        ".phu/run.json",
    ]


def test_build_target_hmms_ignores_leftover_alignment(workspace: Path) -> None:
    proteins_dir = workspace / "out" / "target-proteins"
    proteins_dir.mkdir(parents=True)
    (proteins_dir / "famA.faa").write_text(f">p1\n{SEQ_A}\n")
    (proteins_dir / "famA-aligned.faa").write_text(f">p1\n{SEQ_A}\n")

    screen._build_target_hmms(proteins_dir, workspace / "out")

    assert _listing(workspace / "out" / "target-hmms") == ["famA.hmm"]


def test_jack_output_names_single_seed(workspace: Path) -> None:
    seeds = workspace / "seed.faa"
    seeds.write_text(f">seedA\n{SEQ_A}\n")
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna",
            seed_marker=seeds,
            outdir=out,
            keep_proteins=True,
            save_hmm=True,
        )
    )

    assert _listing(out) == [
        ".phu/prediction.json",
        "final-iteration.hmm",
        "jackhmmer-hits.tsv",
        "jackhmmer-iterations.tsv",
        "matched-contig-ids.txt",
        "matched-contigs.fna",
        "predicted-proteins.faa",
    ]


def test_jack_output_names_multiple_seeds(workspace: Path) -> None:
    seeds = workspace / "seeds.faa"
    seeds.write_text(f">seedA\n{SEQ_A}\n>seedB\n{SEQ_B}\n")
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna",
            seed_marker=seeds,
            outdir=out,
            save_hmm=True,
        )
    )

    assert _listing(out) == [
        "final-iteration-hmms/seedA.hmm",
        "final-iteration-hmms/seedB.hmm",
        "jackhmmer-hits.tsv",
        "jackhmmer-iterations.tsv",
        "matched-contig-ids.txt",
        "matched-contigs.fna",
    ]


def _fake_cluster_run(cmd, **kwargs) -> None:
    """Write what vclust and seqkit would write for each step."""
    args = [str(part) for part in cmd]
    output = Path(args[args.index("-o") + 1])
    if args[1] == "prefilter":
        output.write_text("pairs\n")
    elif args[1] == "align":
        output.write_text("ani\n")
        # vclust derives this name itself; there is no option for it.
        output.with_suffix(".ids.tsv").write_text("ids\n")
    elif args[1] == "cluster":
        assert Path(args[args.index("--ids") + 1]).exists()
        output.write_text("object\tcluster\nc1\tc1\nc2\tc1\n")
    else:
        output.write_text(">c1\nACGT\n")


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        (
            "dereplication",
            ["clusters.tsv", "dereplicated-representatives.fna"],
        ),
        ("votu", ["cluster-representatives.fna", "clusters.tsv"]),
        ("species", ["cluster-representatives.fna", "species-clusters.tsv"]),
    ],
)
def test_cluster_output_names(tmp_path: Path, monkeypatch, mode, expected) -> None:
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c1\nACGT\n>c2\nACGT\n")
    monkeypatch.setattr(cluster, "_binaries", lambda: ("vclust", "seqkit"))
    monkeypatch.setattr(cluster, "run", _fake_cluster_run)
    out = tmp_path / "out"

    cluster._cluster(
        cluster.ClusterConfig(
            mode=cluster.Mode(mode), input_contigs=contigs, output_folder=out
        )
    )

    assert _listing(out) == sorted(
        [
            "ani-pairs.tsv",
            "ani-sequence-ids.tsv",
            "cluster-representative-ids.txt",
            "prefilter-pairs.txt",
            *expected,
        ]
    )


def test_cluster_default_output_folder() -> None:
    config = cluster.ClusterConfig(
        mode=cluster.Mode.votu, input_contigs=Path("contigs.fna")
    )
    assert config.output_folder == Path("phu-cluster")
