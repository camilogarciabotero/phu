"""PHROG families as seeds in `phu jack` (--phrog)."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pyhmmer
import pytest
from test_output_names import SEQ_A, SEQ_B, workspace  # noqa: F401  (fixture)
from typer.testing import CliRunner

import phu.cli as cli_module
from phu.cli import app
from phu.commands import jack
from phu.db.phrog import PhrogAnnotation

runner = CliRunner()


def _build_hmm(name: str, sequence: str) -> pyhmmer.plan7.HMM:
    alphabet = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name=name.encode(), sequence=sequence).digitize(
        alphabet
    )
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        seq, pyhmmer.plan7.Background(alphabet)
    )
    return hmm


PHROGS = {
    "phrog_1": PhrogAnnotation(
        "phrog_1", "terminase large subunit", "head and packaging"
    ),
    "phrog_2": PhrogAnnotation("phrog_2", "portal protein", "head and packaging"),
    "phrog_3": PhrogAnnotation("phrog_3", "holin", "lysis"),
    # A second terminase family: alias:terl selects phrog_1 and phrog_4.
    "phrog_4": PhrogAnnotation(
        "phrog_4", "terminase large subunit", "head and packaging"
    ),
}
MODELS = {
    "phrog_1": SEQ_A,
    "phrog_2": SEQ_B,
    "phrog_3": "MKKLLW" * 6,
    "phrog_4": SEQ_A[2:],
}


@pytest.fixture
def fake_phrogs(monkeypatch):
    monkeypatch.setattr(jack, "phrogs_with_models", lambda: PHROGS)
    monkeypatch.setattr(
        jack,
        "_load_phrog_hmms",
        lambda ids: [(i, _build_hmm(i, MODELS[i])) for i in ids],
    )
    monkeypatch.setattr(jack, "bundle_provenance", lambda: {"bundle_version": "test"})


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---------- configuration ----------


def test_a_seed_is_required(tmp_path):
    with pytest.raises(ValueError, match="seed FASTA, --phrog"):
        jack.JackConfig(input_contigs=tmp_path / "c.fna")


def test_cap_stops_large_selections(fake_phrogs):
    with pytest.raises(
        ValueError, match="selects 4 PHROGs, more than --max-phrog-seeds 3"
    ):
        jack.resolve_phrog_seeds(("category:head", "alias:holin"), cap=3)
    selected, _ = jack.resolve_phrog_seeds(("category:head",), cap=3)
    assert [r.phrog.phrog_id for r in selected] == ["phrog_1", "phrog_2", "phrog_4"]


# ---------- the search from a profile ----------


def test_hmm_seed_iterates_from_the_profile():
    alphabet = pyhmmer.easel.Alphabet.amino()
    targets = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(name=b"c1|gene1", sequence=SEQ_A).digitize(
                alphabet
            ),
            pyhmmer.easel.TextSequence(name=b"c2|gene1", sequence=SEQ_B).digitize(
                alphabet
            ),
        ],
    )
    hits, summary = jack._run_jackhmmer(
        query=_build_hmm("phrog_1", SEQ_A),
        targets=targets,
        seed_id="phrog_1",
        iterations=5,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
    )
    assert [(h.contig, h.prot_id, h.model) for h in hits] == [
        ("c1", "c1|gene1", "phrog_1")
    ]
    assert summary[0]["iteration"] == 1
    assert summary[-1]["converged"] is True
    assert len(summary) <= 5


def test_iterations_cap_the_rounds():
    alphabet = pyhmmer.easel.Alphabet.amino()
    targets = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(name=b"c1|gene1", sequence=SEQ_A).digitize(
                alphabet
            )
        ],
    )
    _, summary = jack._run_jackhmmer(
        _build_hmm("phrog_1", SEQ_A), targets, "phrog_1", 1, 1e-3, 1e-5, 1
    )
    assert len(summary) == 1


def test_seed_that_finds_nothing_ends_cleanly():
    """pyhmmer raises when a round includes nothing; jack reports one round."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    targets = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(name=b"c2|gene1", sequence=SEQ_B).digitize(
                alphabet
            )
        ],
    )
    hits, summary = jack._run_jackhmmer(
        _build_hmm("phrog_1", SEQ_A), targets, "phrog_1", 5, 1e-3, 1e-5, 1
    )
    assert hits == []
    assert summary == [
        {
            "iteration": 1,
            "n_hits": summary[0]["n_hits"],
            "n_included": 0,
            "converged": True,
        }
    ]


# ---------- seed groups ----------


def _hit(contig, protein, seed, score):
    return jack.Hit(
        contig=contig, prot_id=protein, model=seed, bitscore=score, evalue=1e-10
    )


def test_best_hit_per_protein_per_group():
    groups = {"phrog_1": "alias:terl", "phrog_4": "alias:terl", "seedA": "seedA"}
    hits = [
        _hit("c1", "c1|gene1", "phrog_1", 50.0),
        _hit("c1", "c1|gene1", "phrog_4", 80.0),  # same group, same protein: best wins
        _hit("c1", "c1|gene1", "seedA", 40.0),  # other group: kept
        _hit("c1", "c1|gene2", "phrog_1", 30.0),  # other protein: kept
    ]
    kept = jack.best_hit_per_protein_per_group(hits, groups)
    assert [(h.prot_id, h.model, h.bitscore) for h in kept] == [
        ("c1|gene1", "phrog_4", 80.0),
        ("c1|gene1", "seedA", 40.0),
        ("c1|gene2", "phrog_1", 30.0),
    ]


def test_combine_counts_groups_not_seeds():
    groups = {
        "phrog_1": "alias:terl",
        "phrog_4": "alias:terl",
        "phrog_2": "alias:portal",
    }
    # c1 is hit by two terminase PHROGs only; c2 by a terminase and a portal.
    hits = [
        _hit("c1", "c1|gene1", "phrog_1", 50.0),
        _hit("c1", "c1|gene2", "phrog_4", 60.0),
        _hit("c2", "c2|gene1", "phrog_1", 50.0),
        _hit("c2", "c2|gene2", "phrog_2", 70.0),
    ]
    _, all_contigs = jack._choose_top_hits_per_contig(hits, 1, "all", 1, 2, groups)
    assert all_contigs == ["c2"]
    _, two_groups = jack._choose_top_hits_per_contig(hits, 1, "threshold", 2, 2, groups)
    assert two_groups == ["c2"]
    # Without groups every seed counts, and c1 would pass threshold 2.
    _, by_seed = jack._choose_top_hits_per_contig(hits, 1, "threshold", 2, 3)
    assert by_seed == ["c1", "c2"]


# ---------- whole runs ----------


def test_phrog_only_run(workspace, fake_phrogs):  # noqa: F811
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna",
            phrog_queries=("alias:terl", "phrog_2"),
            outdir=out,
        )
    )
    hits = _rows(out / "jackhmmer-hits.tsv")
    # alias:terl (phrog_1, phrog_4, both from SEQ_A) finds c1 and c2;
    # phrog_2 (SEQ_B) finds c1.
    # One hit per contig is kept (--top-per-contig 1).
    assert sorted(h["contig"] for h in hits) == ["c1", "c2"]
    assert {h["seed_source"] for h in hits} == {"phrog"}
    assert next(h for h in hits if h["contig"] == "c2")["seed_group"] == "alias:terl"
    assert {h["seed_id"] for h in _rows(out / "jackhmmer-iterations.tsv")} == {
        "phrog_1",
        "phrog_4",
        "phrog_2",
    }
    record = json.loads((out / ".phu" / "query.json").read_text())
    assert record["seed_fasta"] is None
    assert record["phrog_queries"] == ["alias:terl", "phrog_2"]
    assert record["phrog_seeds"] == [
        {
            "id": "phrog_1",
            "annotation": "terminase large subunit",
            "category": "head and packaging",
            "selected_by": "alias:terl",
        },
        {
            "id": "phrog_4",
            "annotation": "terminase large subunit",
            "category": "head and packaging",
            "selected_by": "alias:terl",
        },
        {
            "id": "phrog_2",
            "annotation": "portal protein",
            "category": "head and packaging",
            "selected_by": "phrog_2",
        },
    ]
    assert record["phu_db"] == {"bundle_version": "test"}


def test_mixed_seeds_count_together(workspace, fake_phrogs):  # noqa: F811
    seeds = workspace / "seed.faa"
    seeds.write_text(f">seedA\n{SEQ_A}\n")
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna",
            seed_marker=seeds,
            phrog_queries=("phrog_1",),
            outdir=out,
            combine_mode="all",
        )
    )
    hits = _rows(out / "jackhmmer-hits.tsv")
    # "all" needs both the FASTA seed and the PHROG seed on a contig.
    assert {h["contig"] for h in hits} == {"c1", "c2"}
    assert {(h["seed_id"], h["seed_source"]) for h in hits} == {
        ("seedA", "fasta"),
        ("phrog_1", "phrog"),
    }


def test_alias_is_one_group_in_a_run(workspace, fake_phrogs):  # noqa: F811
    seeds = workspace / "seed.faa"
    seeds.write_text(f">seedB\n{SEQ_B}\n")
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna",
            seed_marker=seeds,
            phrog_queries=("alias:terl",),
            outdir=out,
            combine_mode="all",
        )
    )
    # Two groups (seedB, alias:terl): c1 carries both SEQ_A and SEQ_B, so it
    # passes "all" although only one of the two terminase PHROGs is counted.
    hits = _rows(out / "jackhmmer-hits.tsv")
    assert {h["contig"] for h in hits} == {"c1"}
    assert {h["seed_group"] for h in hits} == {"seedB", "alias:terl"}
    # Each group adds c1|gene1 once, by its best hit.
    terl = [h for h in hits if h["seed_group"] == "alias:terl"]
    assert [h["protein"] for h in terl] == ["c1|gene1"]
    iterations = _rows(out / "jackhmmer-iterations.tsv")
    assert {(r["seed_id"], r["seed_group"]) for r in iterations} == {
        ("seedB", "seedB"),
        ("phrog_1", "alias:terl"),
        ("phrog_4", "alias:terl"),
    }


def test_threshold_is_checked_against_groups(workspace, fake_phrogs):  # noqa: F811
    with pytest.raises(ValueError, match="number of seed groups \\(1"):
        jack._jack(
            jack.JackConfig(
                input_contigs=workspace / "contigs.fna",
                phrog_queries=("alias:terl",),
                outdir=workspace / "out",
                combine_mode="threshold",
                min_seed_hits=2,
            )
        )


def test_seed_id_clash_is_refused(workspace, fake_phrogs):  # noqa: F811
    seeds = workspace / "seed.faa"
    seeds.write_text(f">phrog_1\n{SEQ_A}\n")
    with pytest.raises(ValueError, match="named like PHROG seeds: phrog_1"):
        jack._jack(
            jack.JackConfig(
                input_contigs=workspace / "contigs.fna",
                seed_marker=seeds,
                phrog_queries=("phrog_1",),
                outdir=workspace / "out",
            )
        )


def test_fasta_only_tables_say_fasta(workspace):  # noqa: F811
    seeds = workspace / "seed.faa"
    seeds.write_text(f">seedA\n{SEQ_A}\n")
    out = workspace / "out"
    jack._jack(
        jack.JackConfig(
            input_contigs=workspace / "contigs.fna", seed_marker=seeds, outdir=out
        )
    )
    assert {h["seed_source"] for h in _rows(out / "jackhmmer-hits.tsv")} == {"fasta"}
    record = json.loads((out / ".phu" / "query.json").read_text())
    assert record["phrog_seeds"] == [] and "phu_db" not in record


# ---------- CLI ----------


def _plain(text: str) -> str:
    return re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)


def test_cli_help_shows_phrog_option():
    result = runner.invoke(app, ["jack", "--help"], terminal_width=200)
    assert result.exit_code == 0
    assert "--phrog" in _plain(result.stdout)
    assert "--max-phrog-seeds" in _plain(result.stdout)


def test_cli_without_any_seed_fails(tmp_path):
    contigs = tmp_path / "c.fna"
    contigs.write_text(">c1\nACGT\n")
    result = runner.invoke(app, ["jack", "-i", str(contigs)])
    assert result.exit_code == 1
    assert "seed FASTA, --phrog" in result.stderr


def test_cli_passes_phrog_queries(tmp_path, monkeypatch):
    contigs = tmp_path / "c.fna"
    contigs.write_text(">c1\nACGT\n")
    seen = {}
    monkeypatch.setattr(cli_module, "_jack", lambda cfg: seen.setdefault("cfg", cfg))
    result = runner.invoke(
        app,
        [
            "jack",
            "-i",
            str(contigs),
            "--phrog",
            "alias:terl",
            "--phrog",
            "phrog_9",
            "--max-phrog-seeds",
            "10",
        ],
    )
    assert result.exit_code == 0, result.output
    assert seen["cfg"].phrog_queries == ("alias:terl", "phrog_9")
    assert seen["cfg"].seed_marker is None and seen["cfg"].max_phrog_seeds == 10


def test_verbose_log_shows_what_queries_resolved_to(workspace, fake_phrogs, capsys):  # noqa: F811
    from phu.core.log import set_verbosity

    set_verbosity(verbose=True)
    try:
        jack._jack(
            jack.JackConfig(
                input_contigs=workspace / "contigs.fna",
                phrog_queries=("alias:terl", "phrog_2", "phrog_1"),
                outdir=workspace / "out",
            )
        )
    finally:
        set_verbosity()
    lines = [
        line
        for line in capsys.readouterr().out.splitlines()
        if line.startswith("[cmd] phu.phrogs")
    ]
    assert lines[:3] == [
        "[cmd] phu.phrogs.resolve query='alias:terl' phrogs=2 ids=phrog_1,phrog_4",
        "[cmd] phu.phrogs.resolve query='phrog_2' phrogs=1 ids=phrog_2",
        # phrog_1 was already selected by alias:terl, so this query adds none.
        "[cmd] phu.phrogs.resolve query='phrog_1' phrogs=0 ids=-",
    ]
    assert lines[3] == (
        "[cmd] phu.phrogs.seed id=phrog_1 group='alias:terl' "
        "annotation='terminase large subunit' category='head and packaging'"
    )
    assert len(lines) == 6
