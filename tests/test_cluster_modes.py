"""Effective vclust settings of each phu cluster mode, and a real vclust run."""

from __future__ import annotations

import random
import re
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from phu.cli import app
from phu.commands import cluster
from phu.commands.cluster import ClusterConfig, Mode, parse_vclust_params

runner = CliRunner()


def _plan(mode: str, params: str = "", **fields) -> cluster.ClusterPlan:
    return ClusterConfig(
        mode=Mode(mode),
        input_contigs=Path("g.fna"),
        vclust_params=parse_vclust_params(params),
        **fields,
    ).plan()


def _commands(plan: cluster.ClusterPlan) -> tuple[list[str], list[str], list[str]]:
    plan.vclust_bin = "vclust"
    return (
        cluster.prefilter_command(plan),
        cluster.align_command(plan),
        cluster.cluster_command(plan),
    )


def _flag(cmd: list[str], name: str):
    return cmd[cmd.index(name) + 1] if name in cmd else None


def _plain(text: str) -> str:
    return re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)


# ---------- recipes ----------


def test_genus_recipe():
    plan = _plan("genus")
    assert plan.algorithm == "complete"
    assert plan.metric == "tani"
    assert plan.thresholds == {"tani": 0.70}  # no query-coverage requirement
    assert plan.prefilter_min_ident == 0.70
    assert plan.align_filters == {}
    assert plan.output.name == "genus-clusters.tsv"
    assert plan.representatives_fna.name == "cluster-representatives.fna"


def test_species_recipe():
    plan = _plan("species")
    assert (plan.algorithm, plan.metric) == ("complete", "tani")
    assert plan.thresholds == {"tani": 0.95}
    assert plan.prefilter_min_ident == 0.70
    assert plan.align_filters == {}
    assert plan.output.name == "species-clusters.tsv"


@pytest.mark.parametrize(
    ("mode", "algorithm", "representatives"),
    [
        ("dereplication", "cd-hit", "dereplicated-representatives.fna"),
        ("votu", "leiden", "cluster-representatives.fna"),
    ],
)
def test_ani_modes_keep_their_defaults(mode, algorithm, representatives):
    plan = _plan(mode)
    prefilter, align, clust = _commands(plan)
    assert plan.algorithm == algorithm and plan.metric == "ani"
    assert plan.thresholds == {"ani": 0.95, "qcov": 0.85}
    assert _flag(prefilter, "--min-ident") == "0.95"
    assert _flag(align, "--out-ani") == "0.95"
    assert _flag(align, "--out-qcov") == "0.85"
    assert _flag(clust, "--ani") == "0.95" and _flag(clust, "--qcov") == "0.85"
    assert plan.output.name == "clusters.tsv"
    assert plan.representatives_fna.name == representatives


@pytest.mark.parametrize(("mode", "tani"), [("genus", "0.7"), ("species", "0.95")])
def test_tani_commands_follow_upstream_recipe(mode, tani):
    prefilter, align, clust = _commands(_plan(mode))
    assert _flag(prefilter, "--min-ident") == "0.7"
    # Every aligned pair is written: no ANI or coverage filter.
    assert not [part for part in align if part.startswith("--out-")]
    assert _flag(clust, "--algorithm") == "complete"
    assert _flag(clust, "--metric") == "tani"
    assert _flag(clust, "--tani") == tani
    assert "--ani" not in clust and "--qcov" not in clust
    assert "--out-repr" in clust


# ---------- overrides ----------


def test_threshold_override_reaches_cluster_and_align_filters():
    prefilter, align, clust = _commands(_plan("votu", "--ani 0.97 --qcov 0.9"))
    assert _flag(clust, "--ani") == "0.97" and _flag(clust, "--qcov") == "0.9"
    assert _flag(align, "--out-ani") == "0.97" and _flag(align, "--out-qcov") == "0.9"
    assert _flag(prefilter, "--min-ident") == "0.95"


def test_overrides_are_not_duplicated():
    prefilter, _, clust = _commands(
        _plan(
            "genus",
            "--min-ident 0.6 --tani 0.75 --algorithm single --leiden-resolution 0.9",
        )
    )
    assert prefilter.count("--min-ident") == 1 and _flag(prefilter, "--min-ident") == "0.6"
    assert clust.count("--tani") == 1 and _flag(clust, "--tani") == "0.75"
    assert clust.count("--algorithm") == 1 and _flag(clust, "--algorithm") == "single"
    assert _flag(clust, "--leiden-resolution") == "0.9"


def test_config_fields_apply_to_the_effective_metric():
    plan = _plan("genus", ani_cutoff=0.8)
    assert plan.thresholds == {"tani": 0.8}
    # --vclust-params take precedence over config fields.
    assert _plan("genus", "--tani 0.75", ani_cutoff=0.8).thresholds == {"tani": 0.75}


def test_metric_change_with_its_threshold():
    plan = _plan("species", "--metric gani --gani 0.9")
    assert plan.metric == "gani" and plan.thresholds == {"gani": 0.9}


def test_safe_align_filter_is_kept():
    _, align, _ = _commands(_plan("genus", "--out-tani 0.7"))
    assert _flag(align, "--out-tani") == "0.7"


@pytest.mark.parametrize(
    ("mode", "params", "message"),
    [
        # Would drop pairs with ANI < 0.95 that genus clustering needs.
        ("genus", "--out-ani 0.95", "--out-ani"),
        ("genus", "--out-qcov 0.5", "--out-qcov"),
        ("species", "--out-tani 0.97", "--out-tani"),
        ("votu", "--out-ani 0.97", "--out-ani"),
        # A tANI threshold must not silently become an ANI threshold.
        ("genus", "--metric ani", "--metric ani needs its own threshold"),
        # Prefilter at 0.95 would drop pairs between 0.90 and 0.95.
        ("votu", "--ani 0.90", "--min-ident"),
        ("species", "--tani 0.6", "--min-ident"),
    ],
)
def test_contradictory_settings_are_rejected(mode, params, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        _plan(mode, params)


def test_lowered_threshold_with_explicit_min_ident_is_accepted():
    plan = _plan("votu", "--ani 0.90 --min-ident 0.9")
    assert plan.prefilter_min_ident == 0.9 and plan.thresholds["ani"] == 0.9


# ---------- CLI ----------


def test_cli_advertises_genus():
    result = runner.invoke(app, ["cluster", "--help"], terminal_width=200)
    output = _plain(result.stdout)
    assert result.exit_code == 0
    assert "dereplication | votu | genus | species" in output
    assert "--mode genus" in output


def test_cli_rejects_contradiction_before_running(tmp_path, monkeypatch):
    contigs = tmp_path / "g.fna"
    contigs.write_text(">a\nACGT\n")
    monkeypatch.setattr(
        cluster, "run", lambda *a, **k: pytest.fail("vclust must not run")
    )
    result = runner.invoke(
        app, ["cluster", "--mode", "genus", "-i", str(contigs), "-p", "--out-ani 0.95"]
    )
    assert result.exit_code == 1
    assert "--out-ani" in result.stderr


def test_cli_rejects_unknown_mode(tmp_path):
    contigs = tmp_path / "g.fna"
    contigs.write_text(">a\nACGT\n")
    result = runner.invoke(app, ["cluster", "--mode", "family", "-i", str(contigs)])
    assert result.exit_code == 1
    assert "genus" in result.stderr


# ---------- real vclust ----------


def _write_genomes(path: Path) -> None:
    """A; A2 = A with 1% substitutions; B = 80% of A + 20% new; C unrelated."""
    rng = random.Random(7)

    def random_dna(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mutate(seq: str, fraction: float) -> str:
        bases = list(seq)
        for i in rng.sample(range(len(bases)), int(len(bases) * fraction)):
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
        return "".join(bases)

    a = random_dna(30_000)
    genomes = {
        "A": a,
        "A2": mutate(a, 0.01),
        "B": a[:24_000] + random_dna(6_000),
        "C": random_dna(30_000),
    }
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in genomes.items()))


def _groups(table: Path) -> set[frozenset[str]]:
    by_cluster: dict[str, set[str]] = {}
    for line in table.read_text().splitlines()[1:]:
        obj, rep = line.split("\t")[:2]
        by_cluster.setdefault(rep, set()).add(obj)
    return {frozenset(members) for members in by_cluster.values()}


def _tani(ani_table: Path) -> dict[frozenset[str], float]:
    lines = ani_table.read_text().splitlines()
    header = lines[0].split("\t")
    q, r, t = header.index("query"), header.index("reference"), header.index("tani")
    return {
        frozenset((row[q], row[r])): float(row[t])
        for row in (line.split("\t") for line in lines[1:])
    }


@pytest.mark.skipif(
    shutil.which("vclust") is None and shutil.which("vclust.py") is None,
    reason="vclust is not installed",
)
def test_real_vclust_separates_genus_from_species(tmp_path):
    genomes = tmp_path / "g.fna"
    _write_genomes(genomes)
    plans = {}
    for mode in ("genus", "species"):
        # One thread: vclust 1.3.1 'align --filter' loses pairs with more
        # threads (see the cluster docs).
        plans[mode] = cluster._cluster(
            ClusterConfig(
                mode=Mode(mode),
                input_contigs=genomes,
                output_folder=tmp_path / mode,
                threads=1,
            )
        )

    tani = _tani(plans["genus"].ani)
    # Expectations come from the computed tANI, not from the construction.
    assert tani[frozenset(("A", "A2"))] >= 0.95
    for other in ("A", "A2"):
        assert 0.70 <= tani[frozenset(("B", other))] < 0.95
    assert all(value < 0.70 for pair, value in tani.items() if "C" in pair)

    assert _groups(plans["genus"].output) == {
        frozenset({"A", "A2", "B"}),
        frozenset({"C"}),
    }
    assert _groups(plans["species"].output) == {
        frozenset({"A", "A2"}),
        frozenset({"B"}),
        frozenset({"C"}),
    }
    reps = plans["genus"].representatives_fna.read_text().count(">")
    assert reps == 2
