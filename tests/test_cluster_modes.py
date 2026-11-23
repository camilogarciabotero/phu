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


def test_cli_passes_genus_to_the_run(tmp_path, monkeypatch):
    import phu.cli as cli_module

    contigs = tmp_path / "g.fna"
    contigs.write_text(">a\nACGT\n")
    seen = {}
    monkeypatch.setattr(cli_module, "_cluster", lambda cfg: seen.setdefault("cfg", cfg))
    result = runner.invoke(
        app, ["cluster", "--mode", "genus", "-i", str(contigs), "-o", str(tmp_path / "o")]
    )
    assert result.exit_code == 0, result.output
    assert seen["cfg"].mode is Mode.genus
    assert seen["cfg"].plan().output.name == "genus-clusters.tsv"


def test_cli_rejects_unknown_mode(tmp_path):
    contigs = tmp_path / "g.fna"
    contigs.write_text(">a\nACGT\n")
    result = runner.invoke(app, ["cluster", "--mode", "family", "-i", str(contigs)])
    assert result.exit_code == 1
    assert "genus" in result.stderr


# ---------- real vclust ----------


def _write_genomes(path: Path) -> None:
    """gA; gA2 = gA with 1% substitutions; gB = 80% of gA + 20% new; gC unrelated."""
    rng = random.Random(7)

    def random_dna(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mutate(seq: str, fraction: float) -> str:
        bases = list(seq)
        for i in rng.sample(range(len(bases)), int(len(bases) * fraction)):
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
        return "".join(bases)

    a = random_dna(30_000)
    # Names of two or more characters: vclust 1.3.1 loses pairs between
    # single-character names when aligning on several threads.
    genomes = {
        "gA": a,
        "gA2": mutate(a, 0.01),
        "gB": a[:24_000] + random_dna(6_000),
        "gC": random_dna(30_000),
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


def _vclust_installed() -> bool:
    from phu.core.exec import CmdNotFound, _executable

    try:
        _executable(["vclust", "vclust.py"])
    except CmdNotFound:
        return False
    return True


@pytest.mark.skipif(not _vclust_installed(), reason="vclust is not installed")
@pytest.mark.parametrize("threads", [1, 0])  # 0 = all cores, the CLI default
def test_real_vclust_separates_genus_from_species(tmp_path, threads):
    genomes = tmp_path / "g.fna"
    _write_genomes(genomes)
    plans = {}
    for mode in ("genus", "species"):
        plans[mode] = cluster._cluster(
            ClusterConfig(
                mode=Mode(mode),
                input_contigs=genomes,
                output_folder=tmp_path / mode,
                threads=threads,
            )
        )

    tani = _tani(plans["genus"].ani)
    # Expectations come from the computed tANI, not from the construction.
    assert tani[frozenset(("gA", "gA2"))] >= 0.95
    for other in ("gA", "gA2"):
        assert 0.70 <= tani[frozenset(("gB", other))] < 0.95
    assert all(value < 0.70 for pair, value in tani.items() if "gC" in pair)

    assert _groups(plans["genus"].output) == {
        frozenset({"gA", "gA2", "gB"}),
        frozenset({"gC"}),
    }
    assert _groups(plans["species"].output) == {
        frozenset({"gA", "gA2"}),
        frozenset({"gB"}),
        frozenset({"gC"}),
    }
    reps = plans["genus"].representatives_fna.read_text().count(">")
    assert reps == 2


def test_single_character_ids_are_detected(tmp_path):
    fasta = tmp_path / "g.fna"
    fasta.write_text(">ctg1 desc\nACGT\n>A\nACGT\n")
    assert cluster.has_single_character_ids(fasta)
    fasta.write_text(">ctg1\nACGT\n>AB x\nACGT\n")
    assert not cluster.has_single_character_ids(fasta)


def _vclust_installed() -> bool:
    from phu.core.exec import CmdNotFound, _executable

    try:
        _executable(["vclust", "vclust.py"])
    except CmdNotFound:
        return False
    return True


@pytest.mark.skipif(not _vclust_installed(), reason="vclust is not installed")
def test_real_vclust_single_character_names_keep_every_pair(tmp_path):
    """vclust 1.3.1 drops these pairs on several threads; phu aligns on one."""
    genomes = tmp_path / "g.fna"
    _write_genomes(genomes)
    renamed = tmp_path / "short.fna"
    renamed.write_text(
        genomes.read_text()
        .replace(">gA2\n", ">D\n")
        .replace(">gA\n", ">A\n")
        .replace(">gB\n", ">B\n")
        .replace(">gC\n", ">C\n")
    )
    plan = cluster._cluster(
        ClusterConfig(
            mode=Mode.genus, input_contigs=renamed, output_folder=tmp_path / "o"
        )
    )
    assert plan.align_threads == 1
    assert _groups(plan.output) == {frozenset({"A", "D", "B"}), frozenset({"C"})}


def test_representative_ids_are_listed_once_per_cluster(tmp_path):
    from phu.commands.cluster import _extract_cluster_ids

    tsv = tmp_path / "clusters.tsv"
    tsv.write_text("object\tcluster\na\ta\na_copy\ta\nb\tb\nc\ta\n")
    out = tmp_path / "ids.txt"

    _extract_cluster_ids(tsv, out)

    assert out.read_text() == "a\nb\n"
