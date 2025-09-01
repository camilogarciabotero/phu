from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from ._exec import run, _executable, CmdNotFound


class Mode(str, Enum):
    dereplication = "dereplication"
    votu = "votu-clustering"
    spp = "spp-clustering"


@dataclass
class SeqClustConfig:
    mode: Mode
    input_contigs: Path
    output_folder: Path = Path("clustered-contigs")
    threads: int = 0  # 0 => all cores
    # defaults reflecting original bash behavior
    ani_cutoff: float = 0.95
    qcov_cutoff: float = 0.85
    metric: str = "ani"  # 'tani' for spp
    algorithm: Optional[str] = None  # set by mode

    def plan(self) -> "SeqClustPlan":
        algorithm = self.algorithm
        metric = self.metric
        ani_cutoff = self.ani_cutoff
        qcov_cutoff: Optional[float] = self.qcov_cutoff
        output_name = "clusters.tsv"

        if self.mode == Mode.dereplication:
            algorithm = "cd-hit"
        elif self.mode == Mode.votu:
            algorithm = "leiden"
        elif self.mode == Mode.spp:
            algorithm = "complete"
            metric = "tani"
            ani_cutoff = 0.95
            qcov_cutoff = None
            output_name = "species.tsv"

        out = self.output_folder
        return SeqClustPlan(
            vclust_bin="",
            seqkit_bin="",
            input_contigs=self.input_contigs,
            out_dir=out,
            fltr=out / "fltr.txt",
            ani=out / "ani.tsv",
            ids=out / "ani.ids.tsv",
            output=out / output_name,
            representatives_ids=out / "cluster_representatives_ids.txt",
            representatives_fna=out
            / (
                "dereplicated_representatives.fna"
                if self.mode == Mode.dereplication
                else "representatives.fna"
            ),
            algorithm=algorithm or "",
            metric=metric,
            ani_cutoff=ani_cutoff,
            qcov_cutoff=qcov_cutoff,
            threads=self.threads,
            mode=self.mode,
        )


@dataclass
class SeqClustPlan:
    vclust_bin: str
    seqkit_bin: str
    input_contigs: Path
    out_dir: Path
    fltr: Path
    ani: Path
    ids: Path
    output: Path
    representatives_ids: Path
    representatives_fna: Path
    algorithm: str
    metric: str
    ani_cutoff: float
    qcov_cutoff: Optional[float]
    threads: int
    mode: Mode


def _threads(n: int) -> int:
    import os

    return max(1, n or (os.cpu_count() or 1))


def _binaries() -> tuple[str, str]:
    """
    Discover required binaries for sequencing clustering.
    """
    vclust = _executable(["vclust", "vclust.py"])
    seqkit = _executable(["seqkit"])
    return vclust, seqkit


def _seqclust(cfg: SeqClustConfig) -> SeqClustPlan:
    """
    Run sequence clustering with vclust.
    """
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")

    vclust_bin, seqkit_bin = _binaries()
    plan = cfg.plan()
    plan.vclust_bin = vclust_bin
    plan.seqkit_bin = seqkit_bin

    plan.out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: prefilter
    min_ident = 0.7 if plan.mode == Mode.spp else 0.95
    print(f"Step 1: Creating pre-alignment filter (min-ident={min_ident})…")
    run(
        [
            plan.vclust_bin,
            "prefilter",
            "-i",
            plan.input_contigs,
            "-o",
            plan.fltr,
            "--min-ident",
            f"{min_ident}",
            "--threads",
            str(_threads(plan.threads)),
        ]
    )

    # Step 2: align
    print("Step 2: Calculating ANI…")
    run(
        [
            plan.vclust_bin,
            "align",
            "-i",
            plan.input_contigs,
            "-o",
            plan.ani,
            "--filter",
            plan.fltr,
            "--out-ani",
            f"{cfg.ani_cutoff}",
            "--out-qcov",
            f"{cfg.qcov_cutoff}",
            "--threads",
            str(_threads(plan.threads)),
        ]
    )

    # Step 3: cluster
    print(f"Step 3: Clustering with {plan.algorithm} (metric={plan.metric})…")
    clustercmd = [
        plan.vclust_bin,
        "cluster",
        "-i",
        plan.ani,
        "-o",
        plan.output,
        "--ids",
        plan.ids,
        "--algorithm",
        plan.algorithm,
        "--metric",
        plan.metric,
        "--out-repr",
    ]
    clustercmd += ["--" + plan.metric, f"{plan.ani_cutoff}"]
    if plan.qcov_cutoff is not None:
        clustercmd += ["--qcov", f"{plan.qcov_cutoff}"]
    run(clustercmd)

    # Summary
    print("Clustering complete. Summary:")
    ngenomes, nclusters = _summarize_tsv(plan.output)
    print(f"Total clusters: {nclusters} from {ngenomes} genomes.")

    # Representatives
    print("Extracting cluster representative IDs…")
    _extract_cluster_ids(plan.output, plan.representatives_ids)

    print("Extracting representative sequences with seqkit…")
    run(
        [
            plan.seqkit_bin,
            "grep",
            "-f",
            plan.representatives_ids,
            str(plan.input_contigs),
            "-o",
            str(plan.representatives_fna),
        ]
    )

    return plan


def _summarize_tsv(tsv: Path) -> tuple[int, int]:
    seen_genomes: set[str] = set()
    seen_clusters: set[str] = set()
    with tsv.open() as fh:
        first = True
        for line in fh:
            if first:
                first = False
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            seen_genomes.add(parts[0])
            seen_clusters.add(parts[1])
    return len(seen_genomes), len(seen_clusters)


def _extract_cluster_ids(tsv: Path, out_ids: Path) -> None:
    with tsv.open() as fh, out_ids.open("w") as out:
        first = True
        for line in fh:
            if first:
                first = False
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out.write(parts[1] + "\n")
