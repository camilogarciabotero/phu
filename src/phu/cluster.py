from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, List, Union
import shlex

from ._exec import run, _executable, CmdNotFound


class Mode(str, Enum):
    dereplication = "dereplication"
    votu = "votu"
    species = "species"


def parse_vclust_params(params_str: str) -> Dict[str, Dict[str, Union[str, int, float, bool]]]:
    """
    Parse vclust parameters from command-line style string.
    
    Example: "--min-kmers 20 --min-ident 0.5 --outfmt lite"
    Returns: {"prefilter": {"min-kmers": 20, "min-ident": 0.5}, "align": {"outfmt": "lite"}}
    """
    if not params_str.strip():
        return {}
    
    # Parse the string into tokens
    try:
        tokens = shlex.split(params_str)
    except ValueError as e:
        raise ValueError(f"Invalid parameter string: {e}")
    
    # Parameter mapping based on vclust wiki
    param_mapping = {
        # Prefilter parameters
        "min-kmers": ("prefilter", int),
        "min-ident": ("prefilter", float),
        "batch-size": ("prefilter", int),
        "kmers-fraction": ("prefilter", float),
        "max-seqs": ("prefilter", int),
        
        # Align parameters
        "outfmt": ("align", str),
        "out-ani": ("align", float),
        "out-qcov": ("align", float),
        
        # Cluster parameters
        "ani": ("cluster", float),
        "tani": ("cluster", float),
        "gani": ("cluster", float),
        "qcov": ("cluster", float),
        "leiden-resolution": ("cluster", float),
        "algorithm": ("cluster", str),
        "metric": ("cluster", str),
    }
    
    result = {"prefilter": {}, "align": {}, "cluster": {}}
    
    i = 0
    while i < len(tokens):
        token = tokens[i]
        
        if not token.startswith("--"):
            raise ValueError(f"Expected parameter name starting with '--', got: {token}")
        
        param_name = token[2:]  # Remove '--'
        
        if param_name not in param_mapping:
            raise ValueError(f"Unknown vclust parameter: --{param_name}")
        
        command, param_type = param_mapping[param_name]
        
        # Check if this is a boolean flag or needs a value
        if param_type == bool:
            result[command][param_name] = True
            i += 1
        else:
            # Need a value
            if i + 1 >= len(tokens):
                raise ValueError(f"Parameter --{param_name} requires a value")
            
            value_str = tokens[i + 1]
            
            try:
                if param_type == int:
                    value = int(value_str)
                elif param_type == float:
                    value = float(value_str)
                else:  # str
                    value = value_str
                
                result[command][param_name] = value
                i += 2
            except ValueError:
                raise ValueError(f"Invalid value for --{param_name}: {value_str} (expected {param_type.__name__})")
    
    # Remove empty sections
    return {k: v for k, v in result.items() if v}


@dataclass
class ClusterConfig:
    mode: Mode
    input_contigs: Path
    output_folder: Path = Path("clustered-contigs")
    threads: int = 0  # 0 => all cores
    # defaults reflecting original bash behavior
    ani_cutoff: float = 0.95
    qcov_cutoff: float = 0.85
    metric: str = "ani"  # 'tani' for species
    algorithm: Optional[str] = None  # set by mode
    # Advanced vclust parameters for customization
    vclust_params: Optional[Dict[str, Dict[str, Union[str, int, float, bool]]]] = field(default_factory=dict)
    """
    Custom vclust parameters organized by command step.
    Example:
    {
        "prefilter": {"min-kmers": 20, "batch-size": 100000, "kmers-fraction": 0.2},
        "align": {"outfmt": "lite", "out-ani": 0.90},
        "cluster": {"leiden-resolution": 0.9}
    }
    """

    def plan(self) -> "ClusterPlan":
        algorithm = self.algorithm
        metric = self.metric
        ani_cutoff = self.ani_cutoff
        qcov_cutoff: Optional[float] = self.qcov_cutoff
        output_name = "clusters.tsv"

        if self.mode == Mode.dereplication:
            algorithm = "cd-hit"
        elif self.mode == Mode.votu:
            algorithm = "leiden"
        elif self.mode == Mode.species:
            algorithm = "complete"
            metric = "tani"
            ani_cutoff = 0.95
            qcov_cutoff = None
            output_name = "species.tsv"

        out = self.output_folder
        return ClusterPlan(
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
            vclust_params=self.vclust_params or {},
        )


@dataclass
class ClusterPlan:
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
    vclust_params: Dict[str, Dict[str, Union[str, int, float, bool]]]


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


def _add_custom_params(cmd: List[str], params: Dict[str, Union[str, int, float, bool]]) -> None:
    """Add custom parameters to a command list."""
    for param, value in params.items():
        if isinstance(value, bool):
            if value:  # Only add flag if True
                cmd.append(f"--{param}")
        else:
            cmd.extend([f"--{param}", str(value)])


def _cluster(cfg: ClusterConfig) -> ClusterPlan:
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
    min_ident = 0.7 if plan.mode == Mode.species else 0.95
    print(f"Step 1: Creating pre-alignment filter (min-ident={min_ident})…")

    prefilter_cmd = [
        plan.vclust_bin,
        "prefilter",
        "-i",
        str(plan.input_contigs),
        "-o",
        str(plan.fltr),
        "--min-ident",
        str(min_ident),
        "--threads",
        str(_threads(plan.threads)),
    ]

    # Add custom prefilter parameters
    if "prefilter" in plan.vclust_params:
        _add_custom_params(prefilter_cmd, plan.vclust_params["prefilter"])

    run(prefilter_cmd)

    # Step 2: align
    print("Step 2: Calculating ANI…")

    align_cmd = [
        plan.vclust_bin,
        "align",
        "-i",
        str(plan.input_contigs),
        "-o",
        str(plan.ani),
        "--filter",
        str(plan.fltr),
        "--threads",
        str(_threads(plan.threads)),
    ]

    # Add default output filtering if not overridden by custom params
    custom_align_params = plan.vclust_params.get("align", {})
    if "out-ani" not in custom_align_params:
        align_cmd.extend(["--out-ani", str(cfg.ani_cutoff)])
    if "out-qcov" not in custom_align_params and cfg.qcov_cutoff is not None:
        align_cmd.extend(["--out-qcov", str(cfg.qcov_cutoff)])

    # Add custom align parameters
    if custom_align_params:
        _add_custom_params(align_cmd, custom_align_params)

    run(align_cmd)

    # Step 3: cluster
    print(f"Step 3: Clustering with {plan.algorithm} (metric={plan.metric})…")

    clustercmd = [
        plan.vclust_bin,
        "cluster",
        "-i",
        str(plan.ani),
        "-o",
        str(plan.output),
        "--ids",
        str(plan.ids),
        "--algorithm",
        plan.algorithm,
        "--metric",
        plan.metric,
        "--out-repr",
    ]

    # Add default metric thresholds if not overridden by custom params
    custom_cluster_params = plan.vclust_params.get("cluster", {})
    if plan.metric not in custom_cluster_params:
        clustercmd.extend([f"--{plan.metric}", str(plan.ani_cutoff)])
    if "qcov" not in custom_cluster_params and plan.qcov_cutoff is not None:
        clustercmd.extend(["--qcov", str(plan.qcov_cutoff)])

    # Add custom cluster parameters
    if custom_cluster_params:
        _add_custom_params(clustercmd, custom_cluster_params)

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