from __future__ import annotations

import os
import shlex
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from phu.core.exec import _executable, run
from phu.core.fasta import _open_text, _record_id, extract_records


class Mode(str, Enum):
    dereplication = "dereplication"
    votu = "votu"
    genus = "genus"
    species = "species"


def parse_vclust_params(
    params_str: str,
) -> dict[str, dict[str, str | int | float | bool]]:
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
        "out-tani": ("align", float),
        "out-gani": ("align", float),
        "out-rcov": ("align", float),
        # Cluster parameters
        "ani": ("cluster", float),
        "tani": ("cluster", float),
        "gani": ("cluster", float),
        "qcov": ("cluster", float),
        "rcov": ("cluster", float),
        "leiden-resolution": ("cluster", float),
        "algorithm": ("cluster", str),
        "metric": ("cluster", str),
    }

    result = {"prefilter": {}, "align": {}, "cluster": {}}

    i = 0
    while i < len(tokens):
        token = tokens[i]

        if not token.startswith("--"):
            raise ValueError(
                f"Expected parameter name starting with '--', got: {token}"
            )

        param_name = token[2:]  # Remove '--'

        if param_name not in param_mapping:
            raise ValueError(f"Unknown vclust parameter: --{param_name}")

        command, param_type = param_mapping[param_name]

        # Check if this is a boolean flag or needs a value
        if param_type is bool:
            result[command][param_name] = True
            i += 1
        else:
            # Need a value
            if i + 1 >= len(tokens):
                raise ValueError(f"Parameter --{param_name} requires a value")

            value_str = tokens[i + 1]

            try:
                if param_type is int:
                    value = int(value_str)
                elif param_type is float:
                    value = float(value_str)
                else:  # str
                    value = value_str

                result[command][param_name] = value
                i += 2
            except ValueError:
                raise ValueError(
                    f"Invalid value for --{param_name}: {value_str} (expected {param_type.__name__})"
                )

    # Remove empty sections
    return {k: v for k, v in result.items() if v}


@dataclass(frozen=True)
class Recipe:
    """Default vclust settings of one mode, after vclust's documented use cases.

    https://github.com/refresh-bio/vclust/wiki/6-Use-cases
    """

    algorithm: str
    metric: str  # ani | tani | gani
    threshold: float  # minimum value of ``metric`` to join a cluster
    qcov: float | None  # minimum query coverage, or None for no requirement
    min_ident: float  # prefilter minimum identity
    # Whether align writes only pairs passing the clustering thresholds
    # (smaller ANI table) or every aligned pair, as upstream's tANI recipe does.
    filter_align: bool
    output_name: str
    representatives_name: str = "cluster-representatives.fna"


RECIPES: dict[Mode, Recipe] = {
    Mode.dereplication: Recipe(
        "cd-hit",
        "ani",
        0.95,
        0.85,
        0.95,
        True,
        "clusters.tsv",
        "dereplicated-representatives.fna",
    ),
    Mode.votu: Recipe("leiden", "ani", 0.95, 0.85, 0.95, True, "clusters.tsv"),
    Mode.genus: Recipe(
        "complete", "tani", 0.70, None, 0.70, False, "genus-clusters.tsv"
    ),
    Mode.species: Recipe(
        "complete", "tani", 0.95, None, 0.70, False, "species-clusters.tsv"
    ),
}

# Similarity measures with both a clustering flag (``--<name>``) and an align
# output filter (``--out-<name>``).
_MEASURES = ("ani", "tani", "gani", "qcov", "rcov")


@dataclass
class ClusterConfig:
    mode: Mode
    input_contigs: Path
    output_folder: Path = Path("phu-cluster")
    threads: int = 0  # 0 => all cores
    # Optional overrides of the mode's recipe (see RECIPES). ``ani_cutoff`` is
    # the threshold of the effective metric (tANI in genus and species modes).
    ani_cutoff: float | None = None
    qcov_cutoff: float | None = None
    metric: str | None = None
    algorithm: str | None = None
    # Advanced vclust parameters by step, from parse_vclust_params; they take
    # precedence over the recipe and the fields above.
    vclust_params: dict[str, dict[str, str | int | float | bool]] | None = field(
        default_factory=dict
    )

    def plan(self) -> ClusterPlan:
        return resolve_plan(self)


@dataclass
class ClusterPlan:
    """Effective settings of every vclust step; commands are built from these only."""

    vclust_bin: str
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
    # Clustering thresholds by measure, e.g. {"tani": 0.7} or {"ani": 0.95, "qcov": 0.85}
    thresholds: dict[str, float]
    prefilter_min_ident: float
    # Align output filters by measure, e.g. {"ani": 0.95, "qcov": 0.85}; empty
    # means every aligned pair is written.
    align_filters: dict[str, float]
    threads: int
    mode: Mode
    # Other user parameters passed through unchanged, by step.
    extra_params: dict[str, dict[str, str | int | float | bool]]
    # Threads for the align step when it must differ from ``threads``.
    align_threads: int | None = None

    @property
    def threshold(self) -> float:
        return self.thresholds[self.metric]


def resolve_plan(cfg: ClusterConfig) -> ClusterPlan:
    """Effective settings: mode recipe, then config fields, then --vclust-params.

    Raises ValueError for combinations that would make a step discard genome
    pairs the clustering needs, or that leave the clustering metric without a
    threshold.
    """
    recipe = RECIPES[cfg.mode]
    params = {step: dict(values) for step, values in (cfg.vclust_params or {}).items()}
    prefilter = params.get("prefilter", {})
    align = params.get("align", {})
    cluster = params.get("cluster", {})

    algorithm = str(cluster.pop("algorithm", None) or cfg.algorithm or recipe.algorithm)
    metric = str(cluster.pop("metric", None) or cfg.metric or recipe.metric)
    if metric not in ("ani", "tani", "gani"):
        raise ValueError(f"Unknown vclust metric {metric!r}: use ani, tani or gani")

    thresholds: dict[str, float] = {}
    if metric == recipe.metric:
        thresholds[metric] = recipe.threshold
    if recipe.qcov is not None:
        thresholds["qcov"] = recipe.qcov
    if cfg.ani_cutoff is not None:
        thresholds[metric] = cfg.ani_cutoff
    if cfg.qcov_cutoff is not None:
        thresholds["qcov"] = cfg.qcov_cutoff
    for measure in _MEASURES:
        if measure in cluster:
            thresholds[measure] = float(cluster.pop(measure))
    if metric not in thresholds:
        raise ValueError(
            f"--metric {metric} needs its own threshold: add --{metric} <value> "
            f"(the {cfg.mode.value} threshold applies to {recipe.metric})"
        )

    min_ident_given = "min-ident" in prefilter
    min_ident = float(prefilter.pop("min-ident", recipe.min_ident))
    if not min_ident_given and thresholds[metric] < recipe.min_ident:
        raise ValueError(
            f"--{metric} {thresholds[metric]} is below the prefilter identity "
            f"{recipe.min_ident}, which would drop pairs before alignment: "
            f"also set --min-ident (for example --min-ident {thresholds[metric]})"
        )

    align_filters: dict[str, float] = {}
    if recipe.filter_align:
        align_filters = {
            measure: value for measure, value in thresholds.items() if value > 0
        }
    for measure in _MEASURES:
        key = f"out-{measure}"
        if key not in align:
            continue
        value = float(align.pop(key))
        if value > thresholds.get(measure, 0.0):
            raise ValueError(
                f"--{key} {value} would remove pairs that --{measure} "
                f"{thresholds.get(measure, 0.0)} keeps for clustering in "
                f"{cfg.mode.value} mode; lower it or remove it"
            )
        align_filters[measure] = value

    out = cfg.output_folder
    return ClusterPlan(
        vclust_bin="",
        input_contigs=cfg.input_contigs,
        out_dir=out,
        fltr=out / "prefilter-pairs.txt",
        ani=out / "ani-pairs.tsv",
        ids=out / "ani-sequence-ids.tsv",
        output=out / recipe.output_name,
        representatives_ids=out / "cluster-representative-ids.txt",
        representatives_fna=out / recipe.representatives_name,
        algorithm=algorithm,
        metric=metric,
        thresholds=thresholds,
        prefilter_min_ident=min_ident,
        align_filters=align_filters,
        threads=cfg.threads,
        mode=cfg.mode,
        extra_params={"prefilter": prefilter, "align": align, "cluster": cluster},
    )


def _threads(n: int) -> int:
    return max(1, n or (os.cpu_count() or 1))


def _binaries() -> str:
    """
    Discover the binary required for sequence clustering.
    """
    return _executable(["vclust", "vclust.py"])


def _add_custom_params(
    cmd: list[str], params: dict[str, str | int | float | bool]
) -> None:
    """Add custom parameters to a command list."""
    for param, value in params.items():
        if isinstance(value, bool):
            if value:  # Only add flag if True
                cmd.append(f"--{param}")
        else:
            cmd.extend([f"--{param}", str(value)])


def prefilter_command(plan: ClusterPlan) -> list[str]:
    cmd = [
        plan.vclust_bin,
        "prefilter",
        "-i",
        str(plan.input_contigs),
        "-o",
        str(plan.fltr),
        "--min-ident",
        str(plan.prefilter_min_ident),
        "--threads",
        str(_threads(plan.threads)),
    ]
    _add_custom_params(cmd, plan.extra_params.get("prefilter", {}))
    return cmd


def align_command(plan: ClusterPlan) -> list[str]:
    cmd = [
        plan.vclust_bin,
        "align",
        "-i",
        str(plan.input_contigs),
        "-o",
        str(plan.ani),
        "--filter",
        str(plan.fltr),
        "--threads",
        str(
            plan.align_threads
            if plan.align_threads is not None
            else _threads(plan.threads)
        ),
    ]
    for measure, value in plan.align_filters.items():
        cmd.extend([f"--out-{measure}", str(value)])
    _add_custom_params(cmd, plan.extra_params.get("align", {}))
    return cmd


def cluster_command(plan: ClusterPlan) -> list[str]:
    cmd = [
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
    ]
    for measure, value in plan.thresholds.items():
        cmd.extend([f"--{measure}", str(value)])
    cmd.append("--out-repr")
    _add_custom_params(cmd, plan.extra_params.get("cluster", {}))
    return cmd


def has_single_character_ids(fasta: Path) -> bool:
    """Whether any record id is one character long.

    vclust 1.3.1 ``align --filter`` on two or more threads skips pairs between
    such sequences, so the align step then runs on one thread.
    """
    with _open_text(fasta) as fh:
        return any(len(_record_id(line)) == 1 for line in fh if line.startswith(">"))


def _cluster(cfg: ClusterConfig) -> ClusterPlan:
    """
    Run sequence clustering with vclust.
    """
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")

    plan = cfg.plan()
    plan.vclust_bin = _binaries()

    plan.out_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"Step 1: Creating pre-alignment filter (min-ident={plan.prefilter_min_ident})…"
    )
    run(prefilter_command(plan))

    if _threads(plan.threads) > 1 and has_single_character_ids(plan.input_contigs):
        plan.align_threads = 1
        print(
            "Note: some sequence names are a single character; vclust 1.3.1 loses "
            "pairs between them on several threads, so the alignment runs on one "
            "thread. Rename the sequences to align in parallel."
        )

    print("Step 2: Calculating ANI…")
    run(align_command(plan))

    # vclust names the sequence-ID table after the ANI table (<ani stem>.ids.tsv)
    # and has no option for it, so it is renamed once the alignment is done.
    vclust_ids = plan.ani.with_suffix(".ids.tsv")
    if vclust_ids.exists():
        vclust_ids.replace(plan.ids)

    print(
        f"Step 3: Clustering with {plan.algorithm} ({plan.metric} >= {plan.threshold})…"
    )
    run(cluster_command(plan))

    # Summary
    print("Clustering complete. Summary:")
    ngenomes, nclusters = _summarize_tsv(plan.output)
    print(f"Total clusters: {nclusters} from {ngenomes} genomes.")

    # Representatives
    print("Extracting cluster representative IDs…")
    _extract_cluster_ids(plan.output, plan.representatives_ids)

    print("Extracting representative sequences…")
    extract_records(
        plan.input_contigs,
        Path(plan.representatives_ids).read_text().split(),
        plan.representatives_fna,
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
