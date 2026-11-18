from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import typer
from typer.core import TyperGroup

from phu import __version__

from phu.core.progress import ProgressReporter, run_click_task
from phu.db.binary import SEARCHED_DATABASES, build_binary_models, has_binary_models
from phu.core.exec import CmdNotFound
from phu.core.log import set_verbosity
from phu.commands.avger import AvgConfig, run_avg
from phu.db.bundle import bundle_provenance, prepare_bundle
from phu.db.fetch import DatabaseError
from phu.commands.cluster import ClusterConfig, Mode, _cluster, parse_vclust_params
from phu.core.gene_prediction import (
    clean_prediction_cache,
)
from phu.commands.clean import find_targets as clean_targets
from phu.commands.clean import human_size as clean_size
from phu.commands.clean import remove as clean_remove
from phu.commands.jack import JackConfig, _jack
from phu.commands.phrogs import app as phrogs_app
from phu.commands.screen import TRIM_METHODS, ScreenConfig, _screen
from phu.commands.normalize_lineage import (
    NormalizeLineageConfig,
    normalize_lineage_file,
)
from phu.commands.simplify_taxa import TaxaConfig, _simplify_taxa

# Order of `phu --help`: panels appear in the order of their first command,
# so Main (the workflows), then Database, then Utils.
COMMAND_ORDER = (
    "screen",
    "jack",
    "avger",
    "cluster",
    "dbs",
    "phrogs",
    "clean",
    "normalize-lineage",
    "simplify-taxa",
)


class _OrderedGroup(TyperGroup):
    def list_commands(self, ctx):
        names = super().list_commands(ctx)
        rank = {name: i for i, name in enumerate(COMMAND_ORDER)}
        return sorted(names, key=lambda name: rank.get(name, len(rank)))


app = typer.Typer(
    cls=_OrderedGroup,
    help="Phage utilities CLI",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,
    no_args_is_help=True,
)
dbs_app = typer.Typer(
    help="Install and check the phu-db database bundle",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)
app.add_typer(dbs_app, name="dbs", rich_help_panel="Database")
app.add_typer(phrogs_app, name="phrogs", rich_help_panel="Utils")

DBS_NOTICE = (
    "`phu dbs {name}` is no longer available: phu keeps its databases as one pinned "
    "bundle. Use `phu dbs prepare` to download, verify and install it ahead of time; "
    "`screen` and `avger` also install it on first use. "
    "See the Databases page of the documentation."
)


@dbs_app.command("prepare")
def dbs_prepare(
    databases: Optional[list[str]] = typer.Argument(
        None,
        help="Accepted for compatibility; the bundle always holds every database",
        hidden=True,
    ),
    bundle: Optional[Path] = typer.Option(
        None,
        "--bundle",
        help="Install from this pre-downloaded phu-db.tar.gz instead of downloading "
        "(same as PHU_DB_BUNDLE)",
        rich_help_panel="Preparation",
    ),
    verify: bool = typer.Option(
        False,
        "--verify",
        help="Re-hash every installed file, and reinstall if any check fails "
        "(same as PHU_DB_VERIFY=1)",
        rich_help_panel="Preparation",
    ),
    connections: Optional[int] = typer.Option(
        None,
        "--connections",
        min=1,
        max=16,
        help="Parts of the bundle to download at once (default 8; same as "
        "PHU_DB_CONNECTIONS). Use 1 for a single connection",
        rich_help_panel="Preparation",
    ),
    all_dbs: bool = typer.Option(False, "--all", hidden=True),
    force_refresh: bool = typer.Option(False, "--force-refresh", hidden=True),
) -> None:
    """Download, verify and install the pinned phu-db bundle.

    The bundle holds Pfam, KOfam, dbCAN, AVG, PHROG and ECOD.

    Optional: `screen` and `avger` install the bundle on first use. Run this
    first on a login node, or before submitting many jobs at once.
    """
    if databases or all_dbs:
        typer.echo(
            "Note: phu-db is one bundle; all databases are installed together.", err=True
        )
    if force_refresh:
        typer.echo(
            "Note: --force-refresh is ignored; the bundle is pinned to this phu release. "
            "Use --verify to re-check the installed files.",
            err=True,
        )
    if connections is not None:
        os.environ["PHU_DB_CONNECTIONS"] = str(connections)
    try:
        root = prepare_bundle(bundle_path=bundle, deep=verify)
        provenance = bundle_provenance()
    except (DatabaseError, OSError) as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    # avger reads binary copies of the model files; build them now, so the
    # first run does not have to and jobs sharing the install find them ready.
    for name in SEARCHED_DATABASES:
        if has_binary_models(name):
            continue
        path = run_click_task(
            f"Building fast-loading models for {name}", build_binary_models, name
        )
        if path.suffix != ".h3m":
            typer.echo(
                f"Note: could not write fast-loading models for {name}; "
                "phu will read the text file instead.",
                err=True,
            )
    typer.echo(f"{provenance['bundle_version']} is ready: {root}")
    for name, release in provenance["upstream_versions"].items():
        typer.echo(f"  {name}: {release}")


def _dbs_removed(name: str) -> None:
    @dbs_app.command(
        name,
        hidden=True,
        context_settings={
            "allow_extra_args": True,
            "ignore_unknown_options": True,
            "help_option_names": [],
        },
    )
    def removed(ctx: typer.Context) -> None:
        typer.echo(DBS_NOTICE.format(name=name), err=True)


for _name in ("list", "status", "refresh", "remove"):
    _dbs_removed(_name)


@app.callback(invoke_without_command=True)
def _root(
    ctx: typer.Context,
    version: bool = typer.Option(
        False, "-v", "--version", is_eager=True, help="Show version and exit."
    ),
    clean_cache: bool = typer.Option(
        False,
        "--clean-cache",
        is_eager=True,
        hidden=True,
        help="Deprecated: use `phu clean --cache`.",
    ),
) -> None:
    # Each invocation starts from, and leaves behind, the default verbosity.
    set_verbosity()
    ctx.call_on_close(set_verbosity)
    if clean_cache:
        typer.echo("`phu --clean-cache` is deprecated; use `phu clean --cache`.", err=True)
        try:
            cache_dir, existed = clean_prediction_cache()
        except OSError as exc:
            typer.secho(f"Failed to remove cache: {exc}", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)
        if existed:
            typer.echo(f"Removed cache directory: {cache_dir}")
        else:
            typer.echo(f"Cache directory does not exist: {cache_dir}")
        raise typer.Exit(0)

    if version:
        typer.echo(f"phu {__version__}")
        raise typer.Exit(0)
    if ctx.invoked_subcommand is None and not ctx.resilient_parsing:
        typer.echo(ctx.get_help())
        raise typer.Exit(0)  # exit code 0 when no subcommand is given


@app.command("clean", rich_help_panel="Utils")
def clean(
    cache: bool = typer.Option(
        False, "--cache", help="Remove cached protein predictions and derived models."
    ),
    dbs: bool = typer.Option(
        False, "--dbs", help="Remove the installed phu-db bundle and partial downloads."
    ),
    all_: bool = typer.Option(False, "--all", help="Both --cache and --dbs."),
    yes: bool = typer.Option(
        False, "--yes", "-y", help="Do not ask before removing the databases."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be removed; remove nothing."
    ),
) -> None:
    """Remove phu's cache and/or its installed databases.

    Without options, lists what phu keeps and its size, and removes nothing.
    The cache is rebuilt by the next run that needs it; the databases must
    be downloaded again (`phu dbs prepare`), so --dbs asks first.
    """
    kinds = tuple(
        kind
        for kind, wanted in (("cache", cache or all_), ("dbs", dbs or all_))
        if wanted
    )
    listing = not kinds
    targets = clean_targets(kinds or ("cache", "dbs"))
    if not targets:
        typer.echo("Nothing to remove.")
        return
    for target in targets:
        typer.echo(f"{target.kind:<6} {clean_size(target.size):>9}  {target.path}")
    if listing:
        typer.echo("Nothing removed. Use --cache, --dbs or --all to remove.")
        return
    if dry_run:
        typer.echo("Dry run: nothing removed.")
        return
    db_targets = [t for t in targets if t.kind == "dbs"]
    if db_targets and not yes:
        total = clean_size(sum(t.size for t in db_targets))
        if not typer.confirm(
            f"Remove the installed databases ({total})? They must be downloaded again",
            default=False,
        ):
            typer.echo("Databases kept.")
            targets = [t for t in targets if t.kind != "dbs"]
    for target in targets:
        try:
            clean_remove(target)
        except OSError as exc:
            typer.secho(f"Could not remove {target.path}: {exc}", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)
        typer.echo(f"Removed {target.path}")


@app.command("cluster", rich_help_panel="Main")
def cluster(
    mode: str = typer.Option(
        ...,
        "--mode",
        help=" | ".join(m.value for m in Mode),
        rich_help_panel="Clustering",
    ),
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", "-i", exists=True, readable=True, help="Input FASTA", rich_help_panel="Input"
    ),
    output_folder: Path = typer.Option(
        Path("phu-cluster"), "--output-folder", "-o", help="Output directory", rich_help_panel="Output"
    ),
    threads: int = typer.Option(
        0, "--threads", "-t", min=0, help="0=all cores; otherwise N threads", rich_help_panel="Runtime"
    ),
    vclust_params: Optional[str] = typer.Option(
        None,
        "--vclust-params",
        "-p",
        help='Custom vclust parameters: "--min-kmers 20 --outfmt lite --ani 0.97"',
        rich_help_panel="Clustering",
    ),
    quiet: bool = typer.Option(
        False, "--quiet", help="Suppress routine progress output.", rich_help_panel="Runtime"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", help="Show additional progress details.", rich_help_panel="Runtime"
    ),
):
    """
    Sequence clustering wrapper around external 'vclust'.

    Modes: dereplication and votu (ANI >= 0.95, coverage >= 0.85), genus
    (tANI >= 0.70) and species (tANI >= 0.95), both with complete linkage.
    Genus and species groups are operational, not ICTV names; use complete or
    near-complete genomes.

    For advanced usage, provide custom vclust parameters as a quoted string.
    See the vclust wiki for parameter details: https://github.com/refresh-bio/vclust/wiki

    Example:
        phu cluster --mode genus -i genomes.fna -o phu-genus
        phu cluster --mode votu --input-contigs genomes.fna --vclust-params="--min-kmers 20 --outfmt lite"
    """
    set_verbosity(quiet=quiet, verbose=verbose)
    # Convert CLI string to enum for compatibility.
    try:
        mode_enum = Mode(mode)
    except ValueError:
        typer.secho(
            "Invalid mode. Use one of: " + ", ".join(m.value for m in Mode),
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)

    # Parse vclust_params
    parsed_params = {}
    if vclust_params:
        try:
            parsed_params = parse_vclust_params(vclust_params)
            typer.echo(f"Using custom vclust parameters: {vclust_params}")
        except ValueError as e:
            typer.secho(
                f"Error parsing vclust parameters: {e}", fg=typer.colors.RED, err=True
            )
            raise typer.Exit(1)

    # Build config
    cfg = ClusterConfig(
        mode=mode_enum,
        input_contigs=input_contigs,
        output_folder=output_folder,
        threads=threads,
        vclust_params=parsed_params,
    )

    try:
        cfg.plan()  # reject contradictory settings before running vclust
    except ValueError as e:
        typer.secho(f"Invalid clustering settings: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    try:
        run_click_task("Clustering contigs", _cluster, cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except CmdNotFound as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        typer.echo(
            "Required executable on PATH: 'vclust' (or 'vclust.py')"
        )
        raise typer.Exit(1)


@app.command("normalize-lineage", rich_help_panel="Utils")
def normalize_lineage(
    input_file: Path = typer.Option(
        ..., "--input-file", "-i", exists=True, readable=True, help="Input delimited table", rich_help_panel="Input"
    ),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output normalized table", rich_help_panel="Output"),
    add_lineage: bool = typer.Option(False, "--add-lineage", "-a", help="Append the compact lineage column", rich_help_panel="Output"),
    lineage_col: str = typer.Option("compact_lineage", "--lineage-col", "-l", help="Name of the lineage column", rich_help_panel="Output"),
    sep: Optional[str] = typer.Option(None, "--sep", "-s", help="Explicit input/output separator", rich_help_panel="Input"),
    strict: bool = typer.Option(False, "--strict", help="Fail on unparsed or rank-mismatched values", rich_help_panel="Quality checks"),
    quiet: bool = typer.Option(False, "--quiet", help="Suppress the QA summary", rich_help_panel="Runtime"),
    verbose: bool = typer.Option(
        False, "--verbose", help="Show additional progress details.", rich_help_panel="Runtime"
    ),
):
    """Normalize vContact lineage predictions without changing column names."""
    set_verbosity(quiet=quiet, verbose=verbose)
    config = NormalizeLineageConfig(
        input_file=input_file,
        output_file=output_file,
        add_lineage=add_lineage,
        lineage_col=lineage_col,
        sep=sep,
        strict=strict,
    )
    try:
        summary = run_click_task(
            "Normalizing lineage",
            normalize_lineage_file,
            config,
        )
        if not quiet:
            typer.echo(
                "QA Summary: "
                + ", ".join(
                    f"{name}={getattr(summary, name)}"
                    for name in (
                        "cells_examined",
                        "normalized",
                        "unchanged",
                        "missing",
                        "unparsed",
                        "rank_mismatch",
                        "multi_candidate",
                        "lineage_skipped_unparsed",
                    )
                ),
                err=True,
            )
            typer.echo(f"Normalized lineage data written to {output_file}", err=True)
    except ValueError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except OSError as exc:
        typer.secho(f"Error processing {input_file}: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


@app.command("simplify-taxa", rich_help_panel="Utils", hidden=True)
def simplify_taxa(
    input_file: Path = typer.Option(
        ...,
        "--input-file",
        "-i",
        exists=True,
        readable=True,
        help="Input vContact final_assignments.csv",
        rich_help_panel="Input",
    ),
    output_file: Path = typer.Option(
        ..., "--output-file", "-o", help="Output file path (.csv or .tsv)", rich_help_panel="Output"
    ),
    add_lineage: bool = typer.Option(
        False,
        "--add-lineage",
        "-a",
        help="Append compact_lineage column from deepest simplified rank",
        rich_help_panel="Output",
    ),
    lineage_col: str = typer.Option(
        "compact_lineage", "--lineage-col", "-l", help="Name of the lineage column", rich_help_panel="Output"
    ),
    sep: Optional[str] = typer.Option(
        None,
        "--sep",
        "-s",
        help="Override delimiter: ',' or '\\t'. Auto-detected from extension if not set",
        rich_help_panel="Input",
    ),
    quiet: bool = typer.Option(
        False, "--quiet", help="Suppress routine progress output.", rich_help_panel="Runtime"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", help="Show additional progress details.", rich_help_panel="Runtime"
    ),
):
    """
    Deprecated: simplify vContact taxonomy prediction columns into compact lineage codes.

    This command is deprecated and will be removed in a future release. It is
    hidden from `phu --help` but still works as before.

    Transforms verbose vContact taxonomy strings like 'novel_genus_1_of_novel_family_2_of_Caudoviricetes'
    into compact codes like 'Caudoviricetes:NF2:NG1'.

    Example:
        phu simplify-taxa -i final_assignments.csv -o simplified.csv --add-lineage
    """
    typer.echo(
        "Warning: 'phu simplify-taxa' is deprecated and will be removed in a future release.",
        err=True,
    )
    set_verbosity(quiet=quiet, verbose=verbose)
    # Build config
    cfg = TaxaConfig(
        input_file=input_file,
        output_file=output_file,
        add_lineage=add_lineage,
        lineage_col=lineage_col,
        sep=sep,
    )

    try:
        run_click_task("Simplifying taxonomy", _simplify_taxa, cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except RuntimeError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


@app.command("avger", rich_help_panel="Main")
def avger(
    input_contigs: Path = typer.Option(
        ...,
        "--input-contigs",
        "-i",
        exists=True,
        readable=True,
        help="Trusted viral contigs FASTA",
        rich_help_panel="Input",
    ),
    output_folder: Path = typer.Option(
        Path("phu-avger"), "--output-folder", "-o", help="Output directory", rich_help_panel="Output"
    ),
    threads: int = typer.Option(1, "--threads", "-t", min=1, help="Threads for prediction and annotation", rich_help_panel="Runtime"),
    mode: str = typer.Option(
        "meta", "--mode", "-m", help="pyrodigal mode: meta|single", rich_help_panel="Prediction"
    ),
    min_gene_len: int = typer.Option(90, "--min-gene-len", min=1, help="Minimum gene length for pyrodigal (nt)", rich_help_panel="Prediction"),
    min_protein_len_aa: int = typer.Option(30, "--min-protein-len-aa", min=1, help="Minimum translated protein length to keep (aa)", rich_help_panel="Prediction"),
    translation_table: Optional[int] = typer.Option(
        None,
        "--ttable",
        "-T",
        min=1,
        help="NCBI translation table; default uses each contig's predicted table",
        rich_help_panel="Prediction",
    ),
    min_amg_weight: float = typer.Option(0.6, "--min-amg-weight", min=0.0, max=1.0, help="Minimum AMG evidence weight", rich_help_panel="Scoring"),
    filter_mode: str = typer.Option("standard", "--filter-mode", help="Candidate filter mode: standard|strict|none", rich_help_panel="Scoring"),
    keep_hits: bool = typer.Option(False, "--keep-hits/--no-keep-hits", help="Keep intermediate annotation hits", rich_help_panel="Output"),
    scaffold_avl_cutoff: float = typer.Option(3.0, "--scaffold-avl-cutoff", min=0.0, help="Minimum scaffold AVL score (exclusive)", rich_help_panel="Scoring"),
    gene_vl_cutoff: float = typer.Option(3.0, "--gene-vl-cutoff", min=0.0, help="Maximum gene VL-score (exclusive)", rich_help_panel="Scoring"),
    gene_v_cutoff: float = typer.Option(10.0, "--gene-v-cutoff", min=0.0, help="Maximum gene V-score (exclusive)", rich_help_panel="Scoring"),
    scoring_evalue: float = typer.Option(1e-5, "--scoring-evalue", min=0.0, help="E-value cutoff for scoring searches", rich_help_panel="Scoring"),
    phrog: bool = typer.Option(True, "--phrog/--no-phrog", help="Search PHROG as a third database, beside Pfam and KOfam", rich_help_panel="Scoring"),
    min_coverage: float = typer.Option(0.5, "--min-coverage", min=0.0, max=1.0, help="Fraction of model and protein a strict hit must align, every database (0 turns it off)", rich_help_panel="Scoring"),
    quiet: bool = typer.Option(
        False, "--quiet", help="Suppress routine progress output.", rich_help_panel="Runtime"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", help="Show additional progress details.", rich_help_panel="Runtime"
    ),
) -> None:
    """Predict and curate putative auxiliary viral genes."""
    set_verbosity(quiet=quiet, verbose=verbose)
    try:
        config = AvgConfig(
            input_contigs=input_contigs,
            output_folder=output_folder,
            threads=threads,
            mode=mode,
            min_gene_len=min_gene_len,
            min_protein_len_aa=min_protein_len_aa,
            translation_table=translation_table,
            min_amg_weight=min_amg_weight,
            filter_mode=filter_mode,
            keep_hits=keep_hits,
            scaffold_avl_cutoff=scaffold_avl_cutoff,
            gene_vl_cutoff=gene_vl_cutoff,
            gene_v_cutoff=gene_v_cutoff,
            scoring_evalue=scoring_evalue,
            phrog=phrog,
            min_coverage=min_coverage,
        )
        if quiet or verbose:
            result = run_avg(
                config, reporter=ProgressReporter(quiet=quiet, verbose=verbose)
            )
        else:
            result = run_avg(config)
    except (FileNotFoundError, ValueError, OSError, DatabaseError) as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    typer.echo(f"Cache hit: {result.prediction.cache_hit}")
    for output in result.outputs:
        typer.echo(f"Results: {output}")


@app.command("screen", rich_help_panel="Main")
def screen(
    input_contigs: Path = typer.Option(
        ...,
        "--input-contigs",
        "-i",
        exists=True,
        readable=True,
        help="Contigs to screen (FASTA).",
    ),
    hmms: list[Path] = typer.Argument(
        None,
        help="HMM files, or database names: PF00001, K00001, phrog_216, ecod:2004.1.1.5, GH128, PUL0621.",
        metavar="[QUERY]...",
        rich_help_panel="Query",
    ),
    all_cazymes: bool = typer.Option(
        False,
        "--all-cazymes",
        help="Search all dbCAN CAZyme families instead of QUERY.",
        rich_help_panel="Other ways to query",
    ),
    all_puls: bool = typer.Option(
        False,
        "--all-puls",
        help="Search for every dbCAN PUL signature instead of QUERY.",
        rich_help_panel="Other ways to query",
    ),
    hmm_mode: str = typer.Option(
        "pure",
        "--hmm-mode",
        "-M",
        help="'pure': each HMM file is one query, even with several models in it; 'mixed': each model is one query.",
        rich_help_panel="Other ways to query",
    ),
    max_evalue: Optional[float] = typer.Option(
        1e-6,
        "--max-evalue",
        "-e",
        help="Highest E-value a hit may have.",
        rich_help_panel="Hit thresholds",
    ),
    min_bitscore: Optional[float] = typer.Option(
        None,
        "--min-bitscore",
        "-b",
        help="Lowest bit score a hit may have.",
        rich_help_panel="Hit thresholds",
    ),
    min_hmm_coverage: float = typer.Option(
        0.35,
        "--min-hmm-coverage",
        min=0.0,
        max=1.0,
        help="Share of your own HMM a hit must cover (0 = off).",
        rich_help_panel="Hit thresholds",
    ),
    min_protein_coverage: Optional[float] = typer.Option(
        None,
        "--min-protein-coverage",
        min=0.0,
        max=1.0,
        help="Share of the protein a hit must cover (off by default).",
        rich_help_panel="Hit thresholds",
    ),
    cut_ga: bool = typer.Option(
        True,
        "--cut-ga/--no-cut-ga",
        help="Use Pfam's curated cutoffs (GA).",
        rich_help_panel="Hit thresholds",
    ),
    use_kofam_thresholds: bool = typer.Option(
        True,
        "--use-kofam-thresholds/--no-use-kofam-thresholds",
        help="Use KOfam's per-KO cutoffs.",
        rich_help_panel="Hit thresholds",
    ),
    use_phrog_thresholds: bool = typer.Option(
        True,
        "--use-phrog-thresholds/--no-use-phrog-thresholds",
        help="For PHROGs, need 60 bits and half of model and protein aligned.",
        rich_help_panel="Hit thresholds",
    ),
    use_ecod_thresholds: bool = typer.Option(
        True,
        "--use-ecod-thresholds/--no-use-ecod-thresholds",
        help="For ECOD families, need half of the model aligned.",
        rich_help_panel="Hit thresholds",
    ),
    combine_mode: str = typer.Option(
        "any",
        "--combine-mode",
        "-c",
        help="Keep contigs hit by 'any' query, 'all' of them, or a 'threshold'.",
        rich_help_panel="Which contigs are kept",
    ),
    min_hmm_hits: int = typer.Option(
        1,
        "--min-hmm-hits",
        "-k",
        help="With --combine-mode threshold: how many queries must hit.",
        rich_help_panel="Which contigs are kept",
    ),
    top_per_contig: int = typer.Option(
        1,
        "--top-per-contig",
        "-n",
        help="Best hits to report per contig.",
        rich_help_panel="Which contigs are kept",
    ),
    output_folder: Path = typer.Option(
        Path("phu-screen"),
        "--output-folder",
        "-o",
        help="Output folder.",
        rich_help_panel="Output",
    ),
    keep_proteins: bool = typer.Option(
        False,
        "--keep-proteins/--no-keep-proteins",
        help="Keep the predicted proteins.",
        rich_help_panel="Output",
    ),
    keep_domtbl: bool = typer.Option(
        True,
        "--keep-domtbl/--no-keep-domtbl",
        help="Keep the hmmsearch domain tables.",
        rich_help_panel="Output",
    ),
    save_target_proteins: bool = typer.Option(
        False,
        "--save-target-proteins/--no-save-target-proteins",
        help="Save the matched proteins of each model (target-proteins/).",
        rich_help_panel="Output",
    ),
    save_target_hmms: bool = typer.Option(
        False,
        "--save-target-hmms/--no-save-target-hmms",
        help="Build an HMM from each model's matched proteins (target-hmms/).",
        rich_help_panel="Output",
    ),
    trim_target_msa: str = typer.Option(
        "none",
        "--trim-target-msa",
        help="Trim alignments before building those HMMs: none, gappyout, strict, strictplus.",
        rich_help_panel="Output",
    ),
    mode: str = typer.Option(
        "meta",
        "--mode",
        "-m",
        help="'meta' for mixed contigs, 'single' for one genome.",
        rich_help_panel="Gene prediction",
    ),
    min_protein_len_aa: int = typer.Option(
        30,
        "--min-protein-len-aa",
        "-g",
        min=1,
        help="Shortest protein to keep (aa).",
        rich_help_panel="Gene prediction",
    ),
    translation_table: Optional[int] = typer.Option(
        None,
        "--ttable",
        "-T",
        min=1,
        help="Genetic code (NCBI table); by default predicted per contig.",
        rich_help_panel="Gene prediction",
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", min=1, help="Threads.", rich_help_panel="Runtime"
    ),
    quiet: bool = typer.Option(
        False, "--quiet", help="Print less.", rich_help_panel="Runtime"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", help="Print each search step.", rich_help_panel="Runtime"
    ),
):
    """
    Find contigs that carry given protein families.

    Proteins are predicted from the contigs and searched with each QUERY: your
    own HMM files, or families from the phu databases by name: Pfam (PF00001),
    KOfam (K00001), PHROG (phrog_216; find them with `phu phrogs`), ECOD
    (ecod:2004.1.1.5 for one family, or a T, H or X group such as ecod:2004.1.1
    for all its families), dbCAN (GH128) or one dbCAN PUL (PUL0621). Each
    database's own cutoffs apply.

    Examples:
        phu screen -i contigs.fa PF00001 PF00589
        phu screen -i contigs.fa my-models/*.hmm
        phu screen -i contigs.fa --combine-mode all K00001 K00002
        phu screen -i contigs.fa --all-cazymes

    Full guide: https://camilogarciabotero.github.io/phu/commands/screen/
    """
    set_verbosity(quiet=quiet, verbose=verbose)
    # Remove duplicates while preserving order
    seen = set()
    unique_hmms = []
    for p in hmms or []:
        if p not in seen:
            unique_hmms.append(p)
            seen.add(p)

    hmm_paths = unique_hmms

    if all_cazymes and all_puls:
        typer.secho(
            "--all-cazymes cannot be combined with --all-puls",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)
    if all_cazymes and hmm_paths:
        typer.secho(
            "--all-cazymes cannot be combined with positional screen queries",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)
    if all_cazymes and combine_mode != "any":
        typer.secho(
            "--combine-mode is not applicable with --all-cazymes",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)
    if all_puls and hmm_paths:
        typer.secho(
            "--all-puls cannot be combined with positional screen queries",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)
    if not all_puls and not all_cazymes and not hmm_paths:
        typer.secho("No HMM files specified", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    if trim_target_msa not in TRIM_METHODS:
        typer.secho(
            f"--trim-target-msa must be one of: {', '.join(TRIM_METHODS)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)
    if trim_target_msa != "none" and not save_target_hmms:
        typer.secho(
            "--trim-target-msa needs --save-target-hmms "
            "(and --save-target-proteins)",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)

    # Build config
    cfg = ScreenConfig(
        input_contigs=input_contigs,
        hmms=hmm_paths,
        all_puls=all_puls,
        all_cazymes=all_cazymes,
        outdir=output_folder,
        mode=mode,
        threads=threads,
        min_bitscore=min_bitscore,
        max_evalue=max_evalue,
        min_hmm_coverage=min_hmm_coverage or None,
        min_protein_coverage=min_protein_coverage or None,
        cut_ga=cut_ga,
        use_kofam_thresholds=use_kofam_thresholds,
        use_phrog_thresholds=use_phrog_thresholds,
        use_ecod_thresholds=use_ecod_thresholds,
        top_per_contig=top_per_contig,
        min_protein_len_aa=min_protein_len_aa,
        translation_table=translation_table,
        keep_proteins=keep_proteins,
        keep_domtbl=keep_domtbl,
        combine_mode=combine_mode,
        min_hmm_hits=min_hmm_hits,
        save_target_proteins=save_target_proteins,
        save_target_hmms=save_target_hmms,
        trim_target_msa=trim_target_msa,
        hmm_mode=hmm_mode,  # New parameter
    )

    try:
        _screen(cfg)
    except (OSError, DatabaseError) as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except ValueError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except CmdNotFound as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        typer.echo("Install the required executable(s) and ensure they are on PATH.")
        raise typer.Exit(1)


@app.command("jack", rich_help_panel="Main")
def jack(
    seed_marker: Optional[Path] = typer.Argument(
        None,
        exists=True,
        readable=True,
        help="Protein FASTA of seed sequences (optional with --phrog).",
        rich_help_panel="Seeds",
    ),
    input_contigs: Path = typer.Option(
        ...,
        "--input-contigs",
        "-i",
        exists=True,
        readable=True,
        help="Contigs to screen (FASTA).",
    ),
    phrog: Optional[list[str]] = typer.Option(
        None,
        "--phrog",
        help="PHROG families as seeds (repeatable): phrog_216, alias:terl, category:tail.",
        rich_help_panel="PHROG seeds",
    ),
    max_phrog_seeds: int = typer.Option(
        50,
        "--max-phrog-seeds",
        min=1,
        help="Stop if --phrog selects more PHROGs than this.",
        rich_help_panel="PHROG seeds",
    ),
    iterations: int = typer.Option(
        5,
        "--iterations",
        "-I",
        min=1,
        help="Most search rounds per seed.",
        rich_help_panel="Search",
    ),
    inc_evalue: float = typer.Option(
        1e-3,
        "--inc-evalue",
        help="E-value a hit needs to join the next round's model.",
        rich_help_panel="Search",
    ),
    max_evalue: Optional[float] = typer.Option(
        1e-5,
        "--max-evalue",
        "-e",
        help="Highest E-value a final hit may have.",
        rich_help_panel="Search",
    ),
    combine_mode: str = typer.Option(
        "any",
        "--combine-mode",
        "-c",
        help="Keep contigs hit by 'any' seed, 'all', or a 'threshold' (one --phrog query = one seed).",
        rich_help_panel="Which contigs are kept",
    ),
    min_seed_hits: int = typer.Option(
        1,
        "--min-seed-hits",
        "-k",
        min=1,
        help="With --combine-mode threshold: how many seeds must hit.",
        rich_help_panel="Which contigs are kept",
    ),
    top_per_contig: int = typer.Option(
        1,
        "--top-per-contig",
        "-n",
        min=1,
        help="Best hits to report per contig.",
        rich_help_panel="Which contigs are kept",
    ),
    output_folder: Path = typer.Option(
        Path("phu-jack"),
        "--output-folder",
        "-o",
        help="Output folder.",
        rich_help_panel="Output",
    ),
    keep_proteins: bool = typer.Option(
        False,
        "--keep-proteins/--no-keep-proteins",
        help="Keep the predicted proteins.",
        rich_help_panel="Output",
    ),
    save_hmm: bool = typer.Option(
        False,
        "--save-hmm/--no-save-hmm",
        help="Save each seed's last-round HMM.",
        rich_help_panel="Output",
    ),
    mode: str = typer.Option(
        "meta",
        "--mode",
        "-m",
        help="'meta' for mixed contigs, 'single' for one genome.",
        rich_help_panel="Gene prediction",
    ),
    min_gene_len: int = typer.Option(
        90,
        "--min-gene-len",
        "-g",
        help="Shortest gene to predict (nt).",
        rich_help_panel="Gene prediction",
    ),
    min_protein_len_aa: int = typer.Option(
        30,
        "--min-protein-len-aa",
        min=1,
        help="Shortest protein to keep (aa).",
        rich_help_panel="Gene prediction",
    ),
    translation_table: Optional[int] = typer.Option(
        None,
        "--ttable",
        "-T",
        min=1,
        help="Genetic code (NCBI table); by default predicted per contig.",
        rich_help_panel="Gene prediction",
    ),
    threads: int = typer.Option(
        1,
        "--threads",
        "-t",
        min=1,
        help="Threads.",
        rich_help_panel="Runtime",
    ),
    quiet: bool = typer.Option(
        False, "--quiet", help="Print less.", rich_help_panel="Runtime"
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Print each search step, and what --phrog queries selected.",
        rich_help_panel="Runtime",
    ),
):
    """
    Find contigs carrying relatives of seed proteins, by iterative search.

    Each seed searches the proteins predicted from the contigs; hits that pass
    --inc-evalue build the model for the next round, so distant family members
    are found round by round. Seeds are protein sequences (a FASTA file),
    PHROG families (--phrog; find them with `phu phrogs`), or both.

    Examples:
        phu jack -i contigs.fa seeds.faa
        phu jack -i contigs.fa --combine-mode all seeds.faa
        phu jack -i contigs.fa --phrog alias:terl --phrog alias:portal
        phu jack -i contigs.fa --iterations 1 --phrog phrog_216

    Full guide: https://camilogarciabotero.github.io/phu/commands/jack/
    """
    set_verbosity(quiet=quiet, verbose=verbose)
    try:
        cfg = JackConfig(
            input_contigs=input_contigs,
            seed_marker=seed_marker,
            phrog_queries=tuple(phrog or ()),
            max_phrog_seeds=max_phrog_seeds,
            outdir=output_folder,
            mode=mode,
            threads=threads,
            iterations=iterations,
            inc_evalue=inc_evalue,
            max_evalue=max_evalue,
            top_per_contig=top_per_contig,
            combine_mode=combine_mode,
            min_seed_hits=min_seed_hits,
            min_gene_len=min_gene_len,
            translation_table=translation_table,
            keep_proteins=keep_proteins,
            save_hmm=save_hmm,
            min_protein_len_aa=min_protein_len_aa,
        )
    except ValueError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    try:
        run_click_task("Iterative marker screening", _jack, cfg)
    except (FileNotFoundError, DatabaseError) as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except ValueError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
