from __future__ import annotations
from pathlib import Path
from typing import List, Optional
import typer

from phu import __version__
from ._exec import CmdNotFound
from .cluster import ClusterConfig, Mode, _cluster, parse_vclust_params
from .gene_prediction_core import clean_prediction_cache
from .jack import JackConfig, _jack
from .pfam_db import (
    get_pfam_database_status,
    prepare_pfam_database,
    refresh_pfam_database,
    remove_pfam_database,
)
from .screen import ScreenConfig, _screen
from .simplify_vcontact_taxa import TaxaConfig, _simplify_taxa

app = typer.Typer(
    help="Phage utilities CLI",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=True,
    no_args_is_help=True
)
dbs_app = typer.Typer(
    help="Manage local phu databases",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
)
app.add_typer(dbs_app, name="dbs", rich_help_panel="Database Management")

SUPPORTED_DBS = ("pfam",)


def _normalize_db_names(databases: List[str], all_dbs: bool) -> List[str]:
    if all_dbs and databases:
        raise ValueError("Use either specific database names or --all, not both")

    selected = list(SUPPORTED_DBS) if all_dbs or not databases else [name.lower() for name in databases]

    unknown = [name for name in selected if name not in SUPPORTED_DBS]
    if unknown:
        raise ValueError(f"Unsupported database(s): {', '.join(sorted(set(unknown)))}")

    # Deduplicate while preserving order.
    return list(dict.fromkeys(selected))


def _db_status_payload(db_name: str) -> dict:
    if db_name == "pfam":
        return get_pfam_database_status()
    raise ValueError(f"Unsupported database: {db_name}")


@dbs_app.command("list")
def dbs_list() -> None:
    """List supported databases and quick readiness state."""
    for db_name in SUPPORTED_DBS:
        status = _db_status_payload(db_name)
        state = "ready" if (status.get("downloaded") and status.get("indexed")) else "not-ready"
        typer.echo(f"{db_name}\t{state}")


@dbs_app.command("status")
def dbs_status(
    databases: Optional[List[str]] = typer.Argument(
        None,
        help="Database names (default: all)",
    ),
    all_dbs: bool = typer.Option(
        False,
        "--all",
        help="Show status for all supported databases",
    ),
) -> None:
    """Show detailed status for one or more databases."""
    try:
        selected = _normalize_db_names(databases or [], all_dbs=all_dbs)
    except ValueError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    for db_name in selected:
        status = _db_status_payload(db_name)
        typer.echo(f"{db_name}:")
        typer.echo(f"  downloaded: {status.get('downloaded')}")
        typer.echo(f"  indexed: {status.get('indexed')}")
        typer.echo(f"  manifest_exists: {status.get('manifest_exists')}")
        typer.echo(f"  model_count: {status.get('model_count')}")
        typer.echo(f"  sparse_cached_count: {status.get('sparse_cached_count')}")
        typer.echo(f"  root: {status.get('root')}")


@dbs_app.command("prepare")
def dbs_prepare(
    databases: Optional[List[str]] = typer.Argument(
        None,
        help="Database names (default: all)",
    ),
    all_dbs: bool = typer.Option(
        False,
        "--all",
        help="Prepare all supported databases",
    ),
    force_refresh: bool = typer.Option(
        False,
        "--force-refresh",
        help="Force re-download/rebuild when supported",
    ),
) -> None:
    """Prepare databases for use by phu commands."""
    try:
        selected = _normalize_db_names(databases or [], all_dbs=all_dbs)
    except ValueError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    try:
        for db_name in selected:
            if db_name == "pfam":
                result = prepare_pfam_database(download=True, index=True, force_refresh=force_refresh)
                typer.echo(f"Prepared {db_name}: {result.get('hmm_path')}")
                typer.echo(f"Index ready: {result.get('offsets_path')}")
    except FileNotFoundError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


@dbs_app.command("refresh")
def dbs_refresh(
    databases: Optional[List[str]] = typer.Argument(
        None,
        help="Database names (default: all)",
    ),
    all_dbs: bool = typer.Option(
        False,
        "--all",
        help="Refresh all supported databases",
    ),
) -> None:
    """Refresh database integrity and repair incomplete state."""
    try:
        selected = _normalize_db_names(databases or [], all_dbs=all_dbs)
    except ValueError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    try:
        for db_name in selected:
            if db_name == "pfam":
                result = refresh_pfam_database()
                typer.echo(f"Refreshed {db_name}: {result.get('hmm_path')}")
                typer.echo(f"Index ready: {result.get('offsets_path')}")
    except FileNotFoundError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


@dbs_app.command("remove")
def dbs_remove(
    databases: Optional[List[str]] = typer.Argument(
        None,
        help="Database names (default: all)",
    ),
    all_dbs: bool = typer.Option(
        False,
        "--all",
        help="Remove all supported databases",
    ),
    yes: bool = typer.Option(
        False,
        "--yes",
        help="Confirm database removal",
    ),
) -> None:
    """Remove local database data."""
    if not yes:
        typer.secho("Refusing to remove databases without --yes", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    try:
        selected = _normalize_db_names(databases or [], all_dbs=all_dbs)
    except ValueError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)

    for db_name in selected:
        if db_name == "pfam":
            removed = remove_pfam_database()
            if removed:
                typer.echo(f"Removed {db_name} database")
            else:
                typer.echo(f"{db_name} database not present")


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
        help="Remove cached protein predictions and exit.",
    ),
) -> None:
    if clean_cache:
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

@app.command("cluster")
def cluster(
    mode: str = typer.Option(
        ..., "--mode", help="dereplication | votu | species"
    ),
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", "-i", exists=True, readable=True, help="Input FASTA"
    ),
    output_folder: Path = typer.Option(
        Path("clustered-contigs"), "--output-folder", "-o", help="Output directory"
    ),
    threads: int = typer.Option(
        0, "--threads","-t", min=0, help="0=all cores; otherwise N threads"
    ),
    vclust_params: Optional[str] = typer.Option(
        None,
        "--vclust-params", "-p",
        help='Custom vclust parameters: "--min-kmers 20 --outfmt lite --ani 0.97"'
    ),
):
    """
    Sequence clustering wrapper around external 'vclust' with three modes.
    
    For advanced usage, provide custom vclust parameters as a quoted string.
    See the vclust wiki for parameter details: https://github.com/refresh-bio/vclust/wiki
    
    Example:
        phu cluster --mode votu --input-contigs genomes.fna --vclust-params="--min-kmers 20 --outfmt lite"
    """
    # Convert CLI string to enum for compatibility.
    try:
        mode_enum = Mode(mode)
    except ValueError:
        typer.secho(
            "Invalid mode. Use one of: dereplication, votu, species",
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
                f"Error parsing vclust parameters: {e}",
                fg=typer.colors.RED,
                err=True
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
        _cluster(cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except CmdNotFound as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        typer.echo(
            "Required executables on PATH: 'vclust' (or 'vclust.py') and 'seqkit'"
        )
        raise typer.Exit(1)

@app.command("simplify-taxa")
def simplify_taxa(
    input_file: Path = typer.Option(
        ..., "--input-file", "-i", exists=True, readable=True, help="Input vContact final_assignments.csv"
    ),
    output_file: Path = typer.Option(
        ..., "--output-file", "-o", help="Output file path (.csv or .tsv)"
    ),
    add_lineage: bool = typer.Option(
        False, "--add-lineage", "-a", help="Append compact_lineage column from deepest simplified rank"
    ),
    lineage_col: str = typer.Option(
        "compact_lineage", "--lineage-col", "-l", help="Name of the lineage column"
    ),
    sep: Optional[str] = typer.Option(
        None, "--sep", "-s", help="Override delimiter: ',' or '\\t'. Auto-detected from extension if not set"
    ),
):
    """
    Simplify vContact taxonomy prediction columns into compact lineage codes.
    
    Transforms verbose vContact taxonomy strings like 'novel_genus_1_of_novel_family_2_of_Caudoviricetes'
    into compact codes like 'Caudoviricetes:NF2:NG1'.
    
    Example:
        phu simplify-taxa -i final_assignments.csv -o simplified.csv --add-lineage
    """
    # Build config
    cfg = TaxaConfig(
        input_file=input_file,
        output_file=output_file,
        add_lineage=add_lineage,
        lineage_col=lineage_col,
        sep=sep,
    )
    
    try:
        _simplify_taxa(cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except RuntimeError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)



@app.command("screen")
def screen(
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", "-i", exists=True, readable=True, help="Input contigs FASTA"
    ),
    hmms: List[Path] = typer.Argument(
        ..., help="HMM files and/or PFAM accessions (e.g. PF00001; supports wildcards like *.hmm)"
    ),
    output_folder: Path = typer.Option(
        Path("phu-screen"), "--output-folder", "-o", help="Output directory"
    ),
    mode: str = typer.Option(
        "meta", "--mode", "-m", help="pyrodigal mode: meta|single"
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", min=1, help="Threads for both pyrodigal and pyhmmer"
    ),
    min_bitscore: Optional[float] = typer.Option(
        None, "--min-bitscore", "-b", help="Minimum bitscore to keep a domain hit"
    ),
    max_evalue: Optional[float] = typer.Option(
        1e-5, "--max-evalue", "-e", help="Maximum independent E-value to keep a domain hit"
    ),
    cut_ga: bool = typer.Option(
        False,
        "--cut-ga/--no-cut-ga",
        help="Use model GA gathering cutoffs during HMM search (PFAM-style thresholding)",
    ),
    top_per_contig: int = typer.Option(
        1, "--top-per-contig", "-n", help="Keep top-N hits per contig (by bitscore)"
    ),
    min_protein_len_aa: int = typer.Option(
        30, "--min-protein-len-aa", "-g", min=1, help="Minimum translated protein length to keep (aa)"
    ),
    translation_table: int = typer.Option(
        11, "--ttable", "-T", help="NCBI translation table for coding sequences"
    ),
    keep_proteins: bool = typer.Option(
        False, "--keep-proteins/--no-keep-proteins", help="Keep the protein FASTA used for searching"
    ),
    keep_domtbl: bool = typer.Option(
        True, "--keep-domtbl/--no-keep-domtbl", help="Keep raw domtblout from hmmsearch"
    ),
    combine_mode: str = typer.Option(
        "any", "--combine-mode", "-c", help="How to combine hits from multiple HMMs: any|all|threshold"
    ),
    min_hmm_hits: int = typer.Option(
        1, "--min-hmm-hits", "-k", help="Minimum number of HMMs that must hit a contig (for threshold mode)"
    ),
    save_target_proteins: bool = typer.Option(
        False, "--save-target-proteins/--no-save-target-proteins", 
        help="Save matched proteins per HMM model in target_proteins/ subfolder"
    ),
    save_target_hmms: bool = typer.Option(
        False, "--save-target-hmms/--no-save-target-hmms", help="Save HMMs built from target proteins in target_hmms/ subfolder"
    ),
    hmm_mode: str = typer.Option(
        "pure", "--hmm-mode", "-M", help="HMM file type: 'pure' (one model per file) or 'mixed' (pressed/concatenated HMMs)"
    )
):
    """
    Screen contigs for protein families using HMMER on predicted CDS.
    
    Supports multiple HMM files with different combination modes:
    - any: Keep contigs matching any HMM (default, most permissive)
    - all: Keep contigs matching all HMMs (most restrictive) 
    - threshold: Keep contigs matching at least --min-hmm-hits HMMs
    
    HMM modes:
    - pure: Each HMM file contains one model (default, most common)
    - mixed: HMM files contain multiple models (pressed/concatenated HMMs)
    
    Examples:
        phu screen -i contigs.fa *.hmm
        phu screen -i contigs.fa PF00001 PF00589
        phu screen -i contigs.fa PF00001 PF00589 --cut-ga
        phu screen -i contigs.fa --combine-mode all file1.hmm file2.hmm file3.hmm
        phu screen -i contigs.fa --combine-mode threshold --min-hmm-hits 5 pfam_database.hmm
        phu screen -i contigs.fa --save-target-proteins *.hmm
    """
    # Remove duplicates while preserving order
    seen = set()
    unique_hmms = []
    for p in hmms:
        if p not in seen:
            unique_hmms.append(p)
            seen.add(p)
    
    hmm_paths = unique_hmms
    
    if not hmm_paths:
        typer.secho("No HMM files specified", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    
    # Build config
    cfg = ScreenConfig(
        input_contigs=input_contigs,
        hmms=hmm_paths,
        outdir=output_folder,
        mode=mode,
        threads=threads,
        min_bitscore=min_bitscore,
        max_evalue=max_evalue,
        cut_ga=cut_ga,
        top_per_contig=top_per_contig,
        min_protein_len_aa=min_protein_len_aa,
        translation_table=translation_table,
        keep_proteins=keep_proteins,
        keep_domtbl=keep_domtbl,
        combine_mode=combine_mode,
        min_hmm_hits=min_hmm_hits,
        save_target_proteins=save_target_proteins,
        save_target_hmms=save_target_hmms,
        hmm_mode=hmm_mode,  # New parameter
    )
    
    try:
        _screen(cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except CmdNotFound as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        typer.echo(
            "Required executables on PATH: 'seqkit'"
        )
        raise typer.Exit(1)


@app.command("jack")
def jack(
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", "-i", exists=True, readable=True, help="Input contigs FASTA"
    ),
    seed_marker: Path = typer.Argument(
        ..., exists=True, readable=True, help="Seed marker protein FASTA (supports one or more sequences)"
    ),
    output_folder: Path = typer.Option(
        Path("phu-jack"), "--output-folder", "-o", help="Output directory"
    ),
    mode: str = typer.Option(
        "meta", "--mode", "-m", help="pyrodigal mode: meta|single"
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", min=1, help="Threads for both pyrodigal and pyhmmer"
    ),
    iterations: int = typer.Option(
        5, "--iterations", "-I", min=1, help="Maximum jackhmmer iterations"
    ),
    inc_evalue: float = typer.Option(
        1e-3, "--inc-evalue", help="Inclusion E-value threshold for iterative jackhmmer"
    ),
    max_evalue: Optional[float] = typer.Option(
        1e-5, "--max-evalue", "-e", help="Maximum independent E-value to keep a final hit"
    ),
    top_per_contig: int = typer.Option(
        1, "--top-per-contig", "-n", min=1, help="Keep top-N hits per contig (by bitscore)"
    ),
    combine_mode: str = typer.Option(
        "any", "--combine-mode", "-c", help="How to combine hits from multiple seed proteins: any|all|threshold"
    ),
    min_seed_hits: int = typer.Option(
        1, "--min-seed-hits", "-k", min=1, help="Minimum number of seeds that must hit a contig (for threshold mode)"
    ),
    min_gene_len: int = typer.Option(
            90, "--min-gene-len", "-g", help="Minimum gene length for pyrodigal (nt)"
        ),
    min_protein_len_aa: int = typer.Option(
        30, "--min-protein-len-aa", min=1, help="Minimum translated protein length to keep (aa)"
    ),
    translation_table: int = typer.Option(
        11, "--ttable", "-T", help="NCBI translation table for coding sequences"
    ),
    keep_proteins: bool = typer.Option(
        False, "--keep-proteins/--no-keep-proteins", help="Keep the protein FASTA used for searching"
    ),
    save_hmm: bool = typer.Option(
        False,
        "--save-hmm/--no-save-hmm",
        help="Save the last jackhmmer iteration HMM as last_iteration.hmm",
    ),
):
    """
    Iteratively screen contigs from one or more seed protein markers with pyhmmer.jackhmmer.

    Combine modes for multi-seed screening:
    - any: keep contigs hit by at least one seed (default)
    - all: keep contigs hit by all seeds
    - threshold: keep contigs hit by at least --min-seed-hits seeds

    Examples:
        phu jack -i contigs.fa marker_seed.faa
        phu jack -i contigs.fa --combine-mode all marker_seeds.faa
        phu jack -i contigs.fa --iterations 7 --inc-evalue 1e-4 marker_seed.faa
    """
    cfg = JackConfig(
        input_contigs=input_contigs,
        seed_marker=seed_marker,
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

    try:
        _jack(cfg)
    except FileNotFoundError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except ValueError as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except CmdNotFound as e:
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        typer.echo(
            "Required executables on PATH: 'seqkit'"
        )
        raise typer.Exit(1)



def main() -> None:
    app()

if __name__ == "__main__":
    main()
