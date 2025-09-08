from __future__ import annotations
from pathlib import Path
from typing import Optional
import typer

from phu import __version__
from .cluster import ClusterConfig, Mode, _cluster, parse_vclust_params
from .simplify_vcontact_taxa import TaxaConfig, _simplify_taxa
from ._exec import CmdNotFound

app = typer.Typer(
    help="Phage utilities CLI",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,
    no_args_is_help=True
)

@app.callback(invoke_without_command=True)
def _root(ctx: typer.Context) -> None:
    # any global init here (env checks, logging setup, etc.)
    if ctx.invoked_subcommand is None and not ctx.resilient_parsing:
        typer.echo(ctx.get_help())
        raise typer.Exit(0)  # exit code 0 when no subcommand is given

@app.command("cluster")
def cluster(
    mode: Mode = typer.Option(
        ..., "--mode", help="dereplication | votu | species"
    ),
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", exists=True, readable=True, help="Input FASTA"
    ),
    output_folder: Path = typer.Option(
        Path("clustered-contigs"), "--output-folder", help="Output directory"
    ),
    threads: int = typer.Option(
        0, "--threads", min=0, help="0=all cores; otherwise N threads"
    ),
    vclust_params: Optional[str] = typer.Option(
        None,
        "--vclust-params",
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
        mode=mode,
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
        ..., "--input", "-i", exists=True, readable=True, help="Input vContact final_assignments.csv"
    ),
    output_file: Path = typer.Option(
        ..., "--output", "-o", help="Output path (.csv or .tsv)"
    ),
    add_lineage: bool = typer.Option(
        False, "--add-lineage", help="Append compact_lineage column from deepest simplified rank"
    ),
    lineage_col: str = typer.Option(
        "compact_lineage", "--lineage-col", help="Name of the lineage column"
    ),
    sep: Optional[str] = typer.Option(
        None, "--sep", help="Override delimiter: ',' or '\\t'. Auto-detected from extension if not set"
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

def main() -> None:
    app()

if __name__ == "__main__":
    main()
