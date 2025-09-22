from __future__ import annotations
from pathlib import Path
from typing import List, Optional
import typer

from phu import __version__
from .cluster import ClusterConfig, Mode, _cluster, parse_vclust_params
from .simplify_vcontact_taxa import TaxaConfig, _simplify_taxa
from .screen import ScreenConfig, _screen
from ._exec import CmdNotFound

app = typer.Typer(
    help="Phage utilities CLI",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=True,
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
        0, "--threads","-t", min=0, help="0=all cores; otherwise N threads"
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
        ..., "--input-file", "-i", exists=True, readable=True, help="Input vContact final_assignments.csv"
    ),
    output_file: Path = typer.Option(
        ..., "--output-file", "-o", help="Output file path (.csv or .tsv)"
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



@app.command("screen")
def screen(
    input_contigs: Path = typer.Option(
        ..., "--input-contigs", "-i", exists=True, readable=True, help="Input contigs FASTA"
    ),
    hmms: List[Path] = typer.Argument(
        ..., help="HMM files (supports wildcards like *.hmm)"
    ),
    output_folder: Path = typer.Option(
        Path("phu-screen"), "--output-folder", "-o", help="Output directory"
    ),
    mode: str = typer.Option(
        "meta", "--mode", help="pyrodigal mode: meta|single"
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", min=1, help="Threads for both pyrodigal and hmmsearch"
    ),
    min_bitscore: Optional[float] = typer.Option(
        None, "--min-bitscore", help="Minimum bitscore to keep a domain hit"
    ),
    max_evalue: Optional[float] = typer.Option(
        1e-5, "--max-evalue", help="Maximum independent E-value to keep a domain hit"
    ),
    top_per_contig: int = typer.Option(
        1, "--top-per-contig", help="Keep top-N hits per contig (by bitscore)"
    ),
    min_gene_len: int = typer.Option(
        90, "--min-gene-len", help="Minimum gene length for pyrodigal (nt)"
    ),
    translation_table: int = typer.Option(
        11, "--ttable", help="NCBI translation table for coding sequences"
    ),
    keep_proteins: bool = typer.Option(
        False, "--keep-proteins/--no-keep-proteins", help="Keep the protein FASTA used for searching"
    ),
    keep_domtbl: bool = typer.Option(
        True, "--keep-domtbl/--no-keep-domtbl", help="Keep raw domtblout from hmmsearch"
    ),
    combine_mode: str = typer.Option(
        "any", "--combine-mode", help="How to combine hits from multiple HMMs: any|all|threshold"
    ),
    min_hmm_hits: int = typer.Option(
        1, "--min-hmm-hits", help="Minimum number of HMMs that must hit a contig (for threshold mode)"
    ),
    save_target_proteins: bool = typer.Option(
        False, "--save-target-proteins/--no-save-target-proteins", 
        help="Save matched proteins per HMM model in target_proteins/ subfolder"
    ),
    save_target_hmms: bool = typer.Option(
        False, "--save-target-hmms/--no-save-target-hmms", help="Save HMMs built from target proteins in target_hmms/ subfolder"
    ),
    hmm_mode: str = typer.Option(
        "pure", "--hmm-mode", help="HMM file type: 'pure' (one model per file) or 'mixed' (pressed/concatenated HMMs)"
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
        top_per_contig=top_per_contig,
        min_gene_len=min_gene_len,
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
            "Required executables on PATH: 'hmmsearch' and 'seqkit'"
        )
        raise typer.Exit(1)



def main() -> None:
    app()

if __name__ == "__main__":
    main()
