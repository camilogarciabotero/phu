"""`phu phrogs`: find PHROG identifiers by alias, category or annotation."""

from __future__ import annotations

from pathlib import Path

import typer

from phu.db import phrog_lookup as lookup
from phu.db.phrog import phrogs_with_models

app = typer.Typer(
    help="Find PHROG identifiers by alias, functional category or annotation, "
    "for use with `phu screen`.",
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)

HEADER = ("phrog_id", "annotation", "category")


def _phrogs():
    return phrogs_with_models()


def _fail(message: str) -> None:
    typer.secho(message, fg=typer.colors.RED, err=True)
    raise typer.Exit(1)


@app.command("categories")
def categories() -> None:
    """List PHROG functional categories and how many PHROGs (with models) each has."""
    counts = lookup.category_counts(_phrogs().values())
    typer.echo("category\tphrogs")
    for name, count in sorted(counts.items(), key=lambda item: -item[1]):
        typer.echo(f"{name}\t{count}")


@app.command("aliases")
def aliases(
    check: bool = typer.Option(
        False, "--check", help="Report alias annotations that match no PHROG."
    ),
) -> None:
    """List the aliases accepted as `alias:<name>`, with their PHROG counts."""
    table = lookup.load_aliases()
    phrogs = _phrogs()
    if check:
        problems = lookup.check_aliases(table, phrogs)
        for name, missing in problems.items():
            typer.echo(f"{name}: no PHROG annotated {', '.join(map(repr, missing))}")
        if problems:
            raise typer.Exit(1)
        typer.echo("Every alias annotation matches at least one PHROG.")
        return
    typer.echo("alias\tphrogs\tdescription")
    for alias in table.aliases.values():
        selected, _ = lookup.resolve([f"alias:{alias.name}"], phrogs, table)
        typer.echo(f"{alias.name}\t{len(selected)}\t{alias.description}")
    for name, target in table.categories.items():
        typer.echo(f"{name}\t(category)\t{target}")


@app.command("search")
def search(
    text: str = typer.Argument(
        ..., help="Words that must all appear in the annotation."
    ),
    category: str | None = typer.Option(
        None, "--category", "-c", help="Only this category (matched forgivingly)."
    ),
    limit: int = typer.Option(
        50, "--limit", "-n", min=0, help="Rows to show; 0 = all."
    ),
) -> None:
    """Explore annotations: PHROGs whose annotation contains every word of TEXT."""
    phrogs = _phrogs()
    table = lookup.load_aliases()
    chosen = None
    if category:
        found = lookup.match_category(
            category, {p.category for p in phrogs.values()}, table.categories
        )
        if len(found) != 1:
            _fail(
                f"--category {category!r} matches "
                + ("; ".join(found) if found else "no category")
                + ". See `phu phrogs categories`."
            )
        chosen = found[0]
    hits = lookup.search(text, phrogs, category=chosen)
    typer.echo("\t".join(HEADER))
    for phrog in hits[: limit or None]:
        typer.echo(f"{phrog.phrog_id}\t{phrog.annotation}\t{phrog.category}")
    shown = min(len(hits), limit) if limit else len(hits)
    typer.echo(f"{len(hits)} PHROG(s) match; {shown} shown.", err=True)


@app.command("resolve")
def resolve(
    queries: list[str] = typer.Argument(
        ...,
        help="phrog_<n>, alias:<name>, category:<name> or 'annot:<annotation>'.",
    ),
    ids_only: bool = typer.Option(
        False, "--ids-only", help="Print only the identifiers, space-separated."
    ),
    output: Path | None = typer.Option(
        None, "--output", "-o", help="Also write the table (TSV) to this file."
    ),
) -> None:
    """Turn queries into PHROG identifiers for `phu screen`.

    Example:
        phu screen -i contigs.fna $(phu phrogs resolve alias:mcp alias:terl --ids-only)
    """
    try:
        selected, notes = lookup.resolve(queries, _phrogs(), lookup.load_aliases())
    except lookup.QueryError as error:
        _fail(str(error))
    for note in notes:
        typer.echo(f"Note: {note}", err=True)

    rows = [
        (r.phrog.phrog_id, r.phrog.annotation, r.phrog.category, r.query)
        for r in selected
    ]
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w") as handle:
            handle.write("\t".join((*HEADER, "query")) + "\n")
            handle.writelines("\t".join(row) + "\n" for row in rows)
    if ids_only:
        typer.echo(" ".join(row[0] for row in rows))
    else:
        typer.echo("\t".join((*HEADER, "query")))
        for row in rows:
            typer.echo("\t".join(row))
    typer.echo(f"{len(rows)} PHROG(s) selected.", err=True)
