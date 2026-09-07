# normalize-lineage

## What does it do?

The `phu normalize-lineage` command normalizes vContact lineage predictions
into compact, rank-aware lineage values without renaming the input columns. It
is intended for tables produced by vContact workflows.

## Synopsis

```bash
phu normalize-lineage -i <INPUT_FILE> -o <OUTPUT_FILE> [OPTIONS]
```

The input delimiter is inferred from the filename unless `--sep` is supplied.
The output table preserves the original columns and is written to the path
provided by `--output-file`.

## Input and output

The command recognizes lineage prediction columns by their taxonomic rank, such
as `genus_prediction`, `family_prediction`, and `order_prediction`. Missing
values remain missing. Values that cannot be parsed are preserved rather than
silently discarded.

For example:

```text
novel_genus_1_of_Viruses
Viruses:NG1
```

The normalized values remain in their original columns. With `--add-lineage`,
a compact lineage column is appended; without it, no new column is added.

## Normalization and quality checks

Existing compact values are left unchanged. The command reports the number of
examined, normalized, unchanged, missing, unparsed, rank-mismatched, and
multi-candidate cells in its QA summary.

The `--strict` option makes quality problems fatal: the command exits non-zero
when values are unparsed or when a value's rank does not match its column. In
non-strict mode, those values are retained and reported. Multiple candidates
separated by `||` are normalized independently.

## Command options

The live CLI groups options into **Input**, **Output**, **Quality checks**, and
**Runtime** panels. Run `phu normalize-lineage --help` for the current terminal
menu.

```text
Usage: phu normalize-lineage [OPTIONS]

 Normalize vContact lineage predictions without changing column names.

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --help  -h        Show this message and exit.                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Input ──────────────────────────────────────────────────────────────────────╮
│ *  --input-file  -i      <path>  Input delimited table [required]            │
│    --sep         -s      <str>   Explicit input/output separator             │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ *  --output-file  -o      <path>  Output normalized table [required]         │
│    --add-lineage  -a              Append the compact lineage column          │
│    --lineage-col  -l      <str>   Name of the lineage column                 │
│                                   [default: compact_lineage]                 │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Quality checks ─────────────────────────────────────────────────────────────╮
│ --strict          Fail on unparsed or rank-mismatched values                 │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --quiet          Suppress the QA summary                                     │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Examples

Normalize a vContact table and append a lineage column:

```bash
phu normalize-lineage \
  --input-file final_assignments.tsv \
  --output-file normalized_assignments.tsv \
  --add-lineage
```

Choose a custom lineage column name:

```bash
phu normalize-lineage \
  -i assignments.csv \
  -o normalized.csv \
  --add-lineage \
  --lineage-col best_taxonomy
```

Require clean, rank-consistent input:

```bash
phu normalize-lineage \
  -i assignments.csv \
  -o normalized.csv \
  --strict
```

## Workflow integration

```bash
vcontact3 --nucleotide viral-genomes.fasta --output-dir vcontact-output
phu normalize-lineage \
  -i vcontact-output/final_assignments.csv \
  -o normalized_assignments.csv \
  --add-lineage
```

The normalized table can then be used for taxonomy filtering, visualization,
and downstream analysis while retaining the source columns.
