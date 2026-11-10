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

Each `novel_<rank>_<n>` node becomes `N` plus the rank letter and its number
(`NR`, `NK`, `NP`, `NC`, `NO`, `NF`, `NSF`, `NG`). The anchor comes first and
the codes follow from the highest rank to the lowest.

### Unplaced ranks and the `default` anchor

vContact3 also writes `unplaced_<rank>` nodes, joined by `_within_`. They carry
no number, so each one becomes `U` plus the same rank letter (`UR`, `UK`, `UP`,
`UC`, `UO`, `UF`, `USF`, `UG`). The placeholder anchor `default` is written as
`DEF`; any other anchor is kept as written.

| Input | Output |
| --- | --- |
| `unplaced_phylum_within_unplaced_kingdom_within_default` | `DEF:UK:UP` |
| `unplaced_class_within_unplaced_phylum_within_unplaced_kingdom_within_default` | `DEF:UK:UP:UC` |
| `unplaced_order_within_Caudoviricetes` | `Caudoviricetes:UO` |
| `novel_genus_1_of_novel_family_2_of_unplaced_order_within_Caudoviricetes` | `Caudoviricetes:UO:NF2:NG1` |

A chain may mix `novel_` and `unplaced_` nodes, and either `_of_` or `_within_`
may follow either kind of node. Every node must be of a strictly lower rank
than the next one toward the anchor. A chain that breaks that order, has no
anchor, or names an unknown rank is reported as unparsed and left as it was.
Compact values that already use `U` codes are left unchanged, so running the
command on its own output changes nothing.

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
