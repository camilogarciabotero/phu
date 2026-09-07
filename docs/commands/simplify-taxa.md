# simplify-taxa

## What does it do?

The `phu simplify-taxa` command converts verbose vContact taxonomy prediction
values into compact lineage codes for downstream analysis. It is designed for
vContact3 `final_assignments.csv` or TSV-style output.

## Synopsis

```bash
phu simplify-taxa -i <INPUT_FILE> -o <OUTPUT_FILE> [OPTIONS]
```

Input and output formats are selected from their file extensions. Use `--sep`
when the input delimiter cannot be inferred from the filename.

## Input and output

The command processes columns ending in `_prediction`, including:

- `realm_prediction`
- `kingdom_prediction`
- `phylum_prediction`
- `class_prediction`
- `order_prediction`
- `family_prediction`
- `subfamily_prediction`
- `genus_prediction`

Other columns and their order are preserved. Values that do not match a
supported vContact pattern are preserved. Multiple candidate values separated
by `||` are processed independently.

CSV output is used unless the output filename ends in `.tsv`; TSV output is
used in that case. With `--add-lineage`, one additional column is appended.

## Transformation logic

A value such as:

```text
novel_genus_1_of_novel_family_2_of_Caudoviricetes
```

becomes:

```text
Caudoviricetes:NF2:NG1
```

Compact codes use these rank prefixes:

- `NK`: novel kingdom
- `NP`: novel phylum
- `NC`: novel class
- `NO`: novel order
- `NF`: novel family
- `NSF`: novel subfamily
- `NG`: novel genus

Some vContact zero-index chains have explicit compatibility handling. Validate
those version-sensitive cases against representative vContact output before
using them as a scientific contract.

## Command options

```text
Usage: phu simplify-taxa [OPTIONS]

 Simplify vContact taxonomy prediction columns into compact lineage codes.

 Transforms verbose vContact taxonomy strings like
 'novel_genus_1_of_novel_family_2_of_Caudoviricetes'
 into compact codes like 'Caudoviricetes:NF2:NG1'.

 Example:
   phu simplify-taxa -i final_assignments.csv -o simplified.csv --add-lineage

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --input-file   -i      <path>  Input vContact final_assignments.csv       │
│                                   [required]                                 │
│ *  --output-file  -o      <path>  Output file path (.csv or .tsv) [required] │
│    --add-lineage  -a              Append compact_lineage column from deepest │
│                                   simplified rank                            │
│    --lineage-col  -l      <str>   Name of the lineage column                 │
│                                   [default: compact_lineage]                 │
│    --sep          -s      <str>   Override delimiter: ',' or '\t'.           │
│                                   Auto-detected from extension if not set    │
│    --quiet                        Suppress routine progress output.          │
│    --verbose                      Show additional progress details.          │
│    --help         -h              Show this message and exit.                │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Examples

Simplify a vContact CSV file:

```bash
phu simplify-taxa \
  --input-file final_assignments.csv \
  --output-file simplified.csv
```

Process TSV input and append the deepest available lineage:

```bash
phu simplify-taxa \
  -i final_assignments.tsv \
  -o simplified.tsv \
  --add-lineage \
  --lineage-col best_taxonomy
```

Override delimiter detection:

```bash
phu simplify-taxa -i assignments.txt -o simplified.csv --sep $'\t'
```

## Lineage column

The `--add-lineage` option appends a column containing the deepest available
simplified rank. The default column name is `compact_lineage`; use
`--lineage-col` to choose another name.

The priority order is:

1. `genus_prediction`
2. `subfamily_prediction`
3. `family_prediction`
4. `order_prediction`
5. `class_prediction`
6. `phylum_prediction`
7. `kingdom_prediction`
8. `realm_prediction`

## Workflow integration

```bash
vcontact3 --nucleotide viral-genomes.fasta --output-dir vcontact-output
phu simplify-taxa \
  -i vcontact-output/final_assignments.csv \
  -o taxonomy_simplified.csv \
  --add-lineage
```

The resulting table can be used for taxonomy filtering, visualization, and
other downstream analyses while retaining the original input fields.
