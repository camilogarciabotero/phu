# normalize-lineage

Normalize vContact lineage predictions into compact, rank-aware lineage codes without renaming input columns.

## Usage

```bash
phu normalize-lineage \
  --input-file final_assignments.tsv \
  --output-file normalized_assignments.tsv \
  --add-lineage
```

The input is read as a delimited table. The separator is inferred by default; use `--sep` when an explicit delimiter is required. The command writes the normalized table to `--output-file` and preserves the original columns.

Lineage values are normalized in columns whose names identify a taxonomic rank, such as `genus_prediction`. For example, `novel_genus_1_of_Viruses` becomes `Viruses:NG1`. Existing compact values remain unchanged. Missing values remain missing, and values that cannot be parsed are preserved.

Use `--lineage-col` to choose the column used for the optional compact lineage output. With `--add-lineage`, that column is appended to the table, or replaced when it already exists. Without this flag, no new column is added.

```bash
phu normalize-lineage -i assignments.csv -o normalized.csv \
  --lineage-col compact_lineage --add-lineage
```

The `--strict` option makes quality problems fatal. The command exits non-zero when values are unparsed or when a value's taxonomic rank does not match its column. Without `--strict`, these values are retained and reported in the QA summary. Use `--quiet` to suppress that summary.

## Options

| Option | Description |
| --- | --- |
| `-i, --input-file` | Input delimited table. |
| `-o, --output-file` | Output normalized table. |
| `-a, --add-lineage` | Append the compact lineage column. |
| `-l, --lineage-col` | Name for the compact lineage column. |
| `-s, --sep` | Explicit input/output separator. |
| `--strict` | Fail on unparsed or rank-mismatched values. |
| `--quiet` | Suppress the QA summary. |
