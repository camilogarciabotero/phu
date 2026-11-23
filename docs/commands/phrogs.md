# phrogs

Find PHROG identifiers by alias, functional category or annotation, and pass
them to [`phu screen`](screen.md), which accepts `phrog_<n>` identifiers as
queries.

```bash
phu phrogs categories                 # categories and how many PHROGs each has
phu phrogs aliases                    # short names such as mcp, terl, portal
phu phrogs search "tail fiber"        # explore annotations
phu phrogs resolve alias:mcp category:lysis --ids-only

# Screen contigs for major capsid and terminase large subunit families
phu screen -i contigs.fna $(phu phrogs resolve alias:mcp alias:terl --ids-only)

# The same queries seed iterative searches in jack
phu jack -i contigs.fna --phrog alias:portal --max-phrog-seeds 60
```

Only the 38,880 PHROG families with a model (`phrog_1` to `phrog_38880`) are
listed, so every identifier printed can be searched. About 5,100 of them have
a known function; the rest are in the category "unknown function". The PHROG
annotation table also numbers entries above 38,880, almost all unannotated
and without an alignment, apparently proteins that formed no family; they
cannot be searched and `phu phrogs` leaves them out. The first use indexes
the PHROG models, which takes a few minutes; `phu screen phrog_<n>` uses the
same index.

## Queries for `resolve`

| Query | Selects |
|---|---|
| `phrog_216` | one PHROG |
| `alias:mcp` | a curated group of annotations (see below) |
| `category:tail` | every PHROG of a functional category |
| `'annot:portal protein'` | PHROGs whose annotation is exactly this text, ignoring letter case |

Queries can be mixed; each PHROG is listed once, under the first query that
selected it. A query that selects nothing stops with an explanation.

`--ids-only` prints the identifiers on one line for the shell; `-o file.tsv`
also writes the table (`phrog_id`, `annotation`, `category`, `query`), a record
of what each name meant in that run.

## Searching annotations

`phu phrogs search` is for exploring. It lists PHROGs whose annotation matches
every word given, in any order, and is forgiving; the `match` column says how
each row matched, and rows are ordered from the most literal match:

| `match` | Example query | Finds |
|---|---|---|
| `exact` | `large terminase`, `baseplat` | each word appears in the annotation |
| `alias:<name>` | `mcp`, `TerL` | the annotations of that alias |
| `synonym:<term>` | `capsid`, `MTase`, `RNR`, `ssb`, `psbA` | PHROG's wording for the term: head protein, methyltransferase, ribonucleotide reductase, single strand DNA binding, photosystem II D1 |
| `variant` | `tail fibre`, `holins`, `tapemeasure`, `antirepressor` | British spelling, plurals, words written together or apart |
| `typo` | `endolysine`, `intergrase` | a close spelling of a word of five letters or more |

Synonyms are listed in `phu/data/phrog_aliases.json` and apply to `search`
only. `resolve` stays literal, apart from category names, so that a run
records exactly what was selected: explore with `search`, then select with
`resolve`.

## Categories

PHROGs v4 has nine functional categories plus "unknown function": tail; head
and packaging; connector; DNA, RNA and nucleotide metabolism; transcription
regulation; lysis; integration and excision; moron, auxiliary metabolic gene
and host takeover; other.

Category names are matched forgivingly, and phu says how it read the name:

| You write | Read as |
|---|---|
| `tail`, `TAIL`, `tails` | tail |
| `head`, `pack` | head and packaging |
| `dna` | DNA, RNA and nucleotide metabolism |
| `moron`, `amg`, `host takeover` | moron, auxiliary metabolic gene and host takeover |
| `lysys` | lysis |

A name that fits several categories (`category:a`) is refused with the
candidates listed.

## Aliases

Aliases are short names for groups of PHROG annotations, stored in
`phu/data/phrog_aliases.json`. Each lists exact PHROG v4 annotation strings,
so `mcp` selects "major head protein", "major capsid protein" and "major head
and protease protein" but not "minor head protein". Run `phu phrogs aliases`
for the current list and counts, and `phu phrogs aliases --check` to confirm
every listed annotation still exists in the installed PHROG table.

| Alias | Annotations |
|---|---|
| `mcp` | major head protein; major capsid protein; major head and protease protein |
| `terl`, `ters` | terminase large / small subunit |
| `portal`, `scaffold`, `head-tail` | portal; head scaffolding; head-tail connector and adaptors |
| `neck`, `collar` | neck and tail-completion proteins; collar proteins |
| `tmp`, `tail-sheath`, `tail-tube`, `baseplate` | tail structure |
| `tail-fiber`, `tail-spike` | tail fibers (not assembly chaperones); tail spikes and receptor-binding tail proteins |
| `holin`, `endolysin`, `spanin` | lysis |
| `integrase`, `excisionase`, `anti-repressor` | integration and regulation |

`moron` and `amg` are shortcuts for the category "moron, auxiliary metabolic
gene and host takeover".

The aliases follow PHROG's own annotation wording; they are a convenience, not
a classification. Annotation text describes a family's best-known member, so
review hits before drawing functional conclusions.

## Source

Terzian P, et al. PHROG: families of prokaryotic virus proteins clustered
using remote homology. NAR Genom Bioinform. 2021;3(3):lqab067.
