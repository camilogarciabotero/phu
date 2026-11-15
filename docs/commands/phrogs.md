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
```

Only PHROGs with a model in the phu-db bundle are listed (38,880 of the
109,404 PHROGs; the others are singletons without an alignment), so every
identifier printed can be searched. Of those, about 5,100 have a known
function; the rest are in the category "unknown function". The first use
indexes the PHROG models, which takes a few minutes; `phu screen phrog_<n>`
uses the same index.

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
