# screen

`phu screen` finds the contigs that carry the protein families you ask for.
It predicts the proteins of your contigs, searches them with profile HMMs, and
keeps the contigs whose proteins match.

Use it to pull viral contigs out of a metagenomic assembly by their marker
proteins, to find contigs that carry a gene of interest, or to collect the
proteins behind a family for a new model.

```bash
phu screen -i contigs.fa PF00589                  # one Pfam family
phu screen -i contigs.fa phrog_216 ecod:2004.1    # a PHROG and an ECOD homology group
phu screen -i contigs.fa my-models/*.hmm          # your own HMMs
```

The matching contigs are written to `phu-screen/matched-contigs.fna`.

!!! note "Input"
    `phu screen` takes **DNA contigs** (FASTA). If you already have proteins,
    search them with pyHMMER or HMMER directly.

## What you can search for

Each query is one of your own HMM files or a family named from the
[phu databases](../databases.md). You can mix kinds in one run.

| Query | Example | What it names | When a hit counts |
| --- | --- | --- | --- |
| Your HMM file | `capsid.hmm` | The model(s) in the file | Aligns to more than 35% of the model (`--min-hmm-coverage`) |
| Pfam | `PF00589`, `PF00589.17` | A Pfam family | Pfam's curated gathering cutoff (`--cut-ga`) |
| KOfam | `K00001` | A KEGG Orthology group | The KO's own score threshold (`--use-kofam-thresholds`) |
| PHROG | `phrog_216` | A PHROG family (`phrog_1` to `phrog_38880`) | 60 bits, and half of both model and protein aligned (`--use-phrog-thresholds`) |
| ECOD | `ecod:2004.1.1.5`, `ecod:2004.1` | An ECOD family, or every family in a larger group | Half of the model aligned (`--use-ecod-thresholds`) |
| dbCAN family | `GH128`, `cbm89` | A CAZyme family (case does not matter) | Domain E-value below 1e-15 and more than 35% of the model aligned |
| dbCAN PUL | `PUL0621` | The CAZyme families of one PUL | Every family found on the same contig |

Every hit must also have an E-value of at most `1e-6` (`--max-evalue`).
`--min-protein-coverage 0.5`, off by default, further requires every hit to
cover half of its protein, whatever the model; this drops single domains
found in long multi-domain proteins. The
first query that names a database family installs the database bundle, if it
is not installed yet (see [Databases](../databases.md)). The full rules are
in [Screen thresholds](../screen-thresholds.md).

### Pfam and KOfam

Pfam accessions may carry a version (`PF00589.17` is read as `PF00589`).
Pfam models keep their curated gathering (GA) cutoffs, also when queried
together with other kinds of model; `--no-cut-ga` turns them off and leaves
only `--max-evalue` and `--min-bitscore`. KOfam models use the per-KO threshold and score type (full sequence or domain) from
KOfam's `ko_list`; `--no-use-kofam-thresholds` turns them off.

### PHROG

PHROGs are phage protein families. PHROG has no curated per-model cutoffs, so
a hit must score at least 60 bits and align to at least half of both the
model and the protein, the rule [`phu avger`](avger.md) uses.
`--no-use-phrog-thresholds` drops that rule. The function and category of
each queried PHROG are printed and recorded in `.phu/query.json`.

To pick PHROGs by name or category instead of number, use
[`phu phrogs`](phrogs.md):

```bash
phu screen -i contigs.fa $(phu phrogs resolve alias:mcp alias:terl --ids-only)
```

### ECOD

ECOD classifies protein **domains** by evolutionary relationship, at four
levels. You can query any of them with `ecod:` and a dotted id:

| Level | Example | Meaning |
| --- | --- | --- |
| F (family) | `ecod:2004.1.1.5` | One family: one model |
| T (topology) | `ecod:2004.1.1` | Every family with this fold topology |
| H (homology) | `ecod:2004.1` | Every family thought to share an ancestor |
| X (possible homology) | `ecod:2004` | Every family in the broadest group |

The bundle has one model per F-group, 32,041 in ECOD v295.2, built from
ECOD's representative domains. Group and family names are in `families.tsv`
of the installed `ecod` database, and hits name their model as
`ecod:<family id>`.

ECOD models describe single domains, often within longer proteins, so a hit
must align to at least half of the model but may cover any part of the
protein. `--no-use-ecod-thresholds` drops that rule. `.phu/query.json` records
what each `ecod:` id expanded to, and the names and Pfam accessions of every
family searched.

!!! tip "Broad groups search many models"
    `ecod:2004.1`, the P-loop homology group, expands to 556 families.
    Expect one domain table per family in the output folder; turn them off
    with `--no-keep-domtbl`.

### dbCAN families and PULs

dbCAN family names are matched against the installed dbCAN models,
whatever their case. A PUL id expands to the CAZyme families of that PUL
(from dbCAN-PUL's `cazymes_predicted_dbcan`), and by default every family must
be found on the same contig; use `--combine-mode threshold` to accept part
of the signature. Only one PUL can be queried per run, and not together with
other queries. To search all CAZymes or all PULs at once, see
[CAZymes and PULs](#cazymes-and-puls).

A PUL match means its CAZyme families occur together on a contig. It does not
show that the contig carries that PUL or predict its substrate.

## Combining queries

With several queries, `--combine-mode` decides which contigs are kept:

| Mode | Keeps a contig when | Use it for |
| --- | --- | --- |
| `any` (default) | At least one query hits | Collecting everything related to any query |
| `all` | Every query hits | Requiring a set of markers together |
| `threshold` | At least `--min-hmm-hits` queries hit | Requiring, for example, 3 of 5 markers |

Each `ecod:` id counts as **one** query, met by a hit to any of
its families: `--combine-mode all ecod:2004.1.1 PF00589` keeps contigs with
some family of that T-group and PF00589.

The same holds for your own HMM files: in the default `--hmm-mode pure`, a
file is one query even when it holds several models. To ask for "two of three
markers" when each marker has several profiles, put each marker's profiles in
one file:

```bash
cat terl/*.hmm > TerL.hmm; cat mcp/*.hmm > MCP.hmm; cat portal/*.hmm > portal.hmm
phu screen -i contigs.fa --combine-mode threshold --min-hmm-hits 2 TerL.hmm MCP.hmm portal.hmm
```

A hit to any profile in `TerL.hmm` counts once, as `TerL`. With
`--hmm-mode mixed`, every model counts as its own query instead.

For each kept contig, the best hit of each matched query is reported
(`--top-per-contig` raises this). `all` and `threshold` check the queries you
gave, not the biology: a contig with every marker you asked for is not proof
of a complete genome, so record the exact queries with your results.

## Output files

| File or folder | Written when | Content |
| --- | --- | --- |
| `matched-contigs.fna` | always | The contigs that passed |
| `matched-contig-ids.txt` | always | Their ids, one per line |
| `{model}-hmmsearch-domains.tbl` | `--keep-domtbl` (default) | hmmsearch per-domain table of each model |
| `.phu/query.json` | always, except `--all-cazymes` / `--all-puls` | The query, how it was expanded and matched, and for database families the bundle version, file hashes and family names |
| `.phu/run.json` | `--all-cazymes` / `--all-puls` | Run manifest with counts, thresholds and files written |
| `query-hmms/` | database queries | The models taken from the databases |
| `predicted-proteins.faa` | `--keep-proteins` | Every predicted protein |
| `.phu/prediction.json` | `--keep-proteins` | Where the cached proteins are and their cache key |
| `target-proteins/{model}.faa` | `--save-target-proteins` | Matched proteins of each model, from kept contigs only |
| `target-hmms/{model}.hmm` | `--save-target-hmms` | An HMM built from each model's matched proteins |

Proteins are named `<contig>|gene<n>`; contig names may themselves contain
`|`.

!!! info "How models are named in output files"
    `{model}` is the model id made safe for file names (unusual characters
    become `_`). For your own HMMs in the default `--hmm-mode pure`, each file
    is one query and its file name (without extension) is the model id, also
    when the file holds several models. If
    two files share a name (`a/foo.hmm`, `b/foo.hmm`), the `NAME` inside each
    HMM tells them apart; if that is identical too, phu warns and pools their
    hits. If two different models would get the same file name (`beta/1` and
    `beta_1`), phu refuses the query rather than mix their results.

    With `--hmm-mode mixed`, each model in a file is its own query (for
    example a concatenated database), and outputs are written per model name
    inside the file.

## Saving matched proteins and building HMMs

`--save-target-proteins` writes the proteins that matched each model, taken
only from contigs that were kept, one file per model.

`--save-target-hmms` (with `--save-target-proteins`) builds a new HMM from
each model's matched proteins: they are aligned with FAMSA and the HMM is
built with pyHMMER, all inside phu. A model with a single protein gets a
single-sequence HMM. Only the proteins of the current run are used.

`--trim-target-msa` trims each alignment with trimAl before the HMM is built:

| Value | What is removed |
| --- | --- |
| `none` | Nothing (default) |
| `gappyout` | Columns judged too gappy from the alignment's own gap distribution; the mildest choice |
| `strict` | Gappy columns and columns with low residue similarity |
| `strictplus` | As `strict`, with a tighter similarity criterion; the most aggressive |

Trimming shortens the model and lowers bit scores, so cutoffs chosen for
untrimmed models do not carry over. The run prints the columns before and
after (`gappyout: 1523 -> 377 columns`). If trimming leaves no columns or
fewer than two sequences, the untrimmed alignment is used, with a warning.

These HMMs are built from whatever matched, so check them before reusing them
as references.

## CAZymes and PULs

Two options search the whole dbCAN collection instead of a query.

**`--all-cazymes`** searches every CAZyme family (AA, CBM, CE, GH, GT and PL,
keeping subfamilies such as `GH5_2`; helper modules such as SLH, cohesin and
dockerin are left out) and keeps contigs with at least one family hit.

```text
phu-screen/
├── matched-contigs.fna
├── matched-contig-ids.txt
├── cazyme-matches.tsv.gz        # one row per contig and family
├── cazyme-family-summary.tsv    # contigs and proteins per family
├── cazyme-class-summary.tsv     # contigs and proteins per class
├── evidence/
│   └── cazyme-hits.tsv.gz       # every hit that passed
└── .phu/
    └── run.json
```

**`--all-puls`** checks every dbCAN-PUL signature: a contig is kept when all
families of at least one PUL occur on it. Each family is searched once. PULs
that name families missing from the installed dbCAN models are skipped and
listed in `.phu/run.json`; that reflects the database, not your contigs.

```text
phu-screen/
├── matched-contigs.fna
├── matched-contig-ids.txt
├── pul-matches.tsv.gz           # one row per contig and matched PUL
├── pul-family-support.tsv.gz    # the family hits behind each match
├── pul-summary.tsv
├── pul-substrate-summary.tsv
├── evidence/
│   └── cazyme-hits.tsv.gz
└── .phu/
    └── run.json
```

The CAZyme evidence is kept even when no PUL matches; the PUL tables then have
only their headers. The per-family HMMs and domain tables of these runs are
temporary and removed at the end.

Counts in the summaries are occurrences, not normalized abundances. The
substrates in `pul-substrate-summary.tsv` are those annotated for the
reference PULs, not predictions for your contigs. This is targeted HMM
screening of CAZyme families, not a full `run_dbcan` annotation.

## How it works

1. **Predict proteins** with pyrodigal-gv, which handles the compact genes
   and alternative genetic codes of phages (`--mode meta`, the default, for
   mixed contigs; `single` for one genome).
   Proteins shorter than `--min-protein-len-aa` are dropped. Predictions are
   cached and reused while the contigs and prediction settings are unchanged
   (see [Cache handling](../cache.md)).
2. **Search** the proteins with each model using pyHMMER.
3. **Filter** the hits with the cutoffs of each model's kind, then keep
   contigs by `--combine-mode`.
4. **Write** the kept contigs and, if asked, their proteins and new HMMs.

Everything runs inside phu: no HMMER, seqkit or aligner program is needed.

## Examples

Keep contigs that carry both a major head protein and a terminase large
subunit (PHROG):

```bash
phu screen -i contigs.fa --combine-mode all phrog_10 phrog_2
```

Keep contigs with at least 3 of 4 marker HMMs, and build new HMMs from what
matched:

```bash
phu screen -i metagenome.fa --combine-mode threshold --min-hmm-hits 3 \
    --save-target-proteins --save-target-hmms \
    marker1.hmm marker2.hmm marker3.hmm marker4.hmm
```

Search a whole ECOD homology group, without the per-family domain tables:

```bash
phu screen -i contigs.fa --no-keep-domtbl ecod:2004.1
```

Be stricter, and use 16 threads:

```bash
phu screen -i contigs.fa --max-evalue 1e-10 --threads 16 protein_family.hmm
```

Turn off Pfam's curated cutoffs and use your own score:

```bash
phu screen -i contigs.fa --no-cut-ga --min-bitscore 30 PF00001
```

Quote paths that contain shell characters:

```bash
phu screen -i 'scaffolds|assembly.fa' protein_family.hmm
```

## Options

```bash
                                                                                
 Usage: phu screen [OPTIONS] [QUERY]...                                         
                                                                                
 Find contigs that carry given protein families.                                
                                                                                
 Proteins are predicted from the contigs and searched with each QUERY: your     
 own HMM files, or families from the phu databases by name: Pfam (PF00001),     
 KOfam (K00001), PHROG (phrog_216; find them with `phu phrogs`), ECOD           
 (ecod:2004.1.1.5 for one family, or a T, H or X group such as ecod:2004.1.1    
 for all its families), dbCAN (GH128) or one dbCAN PUL (PUL0621). Each          
 database's own cutoffs apply.                                                  
                                                                                
 Examples:                                                                      
     phu screen -i contigs.fa PF00001 PF00589                                   
     phu screen -i contigs.fa my-models/*.hmm                                   
     phu screen -i contigs.fa --combine-mode all K00001 K00002                  
     phu screen -i contigs.fa --all-cazymes                                     
                                                                                
 Full guide: https://camilogarciabotero.github.io/phu/commands/screen/          
                                                                                
╭─ Query ──────────────────────────────────────────────────────────────────────╮
│   [QUERY]...      <path>  HMM files, or database names: PF00001, K00001,     │
│                           phrog_216, ecod:2004.1.1.5, GH128, PUL0621.        │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs  -i      <path>  Contigs to screen (FASTA). [required]    │
│    --help           -h              Show this message and exit.              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Other ways to query ────────────────────────────────────────────────────────╮
│ --all-cazymes                 Search all dbCAN CAZyme families instead of    │
│                               QUERY.                                         │
│ --all-puls                    Search for every dbCAN PUL signature instead   │
│                               of QUERY.                                      │
│ --hmm-mode     -M      <str>  'pure': each HMM file is one query, even with  │
│                               several models in it; 'mixed': each model is   │
│                               one query.                                     │
│                               [default: pure]                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Hit thresholds ─────────────────────────────────────────────────────────────╮
│ --max-evalue      -e                      <float>          Highest E-value a │
│                                                            hit may have.     │
│                                                            [default: 1e-06]  │
│ --min-bitscore    -b                      <float>          Lowest bit score  │
│                                                            a hit may have.   │
│ --min-hmm-cover…                          <float range>    Share of your own │
│                                           [0.0<=x<=1.0]    HMM a hit must    │
│                                                            cover (0 = off).  │
│                                                            [default: 0.35]   │
│ --min-protein-c…                          <float range>    Share of the      │
│                                           [0.0<=x<=1.0]    protein a hit     │
│                                                            must cover (off   │
│                                                            by default).      │
│ --cut-ga              --no-cut-ga                          Use Pfam's        │
│                                                            curated cutoffs   │
│                                                            (GA).             │
│                                                            [default: cut-ga] │
│ --use-kofam-thr…      --no-use-kofam-…                     Use KOfam's       │
│                                                            per-KO cutoffs.   │
│                                                            [default:         │
│                                                            use-kofam-thresh… │
│ --use-phrog-thr…      --no-use-phrog-…                     For PHROGs, need  │
│                                                            60 bits and half  │
│                                                            of model and      │
│                                                            protein aligned.  │
│                                                            [default:         │
│                                                            use-phrog-thresh… │
│ --use-ecod-thre…      --no-use-ecod-t…                     For ECOD          │
│                                                            families, need    │
│                                                            half of the model │
│                                                            aligned.          │
│                                                            [default:         │
│                                                            use-ecod-thresho… │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Which contigs are kept ─────────────────────────────────────────────────────╮
│ --combine-mode    -c      <str>  Keep contigs hit by 'any' query, 'all' of   │
│                                  them, or a 'threshold'.                     │
│                                  [default: any]                              │
│ --min-hmm-hits    -k      <int>  With --combine-mode threshold: how many     │
│                                  queries must hit.                           │
│                                  [default: 1]                                │
│ --top-per-contig  -n      <int>  Best hits to report per contig.             │
│                                  [default: 1]                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ --output-folder       -o                         <path>  Output folder.      │
│                                                          [default:           │
│                                                          phu-screen]         │
│ --keep-proteins           --no-keep-proteins             Keep the predicted  │
│                                                          proteins.           │
│                                                          [default:           │
│                                                          no-keep-proteins]   │
│ --keep-domtbl             --no-keep-domtbl               Keep the hmmsearch  │
│                                                          domain tables.      │
│                                                          [default:           │
│                                                          keep-domtbl]        │
│ --save-target-prote…      --no-save-target-p…            Save the matched    │
│                                                          proteins of each    │
│                                                          model               │
│                                                          (target-proteins/). │
│                                                          [default:           │
│                                                          no-save-target-pro… │
│ --save-target-hmms        --no-save-target-h…            Build an HMM from   │
│                                                          each model's        │
│                                                          matched proteins    │
│                                                          (target-hmms/).     │
│                                                          [default:           │
│                                                          no-save-target-hmm… │
│ --trim-target-msa                                <str>   Trim alignments     │
│                                                          before building     │
│                                                          those HMMs: none,   │
│                                                          gappyout, strict,   │
│                                                          strictplus.         │
│                                                          [default: none]     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Gene prediction ────────────────────────────────────────────────────────────╮
│ --mode                -m      <str>               'meta' for mixed contigs,  │
│                                                   'single' for one genome.   │
│                                                   [default: meta]            │
│ --min-protein-len-aa  -g      <int range> [x>=1]  Shortest protein to keep   │
│                                                   (aa).                      │
│                                                   [default: 30]              │
│ --ttable              -T      <int range> [x>=1]  Genetic code (NCBI table); │
│                                                   by default predicted per   │
│                                                   contig.                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --threads  -t      <int range> [x>=1]  Threads. [default: 1]                 │
│ --quiet                                Print less.                           │
│ --verbose                              Print each search step.               │
╰──────────────────────────────────────────────────────────────────────────────╯
```

`--verbose` prints how each step was run (gene prediction, hmmsearch, FAMSA
and extraction settings) on lines starting with `[cmd]`; `--quiet` hides the
timed step lines. See [Run logs](../troubleshooting.md#run-logs).
