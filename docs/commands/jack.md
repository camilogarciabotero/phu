# jack

## What does it do?

The `phu jack` command iteratively screens contigs from one or more seed protein markers using `pyhmmer.hmmer.jackhmmer`. It predicts proteins from the input contigs, searches them against each seed, and uses successive iterations to expand the detected protein family.

This command is useful when a trusted marker protein is available but a profile HMM is not. Interpret expanded hits as sequence-similarity evidence, not as a validated functional or taxonomic assignment.

## Synopsis

```bash
phu jack -i <INPUT_CONTIGS> [OPTIONS] [<SEED_MARKER>] [--phrog QUERY ...]
```

`<SEED_MARKER>` is a protein FASTA file and may contain one or more seed sequences. PHROG families can be seeds too, with `--phrog` (see [PHROG seeds](#phrog-seeds)); give a FASTA, `--phrog`, or both.

**Examples:**

```bash
phu jack -i contigs.fa marker_seed.faa
phu jack -i contigs.fa --combine-mode all marker_seeds.faa
phu jack -i contigs.fa --iterations 7 --inc-evalue 1e-4 marker_seed.faa
phu jack -i contigs.fa --phrog alias:portal --max-phrog-seeds 60
phu jack -i contigs.fa --combine-mode all --phrog phrog_2 marker_seed.faa
```

## Iterative search

Each seed is searched independently with `jackhmmer`. The `--iterations` option sets the maximum number of iterations, while `--inc-evalue` controls which hits are included in the next iteration. Final hits are filtered by `--max-evalue`. The `--top-per-contig` option limits retained hits by bitscore after filtering.

When multiple seed sequences are provided, `--combine-mode` controls contig retention:

- `any`: retain contigs hit by at least one seed.
- `all`: retain contigs hit by every seed.
- `threshold`: retain contigs hit by at least `--min-seed-hits` seeds.

## PHROG seeds

`--phrog QUERY` (repeatable) uses PHROG families as seeds. Queries are those of
[`phu phrogs resolve`](phrogs.md): `phrog_<n>`, `alias:<name>`,
`category:<name>` or `'annot:<annotation>'`. Use `phu phrogs resolve` first to
see what a query selects.

- **How a PHROG seed is searched.** The first iteration searches with the
  PHROG model itself, a profile built from the whole family; later
  iterations use models rebuilt from the included hits, as for sequence
  seeds. `--iterations`, `--inc-evalue` and `--max-evalue` apply to both.
  `--iterations 1` is a plain search with the PHROG model, the same search as
  `phu screen phrog_<n>` but with jack's thresholds. When nothing is included,
  the seed stops after one iteration, reported as converged, as jackhmmer does
  for a sequence seed.
- **Drift.** Models rebuilt from hits can move away from the PHROG family,
  especially when the first iteration includes few or distant proteins.
  Check `jackhmmer-iterations.tsv`, and lower `--iterations` if later
  iterations add unexpected hits.
- **Cap.** Each PHROG is one iterative search, and an alias or category can
  select dozens to thousands of PHROGs (`alias:terl` 54, `category:tail`
  1,238). jack stops before predicting proteins when the selection is larger
  than `--max-phrog-seeds` (default 50). Narrow the query or raise the cap.
- **Seed groups.** `--combine-mode` counts seed groups, not single seeds.
  Each FASTA sequence is a group, and each `--phrog` query is one group
  however many PHROGs it selects: the 54 PHROGs of `alias:terl` are
  alternatives for one marker. So `--combine-mode all --phrog alias:terl
  --phrog alias:portal --max-phrog-seeds 110` keeps contigs with a terminase *and* a portal hit,
  and `--min-seed-hits 2` needs two different queries. A PHROG selected by
  two queries belongs to the first. Within a group, a protein hit by several
  PHROGs is kept once, by its best hit.
- **Mixing.** FASTA and PHROG seeds can be combined. A FASTA record named like
  a selected PHROG (`phrog_2`) is refused, so that every seed ID is unique.
- **What a query meant.** With `--verbose`, jack logs each query with the
  PHROG identifiers it selected, then each PHROG seed with its annotation,
  category and group (lines starting `[cmd] phu.phrogs`). The same is kept in
  `.phu/query.json`.

## Input and output

The input contigs are translated with pyrodigal using `--mode`, `--min-gene-len`, `--min-protein-len-aa`, and optionally `--ttable`. The resulting protein FASTA is cached and reused when the prediction inputs are unchanged. Use `--keep-proteins` to retain the FASTA in the output directory as `predicted-proteins.faa` (with the cache key in `.phu/prediction.json`).

A successful run writes these core outputs:

```text
phu-jack/
├── matched-contigs.fna
├── matched-contig-ids.txt
├── jackhmmer-hits.tsv
└── jackhmmer-iterations.tsv
```

`jackhmmer-hits.tsv` contains final hit records and filtering metadata. `jackhmmer-iterations.tsv` records the per-iteration search summary. Both have a `seed_source` column (`fasta` or `phrog`) and a `seed_group` column (the FASTA seed ID, or the `--phrog` query that selected the PHROG). `.phu/query.json` records the seeds: the FASTA, the `--phrog` queries, and each PHROG seed with its annotation, category and the query that selected it, plus the database bundle's provenance. With `--save-hmm`, the final model is written as `final-iteration.hmm` for a single seed or under `final-iteration-hmms/` for multiple seeds.

## Command options

```text

 Usage: phu jack [OPTIONS] [seed_marker]

 Find contigs carrying relatives of seed proteins, by iterative search.

 Each seed searches the proteins predicted from the contigs; hits that pass
 --inc-evalue build the model for the next round, so distant family members
 are found round by round. Seeds are protein sequences (a FASTA file),
 PHROG families (--phrog; find them with `phu phrogs`), or both.

 Examples:
     phu jack -i contigs.fa seeds.faa
     phu jack -i contigs.fa --combine-mode all seeds.faa
     phu jack -i contigs.fa --phrog alias:terl --max-phrog-seeds 60
     phu jack -i contigs.fa --iterations 1 --phrog phrog_216

 Full guide: https://camilogarciabotero.github.io/phu/commands/jack/

╭─ Seeds ──────────────────────────────────────────────────────────────────────╮
│   seed_marker      <path>  Protein FASTA of seed sequences (optional with    │
│                            --phrog).                                         │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs  -i      <path>  Contigs to screen (FASTA). [required]    │
│    --help           -h              Show this message and exit.              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ PHROG seeds ────────────────────────────────────────────────────────────────╮
│ --phrog                  <str>               PHROG families as seeds         │
│                                              (repeatable): phrog_216,        │
│                                              alias:terl, category:tail.      │
│ --max-phrog-seeds        <int range> [x>=1]  Stop if --phrog selects more    │
│                                              PHROGs than this.               │
│                                              [default: 50]                   │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Search ─────────────────────────────────────────────────────────────────────╮
│ --iterations  -I      <int range> [x>=1]  Most search rounds per seed.       │
│                                           [default: 5]                       │
│ --inc-evalue          <float>             E-value a hit needs to join the    │
│                                           next round's model.                │
│                                           [default: 0.001]                   │
│ --max-evalue  -e      <float>             Highest E-value a final hit may    │
│                                           have.                              │
│                                           [default: 1e-05]                   │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Which contigs are kept ─────────────────────────────────────────────────────╮
│ --combine-mode    -c      <str>               Keep contigs hit by 'any'      │
│                                               seed, 'all', or a 'threshold'  │
│                                               (one --phrog query = one       │
│                                               seed).                         │
│                                               [default: any]                 │
│ --min-seed-hits   -k      <int range> [x>=1]  With --combine-mode threshold: │
│                                               how many seeds must hit.       │
│                                               [default: 1]                   │
│ --top-per-contig  -n      <int range> [x>=1]  Best hits to report per        │
│                                               contig.                        │
│                                               [default: 1]                   │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ --output-folder  -o                        <path>  Output folder.            │
│                                                    [default: phu-jack]       │
│ --keep-proteins      --no-keep-proteins            Keep the predicted        │
│                                                    proteins.                 │
│                                                    [default:                 │
│                                                    no-keep-proteins]         │
│ --save-hmm           --no-save-hmm                 Save each seed's          │
│                                                    last-round HMM.           │
│                                                    [default: no-save-hmm]    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Gene prediction ────────────────────────────────────────────────────────────╮
│ --mode                -m      <str>               'meta' for mixed contigs,  │
│                                                   'single' for one genome.   │
│                                                   [default: meta]            │
│ --min-gene-len        -g      <int>               Shortest gene to predict   │
│                                                   (nt).                      │
│                                                   [default: 90]              │
│ --min-protein-len-aa          <int range> [x>=1]  Shortest protein to keep   │
│                                                   (aa).                      │
│                                                   [default: 30]              │
│ --ttable              -T      <int range> [x>=1]  Genetic code (NCBI table); │
│                                                   by default predicted per   │
│                                                   contig.                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --threads  -t      <int range> [x>=1]  Threads. [default: 1]                 │
│ --quiet                                Print less.                           │
│ --verbose                              Print each search step, and what      │
│                                        --phrog queries selected.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Workflow integration

```bash
# Start with a trusted marker protein.
phu jack \
  --input-contigs viral_assembly.fasta \
  --output-folder phu-jack-capsid \
  --iterations 7 \
  --inc-evalue 1e-4 \
  marker_seed.faa
```

Use `--save-hmm` when the final iterative model is needed for inspection or follow-up searches. Use `--quiet` in scripted workflows and retain the default output tables for provenance.
