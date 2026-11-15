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
phu jack -i contigs.fa --phrog alias:portal
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
- **Mixing.** FASTA and PHROG seeds can be combined; `--combine-mode all` and
  `threshold` count both kinds. A FASTA record named like a selected PHROG
  (`phrog_2`) is refused, so that every seed ID is unique.

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

`jackhmmer-hits.tsv` contains final hit records and filtering metadata. `jackhmmer-iterations.tsv` records the per-iteration search summary. Both have a `seed_source` column, `fasta` or `phrog`. `.phu/query.json` records the seeds: the FASTA, the `--phrog` queries, and each PHROG seed with its annotation, category and the query that selected it, plus the database bundle's provenance. With `--save-hmm`, the final model is written as `final-iteration.hmm` for a single seed or under `final-iteration-hmms/` for multiple seeds.

## Command options

```text
Usage: phu jack [OPTIONS] {seed_marker}

 Iteratively screen contigs from one or more seed protein markers with
 pyhmmer.jackhmmer.

 Combine modes for multi-seed screening:
 - any: keep contigs hit by at least one seed (default)
 - all: keep contigs hit by all seeds
 - threshold: keep contigs hit by at least --min-seed-hits seeds

 Examples:
     phu jack -i contigs.fa marker_seed.faa
     phu jack -i contigs.fa --combine-mode all marker_seeds.faa
     phu jack -i contigs.fa --iterations 7 --inc-evalue 1e-4 marker_seed.faa

╭─ Arguments ───────────────────────────────────────────────────────────────╮
│ *    seed_marker      <path>  Seed marker protein FASTA (supports one or  │
│                               more sequences) [required]                   │
╰───────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs  -i      <path>  Input contigs FASTA [required]           │
│    --help           -h              Show this message and exit.              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ --output-folder  -o                        <path>  Output directory          │
│                                                    [default: phu-jack]       │
│ --keep-proteins      --no-keep-proteins            Keep the protein FASTA    │
│                                                    used for searching        │
│                                                    [default:                 │
│                                                    no-keep-proteins]         │
│ --save-hmm           --no-save-hmm                 Save the last jackhmmer   │
│                                                    iteration HMM as          │
│                                                    final-iteration.hmm       │
│                                                    [default: no-save-hmm]    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Prediction ─────────────────────────────────────────────────────────────────╮
│ --mode                -m      <str>               pyrodigal mode:            │
│                                                   meta|single                │
│                                                   [default: meta]            │
│ --min-gene-len        -g      <int>               Minimum gene length for    │
│                                                   pyrodigal (nt)             │
│                                                   [default: 90]              │
│ --min-protein-len-aa          <int range> [x>=1]  Minimum translated protein │
│                                                   length to keep (aa)        │
│                                                   [default: 30]              │
│ --ttable              -T      <int range> [x>=1]  NCBI translation table;    │
│                                                   default uses each contig's │
│                                                   predicted table            │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --threads  -t      <int range> [x>=1]  Threads for both pyrodigal and        │
│                                        pyhmmer                               │
│                                        [default: 1]                          │
│ --quiet                                Suppress routine progress output.     │
│ --verbose                              Show additional progress details.     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Iterative search ───────────────────────────────────────────────────────────╮
│ --iterations  -I      <int range> [x>=1]  Maximum jackhmmer iterations       │
│                                           [default: 5]                       │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Search thresholds ──────────────────────────────────────────────────────────╮
│ --inc-evalue          <float>  Inclusion E-value threshold for iterative     │
│                                jackhmmer                                     │
│                                [default: 0.001]                              │
│ --max-evalue  -e      <float>  Maximum independent E-value to keep a final   │
│                                hit                                           │
│                                [default: 1e-05]                              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Matching behavior ──────────────────────────────────────────────────────────╮
│ --top-per-contig  -n      <int range> [x>=1]  Keep top-N hits per contig (by │
│                                               bitscore)                      │
│                                               [default: 1]                   │
│ --combine-mode    -c      <str>               How to combine hits from       │
│                                               multiple seed proteins:        │
│                                               any|all|threshold              │
│                                               [default: any]                 │
│ --min-seed-hits   -k      <int range> [x>=1]  Minimum number of seeds that   │
│                                               must hit a contig (for         │
│                                               threshold mode)                │
│                                               [default: 1]                   │
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
