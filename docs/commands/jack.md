# jack

## What does it do?

The `phu jack` command iteratively screens contigs from one or more seed protein markers using `pyhmmer.hmmer.jackhmmer`. It predicts proteins from the input contigs, searches them against each seed, and uses successive iterations to expand the detected protein family.

This command is useful when a trusted marker protein is available but a profile HMM is not. Interpret expanded hits as sequence-similarity evidence, not as a validated functional or taxonomic assignment.

## Synopsis

```bash
phu jack -i <INPUT_CONTIGS> [OPTIONS] <SEED_MARKER>
```

`<SEED_MARKER>` is a protein FASTA file and may contain one or more seed sequences.

**Examples:**

```bash
phu jack -i contigs.fa marker_seed.faa
phu jack -i contigs.fa --combine-mode all marker_seeds.faa
phu jack -i contigs.fa --iterations 7 --inc-evalue 1e-4 marker_seed.faa
```

## Iterative search

Each seed is searched independently with `jackhmmer`. The `--iterations` option sets the maximum number of iterations, while `--inc-evalue` controls which hits are included in the next iteration. Final hits are filtered by `--max-evalue`. The `--top-per-contig` option limits retained hits by bitscore after filtering.

When multiple seed sequences are provided, `--combine-mode` controls contig retention:

- `any`: retain contigs hit by at least one seed.
- `all`: retain contigs hit by every seed.
- `threshold`: retain contigs hit by at least `--min-seed-hits` seeds.

## Input and output

The input contigs are translated with pyrodigal using `--mode`, `--min-gene-len`, `--min-protein-len-aa`, and optionally `--ttable`. The resulting protein FASTA is cached and reused when the prediction inputs are unchanged. Use `--keep-proteins` to retain the FASTA in the output directory.

A successful run writes these core outputs:

```text
phu-jack/
├── screened_contigs.fasta
├── kept_contigs.txt
├── jackhmmer_hits.tsv
├── jackhmmer_iterations.tsv
└── .phu/
    └── run.json
```

`jackhmmer_hits.tsv` contains final hit records and filtering metadata. `jackhmmer_iterations.tsv` records the per-iteration search summary. With `--save-hmm`, the final model is written as `last_iteration.hmm` for a single seed or under `last_iteration_hmms/` for multiple seeds.

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
│                                                    last_iteration.hmm        │
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
