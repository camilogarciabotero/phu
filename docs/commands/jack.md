# jack

## What does it do?

The `phu jack` command screens contigs using an iterative protein-seed search with `pyhmmer.hmmer.jackhmmer`.

It predicts proteins from your contigs, runs iterative searches from one or more seed proteins, and extracts contigs that contain significant included hits.

This implementation supports multi-seed workflows:

- It accepts one seed marker FASTA file containing one or more protein sequences.
- Seed IDs in FASTA headers must be unique.
- You can combine seed evidence with `--combine-mode any|all|threshold`.

## Synopsis

```bash
phu jack --input-contigs [INPUT_CONTIGS] [SEED_MARKER]
```

**Example:**

```bash
phu jack --input-contigs contigs.fa marker_seed.faa
phu jack --input-contigs contigs.fa --combine-mode all marker_seeds.faa
phu jack --input-contigs contigs.fa --combine-mode threshold --min-seed-hits 3 marker_seeds.faa 
```

## How it works

1. Predict proteins from contigs with `pyrodigal`.
    Proteins shorter than `--min-protein-len-aa` are discarded before iterative searching.
2. Run iterative `jackhmmer` from each seed marker against predicted proteins.
3. Keep included final hits and apply final `--max-evalue` filtering.
4. Combine across seeds:
  - `any`: keep contigs hit by at least one seed.
  - `all`: keep contigs hit by all seeds.
  - `threshold`: keep contigs hit by at least `--min-seed-hits` seeds.
5. Keep top hits per contig (`--top-per-contig`) using the selected combine mode.
6. Extract matching contigs to `screened_contigs.fasta`.

## Outputs

By default in `phu-jack/`:

- `screened_contigs.fasta`: extracted contigs that matched iterative search.
- `kept_contigs.txt`: contig IDs kept in final output.
- `jackhmmer_hits.tsv`: final retained protein-level hits with `seed_id`.
- `jackhmmer_iterations.tsv`: iteration-level summary with `seed_id`, `n_hits`, `n_included`, `converged`.

When `--save-hmm` is enabled:

- Single-seed input: writes `last_iteration.hmm`.
- Multi-seed input: writes one HMM per seed in `last_iteration_hmms/*.hmm`.

## Cache handling

Protein prediction is cached and reused when the input contigs and prediction parameters are unchanged. For `phu jack`, changing `--min-gene-len`, `--min-protein-len-aa`, `--mode`, or `--ttable` rebuilds the cached proteins. Changing seed markers, combine mode, or output options does not.

See [cache.md](../cache.md) for the shared cache rules used by both `screen` and `jack`.

## Command Options

```bash
Usage: phu jack [OPTIONS] {seed_marker}

Iteratively screen contigs from one or more seed protein markers with pyhmmer.jackhmmer.

Combine modes for multi-seed screening:
- any: keep contigs hit by at least one seed (default)
- all: keep contigs hit by all seeds
- threshold: keep contigs hit by at least --min-seed-hits seeds

Examples:
    phu jack -i contigs.fa marker_seed.faa
    phu jack -i contigs.fa --combine-mode all marker_seeds.faa
    phu jack -i contigs.fa --iterations 7 --inc-evalue 1e-4 marker_seed.faa

╭─ Arguments ──────────────────────────────────────────────────────────────────────────────────────╮
│ *    seed_marker      <path>  Seed marker protein FASTA (supports one or more sequences)         │
│                               [required]                                                         │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs       -i                        <path>              Input contigs FASTA       │
│                                                                        [required]                │
│    --output-folder       -o                        <path>              Output directory          │
│                                                                        [default: phu-jack]       │
│    --mode                -m                        <str>               pyrodigal mode:           │
│                                                                        meta|single               │
│                                                                        [default: meta]           │
│    --threads             -t                        <int range> [x>=1]  Threads for both          │
│                                                                        pyrodigal and pyhmmer     │
│                                                                        [default: 1]              │
│    --iterations          -I                        <int range> [x>=1]  Maximum jackhmmer         │
│                                                                        iterations                │
│                                                                        [default: 5]              │
│    --inc-evalue                                    <float>             Inclusion E-value         │
│                                                                        threshold for iterative   │
│                                                                        jackhmmer                 │
│                                                                        [default: 0.001]          │
│    --max-evalue          -e                        <float>             Maximum independent       │
│                                                                        E-value to keep a final   │
│                                                                        hit                       │
│                                                                        [default: 1e-05]          │
│    --top-per-contig      -n                        <int range> [x>=1]  Keep top-N hits per       │
│                                                                        contig (by bitscore)      │
│                                                                        [default: 1]              │
│    --combine-mode        -c                        <str>               How to combine hits from  │
│                                                                        multiple seed proteins:   │
│                                                                        any|all|threshold         │
│                                                                        [default: any]            │
│    --min-seed-hits       -k                        <int range> [x>=1]  Minimum number of seeds   │
│                                                                        that must hit a contig    │
│                                                                        (for threshold mode)      │
│                                                                        [default: 1]              │
│    --min-gene-len        -g                        <int>               Minimum gene length for   │
│                                                                        pyrodigal (nt)            │
│                                                                        [default: 90]             │
│    --min-protein-len-aa                            <int range> [x>=1]  Minimum translated        │
│                                                                        protein length to keep    │
│                                                                        (aa)                      │
│                                                                        [default: 30]             │
│    --ttable              -T                        <int>               NCBI translation table    │
│                                                                        for coding sequences      │
│                                                                        [default: 11]             │
│    --keep-proteins           --no-keep-proteins                        Keep the protein FASTA    │
│                                                                        used for searching        │
│                                                                        [default:                 │
│                                                                        no-keep-proteins]         │
│    --save-hmm                --no-save-hmm                             Save the last jackhmmer   │
│                                                                        iteration HMM as          │
│                                                                        last_iteration.hmm        │
│                                                                        [default: no-save-hmm]    │
│    --quiet                                                             Suppress routine progress │
│                                                                        output.                   │
│    --verbose                                                           Show additional progress  │
│                                                                        details.                  │
│    --help                -h                                            Show this message and     │
│                                                                        exit.                     │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
```

## Notes

If you use `--combine-mode threshold`, `--min-seed-hits` cannot exceed the number of seed sequences in `SEED_MARKER`.
