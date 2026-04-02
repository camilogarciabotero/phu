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

## Command options

```bash
Usage: phu jack [OPTIONS] SEED_MARKER

Iteratively screen contigs from one or more seed protein markers with pyhmmer.jackhmmer.

Arguments:
  SEED_MARKER                  Seed marker protein FASTA (supports one or more sequences)

Options:
  -i, --input-contigs PATH     Input contigs FASTA [required]
  -o, --output-folder PATH     Output directory [default: phu-jack]
  -m, --mode TEXT              pyrodigal mode: meta|single [default: meta]
  -t, --threads INTEGER        Threads for pyrodigal and pyhmmer [default: 1]
  -I, --iterations INTEGER     Maximum jackhmmer iterations [default: 5]
      --inc-evalue FLOAT       Inclusion E-value threshold [default: 0.001]
  -e, --max-evalue FLOAT       Max independent E-value for final hits [default: 1e-05]
  -n, --top-per-contig INTEGER Keep top-N hits per contig [default: 1]
  -c, --combine-mode TEXT      Combine seed hits per contig: any|all|threshold [default: any]
  -k, --min-seed-hits INTEGER  Minimum number of seeds required for threshold mode [default: 1]
  -g, --min-gene-len INTEGER   Minimum gene length for pyrodigal [default: 90]
      --min-protein-len-aa INTEGER  Minimum translated protein length to keep [default: 30]
  -T, --ttable INTEGER         NCBI translation table [default: 11]
      --keep-proteins / --no-keep-proteins
                               Keep intermediate proteins FASTA
      --save-hmm / --no-save-hmm
                               Save final jackhmmer HMM output(s)
  -h, --help                   Show help and exit
```

## Notes

If you use `--combine-mode threshold`, `--min-seed-hits` cannot exceed the number of seed sequences in `SEED_MARKER`.
