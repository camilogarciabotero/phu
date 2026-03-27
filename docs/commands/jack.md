# jack

## What does it do?

The `phu jack` command screens contigs using an iterative protein-seed search with `pyhmmer.hmmer.jackhmmer`.

It predicts proteins from your contigs, runs iterative searches from a seed marker protein, and extracts contigs that contain significant included hits.

This first implementation is intentionally strict and simple:

- It accepts one seed marker FASTA file.
- That seed file must contain exactly one protein sequence.
- Extra HMM confirmation is not part of this version.

## Synopsis

```bash
phu jack --input-contigs [INPUT_CONTIGS] [SEED_MARKER]
```

**Example:**

```bash
phu jack --input-contigs contigs.fa marker_seed.faa
```

## How it works

1. Predict proteins from contigs with `pyrodigal`.
2. Run iterative `jackhmmer` from the seed marker against predicted proteins.
3. Keep included final hits and apply final `--max-evalue` filtering.
4. Keep top-N hits per contig (`--top-per-contig`).
5. Extract matching contigs to `screened_contigs.fasta`.

## Outputs

By default in `phu-jack/`:

- `screened_contigs.fasta`: extracted contigs that matched iterative search.
- `kept_contigs.txt`: contig IDs kept in final output.
- `jackhmmer_hits.tsv`: final retained protein-level hits.
- `jackhmmer_iterations.tsv`: iteration-level summary (`n_hits`, `n_included`, `converged`).

## Command options

```bash
Usage: phu jack [OPTIONS] SEED_MARKER

Iteratively screen contigs from a single seed protein marker with pyhmmer.jackhmmer.

Arguments:
  SEED_MARKER                  Seed marker protein FASTA (must contain exactly one sequence)

Options:
  -i, --input-contigs PATH     Input contigs FASTA [required]
  -o, --output-folder PATH     Output directory [default: phu-jack]
  -m, --mode TEXT              pyrodigal mode: meta|single [default: meta]
  -t, --threads INTEGER        Threads for pyrodigal and pyhmmer [default: 1]
  -I, --iterations INTEGER     Maximum jackhmmer iterations [default: 5]
      --inc-evalue FLOAT       Inclusion E-value threshold [default: 0.001]
  -e, --max-evalue FLOAT       Max independent E-value for final hits [default: 1e-05]
  -n, --top-per-contig INTEGER Keep top-N hits per contig [default: 1]
  -g, --min-gene-len INTEGER   Minimum gene length for pyrodigal [default: 90]
  -T, --ttable INTEGER         NCBI translation table [default: 11]
      --keep-proteins / --no-keep-proteins
                               Keep intermediate proteins FASTA
  -h, --help                   Show help and exit
```

## Notes

If your seed FASTA contains multiple sequences, `phu jack` exits with an error. This is by design for the first version to keep behavior explicit and predictable.
