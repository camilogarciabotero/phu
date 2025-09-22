# cluster

Sequence clustering wrapper around the external `vclust` tool with three preconfigured modes for common viral genomics workflows.

## Overview

The `phu cluster` command provides a simplified interface to `vclust`, implementing common use cases (see the [vclust wiki](https://github.com/refresh-bio/vclust/wiki/6-Use-cases)) while allowing advanced customization via `--vclust-params`.

## Synopsis

```bash
phu cluster --mode <MODE> --input-contigs <FASTA_FILE> [OPTIONS]
```

Outputs are placed in an output folder (default `clustered-contigs/`) and include cluster assignment TSVs and representative FASTA files, for example:

```
clustered-contigs/
├── ani.ids.tsv
├── ani.tsv
├── cluster_representatives_ids.txt
├── fltr.txt
├── representatives.fna
└── species.tsv
```

## Modes

- `dereplication` — remove redundant sequences while keeping representatives (cd-hit/ANI-based).

- `votu` — cluster into viral Operational Taxonomic Units (leiden/ANI-based) following MIUViG-style defaults.

- `species` — classify sequences into species (complete/TANI-based) following ICTV-style defaults.

## Default parameters by mode

| Parameter | dereplication | votu | species |
|-----------|---------------|------|---------|
| Algorithm | cd-hit | leiden | complete |
| Metric | ani | ani | tani |
| ANI cutoff | 95% | 95% | 95% |
| Query coverage | 85% | 85% | None |
| Pre-filter min-ident | 95% | 95% | 70% |

## Command Options

```bash 
 Sequence clustering wrapper around external 'vclust' with three modes.      
                                                                             
 For advanced usage, provide custom vclust parameters as a quoted string.    
 See the vclust wiki for parameter details:                                  
 https://github.com/refresh-bio/vclust/wiki                                  
                                                                             
 Example:                                                                    
     phu cluster --mode votu --input-contigs genomes.fna                                               
                                                                             
╭─ Options ─────────────────────────────────────────────────────────────────╮
│ *  --mode                   [dereplication|votu|s  dereplication | votu | │
│                             pecies]                species                │
│                                                    [required]             │
│ *  --input-contigs          PATH                   Input FASTA [required] │
│    --output-folder          PATH                   Output directory       │
│                                                    [default:              │
│                                                    clustered-contigs]     │
│    --threads                INTEGER RANGE [x>=0]   0=all cores; otherwise │
│                                                    N threads              │
│                                                    [default: 0]           │
│    --vclust-params          TEXT                   Custom vclust          │
│                                                    parameters:            │
│                                                    "--min-kmers 20        │
│                                                    --outfmt lite --ani    │
│                                                    0.97"                  │
│    --help           -h                             Show this message and  │
│                                                    exit.                  │
╰───────────────────────────────────────────────────────────────────────────╯
```

## Examples

```bash
# Dereplicate viral contigs (default outputs in clustered-contigs/)
phu cluster --mode dereplication --input-contigs viral_contigs.fna

# Cluster into vOTUs following MIUViG standards
phu cluster --mode votu --input-contigs viral_contigs.fna

# Species classification following ICTV standards
phu cluster --mode species --input-contigs complete_genomes.fna
```

## Advanced: pass custom vclust parameters

```bash {hl_lines="2"}
phu cluster --mode species --input-contigs genomes.fna \
  --vclust-params="--metric tani --tani 0.70" 
```
>Treat species clustering as genus-level by lowering similarity to 70%

## Notes

- `--vclust-params` provides full control over `vclust` parameters; use it to tune behavior for large or highly-redundant datasets.
- For reproducibility, save the full vclust command and parameters used (the tool often logs them in the output folder).

## See also

- vclust wiki: https://github.com/refresh-bio/vclust/wiki/6-Use-cases
