# cluster

Sequence clustering wrapper around external 'vclust' with three preconfigured modes for common viral genomics workflows.

## Overview

The `cluster` command provides a simplified interface to the powerful `vclust` tool, implementing the most common use cases from the [vclust wiki](https://github.com/refresh-bio/vclust/wiki/6-Use-cases) with sensible defaults while allowing advanced customization.

## Basic Usage

```bash
phu cluster --mode <MODE> --input-contigs <FASTA_FILE> [options]
```

Output files include a TSV file with cluster assignments and a FASTA file with representative sequences.

```bash
clustered-contigs/
├── ani.ids.tsv
├── ani.tsv
├── cluster_representatives_ids.txt
├── fltr.txt
├── representatives.fna
└── species.tsv
```


## Modes

### Modes

The `dereplication` mode removes redundant sequences while keeping representatives, using the cd-hit algorithm and ani metric. It is designed to reduce datasets to representative genomes, as outlined in [Section 6.3](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#63-dereplicate-viral-contigs-into-representative-genomes).

The `votu` mode clusters sequences into viral Operational Taxonomic Units, employing the leiden algorithm and ani metric. This mode groups contigs according to MIUViG standards, detailed in [Section 6.2](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#62-assign-viral-contigs-into-votus-following-miuvig-standards).

The `species` mode classifies viruses into species, utilizing the complete algorithm and tani metric. It follows ICTV standards for species classification, as described in [Section 6.1](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#61-classify-viruses-into-species-and-genera-following-ictv-standards).

## Default Parameters by Mode

| Parameter | dereplication | votu | species |
|-----------|---------------|------|---------|
| **Algorithm** | cd-hit | leiden | complete |
| **Metric** | ani | ani | tani |
| **ANI cutoff** | 95% | 95% | 95% |
| **Query coverage** | 85% | 85% | None |
| **Pre-filter min-ident** | 95% | 95% | 70% |

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

### Basic Examples

```bash
# Dereplicate viral contigs
phu cluster --mode dereplication --input-contigs viral_contigs.fna

# Cluster into vOTUs following MIUViG standards
phu cluster --mode votu --input-contigs viral_contigs.fna

# Classify viruses into species following ICTV standards
phu cluster --mode species --input-contigs complete_genomes.fna
```

### Advanced Usage with Custom Parameters

The `--vclust-params` option allows you to customize any vclust parameter while maintaining the convenience of predefined modes. Parameters are automatically routed to the appropriate vclust command (prefilter, align, cluster).

#### Large Dataset Optimization (Wiki Section 6.6)

```bash
# Process large diverse dataset (IMG/VR style)
phu cluster --mode votu --input-contigs large_dataset.fna \
  --vclust-params="--min-kmers 4 --batch-size 2000000 --kmers-fraction 0.2 --outfmt lite"
```

#### Highly Redundant Dataset (Wiki Section 6.7)

```bash
# Process highly redundant sequences
phu cluster --mode votu --input-contigs redundant_genomes.fna \
  --vclust-params="--min-kmers 10 --batch-size 100000 --max-seqs 1000 --outfmt lite --ani 0.97 --qcov 0.95"
```

#### Custom Thresholds

```bash
# More stringent clustering
phu cluster --mode votu --input-contigs genomes.fna \
  --vclust-params="--ani 0.98 --qcov 0.90"

# Species clustering with custom genus threshold
phu cluster --mode species --input-contigs genomes.fna \
  --vclust-params="--tani 0.97"
```

## Comparison with Direct vclust Usage

| Task | phu cluster | Direct vclust |
|------|-------------|---------------|
| **Steps** | Single command | 3 separate commands |
| **Configuration** | Preconfigured modes | Manual parameter setup |
| **Customization** | `--vclust-params` option | Full control |
| **Learning curve** | Minimal | Requires vclust expertise |
| **Use case** | Common workflows | Specialized analyses |