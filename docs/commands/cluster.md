# cluster

Sequence clustering wrapper around external 'vclust' with three preconfigured modes for common viral genomics workflows.

## Overview

The `cluster` command provides a simplified interface to the powerful `vclust` tool, implementing the most common use cases from the [vclust wiki](https://github.com/refresh-bio/vclust/wiki/6-Use-cases) with sensible defaults while allowing advanced customization.

## Basic Usage

```bash
phu cluster --mode <MODE> --input-contigs <FASTA_FILE> [options]
```

## Modes

| Mode | Description | Algorithm | Metric | Use Case |
|------|-------------|-----------|--------|----------|
| `dereplication` | Remove redundant sequences, keep representatives | cd-hit | ani | Reduce dataset to representative genomes (Section 6.3) |
| `votu` | Cluster into viral Operational Taxonomic Units | leiden | ani | Group contigs following MIUViG standards (Section 6.2) |
| `species` | Classify viruses into species | complete | tani | Species classification following ICTV standards (Section 6.1) |

## Default Parameters by Mode

| Parameter | dereplication | votu | species |
|-----------|---------------|------|---------|
| **Algorithm** | cd-hit | leiden | complete |
| **Metric** | ani | ani | tani |
| **ANI cutoff** | 95% | 95% | 95% |
| **Query coverage** | 85% | 85% | None |
| **Pre-filter min-ident** | 95% | 95% | 70% |
| **Output file** | clusters.tsv | clusters.tsv | species.tsv |
| **Representatives file** | dereplicated_representatives.fna | representatives.fna | representatives.fna |

## Command Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--mode` | Choice | **Required** | dereplication \| votu \| species |
| `--input-contigs` | Path | **Required** | Input FASTA file with viral sequences |
| `--output-folder` | Path | clustered-contigs | Output directory for results |
| `--threads` | Integer | 0 | Number of threads (0 = all cores) |
| `--vclust-params` | String | None | Custom vclust parameters (see Advanced Usage) |

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