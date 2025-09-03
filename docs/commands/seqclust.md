# seqclust

Sequence clustering wrapper around external 'vclust' with three preconfigured modes for common viral genomics workflows.

## Overview

The `seqclust` command provides a simplified interface to the powerful `vclust` tool, implementing the most common use cases from the [vclust wiki](https://github.com/refresh-bio/vclust/wiki/6-Use-cases) with sensible defaults while allowing advanced customization.

## Basic Usage

```bash
phu seqclust --mode <MODE> --input-contigs <FASTA> [options]
```

## Modes

### Dereplication Mode
The [`dereplication`](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#63-dereplicate-viral-contigs-into-representative-genomes) mode is designed to remove redundant sequences from your dataset, retaining only representative genomes. It employs the cd-hit algorithm with an ANI (Average Nucleotide Identity) metric to identify and eliminate duplicates, ensuring a streamlined collection of unique sequences. This is particularly useful for reducing large datasets to manageable sizes, as outlined in Section 6.3 of the vclust wiki, making it ideal for preprocessing before further analysis.

### vOTU-Clustering Mode
For clustering into viral Operational Taxonomic Units (vOTUs), the [`votu-clustering`](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#62-assign-viral-contigs-into-votus-following-miuvig-standards) mode uses the leiden algorithm alongside the ANI metric. This approach groups viral contigs according to MIUViG (Minimum Information about an Uncultivated Virus Genome) standards, facilitating the identification of distinct viral populations. As detailed in Section 6.2, it provides a standardized way to delineate vOTUs from metagenomic assemblies, enhancing comparability across studies.

### SPP-Clustering Mode
The [`spp-clustering`](https://github.com/refresh-bio/vclust/wiki/6-Use-cases#61-classify-viruses-into-species-and-genera-following-ictv-standards) mode classifies viruses into species-level groups using the complete algorithm and a TANI (Taxonomic Average Nucleotide Identity) metric. It adheres to ICTV (International Committee on Taxonomy of Viruses) standards for species demarcation, allowing precise taxonomic assignment of viral genomes. This mode, referenced in Section 6.1, is essential for evolutionary and ecological studies requiring accurate species-level resolution.

## Default Parameters by Mode
| Parameter | dereplication | votu-clustering | spp-clustering |
|-----------|---------------|-----------------|----------------|
| **Clustering Algorithm** | cd-hit | leiden | complete |
| **Metric** | ani | ani | tani |
| **ANI cutoff** | 95% | 95% | 95% |
| **Query coverage** | 85% | 85% | None |
| **Pre-filter min-ident** | 95% | 95% | 70% |

!!! tip "Quick run"

    To cluster at the genus level in `spp-clustering` mode, pass `--tani 0.70` via the `--vclust-params` option. This adjusts the taxonomic threshold for genus demarcation according to ICTV standards.


## Command Options

```bash
phu seqclust -h

 Usage: phu seqclust [OPTIONS]                                                            
                                                                                          
 Sequence clustering wrapper around external 'vclust' with three modes.                   
                                                                                          
 For advanced usage, provide custom vclust parameters as a quoted string.                 
 See the vclust wiki for parameter details: https://github.com/refresh-bio/vclust/wiki    
                                                                                          
 Example:                                                                                 
     phu seqclust --mode spp-clustering --input-contigs genomes.fna                       
                                                 
                                                                                          
╭─ Options ──────────────────────────────────────────────────────────────────────────────╮
│ *  --mode                   [dereplication|votu-clusteri  dereplication |              │
│                             ng|spp-clustering]            votu-clustering |            │
│                                                           spp-clustering               │
│                                                           [required]                   │
│ *  --input-contigs          PATH                          Input FASTA [required]       │
│    --output-folder          PATH                          Output directory             │
│                                                           [default: clustered-contigs] │
│    --threads                INTEGER RANGE [x>=0]          0=all cores; otherwise N     │
│                                                           threads                      │
│                                                           [default: 0]                 │
│    --vclust-params          TEXT                          Custom vclust parameters     │
│                                                           like: "--min-kmers 20        │
│                                                           --outfmt lite --ani 0.97"    │
│    --help           -h                                    Show this message and exit.  │
╰────────────────────────────────────────────────────────────────────────────────────────╯
```

## Output Files

For all modes, the following files are generated in the output folder:


```bash
.
├── ani.ids.tsv
├── ani.tsv
├── cluster_representatives_ids.txt
├── fltr.txt
├── representatives.fna
└── species.tsv

1 directory, 6 files
```

| File | Description |
|------|-------------|
| `clusters.tsv` / `species.tsv` | Tab-separated clustering results |
| `representatives.fna` | FASTA file with cluster representative sequences |
| `cluster_representatives_ids.txt` | List of representative sequence IDs |
| `fltr.txt` | Pre-alignment filter results |
| `ani.tsv` | ANI calculation results |
| `ani.ids.tsv` | Sequence ID mappings |

## Examples

### Basic Examples

```bash
# Dereplicate viral contigs
phu seqclust --mode dereplication --input-contigs viral_contigs.fna

# Cluster into vOTUs following MIUViG standards
phu seqclust --mode votu-clustering --input-contigs viral_contigs.fna

# Classify viruses into species following ICTV standards
phu seqclust --mode spp-clustering --input-contigs complete_genomes.fna
```

### Advanced Usage with Custom Parameters

The `--vclust-params` option allows you to customize any vclust parameter while maintaining the convenience of predefined modes. Parameters are automatically routed to the appropriate vclust command (prefilter, align, cluster).

#### Large Dataset Optimization (Wiki Section 6.6)

```bash
# Process large diverse dataset (IMG/VR style)
phu seqclust --mode votu-clustering --input-contigs large_dataset.fna \
  --vclust-params="--min-kmers 4 --batch-size 2000000 --kmers-fraction 0.2 --outfmt lite"
```

#### Highly Redundant Dataset (Wiki Section 6.7)

```bash
# Process highly redundant sequences
phu seqclust --mode votu-clustering --input-contigs redundant_genomes.fna \
  --vclust-params="--min-kmers 10 --batch-size 100000 --max-seqs 1000 --outfmt lite --ani 0.97 --qcov 0.95"
```

#### Custom Thresholds

```bash
# More stringent clustering
phu seqclust --mode votu-clustering --input-contigs genomes.fna \
  --vclust-params="--ani 0.98 --qcov 0.90"

# Species clustering with custom genus threshold
phu seqclust --mode spp-clustering --input-contigs genomes.fna \
  --vclust-params="--tani 0.97"
```

## Supported vclust Parameters

### Prefilter Parameters
- `--min-kmers`: Minimum number of shared k-mers
- `--min-ident`: Minimum sequence identity threshold
- `--batch-size`: Number of sequences processed per batch
- `--kmers-fraction`: Fraction of k-mers to analyze
- `--max-seqs`: Maximum target sequences per query

### Align Parameters
- `--outfmt`: Output format (`lite` for large datasets)
- `--out-ani`: Minimum ANI to report
- `--out-qcov`: Minimum query coverage to report

### Cluster Parameters
- `--ani`, `--tani`, `--gani`: Clustering thresholds for different metrics
- `--qcov`: Query coverage threshold
- `--leiden-resolution`: Leiden algorithm resolution parameter
- `--algorithm`: Override clustering algorithm
- `--metric`: Override distance metric

## Performance Tips

### For Large Datasets (>1M sequences)
```bash
--vclust-params="--batch-size 2000000 --kmers-fraction 0.2 --outfmt lite"
```

### For Highly Redundant Datasets
```bash
--vclust-params="--batch-size 100000 --max-seqs 1000 --kmers-fraction 0.2"
```

### For Memory-Constrained Systems
```bash
--vclust-params="--batch-size 50000 --kmers-fraction 0.1"
```

## Comparison with Direct vclust Usage

| Task | phu seqclust | Direct vclust |
|------|-------------|---------------|
| **Steps** | Single command | 3 separate commands |
| **Configuration** | Preconfigured modes | Manual parameter setup |
| **Customization** | `--vclust-params` option | Full control |
| **Learning curve** | Minimal | Requires vclust expertise |
| **Use case** | Common workflows | Specialized analyses |

## Requirements

- `vclust` (or `vclust.py`) in PATH
- `seqkit` in PATH

## Related

- [vclust documentation](https://github.com/refresh-bio/vclust)
- [vclust wiki use cases](https://github.com/refresh-bio/vclust/wiki/6-Use-cases)
- [MIUViG standards](https://www.nature.com/articles/nbt.4306)
- [ICTV classification](https://ictv.global/)
