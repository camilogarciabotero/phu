# phu - Phage Utilities

<div align="center">
  <strong>A modular toolkit for viral genomics workflows</strong>
</div>

---

## What is phu?

**phu** (phage utilities) is a command-line toolkit designed to streamline viral genomics workflows. It provides intuitive commands that wrap complex bioinformatics utilities behind a consistent interface, making phage and viral sequence analysis more accessible and reproducible.

### Key Features

- **Modular Design**: Clean, focused commands for specific tasks.
- **Easy Installation**: Available through Bioconda, PyPI.
- **Reproducible**: Consistent interface across different utilities.
- **Well Documented**: Comprehensive documentation and examples.

## Quick Start

### Installation

Install `phu` using mamba or conda from the bioconda channel:

```bash
mamba create -n phu bioconda::phu
mamba activate phu
```

### Basic Usage

```bash
phu -h

 Usage: phu [OPTIONS] COMMAND [ARGS]...                                      
                                                                             
 Phage utilities CLI                                                         
                                                                             
╭─ Options ─────────────────────────────────────────────────────────────────╮
│ --help  -h        Show this message and exit.                             │
╰───────────────────────────────────────────────────────────────────────────╯
╭─ Commands ────────────────────────────────────────────────────────────────╮
│ cluster   Sequence clustering wrapper around external 'vclust' with three │
│           modes.                                                          │
╰───────────────────────────────────────────────────────────────────────────╯
```

## Available Commands

### `cluster` - Sequence Clustering

Cluster viral sequences into operational taxonomic units with three specialized modes:

- **`dereplication`** - Remove duplicate sequences
- **`votu`** - Group sequences into viral Operational Taxonomic Units
- **`species`** - Create species-level clusters

**Example:**
```bash
phu cluster --mode votu --input-contigs viral_genomes.fasta --output-folder results/
```

[Learn more about clustering →](commands/cluster.md)


## Use Cases

- **Viral Metagenomics**: Dereplicate and cluster viral contigs from metagenomic assemblies
- **Phage Genomics**: Organize phage genomes into taxonomic groups
- **Comparative Analysis**: Prepare datasets for phylogenetic and comparative genomic studies
- **Database Construction**: Build reference databases of viral sequences

## Contributing

We welcome contributions! Whether it's bug reports, feature requests, or code contributions, please check out our [GitHub repository](https://github.com/camilogarciabotero/phu).

## Citation

If you use phu in your research, please cite:

```
García-Botero, C. (2025). phu: Phage Utilities - A modular toolkit for viral genomics workflows. 
GitHub repository: https://github.com/camilogarciabotero/phu
```

