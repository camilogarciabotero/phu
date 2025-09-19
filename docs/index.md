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
                                                                                              
╭─ Options ──────────────────────────────────────────────────────────────────────────────────╮
│ --help  -h        Show this message and exit.                                              │
╰────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────────────────────────────────────╮
│ screen          Screen contigs for a protein family using HMMER on predicted CDS.          │
│ cluster         Sequence clustering wrapper around external 'vclust' with three modes.     │
│ simplify-taxa   Simplify vContact taxonomy prediction columns into compact lineage codes.  │
╰────────────────────────────────────────────────────────────────────────────────────────────╯
```

## Available Commands

### `screen` - Protein Family Screening

Screen DNA contigs for specific protein families using HMMER on predicted coding sequences. This is particularly useful for identifying viral contigs in metagenomic assemblies or filtering assemblies based on protein content.

**Example:**
```bash
 phu screen --input-contigs assembly.fasta --combine-mode all viral_capsid.hmm portal.hmm
```

[Learn more about protein screening →](commands/screen.md)

### `cluster` - Sequence Clustering


Cluster viral sequences into operational taxonomic units with three specialized modes:

- **`dereplication`** - Remove duplicate sequences
- **`votu`** - Group sequences into viral Operational Taxonomic Units
- **`species`** - Create species-level clusters

**Example:**
```bash
phu cluster --mode votu --input-contigs viral-genomes.fasta
```

[Learn more about clustering →](commands/cluster.md)


### `simplify-taxa` - Taxonomy Simplification

Simplify verbose vContact taxonomy predictions into compact lineage codes for easier analysis and visualization.
**Example:**
```bash
phu simplify-taxa -i final_assignments.csv -o simplified_taxonomy.csv
```

[Learn more about taxonomy simplification →](commands/simplify-taxa.md)

## Use Cases

- **Viral Identification**: Screen metagenomic assemblies for viral contigs using protein markers
- **Multi-marker Analysis**: Find contigs with complete sets of viral proteins (e.g., capsid, portal, primase, terminase)
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

## References

This program uses several key tools and libraries, make sure to acknowledge them when using `phu`:

- [vclust](https://github.com/refresh-bio/vclust): A high-performance clustering tool for viral sequences:
> Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. https://doi.org/10.1038/s41592-025-02701-7

- [seqkit](https://bioinf.shenwei.me/seqkit/): A toolkit for FASTA/Q file manipulation.
> Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.

- [Prodigal](https://github.com/hyattpd/prodigal): A gene prediction tool for prokaryotic genomes.
> Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119

- [pyrodigal](https://pyrodigal.readthedocs.io/en/stable/): A tool for gene prediction in prokaryotic genomes.
> Larralde, M., (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software, 7(72), 4296, https://doi.org/10.21105/joss.04296

- [HMMER](http://hmmer.org/): A suite of tools for sequence analysis using profile hidden Markov models.
> Eddy, S. R. (2011). Accelerated Profile HMM Searches. PLoS Computational Biology, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195

- [pyHMMER](https://pyhmmer.readthedocs.io/en/latest/): Python bindings for HMMER.
> Larralde, M., & Zeller, G. (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad214