# phu - Phage Utilities

<div align="center">
  <img src="assets/phu-logo-black.svg" alt="phu logo" width="180">
</div>

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
phu --version
phu --help
```

```text
Usage: phu [OPTIONS] COMMAND [ARGS]...

 Phage utilities CLI

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --version      -v        Show version and exit.                              │
│ --clean-cache            Remove cached protein predictions and exit.         │
│ --help         -h        Show this message and exit.                         │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Workflow ───────────────────────────────────────────────────────────────────╮
│ cluster            Sequence clustering wrapper around external 'vclust' with │
│                    three modes.                                              │
│ normalize-lineage  Normalize vContact lineage predictions without changing   │
│                    column names.                                             │
│ avger              Predict and curate putative auxiliary viral genes.        │
│ screen             Screen contigs using protein-family or PUL CAZyme         │
│                    signatures.                                               │
│ jack               Iteratively screen contigs from one or more seed protein  │
│                    markers with                                              │
│                    pyhmmer.jackhmmer.                                        │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Database Management ────────────────────────────────────────────────────────╮
│ dbs                Install and check the phu-db database bundle              │
╰──────────────────────────────────────────────────────────────────────────────╯
```


## Available Commands

### `screen` - Protein Family Screening

Screen DNA contigs for specific protein families using HMMER on predicted coding sequences. This is particularly useful for identifying viral contigs in metagenomic assemblies or filtering assemblies based on protein content.

**Example:**
```bash
 phu screen --input-contigs assembly.fasta --combine-mode all viral_capsid.hmm portal.hmm
```

[Learn more about protein screening →](commands/screen.md)

For pass/fail and threshold rules (PFAM and KOfam), see [Screen Thresholds and Decision Logic](screen-thresholds.md).

### `jack` - Seed-Based Iterative Screening

Iteratively screen DNA contigs from one or more seed marker proteins using `pyhmmer.hmmer.jackhmmer`.
When multiple seeds are provided, you can combine evidence with `--combine-mode any|all|threshold`.

**Example:**
```bash
phu jack --input-contigs assembly.fasta marker_seeds.faa --combine-mode all
```

[Learn more about iterative jack screening →](commands/jack.md)

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


### Databases

`screen` and `avger` use Pfam, KOfam, dbCAN and the CheckAMG AVG tables. PHU
ships them as one pinned and checksummed bundle, `phu-db`, which installs
automatically on first use and shows its download progress.

**Example:**
```bash
phu dbs prepare                              # install ahead of time
phu dbs prepare --bundle phu-db.tar.gz       # from a pre-downloaded tarball
```

Set `PHU_DB_FOLDER` to a shared folder so every job reuses one install.

[Learn more about databases and provenance →](databases.md)


### `normalize-lineage` - Normalize vContact Lineage

Normalize vContact lineage predictions while preserving the input table schema.

[Learn more about lineage normalization →](commands/normalize-lineage.md)

## Use Cases

- **Viral Identification**: Screen metagenomic assemblies for viral contigs using protein markers
- **Multi-marker Analysis**: Find contigs with complete sets of viral proteins (e.g., capsid, portal, primase, terminase)
- **Viral Metagenomics**: Dereplicate and cluster viral contigs from metagenomic assemblies
- **Phage Genomics**: Organize phage genomes into taxonomic groups
- **Comparative Analysis**: Prepare datasets for phylogenetic and comparative genomic studies
- **Database Construction**: Build reference databases of viral sequences

## Cache Handling

`phu` caches predicted proteins for both `screen` and `jack`. That means repeated searches can skip the gene-prediction step when the contigs and prediction inputs are unchanged. Changing HMMs, seed markers, combine mode, or output folders does not invalidate the cache.

The cache is rebuilt when you change the contigs, `--mode`, an explicit `--ttable`, or the protein-length filter. Without `--ttable`, meta-mode contigs use the translation table selected by pyrodigal-gv. For `phu screen`, that is `--min-protein-len-aa`. For `phu jack`, both `--min-gene-len` and `--min-protein-len-aa` affect cache reuse.

To clear previous predictions manually, run `phu --clean-cache`.

See the full cache guide in [Cache Handling](https://camilogarciabotero.github.io/phu/cache).

## Contributing

See the repository's [contribution guide](https://github.com/camilogarciabotero/phu/blob/avg/CONTRIBUTING.md),
and [changelog](https://github.com/camilogarciabotero/phu/blob/avg/CHANGELOG.md)
before opening a change.

## Citation

Citation metadata for this development revision is maintained in
[`CITATION.cff`](https://github.com/camilogarciabotero/phu/blob/avg/CITATION.cff).
See [Citation](citation.md) for the distinction between development citations,
the Zenodo concept DOI, and immutable release DOIs.

## References

This program uses several key tools and libraries, make sure to acknowledge them when using `phu`:

- [vclust](https://github.com/refresh-bio/vclust): A high-performance clustering tool for viral sequences:
> Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. https://doi.org/10.1038/s41592-025-02701-7

- [seqkit](https://bioinf.shenwei.me/seqkit/): A toolkit for FASTA/Q file manipulation.
> Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.

- [FAMSA](https://github.com/refresh-bio/FAMSA) and [pyfamsa](https://github.com/althonos/pyfamsa): A multiple sequence aligner for large protein families, and its Python bindings.
> Deorowicz, S., Debudaj-Grabysz, A., & Gudyś, A. (2016). FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific Reports, 6, 33964. https://doi.org/10.1038/srep33964

- [trimAl](https://github.com/inab/trimal) and [pytrimal](https://github.com/althonos/pytrimal): Automated alignment trimming, and its Python bindings.
> Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348

- [Prodigal](https://github.com/hyattpd/prodigal): A gene prediction tool for prokaryotic genomes.
> Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119

- [pyrodigal](https://pyrodigal.readthedocs.io/en/stable/): A tool for gene prediction in prokaryotic genomes.
> Larralde, M., (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software, 7(72), 4296, https://doi.org/10.21105/joss.04296

- [HMMER](http://hmmer.org/): A suite of tools for sequence analysis using profile hidden Markov models.
> Eddy, S. R. (2011). Accelerated Profile HMM Searches. PLoS Computational Biology, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195

- [pyHMMER](https://pyhmmer.readthedocs.io/en/latest/): Python bindings for HMMER.
> Larralde, M., & Zeller, G. (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad214

- [Pfams](https://pfam.xfam.org/): A comprehensive collection of protein families and domains based on hidden Markov models.
> Pfam: The protein families database in 2021: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman, Nucleic Acids Research (2021) doi: 10.1093/nar/gkaa913