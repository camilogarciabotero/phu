
<div align="center">
  <a href="https://anaconda.org/bioconda/phu">
    <img src="https://img.shields.io/conda/vn/bioconda/phu?logo=anaconda&style=flat-square&maxAge=3600" alt="install with bioconda">
  </a>
  <a href="https://anaconda.org/bioconda/phu"> <img src="https://anaconda.org/bioconda/phu/badges/downloads.svg" /> </a>
    <a href="https://github.com/camilogarciabotero/phu/actions/workflows/docs.yaml"><img src="https://github.com/camilogarciabotero/phu/actions/workflows/docs.yaml/badge.svg" alt="docs">
  </a>
  <a href="https://anaconda.org/bioconda/phu"> <img src="https://anaconda.org/bioconda/phu/badges/license.svg" /> </a>
</div>


***
# phu - Phage Utilities

phu (phage utilities) or phutilities, is a modular toolkit for viral genomics workflows. It provides command-line tools to handle common steps in phage bioinformatics pipelines—wrapping complex utilities behind a consistent and intuitive interface.

## Installation

You can install `phu` using `mamba` or `conda` from the `bioconda` channel:

```bash
mamba create -n phu bioconda::phu
```

## Usage

As a command-line tool, `phu` follows a modular structure. You can access different functionalities through subcommands. The general syntax is:

```bash
phu <command> [options]
```

## Commands

- [`cluster`](https://camilogarciabotero.github.io/phu/commands/cluster/): Cluster viral sequences into species or other operational taxonomic units (OTUs).

## Contributing

We welcome contributions to phu! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them.
4. Submit a pull request describing your changes.


## Developers

You can also install the development version of `phu` directly from GitHub:

```bash
git clone https://github.com/camilogarciabotero/phu.git
cd phu
pip install -e .
```

`phu` is also available on PyPI:

```bash
pip install phu
```

## References

This program uses several key tools and libraries, make sure to acknowledge them when using `phu`:

- [vclust](https://github.com/refresh-bio/vclust): A high-performance clustering tool for viral sequences:
> Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. https://doi.org/10.1038/s41592-025-02701-7

- [seqkit](https://bioinf.shenwei.me/seqkit/): A toolkit for FASTA/Q file manipulation.
> Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.