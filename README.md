# phu

<div align="center">
  <a href="http://bioconda.github.io/recipes/phu/README.html">
    <img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="install with bioconda">
  </a>
</div>

phu (phage utilities) or phutilities, is a modular toolkit for viral genomics workflows. It provides command-line tools to handle common steps in phage bioinformatics pipelines—wrapping complex utilities behind a consistent and intuitive interface.

## Installation

> [!WARNING] 
    `phu` is currently in the process of being published on Bioconda. The package may not be immediately available. Please check back soon or follow this repository for updates.

phu will be available through Bioconda. To install (once available), use the following command:

```bash
mamba create -n phu bioconda::phu
```

Meanwhile you can install `phu`tilities using `pip` in an a custom `conda`/`mamba` environment with the current software requirments:

```bash
mamba create -n phu -c bioconda seqkit vclust 
mamba activate phu
python -m pip install phu
```

## Usage

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
