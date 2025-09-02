# phu




phu (phage utilities) or phutilities, is a modular toolkit for viral genomics workflows. It provides command-line tools to handle common steps in phage bioinformatics pipelines—wrapping complex utilities behind a consistent and intuitive interface.

The first implemented module, seqclust, brings flexible clustering strategies for viral contigs and genomes, with additional modules under active development.

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

- `seqclust`: Cluster viral sequences into operational taxonomic units (OTUs) or other current standard (see this [wiki](https://github.com/refresh-bio/vclust/wiki/6-Use-cases))

## Contributing

We welcome contributions to phu! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them.
4. Submit a pull request describing your changes.
