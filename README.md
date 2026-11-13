<div align="center">
  <img src="docs/assets/phu-logo-gray.svg" height="150" width="120"><br/>
  <i> Combating Phage Genomes</i><br/><br/>
</div>


<div align="center">
  <a href="https://anaconda.org/bioconda/phu">
    <img src="https://img.shields.io/conda/vn/bioconda/phu?logo=anaconda&style=flat-square&maxAge=3600" alt="install with bioconda">
  </a>
  <a href="https://anaconda.org/bioconda/phu"> <img src="https://anaconda.org/bioconda/phu/badges/downloads.svg" /> </a>
    <a href="https://github.com/camilogarciabotero/phu/actions/workflows/docs.yaml"><img src="https://github.com/camilogarciabotero/phu/actions/workflows/docs.yaml/badge.svg" alt="docs">
  </a>
  <a href="https://anaconda.org/bioconda/phu"> <img src="https://anaconda.org/bioconda/phu/badges/license.svg" /> </a>
  <a href="https://doi.org/10.5281/zenodo.17180799"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17180799.svg" alt="DOI"></a>
</div>


***
# phu - Phage Utilities

phu (phage utilities) or phutilities, is a modular toolkit for viral genomics workflows. It provides command-line tools to handle common steps in phage bioinformatics pipelines—wrapping complex utilities behind a consistent and intuitive interface.

## Installation

You can install `phu` using `mamba` or `conda` from the `bioconda` channel:

```bash
mamba create -n phu bioconda::phu
```

To try unreleased changes from the `main` branch without installing them using
[uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uvx --from git+https://github.com/camilogarciabotero/phu.git@main phu --help
```


## Usage

As a command-line tool, `phu` follows a modular structure. You can access different functionalities through subcommands. The general syntax is:

```bash
phu <command> [options]
```

For the complete command list and global options, run:

```bash
phu --help
```

The top-level help menu is:

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

## Commands

- [`screen`](https://camilogarciabotero.github.io/phu/commands/screen/): Screen contigs for specific protein families using HMMER on predicted coding sequences.
- [`jack`](https://camilogarciabotero.github.io/phu/commands/jack/): Iteratively screen contigs from one or more seed proteins with jackhmmer and combine seeds hits.
- [`cluster`](https://camilogarciabotero.github.io/phu/commands/cluster/): Cluster viral sequences into species or other operational taxonomic units (OTUs).
- [`avger`](https://camilogarciabotero.github.io/phu/commands/avger/): Predict proteins, annotate them against Pfam and KOfam, and score putative auxiliary viral genes.
- [`normalize-lineage`](https://camilogarciabotero.github.io/phu/commands/normalize-lineage/): Normalize vContact lineages while preserving input columns.
- [`dbs prepare`](https://camilogarciabotero.github.io/phu/databases/#installing-ahead-of-time): Download, verify and install the reference database bundle ahead of time.

Each command page includes the current CLI help menu, input/output behavior,
examples, and workflow guidance.

## Databases

`screen` and `avger` use reference databases: Pfam, KOfam, dbCAN and the CheckAMG AVG tables. `phu` ships them as one pinned and checksummed bundle, [`phu-db`](https://doi.org/10.5281/zenodo.19503347), hosted on Zenodo. The bundle installs automatically the first time a command needs it, or ahead of time with:

```bash
phu dbs prepare
```

The download is about 1.9 GB and shows its progress; the install needs about 12 GB of free disk while it runs and keeps about 9.5 GB. Useful settings:

| Setting | Effect |
| --- | --- |
| `PHU_DB_FOLDER` | Where the databases are installed (default `~/.local/share/phu/db`). Point it at a shared folder so every job reuses one install. |
| `PHU_DB_BUNDLE` or `phu dbs prepare --bundle PATH` | Install from a pre-downloaded `phu-db.tar.gz`, for nodes without internet access. |
| `PHU_DB_CONNECTIONS` or `phu dbs prepare --connections N` | Parts of the bundle downloaded at once (default 8). |
| `PHU_DB_VERIFY=1` or `phu dbs prepare --verify` | Re-hash every installed file and reinstall if any check fails. |

Each run records the bundle version and database releases it used in the `.phu/` folder of its output. See [Databases and Provenance](https://camilogarciabotero.github.io/phu/databases/) for offline and HPC use, what the bundle contains, and citation.

For screening pass/fail and threshold rules, see [`screen thresholds and decision logic`](https://camilogarciabotero.github.io/phu/commands/screen-thresholds/).

## Cache Handling

`phu` caches predicted proteins for both `screen` and `jack` so repeated runs can reuse the same translated proteins when the prediction inputs have not changed. Search settings such as HMM files, seed markers, combine mode, and output folder do not affect the cache.

The cache is rebuilt when you change the contig input, `--mode`, an explicit `--ttable`, or the protein-length filter. Without `--ttable`, meta-mode contigs use the translation table selected by pyrodigal-gv. For `phu screen`, that is `--min-protein-len-aa`. For `phu jack`, both `--min-gene-len` and `--min-protein-len-aa` participate in the cache key.

To remove previously cached predictions, run `phu --clean-cache`.

See the full [Cache Handling guide](https://camilogarciabotero.github.io/phu/cache/), including cache locations, disabling reuse, and cleanup behavior.

## Contributing

We welcome contributions to phu. See [CONTRIBUTING.md](CONTRIBUTING.md) for
testing and documentation requirements. Notable changes are tracked in
[CHANGELOG.md](CHANGELOG.md).

## Developers

You can also install the development version of `phu` directly from GitHub with
[uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
git clone https://github.com/camilogarciabotero/phu.git
cd phu
uv sync
```

`phu` is also available on PyPI and can be installed with
[uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uv tool install phu
```

Run the test suite and documentation build from a synced development checkout:

```bash
uv sync
uv run pytest
uv run mkdocs build --strict
```

The project uses `ruff` for linting when the lint dependency group is installed:

```bash
uv run ruff check src tests
```

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

- [KOFamKOALA](https://www.genome.jp/tools/kofamkoala/): KEGG ortholog assignment based on profile HMM and adaptive score threshold.
> Aramaki T., Blanc-Mathieu R., Endo H., Ohkubo K., Kanehisa M., Goto S., Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2019 Nov 19. pii: btz859. doi: 10.1093/bioinformatics/btz859.

- [V-Score-Search](https://github.com/AnantharamanLab/V-Score-Search): A framework for identifying viral sequences and auxiliary viral genes and [CheckAMG](https://github.com/AnantharamanLab/CheckAMG): Reference data and curation resources for auxiliary metabolic and auxiliary viral gene analysis.

> [Martin et al. (2025)](https://doi.org/10.1038/s41564-025-02095-4): A call for caution in the biological interpretation of viral auxiliary metabolic genes. *Nature Microbiology*, 10, 2122-2129. doi: 10.1038/s41564-025-02095-4

> [Zhou et al. (2026)](https://doi.org/10.1038/s41467-026-72028-0). V- and VL-scores unveil viral signatures and origins of protein families. Nat Commun 17, 5828 (2026) doi: 10.1038/s41467-026-72028-0
