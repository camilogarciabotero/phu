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

## Usage

As a command-line tool, `phu` follows a modular structure. You can access different functionalities through subcommands. The general syntax is:

```bash
phu <command> [options]
```

## Commands

- [`screen`](https://camilogarciabotero.github.io/phu/commands/screen/): Screen contigs for specific protein families using HMMER on predicted coding sequences.
- [`jack`](https://camilogarciabotero.github.io/phu/commands/jack/): Iteratively screen contigs from one or more seed proteins with jackhmmer and combine seeds hits.
- [`cluster`](https://camilogarciabotero.github.io/phu/commands/cluster/): Cluster viral sequences into species or other operational taxonomic units (OTUs).
- [`avger`](https://camilogarciabotero.github.io/phu/commands/avger/): Predict proteins, annotate them against Pfam and KOfam, and score putative auxiliary viral gene candidates.
- [`simplify-taxa`](https://camilogarciabotero.github.io/phu/commands/simplify-taxa/): Simplify vContact taxonomy prediction columns into compact lineage codes.

## Database Management

A special command group, `dbs`, is available to manage local databases used by `phu`. Each database can define its own preparation logic while sharing a common user interface. Built-in backends currently include `pfam`, `kofam`, `dbcan`, `vscore`, and `avg`.

- [`dbs`](https://camilogarciabotero.github.io/phu/commands/dbs/): Manage local databases (list, status, prepare, refresh, remove) for `pfam`, `kofam`, `dbcan`, `vscore`, and `avg`.

For screening pass/fail and threshold rules, see [`screen thresholds and decision logic`](https://camilogarciabotero.github.io/phu/commands/screen-thresholds/).

### AVG workflows

`phu avger` is the annotation and evidence workflow. It predicts proteins from
trusted viral contigs, searches complete Pfam and KOfam databases, adds optional
V-Score-Search annotations, and writes one best passing hit per protein and
database. It reports candidate evidence; it does not establish that a gene is
an AVG by itself.

`phu avger` applies the V/VL/AVL score framework of Zhou et al. and uses
versioned CheckAMG reference tables for AMG, APG, and AReG classification and
curation; see the [avger command guide](https://camilogarciabotero.github.io/phu/commands/avger/)
for attribution and interpretation limits.

The `avg` development work is currently implemented internally and is not a
public `phu` subcommand. Use `phu avger` for the supported annotation and
evidence workflow described in the [avger command guide](https://camilogarciabotero.github.io/phu/commands/avger/).

## Cache Handling

`phu` caches predicted proteins for both `screen` and `jack` so repeated runs can reuse the same translated proteins when the prediction inputs have not changed. Search settings such as HMM files, seed markers, combine mode, and output folder do not affect the cache.

The cache is rebuilt when you change the contig input, `--mode`, `--ttable`, or the protein-length filter. For `phu screen`, that is `--min-protein-len-aa`. For `phu jack`, both `--min-gene-len` and `--min-protein-len-aa` participate in the cache key.

To remove previously cached predictions, run `phu --clean-cache`.

See the full cache guide in [Cache Handling](https://camilogarciabotero.github.io/phu/cache).

## Contributing

We welcome contributions to phu. See [CONTRIBUTING.md](CONTRIBUTING.md) for
testing and documentation requirements.

tracked in [CHANGELOG.md](CHANGELOG.md).

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

- [V-Score-Search](https://github.com/AnantharamanLab/V-Score-Search): A framework for identifying viral sequences and auxiliary viral genes.
> Zhou et al. (2024). V-Score-Search. Preprint. DOI: 10.1101/2024.10.24.619987

- [CheckAMG](https://github.com/AnantharamanLab/CheckAMG): Reference data and curation resources for auxiliary metabolic and auxiliary viral gene analysis.

- [Martin et al. (2025)](https://doi.org/10.1038/s41564-025-02095-4): A call for caution in the biological interpretation of viral auxiliary metabolic genes. *Nature Microbiology*, 10, 2122-2129.