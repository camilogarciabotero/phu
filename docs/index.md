# phu - Phage Utilities

<div align="center">
  <img src="assets/phu-logo-black.svg" alt="phu logo" width="180">
</div>

<div align="center">
  <strong>A modular toolkit for viral genomics workflows</strong>
</div>

---

**phu** (phage utilities) is a command-line toolkit for working with phage
and other viral sequences. Each command answers one common question about
your contigs, from "which of them carry this protein family?" to "which of
them belong to the same species?", and records what it did so the result can
be reproduced.

- **One install, no external programs** for searching: gene prediction, HMM
  searches and alignments run inside phu (only `cluster` calls `vclust`).
- **Reference databases included**: Pfam, KOfam, PHROG, ECOD, dbCAN and the
  CheckAMG tables come as one pinned bundle that installs itself on first
  use.
- **Reproducible**: every run writes the exact query, thresholds, database
  versions and file checksums to a `.phu/` folder beside its results.

## Install

```bash
mamba create -n phu bioconda::phu
mamba activate phu
phu --help
```

The databases download on first use (about 2.9 GB). On a cluster, install
them once on a shared filesystem first; see [Databases](databases.md).

## What do you want to do?

| I want to… | Command | Example |
| --- | --- | --- |
| Find contigs that carry given protein families | [`screen`](commands/screen.md) | `phu screen -i contigs.fa PF00589 phrog_2 ecod:2004.1` |
| Find contigs carrying relatives of proteins I already have | [`jack`](commands/jack.md) | `phu jack -i contigs.fa seeds.faa` |
| Look up PHROG families by name or function | [`phrogs`](commands/phrogs.md) | `phu phrogs search "tail fiber"` |
| Find putative auxiliary viral genes (AVGs) | [`avger`](commands/avger.md) | `phu avger -i viral-contigs.fa` |
| Dereplicate, or group genomes into vOTUs, genera or species | [`cluster`](commands/cluster.md) | `phu cluster --mode votu -i genomes.fa` |
| Tidy vContact lineage predictions | [`normalize-lineage`](commands/normalize-lineage.md) | |
| Install the databases ahead of time | [`dbs prepare`](databases.md#installing-ahead-of-time) | `phu dbs prepare` |
| Free disk space | `clean` | `phu clean --cache` |

## The commands

### `screen`: find contigs by their protein families

Predicts the proteins of your contigs, searches them with profile HMMs, and
keeps the contigs that match. Queries can be your own HMM files or families
named from the databases: Pfam (`PF00589`), KOfam (`K00001`), PHROG
(`phrog_216`), ECOD domains at any level of its hierarchy (`ecod:2004.1`),
dbCAN CAZymes (`GH128`) or a dbCAN PUL (`PUL0621`). Each database's own cutoffs
apply; combine queries with `--combine-mode any|all|threshold`.

```bash
phu screen -i assembly.fasta --combine-mode all capsid.hmm portal.hmm
```

[Screen guide →](commands/screen.md) · [How hits are judged →](screen-thresholds.md)

### `jack`: find relatives of seed proteins

Builds a profile from each seed protein by iterative search
(`pyhmmer.jackhmmer`) and keeps the contigs that carry relatives. Seeds can be
your own proteins or PHROG families.

```bash
phu jack -i assembly.fasta marker_seeds.faa --combine-mode all
phu jack -i assembly.fasta --phrog alias:portal
```

[Jack guide →](commands/jack.md)

### `phrogs`: find PHROG families by name

Turns names into PHROG identifiers for `screen` and `jack`: curated aliases
(`alias:mcp`, `alias:terl`), functional categories (`category:tail`) or a
forgiving search of the annotations.

```bash
phu phrogs search "tail fiber"
phu screen -i assembly.fasta $(phu phrogs resolve alias:mcp alias:terl --ids-only)
```

[PHROG lookup guide →](commands/phrogs.md)

### `avger`: auxiliary viral genes

Annotates the proteins of viral contigs against Pfam, KOfam and PHROG, and
applies CheckAMG-style rules to flag putative auxiliary metabolic,
physiology and regulatory genes.

```bash
phu avger -i viral-contigs.fasta -o phu-avger
```

[AVG guide →](commands/avger.md)

### `cluster`: group sequences

Clusters viral sequences with vclust in four modes: `dereplication`, `votu`,
`genus` (complete genomes, tANI ≥ 70%) and `species` (tANI ≥ 95%).

```bash
phu cluster --mode votu -i viral-genomes.fasta
```

[Clustering guide →](commands/cluster.md)

### `normalize-lineage`: tidy vContact lineages

Normalizes vContact lineage predictions while keeping the input table's
columns. [Lineage guide →](commands/normalize-lineage.md)

### `dbs` and `clean`: databases and disk space

`phu dbs prepare` installs the database bundle ahead of time;
`phu clean` lists what phu keeps between runs (cached protein predictions,
installed databases) and removes it on request.

```bash
phu dbs prepare                          # install the databases now
phu dbs prepare --bundle phu-db.tar.gz   # from a pre-downloaded copy
phu clean                                # show what phu keeps and its size
phu clean --cache                        # remove cached predictions
```

Set `PHU_DB_FOLDER` to a shared folder so every job reuses one install.
[Databases →](databases.md) · [Cache →](cache.md)

## All commands

```text
                                                                                
 Usage: phu [OPTIONS] COMMAND [ARGS]...                                         
                                                                                
 Phage utilities CLI                                                            
                                                                                
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --version  -v        Show version and exit.                                  │
│ --help     -h        Show this message and exit.                             │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Main ───────────────────────────────────────────────────────────────────────╮
│ screen             Find contigs that carry given protein families.           │
│ jack               Find contigs carrying relatives of seed proteins, by      │
│                    iterative search.                                         │
│ avger              Predict and curate putative auxiliary viral genes.        │
│ cluster            Sequence clustering wrapper around external 'vclust'.     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Database ───────────────────────────────────────────────────────────────────╮
│ dbs                Install and check the phu-db database bundle              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Utils ──────────────────────────────────────────────────────────────────────╮
│ phrogs             Find PHROG identifiers by alias, functional category or   │
│                    annotation, for use with `phu screen`.                    │
│ clean              Remove phu's cache and/or its installed databases.        │
│ normalize-lineage  Normalize vContact lineage predictions without changing   │
│                    column names.                                             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Learn more

- [Databases and provenance](databases.md): what the bundle holds, offline
  and HPC installs, and what each run records.
- [Cache handling](cache.md): when predicted proteins are reused.
- [Troubleshooting](troubleshooting.md): install problems and run logs.
- [Citation](citation.md): how to cite phu and the databases.
- [Contributing](https://github.com/camilogarciabotero/phu/blob/main/CONTRIBUTING.md)
  and the [changelog](https://github.com/camilogarciabotero/phu/blob/main/CHANGELOG.md).

## References

phu builds on these tools and databases; please cite the ones your analysis used:

- [vclust](https://github.com/refresh-bio/vclust): A high-performance clustering tool for viral sequences:
> Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. https://doi.org/10.1038/s41592-025-02701-7

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

- [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/): A comprehensive collection of protein families and domains based on hidden Markov models.
> Pfam: The protein families database in 2021: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman, Nucleic Acids Research (2021) doi: 10.1093/nar/gkaa913

- [PHROGs](https://phrogs.lmge.uca.fr/): Prokaryotic virus remote homologous groups.
> Terzian P, et al. (2021). PHROG: families of prokaryotic virus proteins clustered using remote homology. NAR Genomics and Bioinformatics, 3(3), lqab067. https://doi.org/10.1093/nargab/lqab067

- [ECOD](http://prodata.swmed.edu/ecod/): Evolutionary Classification of protein Domains.
> Schaeffer RD, et al. (2025). Nucleic Acids Research, 53(D1), D411–D418. https://doi.org/10.1093/nar/gkae1029
