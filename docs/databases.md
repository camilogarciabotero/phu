# Databases and Provenance

`phu screen` and `phu avger` need reference databases: Pfam, KOfam, dbCAN
and the CheckAMG-derived AVG tables (`phu jack` uses none). PHU ships them as one pinned
bundle, `phu-db`, and installs it automatically the first time a command needs
it. You can also install it ahead of time with
[`phu dbs prepare`](#installing-ahead-of-time).

## Installing ahead of time

```bash
phu dbs prepare
```

This downloads the bundle from Zenodo, verifies it and installs it, then
prints where it is and which upstream releases it holds. It is optional, but
useful on a login node before submitting jobs, so that compute jobs neither
download nor wait for each other. If the bundle is already installed, the
command only checks it and returns.

| Option | Effect |
| --- | --- |
| `--bundle PATH` | Install from a pre-downloaded `phu-db.tar.gz` instead of downloading (same as `PHU_DB_BUNDLE`). |
| `--verify` | Re-hash every installed file and reinstall if any check fails (same as `PHU_DB_VERIFY=1`). |
| `--connections N` | Download `N` parts of the bundle at once, from 1 to 16 (default 8; same as `PHU_DB_CONNECTIONS`). |

The bundle always holds all four databases, so database names given to
`phu dbs prepare` are accepted and ignored.

## What happens on first use

`phu dbs prepare`, or the first database-backed command if you skipped it:

1. Checks free disk space. The install needs about **16 GB** free while it
   runs and keeps about **13 GB** afterwards.
2. Downloads `phu-db.tar.gz` (about 2.2 GB) from Zenodo, or uses a local copy
   if [`PHU_DB_BUNDLE`](#offline-and-hpc-use) is set.
3. Verifies the tarball's exact size and SHA-256 against the values pinned in
   this PHU release. A mismatch stops the install; nothing is used.
4. Unpacks it after checking every member (no links, devices or paths that
   escape the install folder), then verifies every database file's size,
   SHA-256 and record count.
5. Moves the finished tree into place in one step, so an interrupted install
   never leaves a half-installed database behind.

PHU prints what it downloads, how large it is and where it goes. The bundle
is downloaded in several parts at once and an interrupted download continues
where each part stopped. While it downloads, a terminal shows a spinner with a progress bar, the speed and the
time left; a log file (for example a SLURM job's) gets one line per 10%. Later runs
only check file sizes, and line counts of the small tables, which takes a
moment, and do not download again.

If several jobs start at once against the same database folder, the first one
installs and the others wait for it, then reuse its result.

## Where databases are installed

| Setting | Location |
| --- | --- |
| `PHU_DB_FOLDER` set | `$PHU_DB_FOLDER` |
| otherwise, `XDG_DATA_HOME` set | `$XDG_DATA_HOME/phu/db` |
| otherwise | `~/.local/share/phu/db` |

Inside that folder:

```text
phu-db/
  phu-db-0.2.0/
    installed.json        # source, tarball SHA-256, install time
    phu-db.json           # the bundle manifest: upstream URLs, digests, steps
    pfam/  profiles.hmm
    kofam/ profiles.hmm, ko_list.tsv
    dbcan/ profiles.hmm, family_substrate.tsv, pul_rules.tsv
    avg/   v_scores.tsv, positive_{amg,apg,areg}.tsv, filters_{amg,apg,areg}.tsv
    phrog/ profiles.hmm, annotations.tsv
```

When a new PHU release pins a newer bundle, it installs the new version beside
the old one. PHU never deletes other versions, because another PHU release
sharing the folder may still be using them; remove an old `phu-db-<version>/`
folder yourself once nothing needs it.

## Environment variables

| Variable | Effect |
| --- | --- |
| `PHU_DB_FOLDER` | Install root (see above). |
| `PHU_DB_BUNDLE` | Path to a pre-downloaded `phu-db.tar.gz`. It is verified exactly like a download, and PHU never deletes it. |
| `PHU_DB_CONNECTIONS` | How many parts of the bundle to download at once, from 1 to 16 (default 8). Zenodo limits the speed of each connection, so several parts at once download several times faster. Use `1` for a single connection. |
| `PHU_DB_VERIFY=1` | Re-hash every installed file of every database on this run instead of checking sizes only, and discard the indexes and per-model files PHU built from them (they are rebuilt on next use). Reads about 13 GB, so use it when you suspect corruption, not routinely. A failed check reinstalls. |

## Offline and HPC use

Compute nodes often have no internet access, or many jobs would download the
same bundle at once. Download it once and point PHU at it:

```bash
# On a node with internet access: get phu-db.tar.gz for phu-db-0.2.0 from
# https://doi.org/10.5281/zenodo.19503347 (the page lists every version).
export PHU_DB_FOLDER=/shared/project/phu-db      # on a shared filesystem
phu dbs prepare --bundle /shared/project/phu-db.tar.gz   # installs once
phu screen PF00001 -i contigs.fna                 # every job reuses the install
```

Install once on a shared filesystem and every job reuses the same tree. The
install lock is released as soon as a job exits, even if it is killed. On
network filesystems without `flock` support, PHU uses a lock directory
instead. A running install re-stamps it every minute, and PHU clears one that
a killed job left without a stamp for six hours.

PHU builds a small index and per-model files for each database the first time
it is queried, and keeps them beside the database files. The first `phu avger`
run also writes a binary copy of each model file it searches
(`profiles.h3m`, about 5 GB for Pfam, KOfam and PHROG together), which loads
in seconds where the text file takes minutes. The copies are rebuilt
automatically if the database or pyhmmer changes; if they cannot be written,
PHU reads the text files instead. If the install is
read-only for you (for example one a colleague or an administrator installed),
PHU keeps them in your own cache folder instead (`PHU_CACHE_DIR`, default
`~/.cache/phu`), so a shared read-only install works as is.

If the tarball does not match the version this PHU release pins, PHU refuses
it and says so. Download the version named in the error.

## What the bundle contains

`phu-db-0.2.0` contains these pinned upstream releases. PHU records them in
each run's provenance.

| Database | Used by | Upstream release | Files |
| --- | --- | --- | --- |
| `pfam` | `screen`, `avger` | Pfam 37.0 (CC0-1.0) | `profiles.hmm` (21,979 models) |
| `kofam` | `screen`, `avger` | KOfam, genome.jp archive 2026-07-02 | `profiles.hmm`, `ko_list.tsv` (KO thresholds) |
| `dbcan` | `screen` | dbCAN HMMdb V15; dbCAN-PUL substrate table 12-12-2023; dbCAN-PUL workbook snapshot of 2026-09-30 | `profiles.hmm`, `family_substrate.tsv`, `pul_rules.tsv` |
| `avg` | `avger` | CheckAMG v1.1.1 reference tables | V-scores, AMG/APG/AReG positive tables and filter tables |
| `phrog` | `avger` | PHROGs v4 (CC-BY-4.0) | `profiles.hmm` (38,880 models), `annotations.tsv` |

PHROGs is published as HH-suite profiles and MMseqs2 databases, which HMMER
cannot search. The PHROG models in the bundle were built from the official
per-PHROG alignments with pyhmmer, using the default `hmmbuild` settings, one
model per PHROG, named `phrog_1` to `phrog_38880`.

The bundle's own manifest, `phu-db.json`, lists every upstream URL, SHA-256,
byte size, license and record count, and the exact preparation steps that turn
each upstream download into the installed file.

## Provenance

Each run writes the database provenance into the `.phu/` folder of its output
folder: the bundle version, the tarball SHA-256 and each database's upstream
release. `avger`, `screen --all-cazymes` and `screen --all-puls` write it to
`.phu/run.json`; other `screen` queries that name Pfam, KOfam or dbCAN models
write it to `.phu/query.json`, with the SHA-256 of each database file used.
Keep these files with your results.

## Citation

Cite the bundle version you used:

> García-Botero C. *Phu Databases*, version phu-db-0.2.0. Zenodo.
> <https://doi.org/10.5281/zenodo.19503347>

Also cite the upstream databases your analysis relied on, following each
project's own citation guidance. Pfam: Mistry J, et al. *Nucleic Acids Res.*
2021;49(D1):D412–D419. For KOfam, dbCAN and CheckAMG, see their project pages.

## Upgrading from PHU 0.9.x and earlier

Earlier releases downloaded each database separately and managed them with
`phu dbs`. `phu dbs prepare` remains and now installs the single bundle. The
`list`, `status`, `refresh` and `remove` subcommands are gone; they only print
a notice pointing to `phu dbs prepare`.

PHU no longer reads the old per-database folders in `PHU_DB_FOLDER` (`pfam/`,
`kofam/`, `dbcan/`, `avg/` and `vscore/`) and never deletes them. Once the
bundle is installed you can remove them yourself to free about 10 GB; the
bundle lives in `phu-db/` beside them.
