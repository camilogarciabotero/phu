# Databases and Provenance

`phu screen` and `phu avger` need reference databases: Pfam, KOfam, dbCAN
and the CheckAMG-derived AVG tables (`phu jack` uses none). PHU ships them as one pinned
bundle, `phu-db`, and installs it automatically the first time a command needs
it. There is nothing to prepare by hand.

## What happens on first use

The first database-backed command:

1. Checks free disk space. The install needs about **12 GB** free while it
   runs and keeps about **9.5 GB** afterwards.
2. Downloads `phu-db.tar.gz` (about 1.9 GB) from Zenodo, or uses a local copy
   if [`PHU_DB_BUNDLE`](#offline-and-hpc-use) is set.
3. Verifies the tarball's exact size and SHA-256 against the values pinned in
   this PHU release. A mismatch stops the install; nothing is used.
4. Unpacks it after checking every member (no links, devices or paths that
   escape the install folder), then verifies every database file's size,
   SHA-256 and record count.
5. Moves the finished tree into place in one step, so an interrupted install
   never leaves a half-installed database behind.

PHU prints what it downloads, how large it is and where it goes. Later runs
only check file sizes and line counts, which takes a moment, and do not
download again.

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
  phu-db-0.1.1/
    installed.json        # source, tarball SHA-256, install time
    phu-db.json           # the bundle manifest: upstream URLs, digests, steps
    pfam/  profiles.hmm
    kofam/ profiles.hmm, ko_list.tsv
    dbcan/ profiles.hmm, family_substrate.tsv, pul_rules.tsv
    avg/   v_scores.tsv, positive_{amg,apg,areg}.tsv, filters_{amg,apg,areg}.tsv
```

When a new PHU release pins a newer bundle, it installs the new version beside
the old one and removes the old one only after the new install succeeds.

## Environment variables

| Variable | Effect |
| --- | --- |
| `PHU_DB_FOLDER` | Install root (see above). |
| `PHU_DB_BUNDLE` | Path to a pre-downloaded `phu-db.tar.gz`. It is verified exactly like a download, and PHU never deletes it. |
| `PHU_DB_VERIFY=1` | Re-hash every installed file on this run instead of checking sizes only. Reads about 9.5 GB, so use it when you suspect corruption, not routinely. A failed check reinstalls. |

## Offline and HPC use

Compute nodes often have no internet access, or many jobs would download the
same bundle at once. Download it once and point PHU at it:

```bash
# On a node with internet access: get phu-db.tar.gz for phu-db-0.1.1 from
# https://doi.org/10.5281/zenodo.19503347 (the page lists every version).
export PHU_DB_FOLDER=/shared/project/phu-db      # on a shared filesystem
export PHU_DB_BUNDLE=/shared/project/phu-db.tar.gz
phu screen PF00001 -i contigs.fna                 # installs once, then runs
```

Install once on a shared filesystem and every job reuses the same tree. The
install lock is released as soon as a job exits, even if it is killed. On
network filesystems without `flock` support, PHU uses a lock directory
instead, and clears one left by a killed job after six hours.

If the tarball does not match the version this PHU release pins, PHU refuses
it and says so. Download the version named in the error.

## What the bundle contains

`phu-db-0.1.1` contains these pinned upstream releases. PHU records them in
each run's provenance.

| Database | Used by | Upstream release | Files |
| --- | --- | --- | --- |
| `pfam` | `screen`, `avger` | Pfam 37.0 (CC0-1.0) | `profiles.hmm` (21,979 models) |
| `kofam` | `screen`, `avger` | KOfam, genome.jp archive 2026-07-02 | `profiles.hmm`, `ko_list.tsv` (KO thresholds) |
| `dbcan` | `screen` | dbCAN HMMdb V15; dbCAN-PUL substrate table 12-12-2023; dbCAN-PUL workbook snapshot of 2026-09-30 | `profiles.hmm`, `family_substrate.tsv`, `pul_rules.tsv` |
| `avg` | `avger` | CheckAMG v1.1.1 reference tables | V-scores, AMG/APG/AReG positive tables and filter tables |

The bundle's own manifest, `phu-db.json`, lists every upstream URL, SHA-256,
byte size, license and record count, and the exact preparation steps that turn
each upstream download into the installed file.

## Provenance

Each run writes the database provenance to `.phu/run.json` in its output
folder: the bundle version, the tarball SHA-256 and each database's upstream
release. Keep this file with your results.

## Citation

Cite the bundle version you used:

> García-Botero C. *Phu Databases*, version phu-db-0.1.1. Zenodo.
> <https://doi.org/10.5281/zenodo.19503347>

Also cite the upstream databases your analysis relied on, following each
project's own citation guidance. Pfam: Mistry J, et al. *Nucleic Acids Res.*
2021;49(D1):D412–D419. For KOfam, dbCAN and CheckAMG, see their project pages.

## Upgrading from PHU 0.9.x and earlier

Earlier releases downloaded each database separately and managed them with
`phu dbs`. That command group is gone; `phu dbs …` now only prints a notice
that databases install automatically.

After the first successful bundle install, PHU deletes the old per-database
folders in `PHU_DB_FOLDER` (`pfam/`, `kofam/`, `dbcan/`, `avg/` and `vscore/`)
and logs each path and its size. Move anything of your own out of those
folders before upgrading.
