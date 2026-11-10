# Changelog

All notable changes to PHU are recorded here. Development entries describe
behavior in the unreleased branch and may change before a release.

## Unreleased

- Databases now install automatically. The first `screen` or `avger` run that
  needs a database downloads one pinned, checksummed bundle, `phu-db-0.1.1`
  (about 1.9 GB, about 12 GB free disk needed while installing), verifies it
  and installs it under `PHU_DB_FOLDER`. Later runs reuse it.
- Removed the `phu dbs` command group (`list`, `status`, `prepare`, `refresh`,
  `remove`). `phu dbs …` now prints a deprecation notice and exits 0; it will
  be removed in a later release.
- Removed the per-database downloaders for Pfam, KOfam, dbCAN and the AVG/V-score
  tables. Pfam is now pinned to release 37.0 instead of the floating
  `current_release`.
- New environment variables: `PHU_DB_BUNDLE` (install from a pre-downloaded
  tarball, for offline and HPC nodes) and `PHU_DB_VERIFY=1` (re-hash every
  installed file on this run).
- `.phu/run.json` now records the bundle version, tarball SHA-256 and each
  database's upstream release.
- Migration: the old `pfam/`, `kofam/`, `dbcan/`, `avg/` and `vscore/` folders
  in `PHU_DB_FOLDER` are no longer read and can be deleted by hand. See
  [Databases and Provenance](docs/databases.md).

## 0.9.0

- Added the experimental `avger` workflow and AVG-related database backends.
- Added the `dbs` command group for database lifecycle operations.
- Added shared predicted-protein cache handling for `screen` and `jack`.
- Documentation now identifies the development revision explicitly.
- Release work remains for output schemas, provenance manifests, cache
  compatibility guarantees, and stable/development documentation deployments.

## 0.8.0

See the corresponding release notes and immutable release artifacts.
