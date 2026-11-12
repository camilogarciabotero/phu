# Changelog

All notable changes to PHU are recorded here. Development entries describe
behavior in the unreleased branch and may change before a release.

## Unreleased

### Databases

- Databases now install automatically. The first `screen` or `avger` run that
  needs a database downloads one pinned, checksummed bundle, `phu-db-0.1.1`
  (about 1.9 GB, about 12 GB free disk needed while installing), verifies it
  and installs it under `PHU_DB_FOLDER`. Later runs reuse it.
- `phu dbs prepare` now downloads, verifies and installs that single bundle
  ahead of time. New options: `--bundle PATH` (pre-downloaded tarball) and
  `--verify` (re-hash every installed file). Database names, `--all` and
  `--force-refresh` are accepted and ignored.
- Removed `phu dbs list`, `status`, `refresh` and `remove`. They now print a
  notice pointing to `phu dbs prepare` and exit 0; they will be removed in a
  later release.
- Removed the per-database downloaders for Pfam, KOfam, dbCAN and the AVG/V-score
  tables. Pfam is now pinned to release 37.0 instead of the floating
  `current_release`.
- New environment variables: `PHU_DB_BUNDLE` (install from a pre-downloaded
  tarball, for offline and HPC nodes) and `PHU_DB_VERIFY=1` (re-hash every
  installed file on this run).
- `.phu/run.json` now records the bundle version, tarball SHA-256 and each
  database's upstream release. `screen` queries by Pfam, KOfam or dbCAN model
  record the same in `.phu/query.json`, also when no proteins are predicted.
- Installing a newer bundle keeps other `phu-db-<version>/` folders, so two phu
  releases can share one `PHU_DB_FOLDER`.
- Jobs that share an install build each database's index one at a time.
- A read-only shared install works: phu keeps the indexes it builds in the
  user's cache folder when it cannot write beside the database files.
- `PHU_DB_VERIFY=1` (and `phu dbs prepare --verify`) checks every database, not
  only the one a run asks for, and discards the indexes built from them.
- `screen` queries mixing dbCAN families with Pfam or KOfam models now list
  every model in `.phu/query.json`, not only the dbCAN ones.
- Migration: the old `pfam/`, `kofam/`, `dbcan/`, `avg/` and `vscore/` folders
  in `PHU_DB_FOLDER` are no longer read and can be deleted by hand. See
  [Databases and Provenance](docs/databases.md).

### Fixed

- `screen --hmm-mode mixed`: the domain table of an HMM file holding several
  models kept only the hits of its last model. It now holds the hits of every
  model in the file, under one header.
- Step timings were always reported as `0:00:00`. Every step now shows its
  real elapsed time, on a terminal and in redirected output.

### Added

- `normalize-lineage` compacts `unplaced_<rank>` chains. `unplaced_<rank>`
  becomes `U` plus the rank letter (`UR`, `UK`, `UP`, `UC`, `UO`, `UF`, `USF`,
  `UG`), and the anchor `default` becomes `DEF`:
  `unplaced_phylum_within_unplaced_kingdom_within_default` gives `DEF:UK:UP`.
  Chains mixing `novel_` and `unplaced_` nodes are supported, already compact
  values are left alone, and malformed `unplaced_` values (unknown rank, no
  anchor, wrong rank order) are reported as unparsed.
- Run logs for batch jobs: when the output is not a terminal, each step prints
  a start line with a timestamp and an end line with its status and duration,
  for example `[14:04:25] Predicting proteins done in 2m14s`.
- `--verbose` prints one `[cmd]` line per tool: the exact command line of
  external programs (vclust, seqkit, mafft) and the effective parameters of
  in-process steps (pyhmmer, pyrodigal). `screen` gains `--quiet` and
  `--verbose`, and `normalize-lineage` gains `--verbose`, so every command has
  both. See "Run logs" in [Troubleshooting](docs/troubleshooting.md).

### Changed

- `normalize-lineage`: the anchor `default` is written `DEF` in novel chains
  too, so `novel_genus_1_of_default` now gives `DEF:NG1` instead of
  `default:NG1`.
- The three-line Completed / Skipped / Failed summary is no longer printed
  after single steps; each step's own final line carries its status and time.
  `avger`, which reports several steps together, still prints it.
- `screen`: input HMM files that share a file name (`a/foo.hmm`, `b/foo.hmm`)
  no longer fail with "Duplicate HMM model ID". They are told apart by the
  `NAME` inside each HMM, which becomes their model ID and names their domain
  tables and target files. If the internal `NAME` is identical too, `screen`
  warns and pools their hits under one model ID.
- `screen` refuses a query in which two different models would get the same
  output file name once unusual characters are replaced (for example `beta/1`
  and `beta_1`), instead of mixing their domain tables and target proteins.
- `screen --save-target-hmms` builds HMMs only from the protein sets written by
  the current run, and keeps its temporary alignments outside the output
  folder. Files left in `target-proteins/` by an earlier run are not read.
- `normalize-lineage` reports a chain as unparsed when what remains after its
  valid nodes still starts with `novel_` or `unplaced_` (an unknown rank or a
  missing number), instead of treating that remainder as a taxon name.
- With `--verbose`, a `screen` search that is retried without gathering
  cutoffs logs a second `[cmd]` line with the parameters actually used.

### Renamed output files

Output files and folders now follow one kebab-case scheme: the name says what
the file contains, nucleotide FASTA is `.fna`, protein FASTA is `.faa`, ID lists
end in `-ids.txt`, and provenance lives under `.phu/`. These are breaking
changes for scripts that read phu outputs. Table columns and JSON keys are
unchanged.

| Command | Old name | New name |
| --- | --- | --- |
| screen, jack | `screened_contigs.fasta` | `matched-contigs.fna` |
| screen, jack | `kept_contigs.txt` | `matched-contig-ids.txt` |
| screen, jack | `proteins.faa` (`--keep-proteins`) | `predicted-proteins.faa` |
| screen, jack | `.phu_prediction_metadata.json` | `.phu/prediction.json` |
| screen | `hits_<hmm>.domtblout` | `<hmm>-hmmsearch-domains.tbl` |
| screen | `query_manifest.json` | `.phu/query.json` |
| screen | `resolved_pfam_hmms/` | `query-hmms/pfam/` |
| screen | `resolved_kofam_hmms/` | `query-hmms/kofam/` |
| screen | `resolved_dbcan_hmms/` | `query-hmms/dbcan/` |
| screen | `target_proteins/<model>_proteins.mfa` | `target-proteins/<model>.faa` |
| screen | `target_proteins/<model>_proteins.aligned.faa` (temporary) | no longer written to the output folder |
| screen | `target_hmms/` | `target-hmms/` |
| screen | `cazyme_matches.tsv.gz` | `cazyme-matches.tsv.gz` |
| screen | `cazyme_summary.tsv` | `cazyme-family-summary.tsv` |
| screen | `cazyme_class_summary.tsv` | `cazyme-class-summary.tsv` |
| screen | `pul_matches.tsv.gz` | `pul-matches.tsv.gz` |
| screen | `pul_family_support.tsv.gz` | `pul-family-support.tsv.gz` |
| screen | `pul_summary.tsv` | `pul-summary.tsv` |
| screen | `substrate_summary.tsv` | `pul-substrate-summary.tsv` |
| screen | `evidence/cazyme_hits.tsv.gz` | `evidence/cazyme-hits.tsv.gz` |
| jack | `jackhmmer_hits.tsv` | `jackhmmer-hits.tsv` |
| jack | `jackhmmer_iterations.tsv` | `jackhmmer-iterations.tsv` |
| jack | `last_iteration.hmm` | `final-iteration.hmm` |
| jack | `last_iteration_hmms/` | `final-iteration-hmms/` |
| cluster | `clustered-contigs/` (default output folder) | `phu-cluster/` |
| cluster | `fltr.txt` | `prefilter-pairs.txt` |
| cluster | `ani.tsv` | `ani-pairs.tsv` |
| cluster | `ani.ids.tsv` | `ani-sequence-ids.tsv` |
| cluster | `species.tsv` | `species-clusters.tsv` |
| cluster | `cluster_representatives_ids.txt` | `cluster-representative-ids.txt` |
| cluster | `representatives.fna` | `cluster-representatives.fna` |
| cluster | `dereplicated_representatives.fna` | `dereplicated-representatives.fna` |
| avger | `genes.tsv` | `gene-annotations.tsv` |
| avger | `evidence.tsv` | `avg-evidence.tsv` |
| avger | `avg_candidates.tsv` | `avg-candidates.tsv` |
| avger | `avg_predictions.tsv` | `avg-predictions.tsv` |
| avger | `relaxed_hits.tsv.gz` (`--keep-hits`) | `hmm-hits-relaxed.tsv.gz` |
| avger | `strict_hits.tsv.gz` (`--keep-hits`) | `hmm-hits-strict.tsv.gz` |

`clusters.tsv` (cluster, dereplication and votu modes) and `.phu/run.json` keep
their names. The predicted-protein cache layout is unchanged.

## 0.9.0

- Added the experimental `avger` workflow and AVG-related database backends.
- Added the `dbs` command group for database lifecycle operations.
- Added shared predicted-protein cache handling for `screen` and `jack`.
- Documentation now identifies the development revision explicitly.
- Release work remains for output schemas, provenance manifests, cache
  compatibility guarantees, and stable/development documentation deployments.

## 0.8.0

See the corresponding release notes and immutable release artifacts.
