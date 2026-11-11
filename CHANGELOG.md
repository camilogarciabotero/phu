# Changelog

All notable changes to PHU are recorded here. Development entries describe
behavior in the unreleased branch and may change before a release.

## Unreleased

### Databases

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

### Changed

- `screen`: input HMM files that share a file name (`a/foo.hmm`, `b/foo.hmm`)
  no longer fail with "Duplicate HMM model ID". They are told apart by the
  `NAME` inside each HMM, which becomes their model ID and names their domain
  tables and target files. If the internal `NAME` is identical too, `screen`
  warns and pools their hits under one model ID.
- `screen --save-target-hmms` refuses a query holding both a model `X` and a
  model named `X-aligned`, whose proteins would collide with the temporary
  alignment of `X`.

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
| screen | `target_proteins/<model>_proteins.aligned.faa` (temporary) | `target-proteins/<model>-aligned.faa` (temporary) |
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
