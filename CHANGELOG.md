# Changelog

All notable changes to PHU are recorded here. Development entries describe
behavior in the unreleased branch and may change before a release.

## 0.10.0 — 2026-11-18

### CLI

- `phu --help` groups commands into three panels: **Main** (`screen`,
  `jack`, `avger`, `cluster`), **Database** (`dbs`) and **Utils** (`phrogs`,
  `clean`, `normalize-lineage`).
- New `phu clean`: without options it lists what phu keeps (prediction cache,
  installed databases, partial downloads) and its size; `--cache`, `--dbs`
  and `--all` remove them. `--dbs` asks before deleting the databases unless
  `--yes`; `--dry-run` shows what would go. Only phu's own folders are
  removed. `phu --clean-cache` still works but is deprecated and hidden; use
  `phu clean --cache`.

### Screen

- In the default `--hmm-mode pure`, an HMM file may now hold several models
  and counts as one query: concatenate the profiles of one marker into one
  file, and `--combine-mode all` or `threshold` counts markers instead of
  profiles. Hits to any model in the file are reported under the file's
  name, and a note names each such file and its model count. Before, such
  a file was refused; `--hmm-mode mixed` still makes each model its own
  query.
- Fixed: `--combine-mode threshold` kept only the best hits of each contig
  overall, not of each matched query, so `--save-target-proteins` (and
  `--save-target-hmms`) missed the proteins of every query but the
  top-scoring one. It now keeps the best hits of each matched query, as
  `any` and `all` do and as the docs said. Which contigs are kept is
  unchanged.
- Fixed: Pfam's gathering cutoffs were dropped for the whole run when Pfam
  families were queried together with any model without them (your own
  HMMs, KOfam, PHROG, ECOD, dbCAN): the search failed on the missing cutoffs
  and was repeated without any, with only a warning. With `--cut-ga` (the
  default), models that have gathering cutoffs are now searched with them and
  the others without, in the same run.
- `phu screen --help` is shorter and grouped by question: the query, other
  ways to query, hit thresholds, which contigs are kept, output, gene
  prediction and runtime. Option names, short flags and defaults are
  unchanged; the details moved to the online guide.
- `phu screen` searches ECOD families: `ecod:2004.1.1.5` is one F-group, and
  a shorter id stands for every F-group inside a T, H or X group
  (`ecod:2004.1.1`, `ecod:2004.1`, `ecod:2004`). Models come from the new
  `ecod` database of phu-db 0.4.0 (ECOD v295.2, 32,041 families). An ECOD hit
  must align to at least half of the model (`--no-use-ecod-thresholds` drops
  this). With `--combine-mode all` or `threshold`, each `ecod:` id counts as
  one query, met by any of its families. What each id expanded to and the
  families' names are recorded in `.phu/query.json`, and protein files are
  labelled with family names.

### Cluster

- New `phu cluster --mode genus`: genus-level groups of complete genomes,
  complete linkage on total ANI (tANI) ≥ 0.70, following vclust's VIRIDIC-like
  recipe (prefilter `--min-ident 0.7`, every aligned pair kept, no
  query-coverage requirement). Writes `genus-clusters.tsv` and
  `cluster-representatives.fna`. The groups are operational, not ICTV taxa.
- Fixed: `species` mode filtered the alignment output at ANI ≥ 0.95, a
  measure its tANI clustering does not use, and every mode's prefilter
  identity was fixed in the code whatever the threshold. All steps now read
  one effective setting per mode; `genus` and `species` keep every aligned
  pair, as upstream does. Species groups can change in rare cases. The old
  documented way to get genus groups (`--mode species` with `--tani 0.70`)
  kept the ANI 0.95 filter and so did not give genus groups.
- `--vclust-params` flags replace the mode's value instead of being added a
  second time, and `--out-tani`, `--out-gani`, `--out-rcov` and `--rcov` are
  accepted. Settings that would drop pairs the clustering needs, or leave the
  metric without a threshold, stop the run with an explanation before vclust
  starts: an align filter above the clustering threshold, `--metric` without
  its own threshold, or a threshold below the prefilter identity without
  `--min-ident` (for example `--mode votu --vclust-params "--ani 0.90"`,
  which before silently kept the 0.95 prefilter).
- vclust 1.3.1 `align --filter` on two or more threads skips pairs between
  sequences with single-character names, which split genus and species
  groups. When the input has such a name, `phu cluster` now aligns on one
  thread and prints a note.
### Jack

- `phu jack --help` is reorganised like `phu screen --help`: seeds, PHROG seeds,
  search, which contigs are kept, output, gene prediction and runtime, with
  one-line help per option. `--combine-mode` now says that one `--phrog`
  query counts as one seed. Option names, short flags and defaults are
  unchanged.
- `phu jack --phrog QUERY` uses PHROG families as seeds, alone or with a seed
  FASTA (which is now optional). Queries are those of `phu phrogs resolve`
  (`phrog_<n>`, `alias:`, `category:`, `annot:`). The first iteration searches
  with the PHROG model itself; later ones with models rebuilt from the hits.
  A selection larger than `--max-phrog-seeds` (default 50) stops the run
  before proteins are predicted.
- `--combine-mode all` and `threshold` count seed groups: each FASTA seed,
  and each `--phrog` query however many PHROGs it selects. Within a group a
  protein is kept once, by its best hit. For FASTA-only runs nothing changes.
- `jackhmmer-hits.tsv` and `jackhmmer-iterations.tsv` have new
  `seed_source` (`fasta` or `phrog`) and `seed_group` columns after
  `seed_id`. Scripts that read these tables by column position need updating.
- With `--verbose`, jack logs what every `--phrog` query resolved to: the
  PHROG identifiers per query, then each seed's annotation, category and
  group.
- jack now writes `.phu/query.json`, recording its seeds (and, for PHROG
  seeds, their annotations and the database bundle's provenance).

### Phrogs

- New `phu phrogs` command to find PHROG identifiers by name, for
  `phu screen`: `categories`, `aliases`, `search <text>` and
  `resolve <queries>`. Queries are `phrog_<n>`, `alias:<name>` (curated groups
  of exact PHROG annotations such as `mcp`, `terl`, `portal`, `collar`,
  `tail-fiber`, `holin`), `category:<name>` (matched forgivingly: `tails`,
  `head`, `dna`, `moron`/`amg`) and `annot:<annotation>`.
  `--ids-only` prints identifiers for the shell, for example
  `phu screen -i contigs.fna $(phu phrogs resolve alias:mcp --ids-only)`.
  `search` forgives British spelling, plurals, joined words, typos and
  common abbreviations (`capsid`, `TerL`, `MTase`, `RNR`, `ssb`), and says how
  each row matched.

### Avger

- Strict annotations in every database now also need coverage: a hit counts
  only when its domains align over at least half of the model and half of the
  protein (every reported domain counts, overlaps once, as CheckAMG counts
  it). The new `--min-coverage` option sets the fraction; `0` restores the
  old behaviour for Pfam and KOfam. Expect fewer strict Pfam and KOfam
  annotations, most visibly for short domains inside long proteins. The
  relaxed scoring track is unaffected.
- The output of `phu avger` is reorganised around the person reading it.
  `avg-predictions.tsv` is now a short table, one row per predicted gene with
  its class, supporting reference entries and rule values. New files:
  `contig-summary.tsv` (one row per contig), `avg-proteins.faa` (sequences of
  the predicted genes) and `output-info.txt` (the run's counts and a data
  dictionary for the predictions table). The wide per-gene table
  `gene-annotations.tsv` and the audit table `avg-evidence.tsv` move into
  `detailed/`, joined by `filtered-candidates.tsv` (candidates that did not
  become predictions, with the reason) and, with `--keep-hits`, the raw hit
  files. `avg-candidates.tsv` is gone: its rows are the predictions plus the
  filtered candidates. The `zhou_*` columns are renamed `rule_*`. The output
  schema version in `.phu/run.json` is now 3.
- `phu avger` now reads binary copies of the HMM files (about 5 GB), which
  load in seconds instead of minutes. `phu dbs prepare` builds them and shows
  the step; otherwise the first `phu avger` run does, as a visible one-time
  step.
- Writing the avger tables no longer slows down with the size of the
  reference lists: positive and filter references are looked up through an
  index instead of being scanned for every gene.
- Fixed: on inputs with more than 10,000 proteins, `phu avger` searched the
  proteins in batches and each batch computed E-values as if it were the whole
  input, so hits near the E-value cutoff could pass or fail depending on the
  batch they fell in. E-values now always refer to the whole input. Results
  for inputs of up to 10,000 proteins are unchanged.
- `phu avger` searches each database once instead of twice, which roughly
  halves its run time. The relaxed and strict tracks used to run the same
  search separately and differ only in the thresholds applied to the hits;
  one search now feeds both. Results are unchanged.
- `phu avger` now searches PHROG as a third database, beside Pfam and KOfam.
  PHROG scores can make a gene a candidate, PHROG entries support the
  physiology and regulatory classes, and the PHROG filter lists can block a
  class when a gene's PHROG hit is a core viral function. A PHROG annotation
  needs a bit score of at least 60, with at least 50% of the model and 50% of
  the protein aligned. Use
  `--no-phrog` to skip it. Expect longer runs with PHROG on: it has more
  models than Pfam and KOfam.
- New columns in the avger tables: `phrog_id`, `phrog_annotation`,
  `phrog_category`, `phrog_bitscore`, `phrog_evalue`, `phrog_model_coverage`,
  `phrog_protein_coverage`,
  `phrog_scoring_id`, `phrog_v_score`, `phrog_vl_score`, `phrog_scaffold_avl`,
  `phrog_scored_gene_count` and `zhou_phrog_candidate`. The output schema
  version in `.phu/run.json` is now 2.

### Docs

- The home page starts from what you want to do, with a table of tasks and
  commands. The `screen` guide is reorganised around what you can search for
  (with ECOD's levels explained), how queries combine, and what you get; the
  thresholds page opens with one table of the rule for each kind of model.
- Updated for phu-db 0.4.0: ECOD everywhere databases are listed, which
  command uses which database, download and disk sizes, PHROG and ECOD
  citations. The `avger` page shows its current options, links no longer
  point to the old `avg` branch, and the cache page lists `avger` too.

### Databases

- phu now installs phu-db 0.4.0 (doi:10.5281/zenodo.24296611). It adds the
  `ecod` database (ECOD v295.2, F70 representatives, CC-BY-4.0):
  `profiles.hmm` (one model per F-group, 32,041), `msa_a3m.tar.gz` (the
  alignments the models were built from), `sequences.fa`, `families.tsv` and
  `domains.tsv`; every other file is identical to 0.3.0. The bundle manifest
  accepts the ECOD preparation steps and the `fasta` file format.
- phu-db 0.3.0 (doi:10.5281/zenodo.24290671) added
  `phrog/msa_a3m.tar.gz`, the 38,880 PHROG alignments as A3M files for
  HH-suite or MMseqs2; every other file is identical to 0.2.0. The bundle
  manifest accepts the `build-a3m` step and the `tar.gz` file format, checked
  by size and SHA-256.

### Internal

- `phu avger` and `phu screen` check the memory available, honouring a batch
  job's own limit, before loading proteins. `avger` searches all proteins in
  one batch when they fit, so each database's models are read once instead of
  once per 10,000 proteins (results are unchanged: E-values do not depend on
  the batch). `screen` reads the proteins from disk during the search when
  loading them would take more than a tenth of the memory, instead of risking
  running out of it.

- `phu jack` reads the predicted proteins once and reuses them for every
  seed, instead of re-reading the protein file per seed.

- The Python modules are grouped into packages: `phu.commands` (one module per
  command), `phu.avger`, `phu.db` and `phu.core`. The command line is
  unchanged. Code that imported phu modules directly must use the new paths,
  for example `phu.commands.screen` instead of `phu.screen`.

### Requirements

- seqkit is no longer needed. `phu screen`, `phu jack` and `phu cluster` now
  extract matched contigs, target proteins and cluster representatives
  themselves, reading the input once (target proteins were read once per
  model before). The same records are written in the same order; sequence
  lines keep the input's wrapping instead of being rewrapped at 60
  characters as seqkit did. `vclust`, for `phu cluster`, is now the only
  external program.
- `psutil` is now a direct dependency (it already came with pyhmmer).
- PHU now needs Python 3.11 or later. Python 3.10 reached end of life in
  October 2026 and is no longer supported or tested; create a new environment
  with a newer Python to upgrade.

### Screen

- With `--save-target-proteins`, proteins matched by a PHROG query now carry
  the PHROG identifier, function and category in their FASTA headers, as
  KOfam matches carry the KO definition.

- The `hmm_coverage` and `protein_coverage` of a hit now count every reported
  domain (overlaps once) instead of only the best domain, matching `phu
  avger` and CheckAMG. Hits matched in several pieces gain coverage, so the
  dbCAN rule, `--min-hmm-coverage` and `--min-protein-coverage` can keep hits
  they previously dropped. The `hmm_from`/`hmm_to` columns still describe the
  best domain.
- `phu screen` accepts PHROG identifiers as queries, for example
  `phu screen phrog_216 -i contigs.fna`, alone or mixed with Pfam, KOfam and
  dbCAN identifiers and your own HMMs. A PHROG hit needs a bit score of at
  least 60 and at least half of both the model and the protein aligned;
  `--no-use-phrog-thresholds` turns that off. The function and category of
  each queried PHROG are recorded in `.phu/query.json`.
- The default `--max-evalue` of `phu screen` is now `1e-6` (was `1e-5`), for
  every kind of model. Pass `--max-evalue 1e-5` for the previous behaviour.
- Hits to your own HMMs must now cover more than 35% of the model. The new
  `--min-hmm-coverage` option sets the fraction; `0` turns the filter off.
  Pfam, KOfam and dbCAN models are unaffected.
- New `phu screen --min-protein-coverage F` option, off by default: keep a
  hit only when it aligns to at least the fraction `F` of the protein. It
  applies to every kind of model. The CAZyme hit tables gain a
  `protein_coverage` column.

- `phu screen --save-target-hmms` no longer needs `mafft`. Target proteins are
  aligned inside PHU with FAMSA through the `pyfamsa` package, which is now a
  dependency installed with PHU. Models built from several proteins can differ
  slightly from those built from MAFFT alignments. With `--verbose` the step
  is logged as `[cmd] pyfamsa.align ...`.
- New `phu screen --trim-target-msa {none,gappyout,strict,strictplus}` (default
  `none`): trims each target alignment with trimAl, through the new `pytrimal`
  dependency, before its HMM is built. Needs `--save-target-hmms`. Trimmed
  models are shorter and give lower bit scores.
- `pyfamsa` is held to the 0.5 series: its wheels run on glibc 2.17 and later,
  while 0.6 and later need glibc 2.34 or a C++20 compiler to build.

### Databases

- Databases now install automatically. The first `screen` or `avger` run that
  needs a database downloads one pinned, checksummed bundle, `phu-db-0.2.0`
  (about 2.2 GB, about 16 GB free disk needed while installing), verifies it
  and installs it under `PHU_DB_FOLDER`. Later runs reuse it.
- The bundle includes PHROG: HMMER models for the 38,880 PHROGs, built from
  the official alignments, and the PHROG annotation table.
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
- The bundle is downloaded as several parts at once (8 by default; set
  `PHU_DB_CONNECTIONS` or `phu dbs prepare --connections N`), which is several
  times faster from Zenodo. An interrupted download continues every part where
  it stopped; a server without range support falls back to one connection.
- The bundle download shows its progress: on a terminal, a spinner with a
  bar, size, speed and time left; in a log file, one line per 10%.
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

### Deprecated

- `simplify-taxa` is deprecated. It no longer appears in `phu --help` and
  prints a one-line deprecation notice to stderr when run; its behavior and
  output are otherwise unchanged, and `phu simplify-taxa --help` still works.

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
