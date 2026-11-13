# Troubleshooting

## Missing external executable

`cluster` requires the external `vclust` executable. Some workflows may also
use `seqkit`. Confirm that required programs are available on `PATH`:

```bash
command -v vclust
command -v seqkit
```

Install the missing program using the environment or package manager appropriate
for your system, then rerun the command. PHU does not install every external
binary automatically.

## Database install fails

Databases install automatically on first use, or with `phu dbs prepare`; see
[Databases and Provenance](databases.md). Common causes of a failed install:

- **Not enough disk space.** The install needs about 16 GB free in
  `PHU_DB_FOLDER` and reports the exact amount when it stops.
- **No internet access** on the node. Download the bundle elsewhere and set
  `PHU_DB_BUNDLE` to its path.
- **Checksum mismatch.** The download or the `PHU_DB_BUNDLE` file is corrupt or
  is a different bundle version. Download the version named in the error.

- **Certificate errors** (`CERTIFICATE_VERIFY_FAILED`). Some Python installs
  have no certificate bundle. Point Python at the system one, for example
  `export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt`.
- **A slow or stalled download.** The download shows its progress and resumes
  where it stopped when you run the command again. It fetches 8 parts of the
  bundle at once; if your network or proxy objects to that, set
  `PHU_DB_CONNECTIONS=1`.

An interrupted install leaves nothing half-installed; run the command again.
If you suspect an installed file was changed, run `phu dbs prepare --verify`
(or any command once with `PHU_DB_VERIFY=1`): it re-hashes every database
file and reinstalls if a check fails.

A database folder you cannot write to, such as an install shared by a
colleague, works as is: PHU keeps the index files it builds in your own cache
folder. See [Databases and Provenance](databases.md).

## Cache or permissions errors

Set a writable cache location explicitly when the default home cache is not
available:

```bash
PHU_CACHE_DIR="$PWD/.phu-cache" phu --version
```

Use `PHU_CACHE=off` for a diagnostic run without cache reuse. `phu --clean-cache`
removes the resolved cache directory, so verify `PHU_CACHE_DIR` before using it.

## Input and empty results

Check FASTA headers and sequences when a run fails during parsing. A valid FASTA
must contain headers beginning with `>` and nucleotide or protein sequence text
appropriate to the command. No hits is a valid scientific result and is
reported differently from an execution error; inspect the command's stderr and
exit status before treating an empty result as a failure.

## Interrupted runs and stale outputs

Use a new output directory for a clean diagnostic run. Review existing files
before rerunning into a directory that contains results. PHU's commands do not
provide one uniform overwrite contract across all workflows; do not assume that
stale files are removed automatically.

## Debug diagnostics

The workflow commands (`screen`, `jack`, `cluster`, `avger`,
`normalize-lineage`, `simplify-taxa`) accept `--verbose` and `--quiet`.
Preserve the complete command, version, input identifiers, database status,
and error output when reporting a reproducible problem.

## Run logs

When the output is not a terminal (for example a SLURM `.out` file), each step
prints one line when it starts and one when it ends, with the wall-clock time
and the step's real duration, and without colour codes:

```text
[14:02:11] Predicting proteins ...
[14:04:25] Predicting proteins done in 2m14s
[14:04:25] Running HMM searches ...
[14:09:02] Running HMM searches failed after 4m37s: <error message>
```

A step ends with `done in`, `failed after`, `skipped after`, or
`cancelled after`. On a terminal the same steps are shown as live lines with a
running timer.

`--verbose` adds one line per tool invocation, always prefixed with `[cmd]`.
External programs (`vclust`, `seqkit`) are shown as the exact,
shell-quoted command that was executed, including any output redirection.
In-process steps (pyrodigal-gv gene prediction, pyhmmer `hmmsearch` and
`jackhmmer`, pyfamsa alignment, pytrimal trimming) are shown as the function name followed by the `key=value`
parameters that determine the result:

```text
[cmd] pyrodigal_gv.ViralGeneFinder contigs=4 mode=meta meta=True min_gene=90 max_overlap=89 translation_table=auto training_table=none min_protein_len_aa=30 threads=2
[cmd] pyhmmer.hmmsearch hmm_files=1 models=1 sequences=19 cpus=2 parallel=queries bit_cutoffs=gathering hmm_mode=pure ...
[cmd] /usr/bin/seqkit grep -f out/matched-contigs.fna.ids.txt contigs.fna > out/matched-contigs.fna
```

Extract them with `grep '^\[cmd\]' run.log`. Gene prediction prints no
`[cmd]` line when the predicted proteins come from the cache, because nothing
is run. `--quiet` hides the step lines and the `[cmd]` lines; errors are still
printed. Nothing is written to a log file: redirect stdout and stderr to keep
a record.
