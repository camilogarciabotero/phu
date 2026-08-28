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

## Database is not ready

Inspect local database state before running a database-backed command:

```bash
phu dbs status
phu dbs list
```

Prepare only the named backend you need. Preparation may download large files:

```bash
phu dbs prepare pfam
```

An interrupted or incomplete preparation should be reported by `status`; run
`refresh` when supported, or remove and prepare the affected backend again.

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

Use command-specific `--help` to identify supported verbosity controls. Where
available, `--verbose` adds progress details. Preserve the complete command,
version, input identifiers, database status, and error output when reporting a
reproducible problem.
