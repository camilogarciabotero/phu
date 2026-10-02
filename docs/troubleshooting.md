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

Databases install automatically on first use; see
[Databases and Provenance](databases.md). Common causes of a failed install:

- **Not enough disk space.** The install needs about 12 GB free in
  `PHU_DB_FOLDER` and reports the exact amount when it stops.
- **No internet access** on the node. Download the bundle elsewhere and set
  `PHU_DB_BUNDLE` to its path.
- **Checksum mismatch.** The download or the `PHU_DB_BUNDLE` file is corrupt or
  is a different bundle version. Download the version named in the error.

An interrupted install leaves nothing half-installed; run the command again.
If you suspect an installed file was changed, run once with `PHU_DB_VERIFY=1`.

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
