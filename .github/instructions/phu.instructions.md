---
description: "Use when working on phu command implementations, shared prediction/cache logic, CLI interfaces, tests, docs, or dependency metadata. Preserve phu architecture, Typer UX conventions, and cache behavior."
applyTo: "src/phu/**/*.py,tests/**/*.py,docs/**/*.md,README.md,pyproject.toml,.github/instructions/*.md"
---

# phu Development Instructions

## Project Intent

- `phu` is a modular CLI toolkit for viral genomics workflows.
- Keep command behavior explicit, reproducible, and ergonomic for batch bioinformatics usage.
- Prefer extending existing modules over introducing parallel implementations.
- Treat command output formats, exit codes, environment variables, and database
	lifecycle behavior as public interfaces.
- Keep scientific computation separate from filesystem, CLI, and subprocess
	code whenever practical so threshold and boundary behavior can be tested
	without external databases or executables.

## General Development Workflow

- Start from the owning module, nearby tests, and the smallest reproducible
	behavior. Avoid broad refactors while implementing a focused feature.
- Before editing, identify the behavior that should change and one focused test
	that could disconfirm the proposed change.
- After each substantive edit, run the narrowest relevant test or type/lint
	check before expanding the change. Finish with the full test suite when the
	change crosses module boundaries.
- Inspect `git diff` and `git diff --check` before handing off work. Do not
	commit, reset, or discard unrelated user changes.
- Keep generated outputs, local databases, caches, credentials, and temporary
	artifacts out of version control unless they are deliberate small fixtures.

## Data, Reproducibility, and Provenance

- Never hardcode machine-specific paths, credentials, or writable locations.
	Use existing `PHU_*` and XDG configuration conventions and document new
	environment variables.
- Treat downloaded reference data as versioned scientific inputs. Record the
	upstream version, retrieval URL, checksum, normalization rules, and relevant
	threshold provenance in manifests or documentation.
- Use atomic replacement for downloaded/indexed files and validate checksums or
	complete manifests before making data available to a command.
- Make output ordering deterministic. Do not rely on filesystem, hash-map, or
	parallel completion order when writing user-facing tables or manifests.
- Preserve raw identifiers and distinguish database origin from accession when
	joining records; never silently merge equal-looking identifiers from different
	databases.

## API and Compatibility Rules

- Preserve existing public function signatures, CLI flags, output columns, and
	environment-variable semantics unless the task explicitly requires a breaking
	change.
- When a breaking change is necessary, update tests, command documentation,
	README usage examples, and changelog or migration notes together.
- Prefer additive output columns and explicit status values over silently
	changing the meaning of existing fields.
- Keep failure modes actionable: identify the missing input, database, command,
	configuration, or permission and return a non-zero exit status.
- Do not catch broad exceptions around scientific or filesystem operations
	unless the error is converted into a useful domain-specific message and the
	original cause remains available for debugging.

## Repository Architecture

- Preserve separation of concerns:
- CLI entry points and argument parsing in `src/phu/cli.py`.
- Per-command orchestration and config dataclasses in command modules (`screen.py`, `jack.py`, `cluster.py`, `simplify_vcontact_taxa.py`).
- Database lifecycle commands live under the `dbs` Typer command group in `src/phu/cli.py`; keep backend routing there and avoid scattering DB setup logic into feature commands.
- Database backend operations should be implemented as explicit helpers in `src/phu/pfam_db.py` or backend-specific modules with a common contract: `list`, `status`, `prepare`, `refresh`, and `remove`.
- External command discovery/execution helpers in `_exec.py`.
- Shared protein-prediction cache lifecycle in `gene_prediction_core.py`.
- Use `from __future__ import annotations` in Python source files.
- Keep imports, parsing, and configuration side effects out of module import
	time when a command can defer them until execution.
- Use pathlib and structured parsers for paths, TSV/CSV, JSON, and YAML rather
	than ad hoc string manipulation.

## CLI Conventions (Typer)

- Keep `Typer` command style consistent with existing commands.
- Provide long options and short aliases for common inputs (`-i`, `-o`, `-t`, `-m`, `-c`, `-k`) when appropriate.
- Use descriptive help text with practical defaults.
- Validate user-facing parameters early and fail with clear messages.
- Use `typer.secho(..., err=True, fg=typer.colors.RED)` for command errors and `typer.Exit(1)` for failure exits.
- Preserve root-level eager options (`--version`, `--clean-cache`) behavior.
- Group commands by concern when the surface grows: keep workflow commands together and place database lifecycle commands under `phu dbs` with its own help panel.

## Configuration and Validation Patterns

- Use dataclasses for command configuration objects.
- Enforce validation in `__post_init__` with specific `ValueError` messages.
- Keep config objects as the handoff boundary between CLI parsing and execution functions.
- Preserve deterministic behavior for enum-like options (for example, combine modes and command modes).

## Prediction Cache and Reuse Semantics

- `screen` and `jack` both rely on shared prediction caching in `gene_prediction_core.py`; keep them aligned.
- Cache keys must remain deterministic and depend only on prediction inputs (contigs identity + prediction parameters), not search-only parameters.
- Preserve support for `PHU_CACHE_DIR`, `XDG_CACHE_HOME`, and `PHU_CACHE` behavior.
- Keep atomic write patterns (`.partial`, temp manifest + replace) and lock handling for crash/process safety.
- If cache behavior changes, update:
- CLI/help text.
- `docs/cache.md`.
- Tests that assert cache behavior.

## External Tool Wrappers

- Keep wrappers explicit about required external executables and failure modes.
- Continue raising domain-specific errors such as `CmdNotFound` from command execution layers.
- Error messages must state required executables (for example `vclust`, `seqkit`) and expected availability in `PATH`.

## Testing Rules

- Add or update tests for every user-visible behavior change.
- Prefer `pytest` patterns already used in the repository:
- `CliRunner` for command-level behavior.
- `tmp_path` for filesystem workflows.
- `monkeypatch` for environment-dependent behavior.
- `pytest.raises(..., match=...)` for validation branches.
- Keep assertions concrete: exit codes, key help flags, output files, and error messages.
- Test scientific predicates at exact boundaries, including missing, conflicting,
	and below-threshold evidence where those states are part of the contract.
- Use small checked-in fixtures for parser and reference-data tests; do not make
	tests depend on network access, cluster schedulers, or locally installed large
	databases.
- Test cache hits and misses separately from annotation or scoring results.
- For CLI workflows, test both a successful run and the most useful failure
	path for each new validation rule.

## Documentation and Sync Requirements

- Keep docs synchronized with implementation and validation rules.
- When adding/changing flags or command behavior, update:
- `README.md`.
- Matching page in `docs/commands/`.
- `docs/cache.md` for prediction/cache impacts.
- When changing database lifecycle behavior, update the `dbs` command docs and any backend-specific notes, especially for PFAM preparation/indexing semantics.
- Avoid documenting behavior that code does not enforce.

## Dependency and Packaging Guidelines

- Keep runtime dependencies under `[project.dependencies]` in `pyproject.toml`.
- Keep tooling dependencies under `[dependency-groups]` (for example `dev`, `lint`).
- Preserve Python compatibility (`requires-python = ">=3.10"`) unless a deliberate project-wide bump is requested.
- Keep CLI entry point under `[project.scripts]` as `phu = "phu.cli:main"` unless the command surface is intentionally redesigned.
- Prefer the standard library and existing project dependencies. Add a runtime
	dependency only when it removes substantial complexity and its packaging,
	license, and offline-test implications are understood.
- Keep optional or development-only tools out of the runtime import path.

## Change Scope Discipline

- Prefer minimal, focused patches that do not alter unrelated command behavior.
- Keep naming and module organization consistent with the existing codebase.
- If a change introduces new command semantics, include migration notes in docs and tests in the same change set.

## Documentation Validation run

```text
.venv/bin/pytest -q
.venv/bin/mkdocs build --strict --config-file mkdocs.yaml --site-dir /tmp/phu-site
git diff --check
```

# Release Checklist

- [ ] Confirm the release scope and scientific limitations.
- [ ] Run tests, formatting checks, and a strict MkDocs build.
- [ ] Build and install the wheel and sdist in clean environments.
- [ ] Verify package, CLI, docs, changelog, and `CITATION.cff` versions agree.
- [ ] Capture external-tool versions and pinned database releases/checksums.
- [ ] Verify documented output headers and runnable examples.
- [ ] Inspect the sdist for README assets, license, changelog, and citation
      metadata.
- [ ] Update release notes with CLI, default, schema, cache, database, and
      scientific-method changes.
- [ ] Tag the immutable release and record its release DOI.
- [ ] Build stable documentation from the release tag and development
      documentation from the development branch.
- [ ] Confirm the published package, documentation, citation metadata, and
      repository release point to the same version.
