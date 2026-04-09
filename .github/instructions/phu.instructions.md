---
description: "Use when working on phu command implementations, shared prediction/cache logic, CLI interfaces, tests, docs, or dependency metadata. Preserve phu architecture, Typer UX conventions, and cache behavior."
applyTo: "src/phu/**/*.py,tests/**/*.py,docs/**/*.md,README.md,pyproject.toml"
---

# phu Development Instructions

## Project Intent

- `phu` is a modular CLI toolkit for viral genomics workflows.
- Keep command behavior explicit, reproducible, and ergonomic for batch bioinformatics usage.
- Prefer extending existing modules over introducing parallel implementations.

## Repository Architecture

- Preserve separation of concerns:
- CLI entry points and argument parsing in `src/phu/cli.py`.
- Per-command orchestration and config dataclasses in command modules (`screen.py`, `jack.py`, `cluster.py`, `simplify_vcontact_taxa.py`).
- External command discovery/execution helpers in `_exec.py`.
- Shared protein-prediction cache lifecycle in `gene_prediction_core.py`.
- Use `from __future__ import annotations` in Python source files.

## CLI Conventions (Typer)

- Keep `Typer` command style consistent with existing commands.
- Provide long options and short aliases for common inputs (`-i`, `-o`, `-t`, `-m`, `-c`, `-k`) when appropriate.
- Use descriptive help text with practical defaults.
- Validate user-facing parameters early and fail with clear messages.
- Use `typer.secho(..., err=True, fg=typer.colors.RED)` for command errors and `typer.Exit(1)` for failure exits.
- Preserve root-level eager options (`--version`, `--clean-cache`) behavior.

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

## Documentation and Sync Requirements

- Keep docs synchronized with implementation and validation rules.
- When adding/changing flags or command behavior, update:
- `README.md`.
- Matching page in `docs/commands/`.
- `docs/cache.md` for prediction/cache impacts.
- Avoid documenting behavior that code does not enforce.

## Dependency and Packaging Guidelines

- Keep runtime dependencies under `[project.dependencies]` in `pyproject.toml`.
- Keep tooling dependencies under `[dependency-groups]` (for example `dev`, `lint`).
- Preserve Python compatibility (`requires-python = ">=3.10"`) unless a deliberate project-wide bump is requested.
- Keep CLI entry point under `[project.scripts]` as `phu = "phu.cli:main"` unless the command surface is intentionally redesigned.

## Change Scope Discipline

- Prefer minimal, focused patches that do not alter unrelated command behavior.
- Keep naming and module organization consistent with the existing codebase.
- If a change introduces new command semantics, include migration notes in docs and tests in the same change set.
