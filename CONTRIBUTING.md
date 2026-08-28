# Contributing

PHU contributions should preserve scientific correctness, reproducibility, and
stable command interfaces.

## Before opening a pull request

- Read the repository instructions and identify the owning command or module.
- Add or update focused tests for user-visible behavior.
- Update the README and matching command documentation when CLI behavior
  changes.
- Run `uv run pytest -q`, `uv run ruff format --check .`, and
  `uv run mkdocs build --strict`.
- Do not commit databases, caches, credentials, generated sites, or machine-
  specific paths.

Document limitations explicitly when an implementation is experimental or when
scientific validation is incomplete. Keep changes focused and explain any
change to outputs, thresholds, defaults, cache behavior, or dependencies.
