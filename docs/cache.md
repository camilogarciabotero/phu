# Cache Handling

Predicting proteins is often the slowest step of a run. `phu screen`, `phu jack` and `phu avger` keep the predicted proteins in a cache, so a second run on the same contigs, with different queries or thresholds, starts straight at the search.

## What the cache stores

The cache stores the translated protein FASTA generated from the input contigs. Search-specific inputs such as HMM files, seed markers, combine mode, and output folder are not part of the cache key.

## When the cache is reused

The cached proteins are reused when all prediction inputs are the same:

- the same contigs (same file content);
- the same `--mode`, `--min-protein-len-aa` and, in `phu jack`, `--min-gene-len`;
- the same `--ttable`, when you give one.

## When the cache is rebuilt

The cache is invalidated and rebuilt when any of these prediction inputs change:

- `--mode`
- `--min-gene-len` in `phu jack`
- `--min-protein-len-aa`
- an explicit `--ttable`

Without `--ttable`, pyrodigal-gv picks the genetic code of each contig in meta mode, and that choice is part of the cached result.

## How to control it

Cache behavior is controlled by environment variables:

- `PHU_CACHE=off` disables caching and always recomputes proteins.
- `PHU_CACHE_DIR` overrides the cache location.
- If `PHU_CACHE_DIR` is not set, `phu` uses an XDG-style cache directory under your home cache path.

You can also remove all cached predictions explicitly:

```bash
phu clean            # list what phu keeps and its size
phu clean --cache    # remove the prediction cache
```

`--cache` removes the full cache directory resolved by `PHU_CACHE_DIR` (or the default cache path). `phu clean --dbs` removes the installed database bundle and partial downloads, after asking, since they must be downloaded again; `--all` does both, and `--dry-run` only lists.

## Database index files

The cache directory can also hold a `phu-db-derived/` folder. PHU builds a
small index and per-model files for each reference database the first time it
is queried. They normally sit beside the database files; when the database
install is read-only for you, PHU keeps them here instead. `phu clean --cache`
removes them too, and they are rebuilt on the next database query. See
[Databases and Provenance](databases.md).

## Practical cases

- Changing only HMM files in `screen` reuses the cached proteins.
- Changing only seed markers in `jack` reuses the cached proteins.
- Changing the protein-length filter rebuilds the cache, which is expected because it changes the predicted protein set.

## Related docs

- [screen command](commands/screen.md)
- [jack command](commands/jack.md)