# Cache Handling

`phu` caches predicted proteins for both `screen` and `jack` so repeated runs can skip the gene-prediction step when the contigs and prediction inputs are unchanged.

## What the cache stores

The cache stores the translated protein FASTA generated from the input contigs. Search-specific inputs such as HMM files, seed markers, combine mode, and output folder are not part of the cache key.

## When the cache is reused

The cached proteins are reused when all prediction inputs are the same:

- the input contigs file has not changed
- the prediction mode is unchanged
- the protein-length filter is unchanged
- the gene-length filter is unchanged
- the translation table is unchanged

## When the cache is rebuilt

The cache is invalidated and rebuilt when any of these prediction inputs change:

- `--mode`
- `--min-gene-len` in `phu jack`
- `--min-protein-len-aa`
- `--ttable`

For `phu screen`, changing `--min-protein-len-aa`, `--mode`, or `--ttable` forces a rebuild. For `phu jack`, changing either `--min-gene-len` or `--min-protein-len-aa` forces a rebuild because both values feed the shared protein-prediction cache key.

## How to control it

Cache behavior is controlled by environment variables:

- `PHU_CACHE=off` disables caching and always recomputes proteins.
- `PHU_CACHE_DIR` overrides the cache location.
- If `PHU_CACHE_DIR` is not set, `phu` uses an XDG-style cache directory under your home cache path.

You can also remove all cached predictions explicitly:

```bash
phu --clean-cache
```

This removes the full cache directory resolved by `PHU_CACHE_DIR` (or the default cache path) and exits.

## Practical cases

- Changing only HMM files in `screen` reuses the cached proteins.
- Changing only seed markers in `jack` reuses the cached proteins.
- Changing the protein-length filter rebuilds the cache, which is expected because it changes the predicted protein set.

## Related docs

- [screen command](commands/screen.md)
- [jack command](commands/jack.md)