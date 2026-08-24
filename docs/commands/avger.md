# avger

`phu avger` predicts proteins and annotates them against the complete Pfam and
KOfam databases. It writes one deterministic best-hit row per protein and
database to `best_hits.tsv`.

## Synopsis

```bash
phu avger -i contigs.fa -o phu-avger
```

The command reuses the shared prediction cache. KOfam hits use the thresholds
and score types from `ko_list`; Pfam hits use model GA cutoffs. V-score fields
are added for KOfam accessions when the local V-score table is prepared.

## Outputs

- `best_hits.tsv`: best passing Pfam/KOfam hit per protein, with score, E-value,
  coordinates, thresholds, and optional V-score fields.
- `all_hits.tsv.gz`: all passing hits when `--keep-all-hits` is supplied.

Prepare the V-score table before running the default V-score-enriched workflow:

```bash
phu dbs prepare vscore
phu avger -i contigs.fa -o phu-avger
```

Use `--no-use-vscore` to run without V-score lookup.