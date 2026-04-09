# dbs

## What does it do?

The `phu dbs` command group manages local databases used by `phu`. It provides a scalable contract so each database can define its own preparation logic while sharing a common user interface.

Right now, `pfam` is wired first.

## Synopsis

```bash
phu dbs [COMMAND] [OPTIONS] [DATABASES...]
```

## Supported commands

- `phu dbs list`
- `phu dbs status [DATABASES...] [--all]`
- `phu dbs prepare [DATABASES...] [--all] [--force-refresh]`
- `phu dbs refresh [DATABASES...] [--all]`
- `phu dbs remove [DATABASES...] [--all] --yes`

## Behavior contract

For each database, operations are interpreted by that backend:

- **prepare**: make the database ready for runtime use.
- **refresh**: validate local state and repair incomplete or stale data.
- **remove**: delete local data for selected databases.
- **status**: report readiness and local metadata.

For `pfam`, preparation includes:

1. Ensuring `Pfam-A.hmm` is present locally.
2. Building the byte-offset index used for accession lookup.

## Examples

Prepare only PFAM:

```bash
phu dbs prepare pfam
```

Prepare all supported databases:

```bash
phu dbs prepare --all
```

Check status:

```bash
phu dbs status pfam
```

Refresh integrity:

```bash
phu dbs refresh pfam
```

Remove PFAM data:

```bash
phu dbs remove pfam --yes
```

List supported databases:

```bash
phu dbs list
```
