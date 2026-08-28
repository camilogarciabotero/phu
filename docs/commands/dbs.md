# dbs

## What does it do?

The `phu dbs` command group manages local databases used by `phu`. It provides a scalable contract so each database can define its own preparation logic while sharing a common user interface.

Current built-in backends are `pfam`, `kofam`, `vscore`, and `avg`.

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

For `kofam`, preparation includes:

1. Ensuring `kofam.hmm` and `ko_list` are present locally.
2. Building KO metadata index from `ko_list` (`threshold`, `score_type`, etc.).
3. Building the byte-offset index for sparse KO model extraction.

For `vscore`, preparation downloads and validates the normalized V-Score CSV
from the AnantharamanLab repository and stores it locally with a manifest and
checksum. The parser requires the columns `Accession`, `Protein Function`,
`V-Score`, `Normalized VL-score`, `Log10[Hit Number]`, and `Database Origin`.
The normalized VL-score column is required by the Phase 4 AVL calculation.

KOfam metadata drives threshold behavior in `phu screen`; see [screen thresholds and decision logic](screen-thresholds.md).

KOfam support follows the KofamKOALA method for KO assignment and thresholding:

> Aramaki T., Blanc-Mathieu R., Endo H., Ohkubo K., Kanehisa M., Goto S., Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2019 Nov 19. pii: btz859. doi: 10.1093/bioinformatics/btz859.

## Examples

Prepare only PFAM:

```bash
phu dbs prepare pfam
```

Prepare all supported databases:

```bash
phu dbs prepare --all
```

Prepare only KOfam:

```bash
phu dbs prepare kofam
```

Check status:

```bash
phu dbs status pfam kofam
```

Prepare the V-score table:

```bash
phu dbs prepare vscore
```

Refresh integrity:

```bash
phu dbs refresh pfam kofam
```

Remove PFAM data:

```bash
phu dbs remove pfam --yes
```

Remove KOfam data:

```bash
phu dbs remove kofam --yes
```

List supported databases:

```bash
phu dbs list
```
