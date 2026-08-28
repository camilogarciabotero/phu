# Databases and Provenance

PHU database downloads are scientific inputs. Record the backend, local
manifest, release identifier, and checksum with every analysis. The development
registry currently exposes `pfam`, `kofam`, `vscore`, and `avg`.

## Backend matrix

| Backend | Purpose | Upstream source in current build | Local files/indexes | Readiness | Provenance status |
| --- | --- | --- | --- | --- | --- |
| `pfam` | Resolve PFAM accessions and search protein-family HMMs | [Pfam current release](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz) | `Pfam-A.hmm.gz`, `Pfam-A.hmm`, `offsets.json`, optional `models/`, `manifest.json` | HMM and byte-offset index available | The source is a floating `current_release`; pinning, upstream release ID, license, size, and expected checksum are not yet recorded as a stable PHU contract. |
| `kofam` | Resolve KO models and apply KOfam thresholds | Zenodo record `19503464`: [`ko_list.gz`](https://zenodo.org/records/19503464/files/ko_list.gz), [`kofam.hmm.xz`](https://zenodo.org/records/19503464/files/kofam.hmm.xz) | `ko_list.gz`, `ko_list`, `kofam.hmm.xz`, `kofam.hmm`, `ko_metadata.json`, `offsets.json`, optional `models/`, `manifest.json` | HMM, KO metadata, and byte-offset index available | The Zenodo record is pinned in code, but per-file checksums, license, size, and citation metadata are not yet exposed as a stable PHU contract. |
| `vscore` | Provide V-score and normalized VL-score annotations | AnantharamanLab V-Score-Search source configured by the backend | normalized V-score TSV/CSV and `manifest.json` | Validated table and manifest available | Record release URL and checksum from the local manifest. Immutable release, license, size, and update policy require backend metadata work. |
| `avg` | Provide AVG/APG/AReG reference tables and scores | CheckAMG-derived reference source configured by the backend | normalized score, positive, and filter tables plus `manifest.json` | Reference tables and manifest available | The current code reports release `v1.1.1`; upstream URL, checksum policy, license, size, and immutable source identity require a stable provenance contract. |

## Lifecycle

Inspect status before using a backend:

```bash
phu dbs list
phu dbs status pfam kofam vscore avg
```

Prepare an explicitly named backend. Omitting names currently selects all
supported backends, which may download substantial data:

```bash
phu dbs prepare pfam
```

Use `refresh` to validate or repair local state and `remove --yes` to delete a
backend. An `indexed` state is not applicable to every backend; for table-based
backends, readiness means that the validated normalized table and manifest are
available.

## Reporting requirements

A reproducible run should record the backend key, local root, manifest path,
release or snapshot identifier, upstream URL, SHA-256 values, database files,
and preparation timestamp. Do not describe `pfam`'s floating `current_release`
source as immutable. Do not treat a local readiness result as evidence that the
upstream scientific release is permanently identified.

The current development build does not yet guarantee a complete provenance
record for every backend.
