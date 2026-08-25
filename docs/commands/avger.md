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

## Scientific interpretation

`avger` separates three kinds of evidence:

1. **Annotation evidence:** a predicted protein has a Pfam or KOfam profile
    match that passes the database-specific score threshold and independent
    E-value cutoff (`1e-5` by default).
2. **Viral-family evidence:** the matched protein family has a V-score in the
    V-Score-Search table. Higher V-scores indicate stronger virus association
    for the family; they are not proof that every individual protein is viral.
3. **AVG interpretation:** the protein is an auxiliary viral gene candidate
    only when its annotation is evaluated in viral sequence context and a
    curated internal rule supports that interpretation. A metabolic or
    host-like function alone is not sufficient evidence of an AVG.

The V-Score-Search workflow recommends predicting ORFs with `pyrodigal-gv`,
searching Pfam and KEGG/KOfam profiles at maximum E-value `1e-5`, and using
genome-level average V-scores (AV-scores) and average Vᴸ-scores (AVᴸ-scores)
with fragment-size-specific cutoffs for **viral sequence identification**.
It also recommends filtering candidate viral sequences using CheckV
completeness. These are viral-classification criteria, not AVG-classification
criteria. `avger` does not silently apply them to individual proteins; the
input contigs should therefore already have a documented viral-sequence
assessment.

Thus, a row in `best_hits.tsv` is an **AVG candidate record**, not a confirmed
AVG call, when it has passing annotation evidence but no matching curated rule.
Such rows are labeled `unclassified_avg_candidate`. A named classification
requires an internal AVG rule that combines the relevant KOfam/Pfam evidence,
V-score support where appropriate, and the already-established viral context.
Passing a viral AV/AVᴸ cutoff does not by itself make a gene an AVG.

### Evidence levels

| Output situation | Scientific meaning |
| --- | --- |
| No output row for a protein | No Pfam/KOfam annotation passed the configured filters |
| Passing Pfam/KOfam row with V-score | Profile evidence plus family-level virus-association evidence |
| `unclassified_avg_candidate` | Candidate evidence exists, but no curated named class matched |
| Named classification with rule ID and version | A reviewed internal rule matched; inspect its provenance before biological use |

V-scores should not be treated as universal protein-level probability cutoffs.
Their interpretation depends on the database family and the V-Score-Search
methodology. See the [V-Score-Search documentation](https://github.com/AnantharamanLab/V-Score-Search)
and Zhou et al., *V-Score-Search* (preprint,
[doi:10.1101/2024.10.24.619987](https://doi.org/10.1101/2024.10.24.619987)).

### Threshold boundary

The published Table 1 AV-score and AVᴸ-score cutoffs belong to the upstream
viral-sequence classification step (Phase 4). They answer:

> Is this genome or fragment sufficiently virus-like?

AVG classification is a separate downstream step (Phase 5). It answers:

> Given a viral sequence, does this particular gene have evidence consistent
> with an auxiliary viral function?

Therefore, the Table 1 values must not be used as per-protein AVG thresholds.

## AVG candidate criteria

For each contig and database independently, `avger` uses only the best
significant hit for each protein and calculates:

$$
AVL_{db} = \frac{\sum VL\text{-scores of proteins with significant best hits in } db}
{\text{number of proteins with significant best hits in } db}
$$

The denominator is the number of scored significant hits in that database, not
the total number of predicted proteins. Pfam evidence is never combined with a
KOfam contig score, and vice versa.

A database-specific candidate requires all three strict inequalities:

- gene Vᴸ-score `< 3`;
- gene V-score `< 10`;
- same-database contig AVL-score `> 3`.

`pfam_candidate` and `kofam_candidate` are evaluated independently. The
combined `avg_candidate` is true when either database-specific candidate is
true. Boundary values are therefore excluded: Vᴸ-score `3`, V-score `10`, or
AVL-score `3` fails the corresponding condition.

## Flank evidence

For every database-specific candidate, `avger` reports whether another gene
with the same database V-score exactly `10` occurs upstream and downstream
within `10,000 bp`. It reports both nearest distances, edge-incomplete status,
and reason codes. Pfam candidates use Pfam flank evidence; KOfam candidates use
KOfam flank evidence.

Flank support is calculated by default but is not required because the input is
assumed to be a trusted vOTU. Use
`--require-flank-support` when both flank directions must be present for the
candidate to pass. The evidence states are:

- `not_assessable`: required score or coordinate evidence is unavailable;
- `not_candidate`: the strict AVL/V/Vᴸ criteria failed;
- `avg_candidate`: criteria passed, but flank support is absent or required and incomplete;
- `context_supported_avg_candidate`: criteria and both flank directions passed.

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

## Classifications

Classifications use the bundled, versioned rule set maintained with the package.
Rules are tested in declared order; the first matching rule supplies the
classification. Records that do not match remain
`unclassified_avg_candidate`, and the rule-set version is retained in the
output for provenance.