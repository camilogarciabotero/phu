# Screen thresholds

This page is the reference for how `phu screen` decides which hits count and
which contigs are kept. For a task-oriented overview, see the
[screen guide](commands/screen.md#what-you-can-search-for).

These are acceptance rules for profile hits. They do not prove that a contig
is viral or that a gene is an auxiliary viral gene; for that, see the
[avger scientific interpretation](commands/avger.md#scientific-interpretation).

## The decision in brief

For every hit pyHMMER reports as included:

1. Apply the rule of the model's kind (table below).
2. Apply `--min-protein-coverage`, if set.
3. Compare the hit's score with the effective minimum bit score.
4. Keep the hit only if its E-value is at most `--max-evalue` (default `1e-6`).

Contigs with at least one remaining hit are then kept or dropped by
`--combine-mode` (see [Combining hits into contigs](#combining-hits-into-contigs)).

| Model kind | Rule | Turn it off with |
| --- | --- | --- |
| Your own HMM | Aligns to more than `0.35` of the model | `--min-hmm-coverage 0` |
| Pfam | Pfam's gathering (GA) cutoff, applied during the search | `--no-cut-ga` |
| KOfam | The KO's threshold from `ko_list`, on the KO's score type | `--no-use-kofam-thresholds` |
| PHROG | At least 60 bits, and at least `0.5` of model and of protein aligned | `--no-use-phrog-thresholds` |
| ECOD | At least `0.5` of the model aligned | `--no-use-ecod-thresholds` |
| dbCAN | Domain (independent) E-value below `1e-15` and more than `0.35` of the model aligned | always on |

A hit without a usable alignment fails any rule that needs coverage.

## Coverage

Coverage counts every reported domain of a hit, with overlapping domains
counted once, on the model side (`hmm_coverage`) and on the protein side
(`protein_coverage`). This is how `phu avger` and CheckAMG compute it. Both
values are written to the result tables.

`--min-protein-coverage F` is off by default. When set, a hit passes only if
its domains cover at least the fraction `F` of the protein, for every kind of
model and in addition to that kind's rule.

## Your own HMMs

A hit to one of your HMMs must align to more than `--min-hmm-coverage` of the
model (default `0.35`; `0` turns the filter off). This removes matches to a
small piece of the family. The filter does not apply to database models,
which have their own rules. An HMM file of yours is treated as your own even
if its name looks like a database id.

## Pfam

Pfam accessions are resolved to the bundle's Pfam models. With `--cut-ga`
(the default), pyHMMER applies each model's curated gathering cutoff during
the search. With `--no-cut-ga`, only `--min-bitscore` and `--max-evalue`
apply.

!!! warning "Pfam cutoffs are lost when Pfam is mixed with other models"
    Gathering cutoffs exist only in Pfam models. If any model in the run has
    none (your own HMMs, KOfam, PHROG and ECOD models), the whole search is
    repeated without them, Pfam models included, and a warning is printed.
    To keep Pfam's curated cutoffs, query Pfam families in a run of their
    own, or set `--min-bitscore` for the mixed run.

## KOfam

KOfam models come with a threshold and a score type (`full` or `domain`) for
each KO, from KOfam's `ko_list`. With `--use-kofam-thresholds` (the default):

- the score compared is the full-sequence score, or the best domain's score
  when the KO's score type is `domain`;
- the minimum is the KO's threshold, or `--min-bitscore` if that is higher.

So `--min-bitscore` can only make KOfam stricter. The E-value test always uses
the full-sequence E-value, even for KOs scored by domain.

For example, a KO with threshold `136.43` and score type `full` rejects hits
scoring 15 to 20 bits, however good their domain scores look in the domain
table.

## PHROG

PHROG models have no curated cutoffs. A hit must score at least `60` bits and
align to at least `0.5` of the model and `0.5` of the protein, the rule
`phu avger` uses.

## ECOD

ECOD models have no curated cutoffs either. They describe single domains,
often within longer proteins, so a hit must align to at least `0.5` of the
model but may cover any part of the protein. Only the E-value applies with
`--no-use-ecod-thresholds`.

## dbCAN

A dbCAN hit must have a best-domain independent E-value below `1e-15` and
align to more than `0.35` of the model, the cutoffs dbCAN recommends for its
HMMs. These rules are always on.

With `--all-cazymes`, every CAZyme family is searched once and any family hit
keeps its contig. With `--all-puls`, a PUL matches a contig when all its
families hit there, and a contig is kept when at least one PUL matches. PULs
that name families missing from the installed dbCAN models are skipped and
listed in `.phu/run.json`.

## Bit score and E-value

Each hit has a full-sequence bit score and the score of its best domain. The
full-sequence score is used, except for KOs scored by domain (see
[KOfam](#kofam)).

`--min-bitscore`, if set, applies to every model. `--max-evalue` applies to
the hit's full-sequence E-value, never to the per-domain E-values in the
domain table.

## Combining hits into contigs

After filtering, contigs are kept by `--combine-mode`:

- `any`: at least one query hits the contig.
- `all`: every query hits the contig.
- `threshold`: at least `--min-hmm-hits` different queries hit the contig.

Each `ecod:` id is one query, met by a hit to any of its families. For each
kept contig, the best hits of each matched query are reported, up to
`--top-per-contig` (default 1).

## Domain tables

`--keep-domtbl` (the default) keeps hmmsearch's per-domain table of each model
(`<model>-hmmsearch-domains.tbl`). The decisions above are made from pyHMMER's
results directly, not by reading these files back; keep them to audit why a
hit passed or failed.
