# screen thresholds and decision logic

This page documents how `phu screen` decides whether a hit, protein, or contig passes filtering.

It is meant as an implementation-aligned reference for PFAM and KOfam behavior.

These are computational acceptance thresholds for profile annotation, not
biological proof of a viral sequence or an auxiliary viral gene. The published
AV/AVᴸ thresholds classify viral genomes or fragments and are a separate
upstream step. AVG interpretation requires additional viral sequence context
and curated scientific rules; see the
[avger scientific interpretation](avger.md#scientific-interpretation).

## Scope

These rules apply to the screening workflow in `phu screen` after protein prediction and HMM search.

## Core pass/fail flow

For each hit emitted by pyHMMER:

1. Start from hits marked as included by pyHMMER.
2. Compute an effective score.
3. Compute an effective minimum bitscore threshold.
4. Apply score and E-value filters.
5. Group remaining hits by contig and apply `--combine-mode` rules.

Only contigs with at least one remaining hit after all filters can be kept.

## Which score is used

Each hit has:

- Full-sequence bitscore (HMMER "full sequence score").
- Domain bitscore derived as the maximum score among included domains for that hit.

Effective score selection:

- Default: use full-sequence bitscore.
- KOfam model with `score_type = domain`: use domain bitscore when available.
- KOfam model with `score_type = full`: use full-sequence bitscore.

## Threshold precedence

Let:

- `min_bitscore` be the CLI value from `--min-bitscore` (can be unset).
- `ko_threshold` be the KOfam threshold from `ko_list` for that KO (can be missing).

If `--use-kofam-thresholds` is enabled and `ko_threshold` exists:

- If `min_bitscore` is unset: effective minimum bitscore is `ko_threshold`.
- If `min_bitscore` is set: effective minimum bitscore is `max(min_bitscore, ko_threshold)`.

If KOfam thresholds are disabled (or no KO threshold exists), effective minimum bitscore is just `min_bitscore`.

This means user thresholds can only make filtering stricter when KOfam thresholds are active.

## E-value behavior

`--max-evalue` is always applied using the hit independent E-value from the top-level hit.

Important: even when KOfam `score_type` is `domain`, the E-value filter is still based on the hit-level E-value, not domain i-Evalue or c-Evalue from domtblout rows.

## PFAM behavior

PFAM accessions are resolved to local models, then screened like any other HMM model.

Threshold behavior for PFAM depends on CLI options:

- `--cut-ga` on (default): pyHMMER applies profile GA gathering cutoffs during search.
- `--no-cut-ga`: no model GA cutoff is forced by pyHMMER; filtering relies on `--min-bitscore` and `--max-evalue`.

PFAM does not use KOfam `ko_list` thresholds.

## KOfam behavior

KOfam models are resolved by KO ID and enriched with metadata parsed from `ko_list`, including:

- `threshold`
- `score_type` (`full` or `domain`)

When `--use-kofam-thresholds` is enabled (default), KOfam thresholding is applied per KO using the KO `score_type` logic above.

## domtblout interpretation

`--keep-domtbl` keeps raw `domtblout` files for inspection and audit.

In current implementation, pass/fail filtering does not re-parse domtblout text. The selection is performed from in-memory hit objects produced by pyHMMER.

Use domtblout as an audit artifact to interpret why a hit likely passed or failed.

## combine mode after filtering

After score/E-value filtering, contigs are retained by combine mode:

- `any`: keep contigs with at least one passing model; keep top hits per model per contig.
- `all`: keep contigs that match all models.
- `threshold`: keep contigs with at least `--min-hmm-hits` distinct matching models.

## Worked interpretation example

If KO metadata says:

- `threshold = 136.43`
- `score_type = full`

Then hits with full-sequence scores around 15-20 fail thresholding, even if domain scores look reasonable in domtblout.

If `--max-evalue` remains default (`1e-5`), many such hits also fail the E-value filter.
