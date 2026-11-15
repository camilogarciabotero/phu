# cluster

Sequence clustering wrapper around the external
[`vclust`](https://github.com/refresh-bio/vclust) tool, with four preconfigured
modes for common viral genomics workflows. Each mode follows one of vclust's
documented [use cases](https://github.com/refresh-bio/vclust/wiki/6-Use-cases);
`--vclust-params` adjusts them.

## Synopsis

```bash
phu cluster --mode <MODE> --input-contigs <FASTA_FILE> [OPTIONS]
```

## Modes

- `dereplication` — remove redundant sequences, keeping one representative per
  group (greedy CD-HIT-like clustering on ANI and coverage).
- `votu` — viral Operational Taxonomic Units with MIUViG thresholds (Leiden
  clustering on ANI ≥ 95% over ≥ 85% of the shorter sequence).
- `genus` — genus-level groups of complete genomes: complete linkage on total
  ANI (tANI) ≥ 70%.
- `species` — species-level groups of complete genomes: complete linkage on
  tANI ≥ 95%.

`genus` and `species` follow vclust's VIRIDIC-like recipe ("Classify viruses
into species and genera"), which uses the 70% and 95% genome-wide identity
thresholds applied to bacteriophages.

```bash
phu cluster --mode genus -i genomes.fna -o phu-genus
phu cluster --mode species -i genomes.fna -o phu-species
```

## Default settings by mode

| Setting | dereplication | votu | genus | species |
|---|---|---|---|---|
| Algorithm | cd-hit | leiden | complete | complete |
| Metric and threshold | ANI ≥ 0.95 | ANI ≥ 0.95 | tANI ≥ 0.70 | tANI ≥ 0.95 |
| Query coverage | ≥ 0.85 | ≥ 0.85 | none | none |
| Prefilter `--min-ident` | 0.95 | 0.95 | 0.70 | 0.70 |
| Align output filter | ANI ≥ 0.95, qcov ≥ 0.85 | ANI ≥ 0.95, qcov ≥ 0.85 | none (all aligned pairs) | none (all aligned pairs) |
| Cluster table | `clusters.tsv` | `clusters.tsv` | `genus-clusters.tsv` | `species-clusters.tsv` |
| Representatives | `dereplicated-representatives.fna` | `cluster-representatives.fna` | `cluster-representatives.fna` | `cluster-representatives.fna` |

The align output filter of `dereplication` and `votu` always equals their
clustering thresholds, so it only shrinks the ANI table and never removes a
pair the clustering would use.

## ANI, tANI and the clustering algorithm

- **ANI** is the identity of the aligned regions only. Two genomes that share
  one short, nearly identical region have a high ANI, so ANI modes also
  require a query coverage (the aligned fraction of the shorter sequence).
- **tANI** (total ANI) counts identical bases over the full length of both
  genomes, so unaligned sequence lowers it; it combines identity and coverage
  in one number, as VIRIDIC does. A genome sharing 80% of its length with
  another at 99.6% identity has an ANI near 0.99 but a tANI near 0.80: same
  genus, different species.
- **Complete linkage** (`genus`, `species`) puts two genomes in the same group
  only when every pair in the group passes the threshold. Single linkage
  would chain genomes that are each similar to a neighbour but not to each
  other; Leiden (`votu`) and CD-HIT (`dereplication`) are graph and
  greedy-representative methods suited to large sets of contigs.

## Use with complete genomes

Genus and species groups are **operational**: they are numerical clusters, not
ICTV taxa, and phu does not assign taxon names. Use them with complete or
near-complete genomes. tANI is computed over each sequence's full length, so a
fragment of a genome has a low tANI to the complete genome it came from and
can end up in its own genus or species group. For metagenomic contigs of
mixed completeness, use `votu` (ANI with a coverage requirement) instead, and
treat genus and species groups of incomplete sequences as provisional.

## Outputs

```
phu-genus/
├── prefilter-pairs.txt             # vclust prefilter: sequence pairs worth aligning
├── ani-pairs.tsv                   # vclust align: ANI, tANI and coverage per sequence pair
├── ani-sequence-ids.tsv            # vclust align: sequence IDs and lengths used by clustering
├── genus-clusters.tsv              # cluster assignments (see the table above for each mode)
├── cluster-representative-ids.txt  # one representative ID per input sequence's cluster
└── cluster-representatives.fna     # representative sequences
```

Each cluster's representative is its longest sequence (vclust `--out-repr`).

## Advanced: custom vclust parameters

`--vclust-params` takes vclust flags as one quoted string and passes each to
its step (prefilter, align or cluster). Settings are resolved in this order,
later ones winning: the mode's defaults, then `--vclust-params`. A flag given
there replaces the default instead of being added a second time.

```bash
# Stricter genus groups
phu cluster --mode genus -i genomes.fna --vclust-params="--tani 0.75"

# vOTUs at 97% ANI; the align filter follows the new threshold
phu cluster --mode votu -i contigs.fna --vclust-params="--ani 0.97"

# Larger datasets: fewer k-mers and a compact ANI table
phu cluster --mode votu -i contigs.fna --vclust-params="--min-kmers 20 --outfmt lite"
```

phu stops before running vclust when the parameters contradict each other:

| Combination | Why it is rejected |
|---|---|
| An align filter above the clustering threshold of the same measure, e.g. `--mode genus` with `--out-ani 0.95` | Removes pairs the clustering needs. |
| `--metric` without that metric's threshold, e.g. `--mode genus` with `--metric ani` | A tANI threshold would otherwise be read as an ANI one, or none used. |
| A threshold below the mode's prefilter identity without `--min-ident`, e.g. `--mode votu` with `--ani 0.90` | The prefilter at 0.95 would drop pairs before alignment; add `--min-ident 0.9`. |

## Known issue: vclust 1.3.1 with several threads

With vclust 1.3.1, `vclust align --filter` run on two or more threads can skip
genome pairs that the prefilter passed. In a three-genome test it wrote one
of three pairs with `--threads 2` and all three with `--threads 1`. Until this
is fixed upstream, `-t 1` gives complete results for small inputs; compare a
run with `-t 1` when the groups matter.

## Command options

The live CLI groups options into **Clustering**, **Input**, **Output**, and
**Runtime** panels. Run `phu cluster --help` for the current terminal menu.

```

 Sequence clustering wrapper around external 'vclust'.

 Modes: dereplication and votu (ANI >= 0.95, coverage >= 0.85), genus
 (tANI >= 0.70) and species (tANI >= 0.95), both with complete linkage.
 Genus and species groups are operational, not ICTV names; use complete or
 near-complete genomes.

 For advanced usage, provide custom vclust parameters as a quoted string.
 See the vclust wiki for parameter details:
 https://github.com/refresh-bio/vclust/wiki

 Example:
     phu cluster --mode genus -i genomes.fna -o phu-genus
     phu cluster --mode votu --input-contigs genomes.fna
 --vclust-params="--min-kmers 20 --outfmt lite"

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --help  -h        Show this message and exit.                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Clustering ─────────────────────────────────────────────────────────────────╮
│ *  --mode                   <str>  dereplication | votu | genus | species    │
│                                    [required]                                │
│    --vclust-params  -p      <str>  Custom vclust parameters: "--min-kmers 20 │
│                                    --outfmt lite --ani 0.97"                 │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Input ──────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs  -i      <path>  Input FASTA [required]                   │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ --output-folder  -o      <path>  Output directory [default: phu-cluster]     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --threads  -t      <int range> [x>=0]  0=all cores; otherwise N threads      │
│                                        [default: 0]                          │
│ --quiet                                Suppress routine progress output.     │
│ --verbose                              Show additional progress details.     │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## See also

- vclust use cases: https://github.com/refresh-bio/vclust/wiki/6-Use-cases
- Zielezinski A, et al. Ultrafast and accurate sequence alignment and
  clustering of viral genomes. Nat Methods (2025).
  https://doi.org/10.1038/s41592-025-02701-7
