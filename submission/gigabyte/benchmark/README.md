# Benchmark: phu avger against CheckAMG annotate

Stage 1 compares a `phu avger` run with the example outputs that CheckAMG
v1.2 ships, on CheckAMG's own public example data (113 single-contig viruses
and 5 multi-contig vMAGs, 163 contigs, 3,840 proteins). CheckAMG is the
reference for agreement, not ground truth; the two tools are different
methods (see Table 2 of the manuscript).

## Reproduce

```bash
# Inputs and expected outputs, from the CheckAMG repository at tag v1.2:
#   examples/example_data/single_contig_viruses.fasta
#   examples/example_data/multi_contig_vMAGs/*.fasta   (concatenated)
#   examples/example_outputs/annotate/
phu avger -i all_contigs.fna -o phu-out -t 16   # phu 636a9e5, phu-db 0.2.0
python3 compare.py phu-out checkamg-expected
```

`stage1-report.txt` is the run and comparison output (phu 636a9e5,
16 CPUs, 106 s). Proteins are matched by contig and gene ordinal
(3,829 of 3,840 matched).

Caveats recorded for the manuscript: phu's reference tables come from
CheckAMG v1.1.1 while the expected outputs are v1.2, so some disagreement
may come from table changes between those versions; and the shipped outputs
were produced by the CheckAMG authors, not rerun here (stage 2 reruns
CheckAMG annotate on the same machine for timing).
