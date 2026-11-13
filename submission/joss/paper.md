---
title: 'phu: a self-contained command-line toolkit for screening and annotating viral genomes'
tags:
  - Python
  - bioinformatics
  - bacteriophages
  - viral genomics
  - metagenomics
  - profile hidden Markov models
authors:
  - name: Camilo García-Botero
    orcid: 0000-0002-0426-7007
    corresponding: true
    affiliation: 1
  - name: Alejandro Reyes
    orcid: 0000-0003-2907-3265
    affiliation: 1
affiliations:
  - name: Universidad de los Andes, Bogotá, Colombia
    index: 1
date: 12 November 2026
bibliography: paper.bib
---

<!--
DRAFT. Items marked TODO need the authors' input before submission.
JOSS limits: 750-1750 words. Required sections: Summary, Statement of need,
State of the field, Software design, Research impact statement,
AI usage disclosure, Acknowledgements, References.
-->

# Summary

Viruses of bacteria and archaea are recovered from metagenomes as thousands of
assembled contigs, and most analyses then ask the same few questions of them:
which contigs carry a protein family of interest, which contigs belong to the
same species or genus, and which genes might alter host physiology. `phu`
(phage utilities) is a Python command-line toolkit that answers these
questions with one consistent interface. `phu screen` predicts proteins on
contigs and keeps the contigs that match profile hidden Markov models (HMMs),
supplied by the user or named from Pfam, KOfam and dbCAN. `phu jack` does the
same starting from seed protein sequences, by iterative search. `phu cluster`
groups genomes into species- and genus-level clusters. `phu avger` ranks
candidate auxiliary viral genes using published virus-association scores and
curated reference tables. Gene prediction, HMM search, multiple sequence
alignment and alignment trimming run inside the Python process, and the
reference databases are installed from one versioned, checksummed archive, so
a run needs few external programs and records what it used.

# Statement of need

Screening viral contigs for marker proteins is usually done by chaining
separate programs: a gene caller, `hmmsearch`, a table parser, and a sequence
extractor, joined by shell scripts that differ between laboratories. Each link
adds an installation requirement and a place where parameters go unrecorded.
Reference databases add a second problem: Pfam, KOfam and dbCAN are updated in
place, so the same command can give different results months apart unless the
user tracks versions by hand.

`phu` targets researchers who work with viral contigs on shared computing
clusters, where installing compiled tools is often the hardest step. It
replaces the script chain with single commands whose behaviour is documented
and tested, and it makes three things explicit that ad hoc pipelines usually
leave implicit:

- **Contig-level decisions.** Hits are combined per contig (`any`, `all`, or
  a minimum number of models), so the output is the set of contigs that
  satisfy a stated rule rather than a table of protein hits to post-process.
- **Database identity.** All reference data come from a single pinned archive,
  `phu-db`, deposited on Zenodo and verified by checksum on installation. Each
  run writes the database version and file hashes to a manifest.
- **Run provenance.** Each output folder contains JSON manifests with the
  query, the parameters that determine the result, and the software version.

Target users are microbiologists and bioinformaticians studying phages and
other prokaryotic viruses in metagenomes. <!-- TODO: one sentence on the
research questions phu has been used for in the authors' group, described
without unpublished results. -->

# State of the field

Several mature tools address neighbouring problems. geNomad
[@camargo2024genomad] and CheckV [@nayfach2021checkv] identify viral sequences
and assess their quality; `phu` assumes that step is done and works on contigs
already judged viral. Pharokka [@bouras2023pharokka], DRAM-v
[@shaffer2020dram] and VIBRANT [@kieft2020vibrant] annotate every gene of a
genome against fixed databases and report functions; they are not designed to
answer "which of these contigs carry these specific models" with user-supplied
HMMs and an explicit combination rule, nor to build new HMMs from the hits.
CheckAMG and the V-score framework [@zhou2026vscores] provide the scores and
curated tables for auxiliary viral gene analysis, following the terminology of
@martin2025caution; CheckAMG is a comprehensive pipeline with a database of
tens of gigabytes, whereas `phu avger` is a deliberately small implementation
of the candidate-selection rules on top of Pfam and KOfam for exploratory use.
vclust [@zielezinski2025vclust] performs genome clustering; `phu cluster`
calls it with settings for the standard species and genus thresholds and
standardised outputs.

`phu` was built as a separate package rather than as contributions to these
tools because its purpose differs from each of them: it is a thin, uniform
layer for targeted screening that composes existing libraries, not a new
annotation database or a new algorithm. <!-- TODO: authors to confirm this
positioning and add any tool they compared against in practice. -->

# Software design

**In-process libraries instead of external programs.** `phu` uses the Cython
bindings pyrodigal and pyrodigal-gv for gene prediction [@larralde2022pyrodigal;
@hyatt2010prodigal], PyHMMER for `hmmsearch` and `jackhmmer`
[@larralde2023pyhmmer; @eddy2011hmmer], pyfamsa for alignment
[@deorowicz2016famsa] and pytrimal for trimming [@capella2009trimal]. This
removes most executables from the installation, avoids writing and re-parsing
intermediate files, and lets `phu` log the exact parameters of each step. The
cost is that `phu` depends on the availability of binary wheels: for example,
the aligner binding is held to a release series whose wheels run on the older
C libraries still common on clusters. Two external programs remain, vclust and
SeqKit [@shen2024seqkit], for clustering and sequence extraction.

**One pinned database archive instead of per-database downloaders.** Earlier
versions downloaded Pfam [@mistry2021pfam], KOfam [@aramaki2020kofam] and dbCAN
separately from their sources. This was replaced by a single archive built by a
separate preparation tool, pinned by checksum in each `phu` release and
downloaded on first use, in parallel ranges with resumption. The trade-off is
deliberate: users cannot pick up a new Pfam release until a new archive is
published, in exchange for results that are reproducible for a given `phu`
version.

**Cached gene prediction.** Predicted proteins are cached under a key derived
from the input and the prediction parameters, because repeated screens of the
same contigs with different models are the common pattern and prediction would
otherwise be repeated.

**Conservative interpretation in `avger`.** The command separates annotation
evidence, virus-association evidence and class assignment into different
output tables, and its documentation states that results are candidates for
prioritisation, not confirmed auxiliary genes.

The package has an automated test suite run on three Python versions in
continuous integration, a documentation site with one page per command, and
is distributed through PyPI and Bioconda.

# Research impact statement

<!--
TODO (authors): JOSS asks for specific, verifiable evidence here, not
aspirations. Replace this block with what applies:
- Bioconda / PyPI download counts and the date they were read.
- Theses, preprints, papers or pipelines that use phu (with references).
- External users or groups (issues opened by people outside the group).
- A reproducible reference analysis on PUBLIC data, with a link, showing a
  result obtained with phu (for example, screening a public phage genome
  collection for a marker family and clustering the hits).
Do not describe unpublished project data here.
-->

`phu` has been publicly developed since September 2025, with tagged releases
archived on Zenodo and distributed through Bioconda. <!-- TODO: evidence. -->

# AI usage disclosure

<!-- TODO (authors): JOSS requires this to be accurate and complete. Check
every statement below against your own records and edit before submission. -->

Generative AI tools were used in the development of `phu`. GitHub Copilot
(including its coding agent and Autofix) and Anthropic's Claude models, used
through the Claude Code command-line agent (Claude Opus 5.5 and Claude Fable
5.2), assisted with writing and refactoring source code, writing tests, and
drafting documentation; commits made with such assistance carry co-author
trailers in the public history. A first draft of this paper was also produced
with Claude and then revised by the authors. <!-- TODO: state any other tool
or model version used, and the approximate period. --> The authors defined the
scope of the software, made the architectural and design decisions described
above, reviewed and tested all AI-assisted contributions, validated results
against real data, and take responsibility for the content of the software
and of this paper. <!-- TODO: authors to confirm this affirmation. -->

# Acknowledgements

<!-- TODO: funding, computing resources, and people who tested the software. -->

# References
