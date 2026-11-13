---
title: 'phu: a self-contained command-line toolkit for screening and annotating viral genomes'
author:
  - Camilo García-Botero (ORCID 0000-0002-0426-7007; corresponding author)
  - Alejandro Reyes (ORCID 0000-0003-2907-3265)
bibliography: paper.bib
---

<!--
DRAFT for GigaByte, article type "Technical Release".
Structure follows https://gigabytejournal.com/technical-release-description
(checked 2026-11-12). GigaByte accepts Word or LaTeX (Overleaf template), not
Markdown: convert this file when the text is settled, for example
`pandoc manuscript-gigabyte.md --citeproc -o manuscript.docx`.
Items marked TODO need the authors' input. Author addresses and e-mails go on
the title page of the converted file.
Affiliation (TODO: confirm, add department): Universidad de los Andes, Bogotá, Colombia.
-->

# Abstract

<!-- Recommended length: 150 words, no references. -->

Viral genomes recovered from metagenomes are routinely screened for protein
families, grouped into species and genera, and searched for genes that may
alter host physiology. These steps are usually chained from separate programs
with shell scripts, which makes installation and reproducibility difficult.
We present `phu`, a Python command-line toolkit that performs them with one
interface. `phu screen` and `phu jack` select contigs whose predicted proteins
match profile hidden Markov models or seed proteins under an explicit
combination rule; `phu cluster` groups genomes at species and genus level;
`phu avger` ranks candidate auxiliary viral genes. Gene prediction, profile
search, alignment and trimming run in-process, and reference databases are
installed from one pinned, checksummed archive, so each run records exactly
what it used.

**Availability & Implementation:** `phu` is implemented in Python (3.11 or
later), released under the MIT licence, and available from Bioconda and PyPI.
Source code: https://github.com/camilogarciabotero/phu. Documentation:
https://camilogarciabotero.github.io/phu.

# Research Area

<!-- TODO: choose ONE research area and several classifications from
https://gigabytejournal.com/research-area (candidates: Software and Workflows;
classifications such as Bioinformatics, Microbiology, Virology, Metagenomics). -->

# Statement of Need

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

## Related software

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

# Implementation

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

## Commands

<!-- TODO: one short paragraph per command with its inputs and outputs, and a
figure (overview diagram of the four workflows). Figures must not be made or
altered with generative AI. -->

# Validation

<!--
TODO (authors). GigaByte requires that "the tool or method should have
undergone some testing, and this should be summarised in the report", with
test data and expected outputs available to reviewers. It does not need to
outperform other tools. Suggested content, all on PUBLIC data:
1. Correctness: `phu screen` contig selections against a plain
   pyrodigal + hmmsearch pipeline on a public phage genome set.
2. Runtime and memory versus threads for screen, jack and avger.
3. Alignment back-end: models built from FAMSA alignments compared with MAFFT
   alignments (recovery of family members, runtime), with and without
   gappyout trimming.
4. `phu avger` candidates compared with CheckAMG calls on a public virome.
5. Reproducibility: same command, same phu-db version, identical outputs.
Do not use unpublished project data here.
-->

`phu` has an automated test suite of more than 700 tests, run in continuous
integration on Python 3.11, 3.12 and 3.13 together with a build and install
check of the distributed package. <!-- TODO: results of the validation above. -->

# Availability of supporting source code and requirements

- Project name: phu
- Project home page: https://github.com/camilogarciabotero/phu
- Operating system(s): Linux and macOS <!-- TODO: confirm macOS is tested -->
- Programming language: Python
- Other requirements: Python 3.11 or later; vclust and SeqKit for clustering
  and sequence extraction; about 12 GB of free disk to install the reference
  database archive
- License: MIT
- RRID: <!-- TODO: register at SciCrunch.org -->
- bio.tools ID: <!-- TODO: register at bio.tools -->
- Version described: <!-- TODO: the release tagged for submission -->

# Data Availability

<!-- TODO: archive the submitted version and the validation data, then cite
them. Recommended wording follows. -->

Code snapshots and the datasets supporting the results of this article are
available in the Zenodo repository <!-- TODO: DOI of the software snapshot -->.
The `phu-db` reference database archive (version 0.1.1) is available from
Zenodo, https://doi.org/10.5281/zenodo.24121209.

# List of abbreviations

AMG: auxiliary metabolic gene; AVG: auxiliary viral gene; HMM: hidden Markov
model; KO: KEGG Orthology; MSA: multiple sequence alignment.

# Declarations

## Ethics approval and consent to participate

Not applicable.

## Competing interests

<!-- TODO: authors to state. -->

## Authors' contributions

<!-- TODO: CRediT roles for each author. -->

## Funding

<!-- TODO: funders and grant codes. -->

## Acknowledgements

<!-- TODO -->

# Disclosure of use of AI-assisted tools including generative AI

<!-- TODO (authors): GigaByte also asks for in-text citations of the AI tools
with a reference-list entry (tool, model version, date of use, user's name),
the AI outputs as a supplementary file in an open repository, and a mention
in the cover letter. Check every statement below against your own records. -->

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

# References
