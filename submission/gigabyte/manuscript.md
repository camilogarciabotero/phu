<!--
DRAFT for GigaByte, article type: Technical Release.
Structure follows https://gigabytejournal.com/technical-release-description
and the published example DOI 10.46471/gigabyte.194.

GigaByte accepts Word or LaTeX (Overleaf template), not Markdown: convert
this file before submission (for example `pandoc manuscript.md -o manuscript.docx`).
Submission layout: double line spacing, line and page numbers, figures as
separate files, tables in the main file.

Items marked TODO need the authors' input.
-->

# phu: a self-contained command-line toolkit for screening and annotating viral genomes

Camilo García-Botero^1,\*^ (ORCID 0000-0002-0426-7007), Alejandro Reyes^1^ (ORCID 0000-0003-2907-3265)

^1^ Universidad de los Andes, Bogotá, Colombia <!-- TODO: department, full postal address, emails -->

\* Corresponding author

## Abstract

<!-- Recommended 150 words. No references. -->

Viruses of bacteria and archaea are recovered from metagenomes as thousands of
contigs, and many analyses start with the same questions: which contigs carry
a protein family of interest, which belong to the same species, and which
genes may alter host physiology. These are usually answered with chains of
separate programs and scripts that are hard to install and to reproduce. We
present phu, a command-line toolkit that screens contigs with profile hidden
Markov models, searches iteratively from seed proteins, clusters genomes, and
prioritises candidate auxiliary viral genes. Gene prediction, profile search,
alignment and trimming run inside one Python process, reference databases come
from a single versioned and checksummed archive, and every run records its
parameters and database versions.

**Availability & Implementation:** phu is written in Python (3.11 or later),
released under the MIT licence, and available from PyPI and Bioconda. Source
code: https://github.com/camilogarciabotero/phu. Documentation:
https://camilogarciabotero.github.io/phu.

**Research Area:** <!-- TODO: choose one from https://gigabytejournal.com/research-area -->

**Classifications:** <!-- TODO: choose several from the same page (for example bioinformatics, microbiology, virology, metagenomics) -->

## Statement of need

Screening viral contigs for marker proteins is commonly done by joining a gene
caller, `hmmsearch`, a table parser and a sequence extractor with shell
scripts that differ between laboratories. Each program is one more thing to
install, which is often the hardest step on a shared computing cluster, and
each script is a place where parameters go unrecorded. Reference databases add
a second problem: Pfam [1], KOfam [2] and dbCAN are updated in place, so the
same command can give different results months apart unless versions are
tracked by hand.

Existing tools address neighbouring tasks. geNomad [3] and CheckV [4] identify
viral sequences and assess their quality. Pharokka [5], DRAM-v [6] and VIBRANT
[7] annotate every gene of a genome against fixed databases. None is designed
to answer "which of these contigs carry these specific models", with
user-supplied profiles and an explicit rule for combining hits, or to build
new profiles from the hits. For auxiliary viral genes, CheckAMG and the
V-score framework [8], following the terminology of Martin et al. [9], provide
the scores and curated tables; CheckAMG is a comprehensive pipeline with a
database of tens of gigabytes.

phu fills this gap for researchers who work with contigs already judged viral.
It replaces the script chain with documented, tested commands, and makes
explicit three things that ad hoc pipelines usually leave implicit: the rule
that turns protein hits into a set of contigs, the identity of the reference
databases, and the provenance of each run.

## Implementation

phu is a Python package with a command-line interface built on Typer. Its
commands share gene prediction, logging and database handling.

**Commands.** `phu screen` predicts proteins on contigs and keeps the contigs
that match profile hidden Markov models (HMMs), given as files or named from
Pfam, KOfam and dbCAN, including dbCAN polysaccharide utilisation loci. Hits
are combined per contig with a stated rule (any model, all models, or a
minimum number). Optionally it writes the matched proteins per model and
builds a new HMM from them. `phu jack` performs the same contig-level
screening starting from seed protein sequences by iterative search.
`phu cluster` groups genomes into species- and genus-level clusters using
vclust [10]. `phu avger` ranks candidate auxiliary viral genes on trusted
viral contigs from Pfam and KOfam annotations, V- and VL-scores [8] and the
reference tables curated by the CheckAMG project; it keeps annotation
evidence, virus-association evidence and class assignment in separate tables
and reports candidates for prioritisation, not confirmed auxiliary genes. Figure 1 summarises this workflow.
`phu normalize-lineage` standardises lineage strings produced by
gene-sharing-network classifiers.

**In-process libraries.** Gene prediction uses pyrodigal and pyrodigal-gv
[11,12], profile search uses PyHMMER [13,14], multiple sequence alignment uses
FAMSA through pyfamsa [15], and alignment trimming uses trimAl through
pytrimal [16]. Running these inside the process removes most executables from
the installation, avoids writing and re-parsing intermediate files, and lets
phu log the exact parameters of each step. The cost is a dependence on binary
wheels: the aligner binding is held to a release series whose wheels run on
the older C libraries still common on clusters. Two external programs remain,
vclust and SeqKit [17].

**Pinned reference databases.** All reference data are distributed as one
archive, phu-db, built by a separate preparation tool and deposited on Zenodo
[18]. Each phu release pins one archive by checksum; it is downloaded on first
use in parallel ranges with resumption, verified, and installed. Users
therefore cannot pick up a new Pfam release until a new archive is published,
in exchange for results that are reproducible for a given phu version.

**Caching and provenance.** Predicted proteins are cached under a key derived
from the input and the prediction parameters, because repeated screens of the
same contigs with different models are the common pattern. Each output folder
contains JSON manifests with the query, the parameters that determine the
result, the database version and file hashes, and the software version. A
verbose mode prints one line per tool invocation.

<!-- Figure 1: figures/avger-workflow.pdf (source: figures/avger-workflow.svg,
drawn by figures/make_avger_figure.py). TODO: it shows PHROG as a third
database, which is planned and not yet in a release; revise when it lands. -->

**Figure 1. The `phu avger` workflow.** (A) Genes are predicted on trusted
viral contigs. (B) Proteins are searched against Pfam, KOfam and PHROG in a
relaxed track used for scoring and a strict track used for functional calls.
(C) Each hit is given the V- and VL-score of its family, and VL-scores are
averaged per contig and per database (AVL). (D) A gene with low V and VL on a
contig with high AVL is a candidate; candidates are matched to curated
positive lists, filtered, and given a final status. Gene values are
illustrative.

**Testing.** An automated test suite (more than 700 tests) runs on Python
3.11, 3.12 and 3.13 in continuous integration, which also builds the package
and checks that the installed command starts.

## Results and discussion

<!--
TODO (authors). GigaByte asks that the tool "should have undergone some
testing, and this should be summarised in the report", with test data and
expected outputs available to reviewers. Everything here must use PUBLIC
data. Suggested content, to be produced and archived on Zenodo with scripts:

1. Use case on a public phage genome collection: screen for a marker family
   with `phu screen`, cluster the hits with `phu cluster`, report counts and
   run time.
2. Agreement with the conventional chain: same contigs kept by `phu screen`
   and by a gene caller + hmmsearch + manual parsing, with run time and the
   number of installed programs for each.
3. `phu avger` candidates compared with CheckAMG calls on the same public
   genomes, with run time and database size.
4. Effect of building HMMs from hits with FAMSA, with and without trimming,
   on public protein families.

One figure (workflow overview) and one table (run time and agreement) would
be enough for a Technical Release.
-->

## Availability of source code and requirements

- Project name: phu
- Project home page: https://github.com/camilogarciabotero/phu
- Documentation: https://camilogarciabotero.github.io/phu
- Operating system(s): Linux, macOS <!-- TODO: confirm macOS is tested -->
- Programming language: Python 3.11 or later
- Other requirements: vclust 1.0 or later and SeqKit for the `cluster` and `screen` commands; about 12 GB of free disk while the reference databases install
- License: MIT
- RRID: <!-- TODO: register at SciCrunch.org -->
- bio.tools ID: <!-- TODO: register at bio.tools -->
- Distribution: PyPI, Bioconda <!-- TODO: Bioconda is at 0.8.0; update the recipe before submission -->
- Support: GitHub issues

<!-- TODO: state the exact version described (tag a release) and cite its
archived snapshot, for example the Zenodo DOI of that release. -->

## Use of generative AI in manuscript and software

<!-- TODO (authors): GigaByte requires this to be declared here and in the
cover letter, with a numbered reference per tool giving version, model, date
of use and the user's name. Check every statement against your own records. -->

Generative AI tools were used in developing phu and in preparing this
manuscript. GitHub Copilot [19], including its coding agent and Autofix, and
Anthropic's Claude models used through the Claude Code agent [20] assisted
with writing and refactoring source code, writing tests, and drafting
documentation; commits made with such assistance carry co-author trailers in
the public repository history. A first draft of this manuscript was produced
with Claude and revised by the authors. Figure 1 was redrawn as a vector graphic by a
script written with Claude, from a diagram designed by the authors. <!-- TODO:
GigaByte restricts AI-made figures; confirm with the journal that a schematic
drawn this way is acceptable, or redraw it by hand from the SVG. --> The authors defined the scope of the software, made its design
decisions, reviewed and tested all AI-assisted contributions, and are fully
accountable for the software and the manuscript. <!-- TODO: confirm. -->

## Data availability

<!-- TODO: complete once the Results analyses exist. -->
Code snapshots and the datasets supporting the results of this article are
available in the Zenodo repository [TODO]. The reference database archive
phu-db is available from Zenodo [18].

## List of abbreviations

AMG: auxiliary metabolic gene; AVG: auxiliary viral gene; HMM: hidden Markov
model; KO: KEGG Orthology; PUL: polysaccharide utilisation locus.

## Declarations

### Ethics approval and consent to participate

Not applicable.

### Competing interests

The authors declare that they have no competing interests. <!-- TODO: confirm -->

### Funding

<!-- TODO -->

### Authors' contributions

<!-- TODO: CRediT roles. Suggested starting point:
C.G.-B.: Conceptualization, Software, Validation, Writing – original draft.
A.R.: Supervision, Funding acquisition, Writing – review & editing. -->

### Acknowledgements

<!-- TODO: computing resources and people who tested the software. -->

## References

<!-- Numbered in order of citation, Vancouver style, six authors then et al.
Entries marked (verify) were written from memory. URLs need an access date. -->

1. Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL, et al. Pfam: the protein families database in 2021. Nucleic Acids Res. 2021;49:D412-9. doi:10.1093/nar/gkaa913.
2. Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, et al. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2020;36:2251-2. doi:10.1093/bioinformatics/btz859.
3. Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, et al. Identification of mobile genetic elements with geNomad. Nat Biotechnol. 2024;42:1303-12. doi:10.1038/s41587-023-01953-y. (verify)
4. Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol. 2021;39:578-85. doi:10.1038/s41587-020-00774-7. (verify)
5. Bouras G, Nepal R, Houtak G, Psaltis AJ, Wormald PJ, Vreugde S. Pharokka: a fast scalable bacteriophage annotation tool. Bioinformatics. 2023;39:btac776. doi:10.1093/bioinformatics/btac776. (verify)
6. Shaffer M, Borton MA, McGivern BB, Zayed AA, La Rosa SL, Solden LM, et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Res. 2020;48:8883-900. doi:10.1093/nar/gkaa621. (verify)
7. Kieft K, Zhou Z, Anantharaman K. VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. Microbiome. 2020;8:90. doi:10.1186/s40168-020-00867-0. (verify)
8. Zhou K, et al. V- and VL-scores unveil viral signatures and origins of protein families. Nat Commun. 2026;17:5828. doi:10.1038/s41467-026-72028-0. (complete the author list)
9. Martin C, et al. A call for caution in the biological interpretation of viral auxiliary metabolic genes. Nat Microbiol. 2025;10:2122-9. doi:10.1038/s41564-025-02095-4. (complete the author list)
10. Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, et al. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. 2025. doi:10.1038/s41592-025-02701-7.
11. Larralde M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. J Open Source Softw. 2022;7:4296. doi:10.21105/joss.04296.
12. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010;11:119. doi:10.1186/1471-2105-11-119.
13. Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics. 2023;39:btad214. doi:10.1093/bioinformatics/btad214.
14. Eddy SR. Accelerated profile HMM searches. PLoS Comput Biol. 2011;7:e1002195. doi:10.1371/journal.pcbi.1002195.
15. Deorowicz S, Debudaj-Grabysz A, Gudyś A. FAMSA: fast and accurate multiple sequence alignment of huge protein families. Sci Rep. 2016;6:33964. doi:10.1038/srep33964.
16. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. 2009;25:1972-3. doi:10.1093/bioinformatics/btp348.
17. Shen W, Sipos B, Zhao L. SeqKit2: a Swiss army knife for sequence and alignment processing. iMeta. 2024;3:e191. doi:10.1002/imt2.191.
18. phu-db: reference database archive for phu, version 0.1.1. Zenodo. 2026. doi:10.5281/zenodo.24121209. <!-- TODO: authors and title as on the Zenodo record -->
19. GitHub. GitHub Copilot [Large language model coding assistant]. <!-- TODO: versions, dates of use, user's name, URL -->
20. Anthropic. Claude (Claude Opus 5.5 and Claude Fable 5.2, used through Claude Code) [Large language model]. <!-- TODO: dates of use, user's name, URL -->
