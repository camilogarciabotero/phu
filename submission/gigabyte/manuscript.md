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
a second problem: Pfam [1], KOfam [2] and dbCAN [21] are updated in place, so the
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
builds a new HMM from them (Figure 1). `phu jack` performs the same contig-level
screening starting from seed protein sequences by iterative search (Figure 2).
`phu cluster` groups genomes into species- and genus-level clusters using
vclust [10]. `phu avger` ranks candidate auxiliary viral genes on trusted
viral contigs from Pfam, KOfam and PHROG [22] annotations, V- and VL-scores
[8] and the reference tables curated by the CheckAMG project; it keeps
annotation evidence, virus-association evidence and class assignment in
separate tables and reports candidates for prioritisation, not confirmed
auxiliary genes (Figure 3; details below).
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

**Auxiliary viral gene candidates.** <!-- TODO: PHROG support in `phu avger`
is being implemented; recheck every statement in this paragraph against the
merged code and its documentation before submission. --> `phu avger` searches
the predicted proteins against Pfam, KOfam and PHROG and keeps two sets of
hits. A relaxed track (E-value at most 1e-5) is used only for scoring: the
best hit of each protein in each database is given the V- and VL-score of its
family, and VL-scores are averaged per contig and per database (AVL). A strict
track is used for functional calls: Pfam gathering thresholds, KOfam adaptive
thresholds, and for PHROG a bit score of at least 60 with at least half of
both the model and the protein aligned. Coverage is computed as in CheckAMG,
over all reported domains of a hit; the cutoff is stricter than CheckAMG's
default of 0.4. Each database is
searched once and each hit is tested against the thresholds of both tracks;
earlier versions ran one search per track, and removing the second search
halved the run time on a test set of 21 viral contigs (5.9 Mb; 485 s to 245 s
on 16 CPUs) with identical output tables. A gene is a candidate when, in at
least one database, its contig has AVL above 3 while the gene itself has VL
below 3 and V below 10; the three databases are scored and averaged
separately and their values are never mixed. Candidates are then matched to
the CheckAMG positive lists through their strict hits. The lists contain
PHROG entries only for auxiliary physiological genes (77 families) and
auxiliary regulatory genes (655 families), so PHROG can support those two
classes but not an auxiliary metabolic gene call, which still requires Pfam or
KOfam evidence and a minimum AMG weight (0.6 by default). A strict hit to a
family on the CheckAMG filter lists, in any of the three databases, removes
the corresponding class. PHROG is used by default and can be switched off.

**Relationship to CheckAMG.** `phu avger` takes its scores, positive lists,
AMG weights and exclusion lists from CheckAMG [23], and CheckAMG is the
inspiration for the command. The two are nevertheless different methods and
should not be expected to give the same predictions (Table 2). `phu avger` is
a rule-based screen for contigs already known to be viral: a gene must first
pass a fixed rule on V- and VL-scores, computed over the whole contig and
separately for each database, and only then is its annotation matched to the
reference lists. CheckAMG is a broader prediction framework. It annotates
with more databases, judges viral origin from the local neighbourhood of each
gene (flanking viral and mobile-element genes and delineated viral regions)
with a trained classifier, and offers a separate annotation-independent
module based on a protein language model. Neither tool is uniformly the more
stringent: `phu avger` is stricter on candidate eligibility and on KOfam
score thresholds, whereas CheckAMG is stricter on alignment coverage for Pfam
and KEGG hits and on genomic context. The main limitation of `phu avger`
follows from its rule: a whole-contig average can hide a locally atypical
region, and a convincing local viral neighbourhood can sit on a contig whose
average fails. `phu avger` is therefore best described as integrating
V/VL-score-based candidate selection with CheckAMG-derived reference tables,
weights and exclusion rules in a reproducible annotation workflow, and not as
a lightweight implementation of CheckAMG.

<!-- TODO: Table 2 is a source-code comparison of phu (with PHROG support,
pull request 32) and CheckAMG main at commit 66339f9. Before submission:
(1) fix the CheckAMG version compared and state it; phu's tables come from
v1.1.1; (2) re-verify the rows marked [verify] against that version's code;
(3) replace or complement this paragraph with the matched-input benchmark
planned for Results. -->

**Table 2. `phu avger` and CheckAMG compared from their source code.** This
is a comparison of methods, not of performance.

| Aspect | `phu avger` | CheckAMG |
|---|---|---|
| Intended input | Contigs already identified as viral | Nucleotide or protein input; single contigs or multi-contig viral genomes; also broader metagenomic input |
| Functional evidence | Pfam, KOfam and PHROG | Pfam, KEGG, PHROG, FOAM, dbCAN, METABOLIC and CAMPER |
| Candidate logic | Mandatory rule on contig-average VL and gene V and VL, before any classification | Reference-supported functional candidates combined with contextual assessment and curation |
| Viral context | Whole-contig mean VL-score, per database | Local windows, flanking viral signals, hallmark and mobile-element genes, viral-region delineation |
| Confidence | Rule outcomes with audit reasons | Viral-origin probability from a gradient-boosted classifier, reported as low, medium or high [verify] |
| Contextual filters | None | Optional exclusion outside strict viral regions; runs of three or more consecutive candidates excluded by default [verify] |
| Annotation thresholds | Pfam gathering thresholds; KOfam adaptive thresholds; PHROG bit score and coverage | Coverage of at least 0.4 of model and protein for all hits; database thresholds, with a heuristic that can keep some sub-threshold KEGG hits |
| Proteins without annotation | Kept as unclassified candidates when the score rule holds; no annotation-independent predictor | Separate de novo module using protein language model embeddings |
| Conflicting classes | Reported as a class conflict | Resolved by priority: metabolic, then physiological, then regulatory [verify] |
| Software stack | One Python package; no machine-learning inference | Snakemake workflow with machine-learning dependencies; the de novo module benefits from a GPU |

**Pinned reference databases.** All reference data are distributed as one
archive, phu-db, built by a separate preparation tool and deposited on Zenodo
[18] (Table 1). Each phu release pins one archive by checksum; it is
downloaded on first use in parallel ranges with resumption, verified, and
installed, and a copy already on disk can be supplied through the
`PHU_DB_BUNDLE` environment variable for machines without internet access.
The archive carries a manifest, `phu-db.json`, that records for every file
the upstream URL, SHA-256, size, record count and the steps used to prepare
it, so that the archive can be rebuilt and compared. Users therefore cannot
pick up a new Pfam release until a new archive is published, in exchange for
results that are reproducible for a given phu version.

PHROGs are distributed as HH-suite profiles and MMseqs2 databases, which
HMMER cannot search. To keep a single search engine, the archive contains
HMMER models built from the official alignment of each PHROG with PyHMMER,
using the default `hmmbuild` settings and one model per PHROG. The build
fixes the date written into each model, so building twice gives the same
bytes and the result can be pinned by checksum like any other file.

**Table 1. Reference databases in the phu-db archive, version 0.2.0.** The
archive is 2.2 GB compressed and about 12.8 GB installed.

| Database | Upstream version | Content | Used by |
|---|---|---|---|
| Pfam [1] | 37.0 | 21,979 profile HMMs | `screen`, `avger` |
| KOfam [2] | archive of 2026-07-02 | 27,816 profile HMMs with per-model score thresholds | `screen`, `avger` |
| dbCAN [21] | HMMdb V15 | 881 profile HMMs; family substrates and polysaccharide utilisation locus rules | `screen` |
| CheckAMG tables [8] | v1.1.1 | V- and VL-scores; positive and filter lists for three classes of auxiliary viral gene | `avger` |
| PHROGs [22] | v4 | 38,880 profile HMMs built for phu-db from the official alignments; annotation table | `avger` |

<!-- TODO: Table 1 and reference 18 describe phu-db 0.2.0, which is staged on
Zenodo but not yet published. Add its version DOI to reference 18 once it is
published, and confirm that the phu release described here pins it. -->

**Caching and provenance.** Predicted proteins are cached under a key derived
from the input and the prediction parameters, because repeated screens of the
same contigs with different models are the common pattern. Each output folder
contains JSON manifests with the query, the parameters that determine the
result, the database version and file hashes, and the software version. A
verbose mode prints one line per tool invocation.

<!-- Figures are in figures/ as PDF, with the SVG source and the script that
draws each one (make_screen_figure.py, make_jack_figure.py,
make_avger_figure.py; shared helpers in figure_kit.py).
Figure 1: screen-workflow.pdf. Figure 2: jack-workflow.pdf.
Figure 3: avger-workflow.pdf. TODO: Figure 3 shows PHROG as a third
database, which is being implemented and is not yet in a release; recheck the
figure and its legend against the merged code. -->

**Figure 1. The `phu screen` workflow.** (A) The query can mix HMM files with
Pfam, KOfam and dbCAN names, which are resolved from the pinned database
archive; genes are predicted on the contigs. (B) Each model is searched
against the predicted proteins. (C) The models that hit a contig decide
whether it is kept, under the chosen combination rule. (D) Optionally, the
proteins hit by each model are aligned, trimmed and built into a new model.
Contigs, genes and hits are illustrative.

**Figure 2. The `phu jack` workflow.** (A) Seed proteins replace profile
models as the query. (B) Each seed is searched iteratively; hits included in
one round extend the profile searched in the next. (C) Final hits are
filtered and combined per contig. (D) Optionally, the final profile is saved
for use as a query in `phu screen`. Contigs, proteins and hits are
illustrative.

**Figure 3. The `phu avger` workflow.** (A) Genes are predicted on trusted
viral contigs. (B) Proteins are searched once against each of Pfam, KOfam and PHROG;
the hits are kept in a relaxed track used for scoring and a strict track used
for functional calls, which differ only in the thresholds applied.
(C) Each hit is given the V- and VL-score of its family, and VL-scores are
averaged per contig and per database (AVL). (D) A gene with low V and VL on a
contig with high AVL in at least one database is a candidate; candidates are
matched to curated positive lists, filtered, and given a final status. Gene values are
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
3. `phu avger` candidates compared with CheckAMG `annotate` calls on the same
   public genomes, with run time and database size. CheckAMG ships example
   inputs (113 single-contig viruses and 5 multi-contig viral genomes) with
   expected outputs, which allow a first comparison without installing it.
   Report agreement per class with CheckAMG as the reference, not as ground
   truth, and attribute each disagreement to a cause from Table 2. Compare
   with the de novo module separately, if at all.
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
- Other requirements: vclust 1.0 or later and SeqKit for the `cluster` and `screen` commands; about 15 GB of free disk while the reference databases install (2.2 GB download, 12.8 GB installed)
- License: MIT
- RRID: <!-- TODO: register at SciCrunch.org -->
- bio.tools ID: <!-- TODO: register at bio.tools -->
- Distribution: PyPI, Bioconda <!-- TODO: the Bioconda recipe for 0.9.0 was corrected in November 2026; the release described here needs its own recipe update (Python 3.11, pyfamsa, pytrimal, no MAFFT) before submission -->
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
with Claude and revised by the authors. Figures 1 to 3 are vector schematics drawn by
scripts written with Claude; Figure 3 follows a diagram designed by the
authors. <!-- TODO:
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

AMG: auxiliary metabolic gene; APG: auxiliary physiological gene; AReG:
auxiliary regulatory gene; AVG: auxiliary viral gene; AVL: average VL-score of
a contig; HMM: hidden Markov model; KO: KEGG Orthology; PHROG: prokaryotic
virus remote homologous group; PUL: polysaccharide utilisation locus.

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
18. García-Botero C. Phu Databases, version 0.2.0. Zenodo. 2026. doi:10.5281/zenodo.19503347. <!-- TODO: this is the concept DOI, which always resolves to the latest version. Replace it with the version DOI of 0.2.0 once that version is published (0.1.1 is doi:10.5281/zenodo.24121209). -->
19. GitHub. GitHub Copilot [Large language model coding assistant]. <!-- TODO: versions, dates of use, user's name, URL -->
20. Anthropic. Claude (Claude Opus 5.5 and Claude Fable 5.2, used through Claude Code) [Large language model]. <!-- TODO: dates of use, user's name, URL -->
21. Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Res. 2023;51:W115-21. doi:10.1093/nar/gkad328. (verify)
22. Terzian P, Olo Ndela E, Galiez C, Lossouarn J, Pérez Bucio RE, Mom R, et al. PHROG: families of prokaryotic virus proteins clustered using remote homology. NAR Genom Bioinform. 2021;3:lqab067. doi:10.1093/nargab/lqab067. (verify)

23. Anantharaman Lab. CheckAMG: automated discovery and curation of auxiliary metabolic, regulatory and physiology genes encoded by viral genomes. GitHub. https://github.com/AnantharamanLab/CheckAMG. <!-- TODO: the CheckAMG README says a manuscript is in preparation and asks to cite the repository and its Zenodo database record (https://zenodo.org/records/18407279) until then; cite the version compared, and check whether reference 8 was meant to cover CheckAMG where it is cited beside "CheckAMG tables". -->

<!-- TODO: references 21 to 23 were added after the first numbering; renumber
all references in order of first citation when converting to Word or LaTeX. -->
