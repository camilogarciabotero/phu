# screen

## What does it do?

The `phu screen` command helps you find DNA contigs that contain specific protein families. It predicts proteins from your contigs, searches those proteins against Hidden Markov Model (HMM) profiles using **pyHMMER** (a fast Python implementation), and then selects contigs based on the matches and your combination/filtering rules. Think of it as a molecular search engine for pulling out contigs that contain the proteins you care about.

This is especially useful when you have metagenomic assemblies and want to pull out contigs that belong to viruses, or when you're looking for contigs that contain specific metabolic pathways. Target-protein extraction and custom HMM generation are experimental features and are not a validated replacement for a biologically aligned profile-HMM workflow.

For implementation-level details on pass/fail decisions, see [screen thresholds and decision logic](screen-thresholds.md).

## Synopsis

```bash
phu screen --input-contigs [INPUT_CONTIGS] [HMMS...]
```

`[HMMS...]` accepts local HMM paths, PFAM/KOfam identifiers, dbCAN families
(for example `GH128` or `CBM89`), and one dbCAN PUL identifier such as
`PUL0621`.

**Example:**
```bash
phu screen --input-contigs your_contigs.fasta your_protein_family.hmm
```

This simple command will find all contigs in `your_contigs.fasta` that contain proteins matching `your_protein_family.hmm` and save them to a new file called `screened_contigs.fasta` in a folder named `phu-screen`.

If an argument looks like a PFAM accession (`PF` + 5 digits, optionally with a version), `phu` resolves it from a local Pfam-A database automatically.

dbCAN family inputs are case-insensitive and are validated against the local
dbCAN HMM index. A PUL input expands to its ordered CAZyme-family signature
from `cazymes_predicted_dbcan`; by default every resolved family must match on
the same contig. Use `--combine-mode threshold` explicitly for partial
signatures. Only one PUL can be queried per run, and PUL queries cannot be
mixed with other database or local-HMM inputs.

This is PUL CAZyme-signature screening. Family co-occurrence alone does not
identify a biological PUL or predict its substrate. Each run writes
`query_manifest.json` with the original query, expansion, matching rule, and
database hashes.

To screen against every resolvable PUL signature:

```bash
phu screen -i contigs.fa --all-puls
```

Each unique required CAZyme profile is searched once. PUL signatures are
evaluated independently: all required families must occur on the same contig,
and a contig is retained when it matches at least one resolvable PUL.
Unresolved rules are reported and skipped. The output `pul_matches.tsv.gz`
contains one row per matched contig-PUL pair. This remains CAZyme-signature
co-occurrence screening, not proof of a biological PUL or substrate activity.

To screen against all canonical CAZyme profiles in the installed dbCAN index:

```bash
phu screen -i contigs.fa --all-cazymes
```

This searches every AA, CBM, CE, GH, GT, and PL profile once. Ancillary
profiles such as SLH, cohesin, and dockerin are excluded, while subfamily
identities such as `GH5_2` are preserved. A contig is retained after at least
one threshold-passing family hit. The `cazyme_matches.tsv.gz` output records the
retained contig, family, protein, bitscore, E-value, and HMM coverage.

The default all-CAZyme outputs are compact and organized by evidence level:

```text
phu-screen/
├── screened_contigs.fasta
├── kept_contigs.txt
├── cazyme_matches.tsv.gz
├── cazyme_summary.tsv
├── cazyme_class_summary.tsv
├── evidence/
│   └── cazyme_hits.tsv.gz
└── .phu/
  └── run.json
```

`evidence/cazyme_hits.tsv.gz` contains one row per qualifying HMM/domain hit.
`cazyme_matches.tsv.gz` consolidates hits into one row per contig-family,
including unique protein counts, hit counts, protein IDs, and one deterministic
best supporting hit. The family summary reports detected-family occurrence
counts; the class summary reports unique contig and protein counts for the six
canonical classes. These are occurrence counts, not normalized abundance.
The outputs describe targeted HMM-based CAZyme-family evidence, not substrate
prediction or full multi-method `run_dbcan` consensus annotation.

For both `--all-cazymes` and `--all-puls`, per-family HMM files and domtblout
files are temporary search artifacts outside the output directory. They are
removed automatically after success or failure. The durable all-CAZyme result
contains compressed evidence, contig-family matches, family summaries, class
summaries, and `.phu/run.json`. The durable all-PUL result contains compressed
CAZyme evidence, PUL matches, per-family support, PUL summaries, substrate
association summaries, and `.phu/run.json`.

For all-PUL runs, CAZyme evidence is retained even when no complete PUL rule
matches. In that case the PUL tables are header-only and retained-contig
outputs are empty. Reference substrates describe the database PUL annotation;
they are not predictions for the query contigs.

The durable all-PUL output tree is:

```text
phu-screen/
├── screened_contigs.fasta
├── kept_contigs.txt
├── pul_matches.tsv.gz
├── pul_family_support.tsv.gz
├── pul_summary.tsv
├── substrate_summary.tsv
├── evidence/
│   └── cazyme_hits.tsv.gz
└── .phu/
  └── run.json
```

`pul_matches.tsv.gz` contains one row per complete contig-PUL match.
`pul_family_support.tsv.gz` records the supporting family evidence, while the
PUL and substrate summaries aggregate matched reference records. Coordinates
and reference substrates describe supporting database evidence; they are not
validated PUL boundaries or substrate predictions for query contigs.

The number of searched profiles depends on the installed dbCAN snapshot. This
is targeted HMM-based CAZyme-family screening, not full `run_dbcan` consensus
annotation using multiple independent methods, and it does not establish
biological activity or PUL identity.

## How it works

The screen command follows four main steps with several optional enhancements:

**First**, it predicts all possible proteins from your DNA contigs using **pyrodigal** (for standard microbial genes) or **pyrodigal-gv** (for viral genes when available). This step translates your DNA sequences into protein sequences, creating names like `contig123|gene1`, `contig123|gene2`, etc. The tool automatically handles complex contig names with multiple `|` separators.

After translation, proteins shorter than `--min-protein-len-aa` are discarded before HMM searching.

**Second**, it searches the predicted proteins against your HMM profiles using **pyHMMER**, a fast in-memory Python implementation that eliminates the need for external HMMER binaries. Each HMM file is processed efficiently with native Python threading, and results maintain compatibility with standard HMMER formats.

**Third**, it decides which contigs to keep based on the search results and your filtering criteria. This is where the "combine mode" logic becomes important if you're using multiple HMMs, distinguishing between HMMER "targets" (protein sequences) and HMM "models" (profiles that matched).

**Finally**, it extracts the matching contigs from your original file and saves them to the output. The tool can also extract target proteins per model and build custom HMM profiles from those proteins for future use.

## Cache handling

Protein prediction is cached and reused when the input contigs and prediction parameters are unchanged. For `phu screen`, changing `--min-protein-len-aa`, `--mode`, or `--ttable` rebuilds the cached proteins. Changing HMM files, combine mode, or output options does not.

See [cache.md](../cache.md) for the shared cache rules used by both `screen` and `jack`.

## PFAM accession handling

Positional arguments in `phu screen` can be mixed:

- Local HMM files (`capsid.hmm`)
- PFAM accessions (`PF00001`, `PF00589.17`)

When PFAM accessions are used, `phu` will:

1. Ensure a local `Pfam-A.hmm` database exists (download on first use).
2. Resolve each accession to the normalized accession (for example, `PF00589.17` -> `PF00589`).
3. Extract those models and run screening normally.

Tip: prepare PFAM ahead of time to avoid first-run setup delays:

```bash
phu dbs prepare pfam
```

Database location:

- `PHU_DB_FOLDER` if set.
- Otherwise `$XDG_DATA_HOME/phu/db` when `XDG_DATA_HOME` is set.
- Otherwise `~/.local/share/phu/db`.

## Using Multiple HMMs

When you provide multiple HMM files, you need to decide how strict you want to be about matches. There are three ways to combine results:

**Any mode** (the default) keeps contigs that match at least one model. Important detail: when a single contig matches multiple models, "any" preserves the best hit per model (rather than selecting only one overall best hit). As a result a contig that matches model A and model B will yield one protein for A and one protein for B (subject to `--top-per-contig`). This is useful when you want a representative protein per matched model from each contig.

**All mode** applies the current multi-model selection rule. It is intended for requiring evidence across a set of models, but currently has a known limitation when the batch produces no global hits. Do not interpret `all` as proof that every biologically expected marker is present; inspect per-model output and record the exact models queried.

**Threshold mode** lets you specify a minimum number of models that must match. This gives you flexibility between "any" and "all". You might require at least 3 out of 5 models to match, for example.

## Target Data Extraction and HMM Building

The tool offers powerful features for analyzing and reusing your screening results:

**Target Protein Extraction** (`--save-target-proteins`) saves the actual protein sequences that matched each HMM model, organized in separate files per model. These proteins come only from contigs that passed your final filtering criteria.

**Custom HMM Building** (`--save-target-hmms`) builds HMM files from matched target proteins and requires `--save-target-proteins`. For model sets containing multiple proteins, PHU runs the external `mafft` aligner before building the HMM. MAFFT is not installed by `uv` or bundled with PHU; install it in the environment or module used to run PHU, for example with `conda install -c bioconda mafft` or `pixi add mafft`. PHU fails clearly when MAFFT is unavailable. Single-protein outputs use a single-sequence HMM and are not a family alignment.

**Viral Mode Support** - When using viral gene prediction (if pyrodigal-gv is available), the tool is optimized for shorter, more compact viral genes and can handle overlapping gene structures common in viral genomes.

## Understanding Your Results

The main output is `screened_contigs.fasta`, which contains all the contigs that passed your filtering criteria. You'll also get `kept_contigs.txt` with just the names of these contigs.

If you used multiple HMMs, pay attention to how the combine mode affects your results. In "any" mode (see above), a contig that matches multiple models will produce one (best) protein per matched model. In "all" mode, kept contigs must have at least one hit for every model, and the tool returns one best hit per model for each kept contig — so the per-model protein counts will be balanced across models.

When you use the `--save-target-proteins` option, you'll get a folder called `target_proteins` with separate files for each model. Note the distinction:

- In "pure" HMM mode (default), each input HMM file is treated as one model and the output files are grouped by HMM filename (one file per input file).
- In "mixed" HMM mode (used for concatenated/pressed HMM files), a single HMM file can contain multiple model names; in that case `--save-target-proteins` will create one output file per model name found inside the domtblout.

All saved protein FASTA files contain only proteins that come from contigs that were kept in the final `screened_contigs.fasta`, and are de-duplicated per model file.

## Pass/fail behavior summary

The screening decision uses pyHMMER hit objects, not a second parsing pass over domtblout text files.

- Effective score is full-sequence score by default.
- For KOfam models, if `score_type=domain`, effective score switches to domain score.
- If `--use-kofam-thresholds` is enabled, KO thresholds from `ko_list` are applied per model.
- If both `--min-bitscore` and KO threshold exist, the stricter value is used (`max` of both).
- `--max-evalue` is applied to hit-level E-value.

See the full rule set and examples in [screen thresholds and decision logic](screen-thresholds.md).

**Target Data Outputs:**
- `target_proteins/{model}_proteins.mfa` - Proteins matching each model (if `--save-target-proteins`)
- `target_hmms/{model}.hmm` - HMMs built from target proteins after MAFFT alignment (if `--save-target-hmms`)

These outputs respect your HMM mode settings and combination logic, ensuring consistency between your screening criteria and extracted data.

## Command Options

```bash
Usage: phu screen [OPTIONS] HMMS...                                                        
                                                                                            
 Screen contigs for protein families using pyHMMER on predicted CDS.                          
                                                                                            
 Supports multiple HMM files with different combination modes:                              
 - any: Keep contigs matching any HMM (default, most permissive)                            
 - all: Keep contigs matching all HMMs (most restrictive)                                   
 - threshold: Keep contigs matching at least --min-hmm-hits HMMs                            
                                                                                            
 HMM modes:                                                                                 
 - pure: Each HMM file contains one model (default, most common)                            
 - mixed: HMM files contain multiple models (pressed/concatenated HMMs)                     
                                                                                            
 Examples:                                                                                  
  phu screen -i contigs.fa --combine-mode any "*.hmm"
   phu screen -i contigs.fa --combine-mode all file1.hmm file2.hmm file3.hmm
   phu screen -i contigs.fa --combine-mode threshold --min-hmm-hits 5 pfam_database.hmm
  phu screen -i contigs.fa --save-target-proteins "*.hmm"
  phu screen -i contigs.fa --save-target-hmms --save-target-proteins "*.hmm"
                                                                                            
╭─ Arguments ──────────────────────────────────────────────────────────────────────────────╮
│ *    hmms      HMMS...  HMM files (supports wildcards like *.hmm) [required]             │
╰──────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs     -i                        PATH                Input contigs FASTA │
│                                                                      [required]          │
│    --output-folder     -o                        PATH                Output directory    │
│                                                                      [default:           │
│                                                                      phu-screen]         │
│    --mode                                        TEXT                pyrodigal mode:     │
│                                                                      meta|single         │
│                                                                      [default: meta]     │
│    --threads           -t                        INTEGER RANGE       Threads for both    │
│                                                  [x>=1]              pyrodigal and       │
│                                                                      pyHMMER             │
│                                                                      [default: 1]        │
│    --min-bitscore                                FLOAT               Minimum bitscore to │
│                                                                      keep a domain hit   │
│    --max-evalue                                  FLOAT               Maximum independent │
│                                                                      E-value to keep a   │
│                                                                      domain hit          │
│                                                                      [default: 1e-05]    │
│    --top-per-contig                              INTEGER             Keep top-N hits per │
│                                                                      contig (by          │
│                                                                      bitscore)           │
│                                                                      [default: 1]        │
│    --min-protein-len-aa    -g                   INTEGER RANGE       Minimum translated   │
│                                                  [x>=1]              protein length to    │
│                                                                      keep (aa)           │
│                                                                      [default: 30]       │
│    --ttable                                      INTEGER             NCBI translation    │
│                                                                      table for coding    │
│                                                                      sequences           │
│                                                                      [default: 11]       │
│    --keep-proteins         --no-keep-proteins                        Keep the protein    │
│                                                                      FASTA used for      │
│                                                                      searching           │
│                                                                      [default:           │
│                                                                      no-keep-proteins]   │
│    --keep-domtbl           --no-keep-domtbl                          Keep raw domtblout  │
│                                                                      from pyHMMER        │
│                                                                      [default:           │
│                                                                      keep-domtbl]        │
│    --combine-mode                                TEXT                How to combine hits │
│                                                                      from multiple HMMs: │
│                                                                      any|all|threshold   │
│                                                                      [default: any]      │
│    --min-hmm-hits                                INTEGER             Minimum number of   │
│                                                                      HMMs that must hit  │
│                                                                      a contig (for       │
│                                                                      threshold mode)     │
│                                                                      [default: 1]        │
│    --save-target-pro…      --no-save-target-…                        Save matched        │
│                                                                      proteins per HMM    │
│                                                                      model in            │
│                                                                      target_proteins/    │
│                                                                      subfolder           │
│                                                                      [default:           │
│                                                                      no-save-target-pro… │
│    --save-target-hmms      --no-save-target-…                        Build and save HMMs │
│                                                                      from target         │
│                                                                      proteins in         │
│                                                                      target_hmms/        │
│                                                                      subfolder           │
│                                                                      [default:           │
│                                                                      no-save-target-hmm… │
│    --hmm-mode                                    TEXT                HMM file type:       │
│                                                                      'pure' (one model   │
│                                                                      per file) or         │
│                                                                      'mixed'             │
│                                                                      (pressed/concatena… │
│                                                                      HMMs)               │
│                                                                      [default: pure]     │
│    --help              -h                                            Show this message   │
│                                                                      and exit.           │
╰──────────────────────────────────────────────────────────────────────────────────────────╯

```

Use `--output-folder` to change where the results are saved. The default is a folder called `phu-screen` in your current directory.

Use `--threads` to speed things up if you have multiple CPU cores available. This affects both the protein prediction step and the HMMER searches.

Use `--min-protein-len-aa` to keep only proteins at or above a given amino-acid length before running HMM searches.

Use `--max-evalue` to make your searches more or less strict. The default is 1e-5, which is reasonably stringent. Lower values (like 1e-10) are more strict, while higher values (like 1e-3) are more permissive.

`--cut-ga` is enabled by default. pyHMMER applies model-specific GA gathering cutoffs embedded in HMM profiles (especially useful with PFAM models).

Use `--no-cut-ga` if you want to disable GA filtering and rely only on explicit score/E-value filters.

Use `--save-target-proteins` if you want to get the actual protein sequences from the contigs that matched each model. The saved proteins are taken only from contigs that passed final filtering and are grouped per-model (see "HMM modes" above).

Use `--save-target-hmms` together with `--save-target-proteins` to build custom HMMs. For multiple target proteins, PHU requires `mafft` on `PATH` and uses it to create a real multiple-sequence alignment. Install MAFFT separately in the runtime environment, for example `conda install -c bioconda mafft` or `pixi add mafft`. PHU does not vendor or install MAFFT through `uv`.

## Examples

Find contigs with any viral protein (default "any" preserves best-per-model hits):
```bash
phu screen --input-contigs assembly.fasta --combine-mode any viral_capsid.hmm viral_polymerase.hmm
```

Use PFAM accessions directly as positional targets:
```bash
phu screen --input-contigs assembly.fasta PF00001 PF00589
```

Use PFAM with GA gathering cutoffs (default behavior):
```bash
phu screen --input-contigs assembly.fasta PF00001 PF00589
```

Disable GA gathering cutoffs explicitly:
```bash
phu screen --input-contigs assembly.fasta PF00001 PF00589 --no-cut-ga
```

Mix local HMM files with PFAM accessions:
```bash
phu screen --input-contigs contigs.fa capsid.hmm PF00589 --combine-mode all
```

Find contigs that have complete viral genomes (all four proteins):
```bash
phu screen --input-contigs contigs.fa --combine-mode all capsid.hmm portal.hmm primase.hmm terminase.hmm
```

Use multiple threads and save the matching proteins (per-model grouping depends on `--hmm-mode`):
```bash
phu screen --input-contigs large_assembly.fasta --threads 16 --save-target-proteins marker.hmm
```

Be more strict about matches:
```bash
phu screen --input-contigs contigs.fa --max-evalue 1e-10 protein_family.hmm
```

Build custom HMMs from viral proteins after installing MAFFT:
```bash
phu screen --input-contigs viral_assembly.fasta --save-target-proteins --save-target-hmms --combine-mode all capsid.hmm polymerase.hmm
```

Complete viral screening workflow with custom HMM generation:
```bash
phu screen --input-contigs metagenome.fa --save-target-proteins --save-target-hmms --combine-mode threshold --min-hmm-hits 3 viral_marker1.hmm viral_marker2.hmm viral_marker3.hmm viral_marker4.hmm
```

Screen with a shell-safe input path:
```bash
phu screen --input-contigs 'scaffolds_with_complex|names|assembly.fa' protein_family.hmm
```

## What to expect

Runtime depends on input size, model count, available resources, and external tools. PHU does not promise a fixed runtime or scaling factor; benchmark representative data on your environment.

The output size depends on how many contigs match your criteria. In "any" mode, you might get quite a few contigs. In "all" mode, you'll typically get fewer but higher-confidence results.

The search engine is pyHMMER, while sequence extraction uses the external `seqkit` executable. Target-HMM generation additionally requires the external `mafft` executable for multi-sequence alignments.

If you don't get any results, try relaxing your E-value threshold or check that your HMM files are in the correct format. If you get too many results, try using "all" mode instead of "any" mode, or make your E-value threshold more strict.

## Requirements

You need to have **pyrodigal**, **pyHMMER**, and **seqkit** installed and available. **pyrodigal-gv** is optional but recommended for viral genome analysis. Your input contigs should be in FASTA format, and your HMM files should be in HMMER3 format.

**Workflow notes:**
- **No HMMER binary required**: pyHMMER performs the HMM search
- **External extraction**: `seqkit` must be available on `PATH`
- **Optional HMM generation**: `mafft` must be available on `PATH` for multi-sequence target sets
- **Explicit identifiers**: quote shell paths containing special characters

The command expects DNA sequences as input, not protein sequences. If you already have predicted proteins, you should use pyHMMER directly rather than this command.
