# screen

## What does it do?

The `phu screen` command helps you find DNA contigs that contain specific protein families. It predicts proteins from your contigs, searches those proteins against Hidden Markov Model (HMM) profiles using **pyHMMER** (a fast Python implementation), and then selects contigs based on the matches and your combination/filtering rules. Think of it as a molecular search engine for pulling out contigs that contain the proteins you care about.

This is especially useful when you have metagenomic assemblies and want to pull out contigs that belong to viruses, or when you're looking for contigs that contain specific metabolic pathways. Target-protein extraction and custom HMM generation are experimental features and are not a validated replacement for a biologically aligned profile-HMM workflow.

For implementation-level details on pass/fail decisions, see [screen thresholds and decision logic](../screen-thresholds.md).

## Synopsis

```bash
phu screen --input-contigs [INPUT_CONTIGS] [HMMS...]
```

`[HMMS...]` accepts local HMM paths, PFAM/KOfam/PHROG identifiers, ECOD groups (`ecod:<id>`), dbCAN families
(for example `GH128` or `CBM89`), and one dbCAN PUL identifier such as
`PUL0621`.

**Example:**
```bash
phu screen --input-contigs your_contigs.fasta your_protein_family.hmm
```

This simple command will find all contigs in `your_contigs.fasta` that contain proteins matching `your_protein_family.hmm` and save them to a new file called `matched-contigs.fna` in a folder named `phu-screen`.

If an argument looks like a PFAM accession (`PF` + 5 digits, optionally with a version), `phu` resolves it from a local Pfam-A database automatically.

dbCAN family inputs are case-insensitive and are validated against the local
dbCAN HMM index. A PUL input expands to its ordered CAZyme-family signature
from `cazymes_predicted_dbcan`; by default every resolved family must match on
the same contig. Use `--combine-mode threshold` explicitly for partial
signatures. Only one PUL can be queried per run, and PUL queries cannot be
mixed with other database or local-HMM inputs.

This is PUL CAZyme-signature screening. Family co-occurrence alone does not
identify a biological PUL or predict its substrate. Each run writes
`.phu/query.json` with the original query, expansion and matching rule. When
the query names Pfam, KOfam or dbCAN models, it also records the `phu-db`
bundle version and the SHA-256 of each database file used.

To screen against every resolvable PUL signature:

```bash
phu screen -i contigs.fa --all-puls
```

Each unique required CAZyme profile is searched once. PUL signatures are
evaluated independently: all required families must occur on the same contig,
and a contig is retained when it matches at least one resolvable PUL.
Unresolved rules are reported and skipped. The output `pul-matches.tsv.gz`
contains one row per matched contig-PUL pair. This remains CAZyme-signature
co-occurrence screening, not proof of a biological PUL or substrate activity.

To screen against all canonical CAZyme profiles in the installed dbCAN index:

```bash
phu screen -i contigs.fa --all-cazymes
```

This searches every AA, CBM, CE, GH, GT, and PL profile once. Ancillary
profiles such as SLH, cohesin, and dockerin are excluded, while subfamily
identities such as `GH5_2` are preserved. A contig is retained after at least
one threshold-passing family hit. The `cazyme-matches.tsv.gz` output records the
retained contig, family, protein, bitscore, E-value, and HMM coverage.

The default all-CAZyme outputs are compact and organized by evidence level:

```text
phu-screen/
├── matched-contigs.fna
├── matched-contig-ids.txt
├── cazyme-matches.tsv.gz
├── cazyme-family-summary.tsv
├── cazyme-class-summary.tsv
├── evidence/
│   └── cazyme-hits.tsv.gz
└── .phu/
  └── run.json
```

`evidence/cazyme-hits.tsv.gz` contains one row per qualifying HMM/domain hit.
`cazyme-matches.tsv.gz` consolidates hits into one row per contig-family,
including unique protein counts, hit counts, protein IDs, and one deterministic
best supporting hit. The family summary reports detected-family occurrence
counts; the class summary reports unique contig and protein counts for the six
canonical classes. These are occurrence counts, not normalized abundance.
The outputs describe targeted HMM-based CAZyme-family evidence, not substrate
prediction or full multi-method `run_dbcan` consensus annotation.

For both `--all-cazymes` and `--all-puls`, per-family HMM files and hmmsearch
domain tables are temporary search artifacts outside the output directory. They are
removed automatically after success or failure. The durable all-CAZyme result
contains compressed evidence, contig-family matches, family summaries, class
summaries, and `.phu/run.json`. The durable all-PUL result contains compressed
CAZyme evidence, PUL matches, per-family support, PUL summaries, substrate
association summaries, and `.phu/run.json`.

For all-PUL runs, CAZyme evidence is retained even when no complete PUL rule
matches. In that case the PUL tables are header-only and retained-contig
outputs are empty. Reference substrates describe the database PUL annotation;
they are not predictions for the query contigs.

PUL rules that reference families unavailable in the installed dbCAN snapshot
are reported as unresolved and skipped. The unresolved-rule count therefore
describes database coverage, not failed HMM searches or evidence that the query
contigs lack those PULs. Rebuild or refresh the dbCAN database when a newer
snapshot is required.

The durable all-PUL output tree is:

```text
phu-screen/
├── matched-contigs.fna
├── matched-contig-ids.txt
├── pul-matches.tsv.gz
├── pul-family-support.tsv.gz
├── pul-summary.tsv
├── pul-substrate-summary.tsv
├── evidence/
│   └── cazyme-hits.tsv.gz
└── .phu/
  └── run.json
```

`pul-matches.tsv.gz` contains one row per complete contig-PUL match.
`pul-family-support.tsv.gz` records the supporting family evidence, while the
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

1. Use the installed Pfam database, installing the database bundle first if this is the first database-backed run.
2. Resolve each accession to the normalized accession (for example, `PF00589.17` -> `PF00589`).
3. Extract those models and run screening normally.

The first database-backed run installs the database bundle (about 2.2 GB
download). For install location, offline use and disk needs, see
[Databases and Provenance](../databases.md).

## Using Multiple HMMs

When you provide multiple HMM files, you need to decide how strict you want to be about matches. There are three ways to combine results:

**Any mode** (the default) keeps contigs that match at least one model. Important detail: when a single contig matches multiple models, "any" preserves the best hit per model (rather than selecting only one overall best hit). As a result a contig that matches model A and model B will yield one protein for A and one protein for B (subject to `--top-per-contig`). This is useful when you want a representative protein per matched model from each contig.

**All mode** applies the current multi-model selection rule. It is intended for requiring evidence across a set of models, but currently has a known limitation when the batch produces no global hits. Do not interpret `all` as proof that every biologically expected marker is present; inspect per-model output and record the exact models queried.

**Threshold mode** lets you specify a minimum number of models that must match. This gives you flexibility between "any" and "all". You might require at least 3 out of 5 models to match, for example.

## Target Data Extraction and HMM Building

The tool offers powerful features for analyzing and reusing your screening results:

**Target Protein Extraction** (`--save-target-proteins`) saves the actual protein sequences that matched each HMM model, organized in separate files per model. These proteins come only from contigs that passed your final filtering criteria.

**Custom HMM Building** (`--save-target-hmms`) builds HMM files from matched target proteins and requires `--save-target-proteins`. For model sets containing multiple proteins, PHU aligns them with FAMSA before building the HMM. FAMSA runs inside PHU through the `pyfamsa` package, which is installed with PHU, so no aligner program is needed. Single-protein outputs use a single-sequence HMM and are not a family alignment.

**Viral Mode Support** - When using viral gene prediction (if pyrodigal-gv is available), the tool is optimized for shorter, more compact viral genes and can handle overlapping gene structures common in viral genomes.

## Understanding Your Results

The main output is `matched-contigs.fna`, which contains all the contigs that passed your filtering criteria. You'll also get `matched-contig-ids.txt` with just the names of these contigs.

If you used multiple HMMs, pay attention to how the combine mode affects your results. In "any" mode (see above), a contig that matches multiple models will produce one (best) protein per matched model. In "all" mode, kept contigs must have at least one hit for every model, and the tool returns one best hit per model for each kept contig — so the per-model protein counts will be balanced across models.

When you use the `--save-target-proteins` option, you'll get a folder called `target-proteins` with separate files for each model. Note the distinction:

- In "pure" HMM mode (default), each input HMM file is treated as one model and the output files are grouped by HMM filename (one file per input file).
- In "mixed" HMM mode (used for concatenated/pressed HMM files), a single HMM file can contain multiple model names; in that case `--save-target-proteins` will create one output file per model name found inside the HMM file.

All saved protein FASTA files contain only proteins that come from contigs that were kept in the final `matched-contigs.fna`, and are de-duplicated per model file.

## Output files

| File or folder | Written when | Content |
| --- | --- | --- |
| `matched-contigs.fna` | always | Nucleotide FASTA of the contigs that passed the screen |
| `matched-contig-ids.txt` | always | IDs of those contigs, one per line |
| `{hmm}-hmmsearch-domains.tbl` | `--keep-domtbl` (default) | hmmsearch per-domain table, one per input HMM file |
| `.phu/query.json` | always, except `--all-cazymes` / `--all-puls` | Original query, expansion and matching rule; for Pfam, KOfam, PHROG, ECOD and dbCAN queries also the `phu-db` bundle version (for PHROG the function and category of each family; for ECOD what each `ecod:` id expanded to and the names of each F-group) and database file hashes |
| `.phu/run.json` | `--all-cazymes` / `--all-puls` | Run manifest with counts, thresholds and `generated_files` |
| `query-hmms/pfam/`, `query-hmms/kofam/`, `query-hmms/dbcan/` | Pfam, KO, dbCAN-family or PUL queries | The HMM files resolved from the local databases |
| `predicted-proteins.faa` | `--keep-proteins` | Protein FASTA of all predicted proteins |
| `.phu/prediction.json` | `--keep-proteins` | Protein-prediction cache key and cache location |
| `target-proteins/{model}.faa` | `--save-target-proteins` | Matched proteins of each model, from kept contigs only |
| `target-hmms/{model}.hmm` | `--save-target-hmms` | HMM built from each model's target proteins |
| `cazyme-matches.tsv.gz`, `cazyme-family-summary.tsv`, `cazyme-class-summary.tsv` | `--all-cazymes` | Contig-family matches and per-family / per-class summaries |
| `pul-matches.tsv.gz`, `pul-family-support.tsv.gz`, `pul-summary.tsv`, `pul-substrate-summary.tsv` | `--all-puls` | PUL signature matches, per-family support and summaries |
| `evidence/cazyme-hits.tsv.gz` | `--all-cazymes` / `--all-puls` | Every qualifying CAZyme HMM hit |

`{hmm}` is the HMM file name without its extension, and in pure mode it is also
the model ID. If two input files share a file name (`a/foo.hmm` and
`b/foo.hmm`), `phu` tells them apart by the `NAME` field inside each HMM: that
name becomes the model ID and names their domain tables and target files. If
both the file name and the internal `NAME` are identical, the files cannot be
told apart: `phu` prints a warning, pools their hits under one model ID and
writes them to one domain table.

Model names are made safe for file names by replacing unusual characters with
`_`. If two different models would end up with the same file name (for example
`beta/1` and `beta_1`), `phu` refuses the query rather than mix their results;
rename one of them.

With `--save-target-hmms`, HMMs are built only from the protein sets written
by the current run; files left in `target-proteins/` by an earlier run are not
read. The alignments used to build them are temporary and are not written into
the output folder.

**Alignment trimming** (`--trim-target-msa`, off by default) removes columns
from each alignment with trimAl (through `pytrimal`, installed with PHU)
before the HMM is built. It needs `--save-target-hmms`.

| Value | What is removed |
| --- | --- |
| `none` | Nothing (default) |
| `gappyout` | Columns judged too gappy from the alignment's own gap distribution; the mildest choice |
| `strict` | Gappy columns and columns with low residue similarity |
| `strictplus` | As `strict`, with a tighter similarity criterion; the most aggressive |

Trimming shortens the model, most of all for protein sets with very uneven
lengths, and it lowers hit bit scores, so cutoffs chosen for untrimmed models
do not carry over. The run prints the column count before and after, for
example `(gappyout: 1523 -> 377 columns)`. If trimming leaves no columns, or
fewer than two sequences, that model is built from the untrimmed alignment and
a warning is printed.

## Pass/fail behavior summary

The screening decision uses pyHMMER hit objects, not a second parsing pass over the domain table text files.

- Effective score is full-sequence score by default.
- For KOfam models, if `score_type=domain`, effective score switches to domain score.
- If `--use-kofam-thresholds` is enabled, KO thresholds from `ko_list` are applied per model.
- If both `--min-bitscore` and KO threshold exist, the stricter value is used (`max` of both).
- `--max-evalue` is applied to hit-level E-value.

See the full rule set and examples in [screen thresholds and decision logic](../screen-thresholds.md).

**Target Data Outputs:**
- `target-proteins/{model}.faa` - Proteins matching each model (if `--save-target-proteins`)
- `target-hmms/{model}.hmm` - HMMs built from target proteins after FAMSA alignment (if `--save-target-hmms`)

These outputs respect your HMM mode settings and combination logic, ensuring consistency between your screening criteria and extracted data.

## Command Options

```bash
                                                                                
 Usage: phu screen [OPTIONS] [QUERY]...                                         
                                                                                
 Find contigs that carry given protein families.                                
                                                                                
 Proteins are predicted from the contigs and searched with each QUERY: your     
 own HMM files, or families from the phu databases by name: Pfam (PF00001),     
 KOfam (K00001), PHROG (phrog_216; find them with `phu phrogs`), ECOD           
 (ecod:2004.1.1.5 for one family, or a T, H or X group such as ecod:2004.1.1    
 for all its families), dbCAN (GH128) or one dbCAN PUL (PUL0621). Each          
 database's own cutoffs apply.                                                  
                                                                                
 Examples:                                                                      
     phu screen -i contigs.fa PF00001 PF00589                                   
     phu screen -i contigs.fa my-models/*.hmm                                   
     phu screen -i contigs.fa --combine-mode all K00001 K00002                  
     phu screen -i contigs.fa --all-cazymes                                     
                                                                                
 Full guide: https://camilogarciabotero.github.io/phu/commands/screen/          
                                                                                
╭─ Query ──────────────────────────────────────────────────────────────────────╮
│   [QUERY]...      <path>  HMM files, or database names: PF00001, K00001,     │
│                           phrog_216, ecod:2004.1.1.5, GH128, PUL0621.        │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --input-contigs  -i      <path>  Contigs to screen (FASTA). [required]    │
│    --help           -h              Show this message and exit.              │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Other ways to query ────────────────────────────────────────────────────────╮
│ --all-cazymes                 Search all dbCAN CAZyme families instead of    │
│                               QUERY.                                         │
│ --all-puls                    Search for every dbCAN PUL signature instead   │
│                               of QUERY.                                      │
│ --hmm-mode     -M      <str>  'pure': one model per HMM file; 'mixed': files │
│                               with many models.                              │
│                               [default: pure]                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Hit thresholds ─────────────────────────────────────────────────────────────╮
│ --max-evalue      -e                      <float>          Highest E-value a │
│                                                            hit may have.     │
│                                                            [default: 1e-06]  │
│ --min-bitscore    -b                      <float>          Lowest bit score  │
│                                                            a hit may have.   │
│ --min-hmm-cover…                          <float range>    Share of your own │
│                                           [0.0<=x<=1.0]    HMM a hit must    │
│                                                            cover (0 = off).  │
│                                                            [default: 0.35]   │
│ --min-protein-c…                          <float range>    Share of the      │
│                                           [0.0<=x<=1.0]    protein a hit     │
│                                                            must cover (off   │
│                                                            by default).      │
│ --cut-ga              --no-cut-ga                          Use Pfam's        │
│                                                            curated cutoffs   │
│                                                            (GA).             │
│                                                            [default: cut-ga] │
│ --use-kofam-thr…      --no-use-kofam-…                     Use KOfam's       │
│                                                            per-KO cutoffs.   │
│                                                            [default:         │
│                                                            use-kofam-thresh… │
│ --use-phrog-thr…      --no-use-phrog-…                     For PHROGs, need  │
│                                                            60 bits and half  │
│                                                            of model and      │
│                                                            protein aligned.  │
│                                                            [default:         │
│                                                            use-phrog-thresh… │
│ --use-ecod-thre…      --no-use-ecod-t…                     For ECOD          │
│                                                            families, need    │
│                                                            half of the model │
│                                                            aligned.          │
│                                                            [default:         │
│                                                            use-ecod-thresho… │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Which contigs are kept ─────────────────────────────────────────────────────╮
│ --combine-mode    -c      <str>  Keep contigs hit by 'any' query, 'all' of   │
│                                  them, or a 'threshold'.                     │
│                                  [default: any]                              │
│ --min-hmm-hits    -k      <int>  With --combine-mode threshold: how many     │
│                                  queries must hit.                           │
│                                  [default: 1]                                │
│ --top-per-contig  -n      <int>  Best hits to report per contig.             │
│                                  [default: 1]                                │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Output ─────────────────────────────────────────────────────────────────────╮
│ --output-folder       -o                         <path>  Output folder.      │
│                                                          [default:           │
│                                                          phu-screen]         │
│ --keep-proteins           --no-keep-proteins             Keep the predicted  │
│                                                          proteins.           │
│                                                          [default:           │
│                                                          no-keep-proteins]   │
│ --keep-domtbl             --no-keep-domtbl               Keep the hmmsearch  │
│                                                          domain tables.      │
│                                                          [default:           │
│                                                          keep-domtbl]        │
│ --save-target-prote…      --no-save-target-p…            Save the matched    │
│                                                          proteins of each    │
│                                                          model               │
│                                                          (target-proteins/). │
│                                                          [default:           │
│                                                          no-save-target-pro… │
│ --save-target-hmms        --no-save-target-h…            Build an HMM from   │
│                                                          each model's        │
│                                                          matched proteins    │
│                                                          (target-hmms/).     │
│                                                          [default:           │
│                                                          no-save-target-hmm… │
│ --trim-target-msa                                <str>   Trim alignments     │
│                                                          before building     │
│                                                          those HMMs: none,   │
│                                                          gappyout, strict,   │
│                                                          strictplus.         │
│                                                          [default: none]     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Gene prediction ────────────────────────────────────────────────────────────╮
│ --mode                -m      <str>               'meta' for mixed contigs,  │
│                                                   'single' for one genome.   │
│                                                   [default: meta]            │
│ --min-protein-len-aa  -g      <int range> [x>=1]  Shortest protein to keep   │
│                                                   (aa).                      │
│                                                   [default: 30]              │
│ --ttable              -T      <int range> [x>=1]  Genetic code (NCBI table); │
│                                                   by default predicted per   │
│                                                   contig.                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Runtime ────────────────────────────────────────────────────────────────────╮
│ --threads  -t      <int range> [x>=1]  Threads. [default: 1]                 │
│ --quiet                                Print less.                           │
│ --verbose                              Print each search step.               │
╰──────────────────────────────────────────────────────────────────────────────╯
```

Use `--verbose` to record how each tool was run: the effective parameters of the in-process gene prediction, `hmmsearch`, FAMSA alignment and sequence extraction steps, are printed on lines that start with `[cmd]`. Use `--quiet` to hide the timed step lines. See [Run logs](../troubleshooting.md#run-logs) for the log format.

Use `--output-folder` to change where the results are saved. The default is a folder called `phu-screen` in your current directory.

Use `--threads` to speed things up if you have multiple CPU cores available. This affects both the protein prediction step and the HMMER searches.

Use `--min-protein-len-aa` to keep only proteins at or above a given amino-acid length before running HMM searches.

Use `--max-evalue` to make your searches more or less strict. The default is 1e-6, which is reasonably stringent, and it applies to every kind of model: your own HMMs, Pfam, KOfam and dbCAN. Lower values (like 1e-10) are more strict, while higher values (like 1e-3) are more permissive.

Use `--min-hmm-coverage` to control how much of one of your own HMMs a hit must cover. The default is 0.35: a hit is kept only when its domains align to more than 35% of the model (every reported domain counts, overlaps once), which removes matches to a small fragment of the family. It applies only to HMMs you supply; Pfam and KOfam models keep their curated score thresholds, and dbCAN models always use the same 0.35 rule. Set it to 0 to turn it off. The coverage of every hit is written to the result tables as `hmm_coverage`.

Use `--min-protein-coverage` to control how much of the protein a hit must cover. It is off by default. With `--min-protein-coverage 0.5`, a hit is kept only when its domains align to at least half of the protein, which rejects a small domain found inside a long protein. Unlike `--min-hmm-coverage`, it applies to every kind of model, including Pfam, KOfam and dbCAN, so expect it to remove single-domain matches in multi-domain proteins.

PHROG families are queried by identifier, for example `phu screen phrog_216 -i contigs.fna`, and can be mixed with other identifiers and your own HMMs. Models exist for `phrog_1` to `phrog_38880`. PHROG has no curated per-model cutoffs, so a PHROG hit must score at least 60 bits and align to at least half of both the model and the protein, the rule `phu avger` uses. `--no-use-phrog-thresholds` drops that rule and leaves only `--max-evalue` and `--min-bitscore`. The function and category of each queried PHROG are printed and recorded in `.phu/query.json`. To pick PHROGs by name or category instead of number, use [`phu phrogs`](phrogs.md), for example `phu screen -i contigs.fna $(phu phrogs resolve alias:mcp --ids-only)`.

ECOD families are queried as `ecod:<id>` at any level of the ECOD hierarchy. A four-part F-group id such as `ecod:2004.1.1.5` is one family; a shorter id stands for every F-group inside it: `ecod:2004.1.1` (a T-group), `ecod:2004.1` (an H-group, a homology group) or `ecod:2004` (an X-group). The bundle has one model per F-group (32,041 in ECOD v295.2), built from ECOD's F70 representatives; names are in `families.tsv` of the installed `ecod` database. Hits report models as `ecod:<f_id>`. ECOD models are single domains, often of longer proteins, and have no curated cutoffs, so an ECOD hit must align to at least half of the model; protein coverage is not required. `--no-use-ecod-thresholds` drops that rule. What each `ecod:` id expanded to, and the F-, T-, H- and X-group names and Pfam accessions of every family, are recorded in `.phu/query.json`.

`--cut-ga` is enabled by default. pyHMMER applies model-specific GA gathering cutoffs embedded in HMM profiles (especially useful with PFAM models).

Use `--no-cut-ga` if you want to disable GA filtering and rely only on explicit score/E-value filters.

Use `--save-target-proteins` if you want to get the actual protein sequences from the contigs that matched each model. The saved proteins are taken only from contigs that passed final filtering and are grouped per-model (see "HMM modes" above).

Use `--save-target-hmms` together with `--save-target-proteins` to build custom HMMs. For multiple target proteins, PHU creates a real multiple-sequence alignment with FAMSA (through `pyfamsa`, installed with PHU). No separate aligner has to be installed.

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

Build custom HMMs from viral proteins (add `--trim-target-msa gappyout` to trim each alignment first):
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

The search engine is pyHMMER, and matched contigs and proteins are extracted in-process; `phu screen` runs no external program. Target-HMM generation aligns multi-sequence sets in-process with FAMSA (`pyfamsa`).

If you don't get any results, try relaxing your E-value threshold or check that your HMM files are in the correct format. If you get too many results, try using "all" mode instead of "any" mode, or make your E-value threshold more strict.

## Requirements

You need to have **pyrodigal** and **pyHMMER** installed (both are installed with PHU). **pyrodigal-gv** is optional but recommended for viral genome analysis. Your input contigs should be in FASTA format, and your HMM files should be in HMMER3 format.

**Workflow notes:**
- **No HMMER binary required**: pyHMMER performs the HMM search
- **No seqkit required**: matched contigs and proteins are extracted in-process
- **Optional HMM generation**: no extra program; multi-sequence target sets are aligned with FAMSA through `pyfamsa`
- **Explicit identifiers**: quote shell paths containing special characters

The command expects DNA sequences as input, not protein sequences. If you already have predicted proteins, you should use pyHMMER directly rather than this command.
