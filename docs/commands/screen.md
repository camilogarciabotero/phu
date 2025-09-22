# screen

## What does it do?

The `phu screen` command helps you find DNA contigs that contain specific protein families. It predicts proteins from your contigs, searches those proteins against Hidden Markov Model (HMM) profiles using **pyHMMER** (a fast Python implementation), and then selects contigs based on the matches and your combination/filtering rules. Think of it as a molecular search engine for pulling out contigs that contain the proteins you care about.

This is especially useful when you have metagenomic assemblies and want to pull out contigs that belong to viruses, or when you're looking for contigs that contain specific metabolic pathways. The tool now includes advanced features like building custom HMMs from your target proteins and specialized viral gene prediction.

## Synopsis

```bash
phu screen --input-contigs [INPUT_CONTIGS] [HMMS...]
```

**Example:**
```bash
phu screen --input-contigs your_contigs.fasta your_protein_family.hmm
```

This simple command will find all contigs in `your_contigs.fasta` that contain proteins matching `your_protein_family.hmm` and save them to a new file called `screened_contigs.fasta` in a folder named `phu-screen`.

## How it works

The screen command follows four main steps with several optional enhancements:

**First**, it predicts all possible proteins from your DNA contigs using **pyrodigal** (for standard microbial genes) or **pyrodigal-gv** (for viral genes when available). This step translates your DNA sequences into protein sequences, creating names like `contig123|gene1`, `contig123|gene2`, etc. The tool automatically handles complex contig names with multiple `|` separators.

**Second**, it searches the predicted proteins against your HMM profiles using **pyHMMER**, a fast in-memory Python implementation that eliminates the need for external HMMER binaries. Each HMM file is processed efficiently with native Python threading, and results maintain compatibility with standard HMMER formats.

**Third**, it decides which contigs to keep based on the search results and your filtering criteria. This is where the "combine mode" logic becomes important if you're using multiple HMMs, distinguishing between HMMER "targets" (protein sequences) and HMM "models" (profiles that matched).

**Finally**, it extracts the matching contigs from your original file and saves them to the output. The tool can also extract target proteins per model and build custom HMM profiles from those proteins for future use.

## Using Multiple HMMs

When you provide multiple HMM files, you need to decide how strict you want to be about matches. There are three ways to combine results:

**Any mode** (the default) keeps contigs that match at least one model. Important detail: when a single contig matches multiple models, "any" preserves the best hit per model (rather than selecting only one overall best hit). As a result a contig that matches model A and model B will yield one protein for A and one protein for B (subject to `--top-per-contig`). This is useful when you want a representative protein per matched model from each contig.

**All mode** only keeps contigs that match every single HMM you provided. This is very strict and useful when you need complete sets of proteins. For instance, if you're looking for complete viral genomes that must have all four proteins (capsid, portal, primase, and terminase), you would use "all" mode.

**Threshold mode** lets you specify a minimum number of models that must match. This gives you flexibility between "any" and "all". You might require at least 3 out of 5 models to match, for example.

## Target Data Extraction and HMM Building

The tool offers powerful features for analyzing and reusing your screening results:

**Target Protein Extraction** (`--save-target-proteins`) saves the actual protein sequences that matched each HMM model, organized in separate files per model. These proteins come only from contigs that passed your final filtering criteria.

**Custom HMM Building** (`--save-target-hmms`) automatically builds new HMM profiles from your target proteins. This works independently of protein saving - if you only want HMMs, the tool will extract proteins temporarily, build the HMMs, then clean up. For single sequences, it builds individual HMMs; for multiple sequences, it creates simple alignments by padding to equal length before building consensus models.

**Viral Mode Support** - When using viral gene prediction (if pyrodigal-gv is available), the tool is optimized for shorter, more compact viral genes and can handle overlapping gene structures common in viral genomes.

## Understanding Your Results

The main output is `screened_contigs.fasta`, which contains all the contigs that passed your filtering criteria. You'll also get `kept_contigs.txt` with just the names of these contigs.

If you used multiple HMMs, pay attention to how the combine mode affects your results. In "any" mode (see above), a contig that matches multiple models will produce one (best) protein per matched model. In "all" mode, kept contigs must have at least one hit for every model, and the tool returns one best hit per model for each kept contig — so the per-model protein counts will be balanced across models.

When you use the `--save-target-proteins` option, you'll get a folder called `target_proteins` with separate files for each model. Note the distinction:

- In "pure" HMM mode (default), each input HMM file is treated as one model and the output files are grouped by HMM filename (one file per input file).
- In "mixed" HMM mode (used for concatenated/pressed HMM files), a single HMM file can contain multiple model names; in that case `--save-target-proteins` will create one output file per model name found inside the domtblout.

All saved protein FASTA files contain only proteins that come from contigs that were kept in the final `screened_contigs.fasta`, and are de-duplicated per model file.

**New Target Data Outputs:**
- `target_proteins/{model}_proteins.mfa` - Proteins matching each model (if `--save-target-proteins`)
- `target_hmms/{model}.hmm` - Custom HMMs built from target proteins (if `--save-target-hmms`)

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
   phu screen -i contigs.fa --combine-mode any *.hmm
   phu screen -i contigs.fa --combine-mode all file1.hmm file2.hmm file3.hmm
   phu screen -i contigs.fa --combine-mode threshold --min-hmm-hits 5 pfam_database.hmm
   phu screen -i contigs.fa --save-target-proteins *.hmm
   phu screen -i contigs.fa --save-target-hmms *.hmm
                                                                                            
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
│    --min-gene-len                                INTEGER             Minimum gene length │
│                                                                      for pyrodigal (nt)  │
│                                                                      [default: 90]       │
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

Use `--max-evalue` to make your searches more or less strict. The default is 1e-5, which is reasonably stringent. Lower values (like 1e-10) are more strict, while higher values (like 1e-3) are more permissive.

Use `--save-target-proteins` if you want to get the actual protein sequences from the contigs that matched each model. The saved proteins are taken only from contigs that passed final filtering and are grouped per-model (see "HMM modes" above).

Use `--save-target-hmms` to build custom HMM profiles from your target proteins. This works independently of `--save-target-proteins` - the tool can extract proteins temporarily just for HMM building if needed.

## Examples

Find contigs with any viral protein (default "any" preserves best-per-model hits):
```bash
phu screen --input-contigs assembly.fasta --combine-mode any viral_capsid.hmm viral_polymerase.hmm
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

Build custom HMMs from viral proteins (works with or without saving proteins):
```bash
phu screen --input-contigs viral_assembly.fasta --save-target-hmms --combine-mode all capsid.hmm polymerase.hmm
```

Complete viral screening workflow with custom HMM generation:
```bash
phu screen --input-contigs metagenome.fa --save-target-proteins --save-target-hmms --combine-mode threshold --min-hmm-hits 3 viral_marker1.hmm viral_marker2.hmm viral_marker3.hmm viral_marker4.hmm
```

Screen with complex contig names (automatically handled):
```bash
phu screen --input-contigs scaffolds_with_complex|names|assembly.fa protein_family.hmm
```

## What to expect

Gene prediction usually takes 1-2 minutes per million base pairs of input. The pyHMMER searches are significantly faster than traditional HMMER due to in-memory processing and native Python threading. Performance scales well with the `--threads` option.

The output size depends on how many contigs match your criteria. In "any" mode, you might get quite a few contigs. In "all" mode, you'll typically get fewer but higher-confidence results.

**New performance benefits:**
- **Faster execution**: pyHMMER eliminates subprocess overhead and file I/O
- **Better memory usage**: In-memory processing of all data
- **No binary dependencies**: Pure Python implementation
- **Robust handling**: Automatically deals with complex contig naming schemes

If you don't get any results, try relaxing your E-value threshold or check that your HMM files are in the correct format. If you get too many results, try using "all" mode instead of "any" mode, or make your E-value threshold more strict.

## Requirements

You need to have **pyrodigal**, **pyHMMER**, and **seqkit** installed and available. **pyrodigal-gv** is optional but recommended for viral genome analysis. Your input contigs should be in FASTA format, and your HMM files should be in HMMER3 format.

**Key improvements:**
- **No HMMER binary required**: pyHMMER provides a pure Python implementation
- **Automatic viral support**: Uses pyrodigal-gv when available for viral gene prediction
- **Robust contig handling**: Automatically handles complex contig naming with multiple separators
- **Flexible HMM building**: Create custom profiles with or without saving intermediate proteins

The command expects DNA sequences as input, not protein sequences. If you already have predicted proteins, you should use pyHMMER directly rather than this command.
