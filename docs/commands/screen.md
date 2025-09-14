# Screen Command

## What does it do?

The `phu screen` command helps you find DNA contigs that contain specific protein families. It works by predicting all the proteins in your contigs, then searching those proteins against Hidden Markov Model (HMM) profiles to find matches. Think of it as a molecular search engine for finding contigs with the proteins you care about.

This is especially useful when you have metagenomic assemblies and want to pull out contigs that belong to viruses, or when you're looking for contigs that contain specific metabolic pathways.

## Basic Usage

```bash
phu screen --input-contigs your_contigs.fasta your_protein_family.hmm
```

This simple command will find all contigs in `your_contigs.fasta` that contain proteins matching `your_protein_family.hmm` and save them to a new file called `screened_contigs.fasta` in a folder named `phu-screen`.

## How it works

The screen command follows four main steps:

First, it predicts all possible proteins from your DNA contigs using a tool called `pyrodigal`. This step also translates your DNA sequences into protein sequences, creating names like `contig123|gene1`, `contig123|gene2`, and so on.

Second, it searches these predicted proteins against your HMM profiles using HMMER. Each HMM file gets searched separately, and the results are saved in individual files.

Third, it decides which contigs to keep based on the search results and your filtering criteria. This is where the "combine mode" logic becomes important if you're using multiple HMMs.

Finally, it extracts the matching contigs from your original file and saves them to the output. If you've requested it, it also saves the matching proteins organized by which HMM they matched.

## Using Multiple HMMs

When you provide multiple HMM files, you need to decide how strict you want to be about matches. There are three ways to combine results:

**Any mode** (the default) keeps contigs that match at least one of your HMMs. This is the most permissive option and good for broad searches. For example, if you're looking for any viral marker, you might use several different viral HMMs and keep contigs that match any of them.

**All mode** only keeps contigs that match every single HMM you provided. This is very strict and useful when you need complete sets of proteins. For instance, if you're looking for complete viral genomes that must have all four proteins (capsid, portal, primase, and terminase), you would use "all" mode.

**Threshold mode** lets you specify a minimum number of HMMs that must match. This gives you flexibility between "any" and "all". You might require at least 3 out of 5 HMMs to match, for example.

## Understanding Your Results

The main output is `screened_contigs.fasta`, which contains all the contigs that passed your filtering criteria. You'll also get `kept_contigs.txt` with just the names of these contigs.

If you used multiple HMMs, pay attention to how the combine mode affects your results. In "any" mode, you might get different numbers of proteins for each HMM because some contigs only matched some of the HMMs. In "all" mode, you'll get exactly the same number of proteins for each HMM because every kept contig had to match all of them.

When you use the `--save-target-proteins` option, you'll get a folder called `target_proteins` with separate files for each HMM. These contain only the proteins from contigs that made it into your final results, not all the proteins that matched during the search.

## Common Options

Use `--outdir` to change where the results are saved. The default is a folder called `phu-screen` in your current directory.

Use `--threads` to speed things up if you have multiple CPU cores available. This affects both the protein prediction step and the HMMER searches.

Use `--max-evalue` to make your searches more or less strict. The default is 1e-5, which is reasonably stringent. Lower values (like 1e-10) are more strict, while higher values (like 1e-3) are more permissive.

Use `--save-target-proteins` if you want to get the actual protein sequences that matched each HMM, not just the contigs.

## Examples

Find contigs with any viral protein:
```bash
phu screen --input-contigs assembly.fasta viral_capsid.hmm viral_polymerase.hmm
```

Find contigs that have complete viral genomes (all four proteins):
```bash
phu screen --input-contigs contigs.fa capsid.hmm portal.hmm primase.hmm terminase.hmm --combine-mode all
```

Use multiple threads and save the matching proteins:
```bash
phu screen --input-contigs large_assembly.fasta marker.hmm --threads 16 --save-target-proteins
```

Be more strict about matches:
```bash
phu screen --input-contigs contigs.fa protein_family.hmm --max-evalue 1e-10
```

## What to expect

Gene prediction usually takes 1-2 minutes per million base pairs of input. The HMMER searches take longer and depend on how many proteins were predicted and how many HMMs you're using. Using more threads helps significantly.

The output size depends on how many contigs match your criteria. In "any" mode, you might get quite a few contigs. In "all" mode, you'll typically get fewer but higher-confidence results.

If you don't get any results, try relaxing your E-value threshold or check that your HMM files are in the correct format. If you get too many results, try using "all" mode instead of "any" mode, or make your E-value threshold more strict.

## Requirements

You need to have pyrodigal, HMMER, and seqkit installed and available in your PATH. Your input contigs should be in FASTA format, and your HMM files should be in HMMER3 format (usually with .hmm extension).

The command expects DNA sequences as input, not protein sequences. If you already have predicted proteins, you should use HMMER directly rather than this command.
