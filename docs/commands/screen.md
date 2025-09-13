# Screen

## Overview

The `phu screen` command screens DNA contigs for specific protein families using HMMER on predicted coding sequences (CDS). It provides a streamlined workflow that:

1. Predicts proteins from input contigs using pyrodigal
2. Searches predicted proteins against HMM profiles using hmmsearch
3. Identifies and filters significant hits
4. Extracts contigs containing matches

## Synopsis

```bash
phu screen [OPTIONS] --input-contigs FASTA_FILE --hmm HMM_FILE
```

## Description

The screen command is designed for high-throughput screening of genomic contigs to identify sequences containing specific protein domains or families. It's particularly useful for:

- Identifying phage contigs in metagenomic assemblies
- Screening for specific protein families across large contig datasets
- Filtering assemblies based on protein domain content
- Quality control and contamination detection

## Required Arguments

- `--input-contigs PATH`: Input contig sequences in FASTA format
- `--hmm PATH`: HMM profile file (single profile or database)

## Optional Arguments

### Output Control
- `--outdir PATH`: Output directory (default: `phu-screen`)
- `--keep-proteins`: Keep intermediate protein FASTA file
- `--keep-domtbl`: Keep HMMER domain table output (default: True)

### Gene Prediction Options
- `--mode {meta,single}`: Pyrodigal mode (default: `meta`)
  - `meta`: Use pre-trained metagenomic models
  - `single`: Train on input sequences (slower, for single genomes)
- `--min-gene-len INT`: Minimum gene length in nucleotides (default: 90)
- `--translation-table INT`: Genetic code table (default: 11)

### Filtering Options
- `--min-bitscore FLOAT`: Minimum bitscore threshold
- `--max-evalue FLOAT`: Maximum E-value threshold (default: 1e-5)
- `--top-per-contig INT`: Maximum hits to keep per contig (default: 1)

### Performance Options
- `--threads INT`: Number of CPU threads for HMMER (default: 1)

## Output Files

The screen command creates several output files in the specified output directory:

### Primary Outputs
- `screened_contigs.fasta`: Contigs containing significant hits
- `kept_contigs.txt`: List of contig IDs that passed filtering
- `kept_hits.json`: Detailed information about significant hits

### Intermediate Files (optional)
- `proteins.faa`: Predicted protein sequences (if `--keep-proteins`)
- `hits.domtblout`: Raw HMMER domain table output (if `--keep-domtbl`)

## Workflow Details

### 1. Gene Prediction
Uses pyrodigal to predict coding sequences from input contigs:
- Creates protein IDs in format: `{contig_id}|gene{number}`
- Handles both metagenomic and single-genome modes
- Filters genes by minimum length

### 2. Protein Search
Runs hmmsearch against predicted proteins:
- Uses proper HMMER command structure: `hmmsearch [options] <hmmfile> <seqfile>`
- Includes threading support with `--cpu` parameter
- Generates domain table output with `--domtblout`

### 3. Hit Processing
Parses HMMER results and applies filtering:
- Extracts contig IDs from protein sequence names
- Applies E-value and bitscore thresholds
- Selects top hits per contig based on bitscore

### 4. Contig Extraction
Uses seqkit to extract matching contigs:
- Creates temporary ID list file
- Extracts sequences matching the IDs
- Preserves original contig names and descriptions

## Examples

### Basic Usage
```bash
# Screen contigs for a specific protein family
phu screen --input-contigs assembly.fasta --hmm pfam_domain.hmm

# Use custom output directory
phu screen --input-contigs contigs.fa --hmm family.hmm --outdir screening_results
```

### Advanced Filtering
```bash
# Apply strict E-value threshold
phu screen --input-contigs assembly.fasta --hmm domain.hmm --max-evalue 1e-10

# Keep multiple hits per contig
phu screen --input-contigs contigs.fa --hmm family.hmm --top-per-contig 3

# Use bitscore threshold
phu screen --input-contigs assembly.fasta --hmm domain.hmm --min-bitscore 50
```

### Performance Optimization
```bash
# Use multiple threads
phu screen --input-contigs large_assembly.fasta --hmm domain.hmm --threads 16

# Single genome mode (more accurate for complete genomes)
phu screen --input-contigs genome.fasta --hmm family.hmm --mode single
```

### File Management
```bash
# Keep intermediate files for inspection
phu screen --input-contigs contigs.fa --hmm domain.hmm \
    --keep-proteins --keep-domtbl

# Clean run (minimal output files)
phu screen --input-contigs assembly.fasta --hmm family.hmm \
    --no-keep-proteins --no-keep-domtbl
```

## Understanding Results

### Hit Quality Assessment
- **E-value**: Statistical significance (lower = better)
- **Bitscore**: Raw match score (higher = better)
- **Domain coverage**: Check alignment coordinates in JSON output

### Common Patterns
- **High E-value with low bitscore**: Likely false positive
- **Multiple hits per contig**: May indicate domain repeats or multidomain proteins
- **Edge effects**: Partial domains at contig ends

### Troubleshooting
- **No proteins predicted**: Check input sequence format and minimum gene length
- **No hits found**: Verify HMM file format and consider relaxing thresholds
- **seqkit errors**: Ensure contig IDs don't contain special characters

## Integration with Other Tools

### Input Preparation
```bash
# From SPAdes assembly
cp contigs.fasta input_contigs.fasta

# Filter by length first
seqkit seq -m 1000 assembly.fasta > filtered_contigs.fasta
```

### Downstream Analysis
```bash
# Count results
wc -l phu-screen/kept_contigs.txt

# Get statistics
seqkit stats phu-screen/screened_contigs.fasta

# Further analysis
blastn -query phu-screen/screened_contigs.fasta -db nt
```

## Dependencies

### Required Software
- **pyrodigal**: Gene prediction (≥3.0)
- **HMMER**: Protein profile searches (≥3.3)
- **seqkit**: Sequence manipulation

### File Formats
- **Input contigs**: FASTA format (DNA sequences)
- **HMM profiles**: HMMER3 format (.hmm files)
- **Output**: FASTA, text, and JSON formats

## Citation

If you use the screen command in your research, please cite:
- HMMER: [Eddy, S.R. "Accelerated profile HMM searches." PLoS Comp. Biol. 7:e1002195, 2011]
- Pyrodigal: [Larralde, M. "Pyrodigal: Python bindings and interface to Prodigal." JOSS 7:4296, 2022]
