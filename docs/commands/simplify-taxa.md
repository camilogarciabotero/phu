# simplify-taxa

Simplify vContact taxonomy prediction columns into compact lineage codes for streamlined downstream analysis.

## Overview

The `simplify-taxa` command processes vContact3 final assignment files, transforming verbose taxonomic prediction strings into compact, standardized lineage codes. This transformation makes taxonomic data more manageable for visualization, filtering, and analysis while preserving hierarchical relationships.

## Basic Usage

```bash
phu simplify-taxa -i <INPUT_FILE> -o <OUTPUT_FILE> [options]
```

Input accepts CSV or TSV files from vContact3's `final_assignments.csv` output. Output format is automatically detected from file extension.

## Input/Output Formats

### Supported Formats
- **Input**: CSV, TSV (auto-detected from extension or `--sep` parameter)
- **Output**: CSV, TSV (auto-detected from file extension)

### Expected Input Columns
The command automatically detects and processes any columns matching the pattern `*_prediction`:
- `kingdom_prediction`
- `phylum_prediction` 
- `class_prediction`
- `order_prediction`
- `family_prediction`
- `subfamily_prediction`
- `genus_prediction`
- `realm_prediction` (if present)

## Transformation Logic

### Before Transformation
```
novel_genus_1_of_novel_family_2_of_novel_order_3_of_Caudoviricetes
```

### After Transformation  
```
Caudoviricetes:NO3:NF2:NG1
```

### Compact Code Format

The transformation uses standardized rank codes:
- `NK` = Novel Kingdom
- `NP` = Novel Phylum  
- `NC` = Novel Class
- `NO` = Novel Order
- `NF` = Novel Family
- `NSF` = Novel Subfamily
- `NG` = Novel Genus

## Command Options

```bash
 Simplify vContact taxonomy prediction columns into compact lineage codes.   
                                                                             
 Transforms verbose vContact taxonomy strings like                           
 'novel_genus_1_of_novel_family_2_of_Caudoviricetes' into compact codes like 
 'Caudoviricetes:NF2:NG1'.                                                   
                                                                             
 Example:                                                                    
     phu simplify-taxa -i final_assignments.csv -o simplified.csv           
     --add-lineage                                                           
                                                                             
╭─ Options ─────────────────────────────────────────────────────────────────╮
│ *  --input           -i  PATH  Input vContact final_assignments.csv       │
│                              [required]                                   │
│ *  --output          -o  PATH  Output path (.csv or .tsv) [required]      │
│    --add-lineage         FLAG  Append compact_lineage column from deepest │
│                              simplified rank                              │
│    --lineage-col         TEXT  Name of the lineage column                 │
│                              [default: compact_lineage]                   │
│    --sep                 TEXT  Override delimiter: ',' or '\t'.           │
│                              Auto-detected from extension if not set      │
│    --help            -h        Show this message and exit.                │
╰───────────────────────────────────────────────────────────────────────────╯
```

## Examples

### Basic Usage

```bash
# Simplify taxonomy predictions in CSV format
phu simplify-taxa -i final_assignments.csv -o simplified_taxonomy.csv

# Process TSV format with automatic detection
phu simplify-taxa -i final_assignments.tsv -o simplified_taxonomy.tsv

# Override input delimiter detection
phu simplify-taxa -i data.txt -o output.csv --sep "\t"
```

### Advanced Usage

```bash
# Add compact lineage column with deepest available classification
phu simplify-taxa -i final_assignments.csv -o simplified.csv --add-lineage

# Customize lineage column name
phu simplify-taxa -i final_assignments.csv -o simplified.csv \
  --add-lineage --lineage-col "best_taxonomy"
```

## Lineage Column Feature

The `--add-lineage` option creates an additional column containing the deepest (most specific) available taxonomic classification for each sequence.

### Priority Order (Most → Least Specific)
1. `genus_prediction`
2. `subfamily_prediction`
3. `family_prediction` 
4. `order_prediction`
5. `class_prediction`
6. `phylum_prediction`
7. `kingdom_prediction`
8. `realm_prediction`

### Example Output
| Sequence | genus_prediction | family_prediction | compact_lineage |
|----------|------------------|-------------------|-----------------|
| seq1 | Caudoviricetes:NF2:NG1 | Caudoviricetes:NF2 | Caudoviricetes:NF2:NG1 |
| seq2 | - | Caudoviricetes:NF5 | Caudoviricetes:NF5 |
| seq3 | - | - | - |

## Special Cases Handled

### Edge Cases for "0" Chains
The tool correctly handles vContact2's special "0" designation patterns:

```bash
# Input
novel_class_0_of_novel_phylum_0_of_novel_kingdom_5_of_Duplodnaviria

# Output  
Duplodnaviria:NK5:NP0:NC0
```

### Multiple Candidates
When vContact2 provides multiple taxonomic candidates (separated by `||`), each is processed independently:

```bash
# Input
Caudoviricetes:NF1:NG2||Caudoviricetes:NF3:NG4

# Output
Caudoviricetes:NF1:NG2||Caudoviricetes:NF3:NG4
```

## Quality Assessment

After processing, the command provides a summary showing remaining `novel_` strings for quality control:

```
QA Summary:
  genus_prediction: 45 remaining 'novel_' strings
  family_prediction: 12 remaining 'novel_' strings
  order_prediction: 3 remaining 'novel_' strings
```

## Workflow Integration

### Typical Bioinformatics Pipeline

```bash
# 1. Run vContact3 (external)
vcontact3 --nucleotide <viral-genome.fasta> --output-dir <vcontact-output>

# 2. Simplify taxonomy predictions
phu simplify-taxa -i vcontact_output/final_assignments.csv \
  -o taxonomy_simplified.csv --add-lineage

# 3. Use simplified taxonomy for downstream analysis
# - Phylogenetic visualization
# - Diversity analysis  
# - Taxonomic filtering
```

## Comparison with Manual Processing

| Task | phu simplify-taxa | Manual Processing |
|------|-------------------|-------------------|
| **Complexity** | Single command | Custom scripts/regex |
| **Edge cases** | Automatically handled | Error-prone |
| **Consistency** | Standardized format | Variable approaches |
| **Speed** | Optimized pandas operations | Slower loops |
| **Maintenance** | Built-in updates | Manual fixes needed |

## Output File Structure

The output file preserves the original structure while transforming taxonomy columns:

```
Original columns + Simplified *_prediction columns [+ compact_lineage column]
```

All non-taxonomy columns remain unchanged, ensuring compatibility with existing workflows.