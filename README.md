# TrimNN Analysis Pipeline

A comprehensive pipeline for analyzing spatial transcriptomics data, including motif extraction, differential expression analysis, GO and pathway analysis, and cell-cell communication analysis.

## Overview

The pipeline performs the following analyses:
1. Motif extraction from spatial data
2. Differential Expression Gene (DEG) analysis
3. Gene Ontology (GO) and pathway analysis
4. CellChat analysis for cell-cell communication
5. Size and Effect Size analysis (Size3ES)

## Dataset Specifications

The pipeline is optimized for analyzing spatial transcriptomics data with the following specifications:

### Sample Information
- **Species**: Mouse (transgenic AD model)
- **Age Groups**: 
  - 8-month-old: control8 and sample8
  - 13-month-old: control13 and sample13
- **Conditions**: AD vs. Wild-type (Control)
- **Replicates**: 2 replicates each for AD and control at both timepoints

### Data Structure
- Each sample contains spatial gene expression data
- Cell type annotations are provided for each cell
- Expression matrices are in sparse format
- Raw count data is used for differential expression analysis

## Installation

### Prerequisites
- Python 3.9+
- R 4.0+
- Conda

### Setup
1. Clone the repository:
```bash
git clone https://github.com/SathvikSai1999/pattern_analysis_pipelie.git
cd pattern_analysis_pipelie
```

2. Make scripts executable:
```bash
chmod +x *.py *.r *.sh
```

3. Install dependencies:
```bash
./install_dependencies.sh
```

## Usage

### Basic Usage
```bash
# Run entire pipeline with default settings
./run_pattern_analysis_pipeline.sh

# Run specific analyses
./run_pattern_analysis_pipeline.sh --run-only deg,go
```

### Complete Pipeline Example
```bash
./run_pattern_analysis_pipeline.sh \
    # Motif Analysis
    --motif-input-dir data/motifs \
    --motif-output-dir output/motifs \
    
    # DEG Analysis
    --deg-input-dir data/DEG \
    --deg-output-dir output/DEG \
    --DEG-method both \
    --DEG-p-adj-cutoff 0.01 \
    
    # GO Analysis
    --go-input-dir data/GO \
    --go-output-dir output/GO \
    --go-pvalue-cutoff 0.005 \
    --pathway-db KEGG \
    --species mouse \
    --go-qvalue-cutoff 0.01 \
    --go-show-category 30 \
    
    # CellChat Analysis
    --cellchat-input-dir data/cellchat \
    --cellchat-output-dir output/cellchat \
    --cellchat-type mean \
    --cellchat-trim 0.2 \
    --cellchat-search "ECM-Receptor" \
    
    # Size3ES Analysis
    --size3es-input-dir data/size3es \
    --size3es-output-dir output/size3es \
    --size3es-months "control8,sample8,control13,sample13" \
    --size3es-replicates "1,2" \
    --size3es-control-group control \
    --size3es-pvalue-cutoff 0.05 \
    --size3es-effect-size-threshold 0.1 \
    --size3es-multiple-testing-correction BH
```

## Parameters

### General Options
- `--run-only`: Specify analyses to run [default: all]
  - Options: all, deg, go, cellchat, size3es
  - Multiple values allowed (comma-separated)

### Analysis-Specific Parameters

#### 1. DEG Analysis
- `--DEG-method`: Analysis method [default: both]
  - Options: deseq2, wilcox, both
  - Multiple values allowed (comma-separated)

  **Example Usage Cases:**

  1. **DESeq2 Method** (Recommended for bulk RNA-seq and large single-cell datasets):
     ```bash
     ./run_pattern_analysis_pipeline.sh \
         --run-only deg \
         --DEG-method deseq2 \
         --deg-input-dir data/DEG \
         --deg-output-dir output/DEG \
         --DEG-p-adj-cutoff 0.01
     ```
     - Uses negative binomial model
     - Accounts for biological variation
     - Best for datasets with multiple replicates
     - Output: CSV files with log2FoldChange, p-value, and adjusted p-value

  2. **Wilcoxon Method** (Recommended for small datasets or when data is not normally distributed):
     ```bash
     ./run_pattern_analysis_pipeline.sh \
         --run-only deg \
         --DEG-method wilcox \
         --deg-input-dir data/DEG \
         --deg-output-dir output/DEG \
         --DEG-p-adj-cutoff 0.01
     ```
     - Non-parametric test
     - No assumption of normal distribution
     - Good for small sample sizes
     - Output: CSV files with fold change and p-values

  3. **Both Methods** (Recommended for validation and comparison):
     ```bash
     ./run_pattern_analysis_pipeline.sh \
         --run-only deg \
         --DEG-method both \
         --deg-input-dir data/DEG \
         --deg-output-dir output/DEG \
         --DEG-p-adj-cutoff 0.01
     ```
     - Runs both DESeq2 and Wilcoxon analyses
     - Allows comparison of results
     - Useful for method validation
     - Output: Separate CSV files for each method

- `--DEG-p-adj-cutoff`: Adjusted p-value cutoff [default: 0.05]
  - Range: 0-1

#### 2. GO Analysis
- `--go-pvalue-cutoff`: P-value cutoff [default: 0.01]
  - Range: 0-1
- `--pathway-db`: Pathway database [default: Reactome]
  - Options: Reactome, KEGG, GO
- `--species`: Organism/species [default: mouse]
  - Options: mouse, human
- `--go-qvalue-cutoff`: Q-value cutoff [default: 0.05]
  - Range: 0-1
- `--go-show-category`: Categories to show [default: 20]
  - Range: 1-100

#### 3. CellChat Analysis
- `--cellchat-type`: Communication probability type [default: truncatedMean]
  - Options: truncatedMean, mean, median
- `--cellchat-trim`: Trimming factor [default: 0.1]
  - Range: 0-1
- `--cellchat-search`: Search type [default: Secreted Signaling]
  - Options: Secreted Signaling, ECM-Receptor, Cell-Cell Contact

#### 4. Size3ES Analysis
The Size3ES analysis examines cell type interactions and their effect sizes across different time points and replicates. It uses a combination of Fisher's exact test and Cramer's V to assess the significance and strength of cell type interactions.

- `--size3es-months`: Months to analyze [default: control8,sample8,control13,sample13]
  - Options: control8, sample8, control13, sample13
- `--size3es-replicates`: Replicates to use [default: 1,2]
  - Options: 1, 2
- `--size3es-control-group`: Control group name [default: control]
- `--size3es-pvalue-cutoff`: P-value cutoff [default: 0.05]
  - Range: 0-1
  - Used for filtering significant interactions
- `--size3es-effect-size-threshold`: Effect size threshold [default: 0.1]
  - Range: 0-1
  - Minimum Cramer's V value to consider an interaction meaningful
- `--size3es-multiple-testing-correction`: Correction method [default: BH]
  - Options: BH, bonferroni, holm, hochberg, none
  - BH (Benjamini-Hochberg) recommended for most cases

### Output Format
The Size3ES analysis generates CSV files with the following columns:
- `tri_celltype`: Cell type triplet (e.g., "Astro&CA1&CTX-Ex")
- `occurrence_num`: Number of occurrences
- `rank`: Interaction rank
- `meta_p`: Combined p-value from Fisher's exact tests
- `effect_size_1`: Effect size for first cell type
- `effect_size_2`: Effect size for second cell type
- `effect_size_3`: Effect size for third cell type

## Input/Output Structure

### Input Files
- `motif_ids.csv`: Motif identifiers
- `matrix_raw.txt`: Raw gene expression matrix
- `cell_types.txt`: Cell type annotations

### Output Directories
Each analysis step has its own input and output directories:
```
project_root/
├── data/                  # Input directories
│   ├── motifs/
│   ├── DEG/
│   ├── GO/
│   ├── cellchat/
│   └── size3es/
├── output/               # Output directories
│   ├── motifs/
│   ├── DEG/
│   ├── GO/
│   ├── cellchat/
│   └── size3es/
└── logs/                # Pipeline execution logs
```
