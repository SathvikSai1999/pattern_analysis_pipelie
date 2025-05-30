# TrimNN Analysis Pipeline

A comprehensive pipeline for analyzing spatial transcriptomics data, including motif extraction, differential expression analysis, GO analysis, pathway analysis, and cell-cell communication analysis.

## Overview

The pipeline performs the following analyses:
1. Motif extraction from spatial data
2. Differential Expression Gene (DEG) analysis
3. Gene Ontology (GO) analysis
4. Pathway analysis
5. CellChat analysis for cell-cell communication
6. Size and Effect Size analysis (Size3ES)

## Dataset Specifications

The pipeline is optimized for analyzing spatial transcriptomics data with the following specifications:

### Sample Information
- **Species**: Mouse (transgenic AD model)
- **Age Groups**: 
  - 8-month-old: control_8 and sample_8
  - 13-month-old: control_13 and sample_13
- **Conditions**: AD/Sample vs. Wild-type (Control)
- **Replicates**: 2 replicates each for AD/Sample and control at both timepoints

### Replicate Structure
|   Time Point   | Condition | Number of Replicates |
|----------------|-----------|---------------------|
| 8 months old   | AD/Sample |           2         |
| 8 months old   |  Control  |           2         |
| 13 months old  | AD/Sample |           2         |
| 13 months old  |  Control  |           2         |

Each replicate is labeled with an underscore and number (e.g., control_8_1, control_8_2, sample_8_1, sample_8_2).

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
./run_pattern_analysis_pipeline.sh --run-only deg,go,pathway
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
    --go-qvalue-cutoff 0.01 \
    --go-show-category 30 \
    
    # Pathway Analysis
    --pathway-input-dir data/pathway \
    --pathway-output-dir output/pathway \
    --pathway-pvalue-cutoff 0.005 \
    --pathway-qvalue-cutoff 0.01 \
    --pathway-show-category 30 \
    --pathway-db KEGG \
    
    # CellChat Analysis
    --cellchat-input-dir data/cellchat \
    --cellchat-output-dir output/cellchat \
    --cellchat-type mean \
    --cellchat-trim 0.2 \
    --cellchat-search "ECM-Receptor" \
    
    # Size3ES Analysis
    --size3es-control- data/size3es/spatial_8months-control-replicate_1.csv \
    --size3es-control- data/size3es/spatial_8months-control-replicate_2.csv \
    --size3es-control- data/size3es/spatial_13months-control-replicate_1.csv \
    --size3es-control- data/size3es/spatial_13months-control-replicate_2.csv \
    --size3es-sample- data/size3es/spatial_8months-disease-replicate_1.csv \
    --size3es-sample- data/size3es/spatial_8months-disease-replicate_2.csv \
    --size3es-sample- data/size3es/spatial_13months-disease-replicate_1.csv \
    --size3es-sample- data/size3es/spatial_13months-disease-replicate_2.csv \
    --size3es-output-dir output/size3es \
    --size3es-pvalue-cutoff 0.05 \
    --size3es-effect-size-threshold 0.1 \
    --size3es-multiple-testing-correction BH
```

## Parameters

### General Options
- `--run-only`: Specify analyses to run [default: all]
  - Options: all, deg, go, pathway, cellchat, size3es
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
- `--go-qvalue-cutoff`: Q-value cutoff [default: 0.05]
  - Range: 0-1
- `--go-show-category`: Categories to show [default: 20]
  - Range: 1-100

#### 3. Pathway Analysis
- `--pathway-pvalue-cutoff`: P-value cutoff [default: 0.01]
  - Range: 0-1
- `--pathway-qvalue-cutoff`: Q-value cutoff [default: 0.05]
  - Range: 0-1
- `--pathway-show-category`: Categories to show [default: 20]
  - Range: 1-100
- `--pathway-db`: Pathway database [default: Reactome]
  - Options: Reactome, KEGG

#### 4. CellChat Analysis
- `--cellchat-type`: Communication probability type [default: truncatedMean]
  - Options: truncatedMean, mean, median
- `--cellchat-trim`: Trimming factor [default: 0.1]
  - Range: 0-1
- `--cellchat-search`: Search type [default: Secreted Signaling]
  - Options: Secreted Signaling, ECM-Receptor, Cell-Cell Contact

#### 5. Size3ES Analysis
The Size3ES analysis examines cell type interactions and their effect sizes between control and sample conditions. It uses a combination of Fisher's exact test and Cramer's V to assess the significance and strength of cell type interactions.

##### Input Data Structure
Expected directory structure:
```
data/
├── control/
│   ├── 8m/
│   │   ├── r1_spatial.csv
│   │   └── r2_spatial.csv
│   └── 13m/
│       ├── r1_spatial.csv
│       └── r2_spatial.csv
└── sample/
    ├── 8m/
    │   ├── r1_spatial.csv
    │   └── r2_spatial.csv
    └── 13m/
        ├── r1_spatial.csv
        └── r2_spatial.csv
```

##### Command Line Arguments
- `--size3es-control-dir`: Directory containing control samples
- `--size3es-sample-dir`: Directory containing sample/AD samples
- `--size3es-timepoints`: Comma-separated list of time points (e.g., "8,13")
- `--size3es-replicates`: Comma-separated list of replicates (e.g., "1,2")

##### Output Files
- Control results (in `output/size3es/`):
  - `control_8m_tri_ES&P.csv`: 8-month results
  - `control_13m_tri_ES&P.csv`: 13-month results

- Sample results (in `output/size3es/`):
  - `sample_8m_tri_ES&P.csv`: 8-month results
  - `sample_13m_tri_ES&P.csv`: 13-month results

- Comparison results (in `output/size3es/`):
  - `comparison_8m_results.csv`: 8-month comparison
  - `comparison_13m_results.csv`: 13-month comparison

##### Statistical Parameters
- `--size3es-pvalue-cutoff`: P-value threshold [default: 0.05]
  - Range: 0-1
  - Used for filtering significant interactions

- `--size3es-effect-size-threshold`: Minimum effect size [default: 0.1]
  - Range: 0-1
  - Minimum Cramer's V value to consider an interaction meaningful

- `--size3es-multiple-testing-correction`: P-value correction method [default: BH]
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
