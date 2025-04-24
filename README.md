# TrimNN Analysis Pipeline

A comprehensive pipeline for analyzing spatial transcriptomics data, including motif extraction, differential expression analysis, GO and pathway analysis, and cell-cell communication analysis.

## Overview

The pipeline utilizes various scripts to extract motifs, perform differential gene expression analysis, and analyze cell interactions using the CellChat package. It performs the following analyses:
1. Motif extraction from spatial data
2. Differential Expression Gene (DEG) analysis
3. Gene Ontology (GO) and pathway analysis
4. CellChat analysis for cell-cell communication
5. Size and Effect Size analysis

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

## Running the Pipeline

### Basic Usage
```bash
./run_pattern_analysis_pipeline.sh
```
This will run the pipeline with all default parameters.

### Advanced Usage with Parameters
```bash
./run_pattern_analysis_pipeline.sh \
  --hop-distance 3 \
  --DEG-p-adj-cutoff 0.05 \
  --go-pvalue-cutoff 0.01 \
  --pathway-db Reactome \
  --go-qvalue-cutoff 0.05 \
  --go-ontology-types "BP,MF,CC" \
  --go-show-category 20 \
  --cellchat-type truncatedMean \
  --cellchat-trim 0.1 \
  --cellchat-search "Secreted Signaling" \
  --size3es-months "8,13" \
  --size3es-replicates "1,2" \
  --size3es-control-group control \
  --size3es-pvalue-cutoff 0.05 \
  --size3es-effect-size-threshold 0.1 \
  --size3es-multiple-testing-correction BH
```

## Table of Contents

### Input Files from TrimNN
- `motif_ids.csv`: Contains motif identifiers
- `matrix_raw.txt`: Raw gene expression matrix
- `cell_types.txt`: Cell type annotations

### Pipeline Scripts

The `run_pattern_analysis_pipeline.sh` executes the following scripts in sequence:

#### 1. extract_motifs.py
**Purpose**: Extracts cellular patterns (motifs) from spatial transcriptomics data.

**Parameters**:
- `--hop-distance`: Hop distance for motif extraction [default: 3, range: 1-10]

**Outputs**:
- `cell_type_matrix.txt`: Matrix of cell types in each pattern
- `gene_expression_matrix.txt`: Gene expression values for each cell
- `mtf&non-motif_label.txt`: Labels indicating if it's a motif (1) or non-motif (0)

#### 2. DEG.r
**Purpose**: Performs differential gene expression analysis between motifs and non-motifs.

**Parameters**:
- `--DEG-p-adj-cutoff`: Adjusted p-value cutoff [default: 0.05, range: 0-1]

**Outputs**:
- `output/DEG/cccm.csv`: CSV file containing statistics for each gene, including log2 fold change and p-values

#### 3. extract_deg.py
**Purpose**: Filters DEG results to identify significantly overexpressed genes.

**Outputs**:
- `output/deg_overexpressed/original/cccm.txt`: List of significant genes in original patterns
- `output/deg_overexpressed/3hop/cccm.txt`: List of significant genes in 3-hop patterns

#### 4. go_and_pathway.r
**Purpose**: Conducts Gene Ontology (GO) and pathway analysis on differentially expressed genes.

**Parameters**:
- `--go-pvalue-cutoff`: P-value cutoff [default: 0.01, range: 0-1]
- `--pathway-db`: Pathway database selection [default: Reactome, options: Reactome, KEGG, GO]
- `--go-qvalue-cutoff`: Q-value cutoff [default: 0.05, range: 0-1]
- `--go-ontology-types`: Comma-separated list of GO ontology types [default: BP,MF,CC, options: BP, MF, CC]
- `--go-show-category`: Number of categories to show [default: 20, range: 1-100]

**Outputs**:
- `output/Go/original/cccm_BP.csv`: GO analysis results for Biological Processes
- `output/Go/original/cccm_MF.csv`: GO analysis results for Molecular Functions
- `output/Go/original/cccm_CC.csv`: GO analysis results for Cellular Components
- `output/Go/3hop/cccm_BP.csv`: GO analysis results for 3-hop patterns
- `output/pathway/original/cccm.csv`: Pathway analysis results
- `output/pathway/3hop/cccm.csv`: Pathway analysis results for 3-hop patterns

#### 5. generate_rank_files.py
**Purpose**: Generates ranked lists of genes and patterns for downstream analysis.

**Outputs**:
- `original_rank.txt`: Ranked list of genes for original patterns
- `3hop_rank.txt`: Ranked list of genes for 3-hop extended patterns
- `original_pattern_rank.txt`: Ranked list of patterns for original analysis
- `3hop_pattern_rank.txt`: Ranked list of patterns for 3-hop analysis

#### 6. size3ES&P.r
**Purpose**: Analyzes effect sizes and statistical significance of cell type interactions.

**Parameters**:
- `--size3es-months`: Comma-separated list of months [default: 8,13, options: 8, 13]
- `--size3es-replicates`: Comma-separated list of replicates [default: 1,2, options: 1, 2]
- `--size3es-control-group`: Control group name [default: control]
- `--size3es-pvalue-cutoff`: P-value cutoff [default: 0.05, range: 0-1]
- `--size3es-effect-size-threshold`: Effect size threshold [default: 0.1, range: 0-1]
- `--size3es-multiple-testing-correction`: Multiple testing correction method [default: BH, options: BH, bonferroni, holm, hochberg, none]

**Outputs**:
- `effect_size/AD_8_control_1_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 8, replicate 1
- `effect_size/AD_8_control_2_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 8, replicate 2
- `effect_size/AD_13_control_1_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 13, replicate 1
- `effect_size/AD_13_control_2_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 13, replicate 2

#### 7. cellchat.r
**Purpose**: Analyzes cell-cell communication using the CellChat package.

**Parameters**:
- `--cellchat-type`: Communication probability type [default: truncatedMean, options: truncatedMean, mean, median]
- `--cellchat-trim`: Trimming factor [default: 0.1, range: 0-1]
- `--cellchat-search`: Search type [default: Secreted Signaling, options: Secreted Signaling, ECM-Receptor, Cell-Cell Contact]

**Outputs**:
- Communication probability matrices and visualizations (specific output files may vary based on the analysis)

