# pattern_analysis_pipeline

## Overview: 
The pipeline utilizes various scripts to extract motifs, perform differential gene expression analysis, and analyze cell interactions using the CellChat package.


# Running the pipeline
1. install_dependencies.sh
2. run_pttern_analysis_pipeline.sh


## Table of Contents
- [Input Files](#input-files) from Trimnn (motif_ids.csv),(matrix_raw.txt), (cell_types.txt)

- [Scripts](#scripts): run_pttern_analysis_pipeline.sh executes the follwoing scripts

1. [extract_motifs.py](#extract_motifspy)
- **Purpose**: Extracts cellular patterns (motifs) from spatial transcriptomics data.
- **Outputs**: 
  - `cell_type_matrix.txt`: Matrix of cell types in each pattern.
  - `gene_expression_matrix.txt`: Gene expression values for each cell.
  - `mtf&non-motif_label.txt`: Labels indicating if it's a motif (1) or non-motif (0).
  
  
2.  [DEG.r](#degr)
- **Purpose**: Performs differential gene expression analysis between motifs and non-motifs.
- **Outputs**: 
  - `output/DEG/cccm.csv`: CSV file containing statistics for each gene, including log2 fold change and p-values.
  
  
3.  [extract_deg.py](#extract_degpy)
- **Purpose**: Filters DEG results to identify significantly overexpressed genes.
- **Outputs**: 
  - `output/deg_overexpressed/original/cccm.txt`: List of significant genes in original patterns.
  - `output/deg_overexpressed/3hop/cccm.txt`: List of significant genes in 3-hop patterns.
  
4.  [go_and_pathway.r](#go_and_pathwayr)
- **Purpose**: Conducts Gene Ontology (GO) and pathway analysis on differentially expressed genes.
- **Outputs**: 
  - `output/Go/original/cccm_BP.csv`: GO analysis results for Biological Processes.
  - `output/Go/original/cccm_MF.csv`: GO analysis results for Molecular Functions.
  - `output/Go/original/cccm_CC.csv`: GO analysis results for Cellular Components.
  - `output/Go/3hop/cccm_BP.csv`: GO analysis results for 3-hop patterns.
  - `output/pathway/original/cccm.csv`: Pathway analysis results.
  - `output/pathway/3hop/cccm.csv`: Pathway analysis results for 3-hop patterns.
  
5.  [generate_rank_files.py] (#works on exteact_motif output files)
- **Purpose**: Generates ranked lists of genes and patterns for downstream analysis.
- **Outputs**: 
  - `original_rank.txt`: Ranked list of genes for original patterns.
  - `3hop_rank.txt`: Ranked list of genes for 3-hop extended patterns.
  - `original_pattern_rank.txt`: Ranked list of patterns for original analysis.
  - `3hop_pattern_rank.txt`: Ranked list of patterns for 3-hop analysis.
  
6.  [size3ES&P.r](#size3espr)
- **Purpose**: Analyzes effect sizes and statistical significance of cell type interactions.
- **Outputs**: 
  - `effect_size/AD_8_control_1_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 8, replicate 1.
  - `effect_size/AD_8_control_2_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 8, replicate 2.
  - `effect_size/AD_13_control_1_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 13, replicate 1.
  - `effect_size/AD_13_control_2_tri_ES&P.csv`: Effect sizes and p-values for interactions in month 13, replicate 2.
  
7.  [cellchat.r](#cellchatr)
- **Purpose**: Analyzes cell-cell communication using the CellChat package.
- **Outputs**: 
  - Communication probability matrices and visualizations (specific output files may vary based on the analysis).
  
  


```

