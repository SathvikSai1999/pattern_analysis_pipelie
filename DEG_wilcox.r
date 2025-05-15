library(Matrix)
library(readr)

# Get input directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

gene_mtx <- file.path(input_dir, '13months-disease-replicate_1_matrix_raw.txt')
gene <- read.csv(file.path(input_dir, 'ad_gene.txt'), header = FALSE)
mtx <- readMM(gene_mtx)

folder_path <- file.path(input_dir, 'tree_search/cell_label/disease/')
files <- list.files(folder_path)
rownames(mtx) <- gene$V1
CCC_dense <- as.matrix(mtx)

for (f in files) {
  output_file <- file.path(output_dir, 'wilcox', paste0(f, '_deg.csv'))
  if (file.exists(output_file)) {
    next
  }
  
  group <- read.csv(file.path(folder_path, f, 'mtf&non-motif_label.txt'), header = FALSE)
  group$V1 <- factor(group$V1)
  colnames(CCC_dense) <- seq_len(ncol(CCC_dense))  # ensure cell columns are indexed
  labels <- group$V1

  # Make sure group labels match number of cells
  if (length(labels) != ncol(CCC_dense)) {
    warning(paste("Label length does not match number of cells in", f))
    next
  }

  res_df <- data.frame(
    gene = rownames(CCC_dense),
    p_val = NA,
    avg_logFC = NA
  )
  print(f)

  # Perform Wilcoxon rank-sum test per gene
  for (i in seq_len(nrow(CCC_dense))) {
    expr <- CCC_dense[i, ]
    group1 <- expr[labels == levels(labels)[1]] #reference condition: non-motif
    group2 <- expr[labels == levels(labels)[2]]

    if (length(unique(expr)) > 1) {
      test <- wilcox.test(group1, group2)
      res_df$p_val[i] <- test$p.value
      res_df$avg_logFC[i] <- log2(mean(group2 + 1) / mean(group1 + 1))
    }
  }

  res_df$p_adj <- p.adjust(res_df$p_val, method = "fdr")

  # Create output directory if it doesn't exist
  output_path <- file.path(output_dir, 'wilcox')
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  write.csv(res_df, file = output_file, row.names = FALSE)
}
