library(Matrix)
library(DESeq2)
library(readr)

# Get input directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

gene_mtx <- file.path(input_dir, '13months-disease-replicate_1_matrix_raw.txt')
gene <- read.csv(file.path(input_dir, 'ad_gene.txt'), header=FALSE)
mtx <- readMM(gene_mtx)

folder_path <- file.path(input_dir, 'tree_search/cell_label/disease/')

files <- list.files(folder_path)
files <- sample(files)

for (f in files){
  output_file <- file.path(output_dir, 'deseq2', paste0(f, '_deg.csv'))
  if (file.exists(output_file)) {
  next
}

  group <- read.csv(file.path(folder_path, f, 'mtf&non-motif_label.txt'), header=FALSE)
group$V1 <- factor(group$V1)
rownames(mtx) <- gene$V1
CCC_dense <- as.matrix(mtx)

print(f)
#DEG
dds <- DESeqDataSetFromMatrix(countData = CCC_dense,
                            colData = group,
                            design = ~ V1)

dds <- estimateSizeFactors(dds, type="poscounts")
library(scran)
scr <- computeSumFactors(dds)
# use scran's sum factors:
sizeFactors(dds) <- sizeFactors(scr)

  dds <- DESeq(dds, test="LRT", reduced=~1, useT=TRUE, minReplicatesForReplace=Inf, minmu=1e-6)
res <- results(dds)

  # Create output directory if it doesn't exist
  output_path <- file.path(output_dir, 'deseq2')
  if (!file.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}
  write.csv(as.data.frame(res), file = output_file)
}

