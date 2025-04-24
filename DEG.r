library(Matrix)
library(DESeq2)
library(parallel)

# Set number of cores for parallel processing
num_cores <- detectCores()
options(mc.cores = num_cores)

mtf = 'cccm'
for (region in c('original','3hop')){
    print(paste("Processing region:", region))
    
    gene_mtx <- paste0('data/matrix_raw.txt')
    gene <- read.csv('data/ad_gene.txt',header=FALSE)
    mtx <- readMM(gene_mtx)
    group <- read.csv(paste0('data/extracted_matrix/',region,'/',mtf,'/mtf&non-motif_label.txt'),header=FALSE)
    group$V1 <- factor(group$V1)
    rownames(mtx) <- gene$V1
    CCC_dense <- as.matrix(mtx)
    CCC_dense <- CCC_dense + 1

    #DEG with parallel processing
    print("Starting DESeq2 analysis...")
    dds <- DESeqDataSetFromMatrix(countData = CCC_dense,
                                colData = group,
                                design = ~ V1)
    dds <- DESeq(dds, parallel = TRUE)
    res <- results(dds)

    path <- paste0("output/DEG/",region,'/')
    if (!file.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    write.csv(as.data.frame(res), file = paste0(path, mtf,'.csv'))
    print(paste("Completed processing region:", region))
}
