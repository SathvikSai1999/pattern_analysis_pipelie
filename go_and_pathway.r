library(ggplot2)
library(scales)
library(ggpubr)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)
library(ComplexHeatmap)
library(ComplexUpset)
library(ggvenn)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
qvalue_cutoff <- 0.05
ontology_types <- c("BP", "MF", "CC")
show_category <- 20

# Parse arguments
for (arg in args) {
    if (grepl("--qvalue-cutoff=", arg)) {
        qvalue_cutoff <- as.numeric(sub("--qvalue-cutoff=", "", arg))
    } else if (grepl("--ontology-types=", arg)) {
        ontology_types <- strsplit(sub("--ontology-types=", "", arg), ",")[[1]]
    } else if (grepl("--show-category=", arg)) {
        show_category <- as.numeric(sub("--show-category=", "", arg))
    }
}

motifs <- c('cccm')

go <- function(motifs){
   for (region in c('3hop','original'))
    for (mtf in motifs)
        for (ont in ontology_types)
        {
            gene <- read.csv(paste0("output/deg_overexpressed/",region,'/',mtf,".txt"),head=FALSE)

            gene_id <- bitr(gene$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

            ego <- enrichGO(gene = gene$V1,
                            OrgDb = org.Mm.eg.db,
                            keyType = 'SYMBOL',
                            ont = ont,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = qvalue_cutoff)
            print(paste("Processing GO analysis for:", mtf, "Region:", region, "Ontology:", ont))
            print(paste("Number of enriched terms found:", nrow(summary(ego))))

            if (nrow(summary(ego)) > 0){
                path <- paste0("output/Go/",region,'/')
                if (!file.exists(path)) {
                    dir.create(path, recursive = TRUE)
                }
                write.csv(ego, file = paste0("output/Go/",region,'/',mtf,"_",ont,".csv"), row.names = FALSE)
                # Save only PNG, no PDF
                png(file=paste0("output/Go/",region,'/',mtf,"_",ont,".png"),width = 10,height = 10,units = "in", res = 600)
                print(dotplot(ego, showCategory=show_category, color="qvalue", label_format = 50))
                dev.off()
            }
        }
}

pathway <- function(motifs){
    for (region in c('3hop','original'))
        for (mtf in motifs)
        {
            gene <- read.csv(paste0("output/deg_overexpressed/",region,'/',mtf,".txt"),head=FALSE)

            gene_id <- bitr(gene$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

            # Pathway analysis
            Gene_IDs <- mapIds(org.Mm.eg.db,
                               keys=gene$V1,
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
            egoKEGG <- enrichPathway(Gene_IDs, organism='mouse')

            print(paste("Processing Pathway analysis for:", mtf, "Region:", region))
            print(paste("Number of enriched pathways found:", nrow(summary(egoKEGG))))

            if (nrow(summary(egoKEGG)) > 0){
                path <- paste0("output/pathway/",region,'/')
                if (!file.exists(path)) {
                    dir.create(path, recursive = TRUE)
                }
                write.csv(egoKEGG, file = paste0("output/pathway/",region,'/',mtf,".csv"), row.names = FALSE)
                # Save only PNG, no PDF
                png(file=paste0("output/pathway/",region,'/',mtf,".png"),width = 10,height = 10,units = "in",res = 600)
                print(dotplot(egoKEGG, showCategory=show_category, color="qvalue", label_format = 50))
                dev.off()
            }
        }
}

# Call functions to execute the analysis
go(motifs)
pathway(motifs)

