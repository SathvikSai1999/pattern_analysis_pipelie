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

motifs <- c('cccm')

go <- function(motifs){
   for (region in c('3hop','original'))
    for (mtf in motifs)
        for (ont in c("BP","MF","CC"))
        {
            gene <- read.csv(paste0("output/deg_overexpressed/",region,'/',mtf,".txt"),head=FALSE)

            gene_id <- bitr(gene$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

            ego <- enrichGO(gene = gene$V1,
                            OrgDb = org.Mm.eg.db,
                            keyType = 'SYMBOL',
                            ont = ont,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05)
            print(paste("Processing GO analysis for:", mtf, "Region:", region, "Ontology:", ont))
            print(paste("Number of enriched terms found:", nrow(summary(ego))))

            if (nrow(summary(ego)) > 0){
                path <- paste0("output/Go/",region,'/')
                if (!file.exists(path)) {
                    dir.create(path, recursive = TRUE)
                }
                write.csv(ego, file = paste0("output/Go/",region,'/',mtf,"_",ont,".csv"), row.names = FALSE)
                GeneEnrich <- dotplot(ego, showCategory=20, color="qvalue", label_format = 50)
                png(file=paste0("output/Go/",region,'/',mtf,"_",ont,".png"),width = 10,height = 10,units = "in", res = 600)
                print(GeneEnrich)
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

                ReacEnrich <- dotplot(egoKEGG, showCategory=20, color="qvalue", label_format = 50)
                png(file=paste0("output/pathway/",region,'/',mtf,".png"),width = 10,height = 10,units = "in",res = 600)
                print(ReacEnrich)
                dev.off()
            }
        }
}

#  Call functions to execute the analysis
go(motifs)
pathway(motifs)

