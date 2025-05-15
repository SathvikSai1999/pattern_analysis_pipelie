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
species <- "mouse"  # Default species
pathway_db <- "Reactome"  # Default pathway database

# Function to check required packages
check_required_packages <- function(species) {
    required_packages <- c("clusterProfiler", "ReactomePA")
    if (species == "mouse") {
        required_packages <- c(required_packages, "org.Mm.eg.db")
    } else if (species == "human") {
        required_packages <- c(required_packages, "org.Hs.eg.db")
    }
    
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("Required package", pkg, "is not installed. Please install it first."))
        }
    }
}

# Function to check input files
check_input_files <- function(region, mtf) {
    file_path <- paste0("output/deg_overexpressed/", region, '/', mtf, ".txt")
    if (!file.exists(file_path)) {
        stop(paste("Input file not found:", file_path))
    }
}

# Function to check gene conversion
check_gene_conversion <- function(gene_ids, original_genes) {
    if (nrow(gene_ids) == 0) {
        stop("No genes could be converted to Entrez IDs. Please check your gene symbols.")
    }
    if (nrow(gene_ids) < length(original_genes)) {
        warning(paste("Only", nrow(gene_ids), "out of", length(original_genes), 
                     "genes were successfully converted to Entrez IDs."))
    }
    return(gene_ids)
}

# Function to get organism code
get_organism_code <- function(species) {
    switch(tolower(species),
           "mouse" = "mmu",
           "human" = "hsa",
           "rat" = "rno",
           "zebrafish" = "dre",
           "fruitfly" = "dme",
           "yeast" = "sce",
           stop("Unsupported species. Please choose from: mouse, human, rat, zebrafish, fruitfly, yeast")
    )
}

# Function to get organism database
get_organism_db <- function(species) {
    switch(tolower(species),
           "mouse" = org.Mm.eg.db,
           "human" = org.Hs.eg.db,
           stop("Unsupported species. Please choose from: mouse, human")
    )
}

# Parse arguments
for (arg in args) {
    if (grepl("--qvalue-cutoff=", arg)) {
        qvalue_cutoff <- as.numeric(sub("--qvalue-cutoff=", "", arg))
    } else if (grepl("--ontology-types=", arg)) {
        ontology_types <- strsplit(sub("--ontology-types=", "", arg), ",")[[1]]
    } else if (grepl("--show-category=", arg)) {
        show_category <- as.numeric(sub("--show-category=", "", arg))
    } else if (grepl("--species=", arg)) {
        species <- sub("--species=", "", arg)
    } else if (grepl("--pathway-db=", arg)) {
        pathway_db <- sub("--pathway-db=", "", arg)
    }
}

# Check required packages
check_required_packages(species)

# Get organism database and code
org_db <- get_organism_db(species)
org_code <- get_organism_code(species)

motifs <- c('cccm')

go <- function(motifs){
   for (region in c('3hop','original'))
    for (mtf in motifs) {
        # Check input file
        check_input_files(region, mtf)
        
        for (ont in ontology_types) {
            gene <- read.csv(paste0("output/deg_overexpressed/",region,'/',mtf,".txt"),head=FALSE)

            # Convert gene symbols to Entrez IDs with error handling
            tryCatch({
                gene_id <- bitr(gene$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
                gene_id <- check_gene_conversion(gene_id, gene$V1)

            ego <- enrichGO(gene = gene$V1,
                              OrgDb = org_db,
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
                png(file=paste0("output/Go/",region,'/',mtf,"_",ont,".png"),width = 10,height = 10,units = "in", res = 600)
                print(dotplot(ego, showCategory=show_category, color="qvalue", label_format = 50))
                dev.off()
                } else {
                    warning(paste("No enriched terms found for", mtf, "Region:", region, "Ontology:", ont))
                }
            }, error = function(e) {
                warning(paste("Error in GO analysis for", mtf, "Region:", region, "Ontology:", ont, ":", e$message))
            })
            }
        }
}

pathway <- function(motifs){
    for (region in c('3hop','original'))
        for (mtf in motifs) {
            # Check input file
            check_input_files(region, mtf)
            
            tryCatch({
            gene <- read.csv(paste0("output/deg_overexpressed/",region,'/',mtf,".txt"),head=FALSE)

                # Convert gene symbols to Entrez IDs with error handling
                gene_id <- bitr(gene$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
                gene_id <- check_gene_conversion(gene_id, gene$V1)

                # Pathway analysis based on selected database
                Gene_IDs <- mapIds(org_db,
                               keys=gene$V1,
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
                
                if (pathway_db == "Reactome") {
                    egoPath <- enrichPathway(Gene_IDs, organism=tolower(species))
                } else if (pathway_db == "KEGG") {
                    egoPath <- enrichKEGG(gene = Gene_IDs,
                                        organism = org_code,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = qvalue_cutoff)
                } else {
                    stop("Unsupported pathway database. Please choose from: Reactome, KEGG")
                }

            print(paste("Processing Pathway analysis for:", mtf, "Region:", region))
                print(paste("Number of enriched pathways found:", nrow(summary(egoPath))))

                if (nrow(summary(egoPath)) > 0){
                path <- paste0("output/pathway/",region,'/')
                if (!file.exists(path)) {
                    dir.create(path, recursive = TRUE)
                }
                    write.csv(egoPath, file = paste0("output/pathway/",region,'/',mtf,".csv"), row.names = FALSE)
                png(file=paste0("output/pathway/",region,'/',mtf,".png"),width = 10,height = 10,units = "in",res = 600)
                    print(dotplot(egoPath, showCategory=show_category, color="qvalue", label_format = 50))
                dev.off()
                } else {
                    warning(paste("No enriched pathways found for", mtf, "Region:", region))
            }
            }, error = function(e) {
                warning(paste("Error in Pathway analysis for", mtf, "Region:", region, ":", e$message))
            })
        }
}

# Call functions to execute the analysis
go(motifs)
pathway(motifs)

