library(CellChat)
library(patchwork)
library(Matrix)
options(stringsAsFactors = FALSE)

mtf <- 'cccm'

CellChat <- function(m,mtf){
    for(mtf in c('cccm'))
    {
        path <- paste0("output/CellChat(3hop)/")
        if (!file.exists(path)) {
            dir.create(path, recursive = TRUE)
        }

        # Read and validate matrix
        matrix_path <- paste0('data/extracted_matrix/3hop/',mtf,'/matrix.txt')
        if (!file.exists(matrix_path)) {
            stop(paste("Matrix file not found:", matrix_path))
        }
        
        tryCatch({
            # First check if the file is a valid MatrixMarket format
            header <- readLines(matrix_path, n = 1)
            if (!grepl("^%%MatrixMarket", header)) {
                stop("File is not in MatrixMarket format")
            }
            
            # Read dimensions from the second line
            dims <- readLines(matrix_path, n = 2)[2]
            dims <- as.numeric(strsplit(dims, " ")[[1]])
            nrows <- dims[1]
            ncols <- dims[2]
            
            print(paste("Matrix dimensions from header:", nrows, "x", ncols))
            
            # Create a temporary file with corrected indices
            temp_file <- tempfile()
            system(paste("awk 'NR<=2 {print} NR>2 && $1<=", nrows, " {print}'", matrix_path, ">", temp_file))
            
            # Read the corrected matrix
            data.input <- readMM(temp_file)
            unlink(temp_file)  # Clean up temporary file
            
            print("Matrix dimensions after reading:")
            print(dim(data.input))
            
            # Validate matrix dimensions
            if (any(dim(data.input) == 0)) {
                stop("Matrix has zero dimensions")
            }
            
            if (nrow(data.input) != nrows || ncol(data.input) != ncols) {
                stop(paste("Matrix dimensions mismatch. Expected:", nrows, "x", ncols, 
                         "Got:", nrow(data.input), "x", ncol(data.input)))
            }
            
            # Read and validate gene names
            gene_name <- read.csv('data/ad_gene.txt', header=FALSE, nrows=nrows)
            print(paste("Read", nrow(gene_name), "gene names"))
            
            if (nrow(gene_name) != nrow(data.input)) {
                stop(paste("Gene names length (", nrow(gene_name), 
                         ") does not match matrix rows (", nrow(data.input), ")"))
            }

            # Read and validate cell types
            celltype_path <- paste0('data/extracted_matrix/3hop/',mtf,'/cell_types.txt')
            if (!file.exists(celltype_path)) {
                stop(paste("Cell types file not found:", celltype_path))
            }
            
            celltype <- readLines(celltype_path)
            print("Cell types length:")
            print(length(celltype))
            
            if (length(celltype) != ncol(data.input)) {
                stop(paste("Cell types length (", length(celltype), 
                         ") does not match matrix columns (", ncol(data.input), ")"))
            }

            meta <- data.frame(labels = celltype, row.names = 1:length(celltype))
            data.input@Dimnames[[1]] <- gene_name$V1
            data.input@Dimnames[[2]] <- 1:data.input@Dim[2]

            cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

            # use CellChatDB.mouse if running on mouse data
            CellChatDB <- CellChatDB.mouse 
            CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
            cellchat@DB <- CellChatDB.use

            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cellchat <- computeCommunProb(cellchat)

            cellchat_truncatedMean <- computeCommunProb(cellchat,type =  "truncatedMean", trim = 0.1)
            cellchat_truncatedMean <- computeCommunProbPathway(cellchat_truncatedMean)
            cellchat_truncatedMean <- aggregateNet(cellchat_truncatedMean)

            write.csv(cellchat_truncatedMean@net$weight,paste0(path,mtf,'_weight_truncatedMean.csv'))
            write.csv(cellchat_truncatedMean@net$count,paste0(path,mtf,'_count_truncatedMean.csv'))
            
        }, error = function(e) {
            print(paste("Error in processing matrix:", e$message))
            print("Matrix file header:")
            system(paste("head -n 2", matrix_path))
            stop(e)
        })
    }
}

CellChat(mtf)
