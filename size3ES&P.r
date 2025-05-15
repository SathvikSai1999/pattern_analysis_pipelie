library(gtools)
library(metap)
library(rcompanion)
library(parallel)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
months <- c(8,13)
replicates <- c(1,2)
control_group <- c('control')
p_value_cutoff <- 0.05
effect_size_threshold <- 0.1
multiple_testing_correction <- "BH"
output_dir <- "output/size3es"

# Parse arguments
for (arg in args) {
    if (grepl("--months=", arg)) {
        months <- as.numeric(strsplit(sub("--months=", "", arg), ",")[[1]])
    } else if (grepl("--replicates=", arg)) {
        replicates <- as.numeric(strsplit(sub("--replicates=", "", arg), ",")[[1]])
    } else if (grepl("--control-group=", arg)) {
        control_group <- sub("--control-group=", "", arg)
    } else if (grepl("--p-value-cutoff=", arg)) {
        p_value_cutoff <- as.numeric(sub("--p-value-cutoff=", "", arg))
    } else if (grepl("--effect-size-threshold=", arg)) {
        effect_size_threshold <- as.numeric(sub("--effect-size-threshold=", "", arg))
    } else if (grepl("--multiple-testing-correction=", arg)) {
        multiple_testing_correction <- sub("--multiple-testing-correction=", "", arg)
    } else if (grepl("--output-dir=", arg)) {
        output_dir <- sub("--output-dir=", "", arg)
    }
}

# Set number of cores for parallel processing
num_cores <- detectCores()
options(mc.cores = num_cores)

# Create result_rank directory if it doesn't exist
if (!dir.exists('result_rank')) {
    dir.create('result_rank')
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

process_replicate <- function(m, c, r) {
    print(sprintf("Processing month %d, replicate %d", m, r))
    
    # Check if input files exist
    node_rank_file <- sprintf('result_rank/AD_%dm_r%d_control_node_rank.csv', m, r)
    tri_rank_file <- sprintf('result_rank/AD_%dm_r%d_control_tri_rank.csv', m, r)
    
    if (!file.exists(node_rank_file)) {
        warning(sprintf("Node rank file not found: %s", node_rank_file))
        return(NULL)
    }
    
    if (!file.exists(tri_rank_file)) {
        warning(sprintf("Tri rank file not found: %s", tri_rank_file))
        return(NULL)
    }
    
    tryCatch({
        data <- read.csv(tri_rank_file, stringsAsFactors = FALSE)
        data <- data[order(data$tri_celltype),]
        
        # Initialize dictionary with ranks
        dictionary <- setNames(data$rank, data$tri_celltype)
        
        # Process data in parallel
        results <- mclapply(data$tri_celltype, function(perm) {
            split_string <- strsplit(perm, "&", fixed = TRUE)[[1]]
            effect_size <- list(NULL, NULL, NULL)
            pvalue <- NULL
            
            for (i in seq_len(length(split_string))) {
                cell = split_string[i]
                if (i==1) {idx <- c(2,3)}
                if (i==2) {idx <- c(3,1)}
                if (i==3) {idx <- c(1,2)}
                cell_1 = split_string[idx[1]]
                cell_2 = split_string[idx[2]]
                
                table_test <- matrix(0, nrow = 2, ncol = 2)
                
                for (tri in data$tri_celltype) {
                    tri_cells <- strsplit(tri, "&", fixed = TRUE)[[1]]
                    
                    # Check if cell is in the triplet
                    if (cell %in% tri_cells) {
                        # Check if cell_1 and cell_2 are in the remaining cells
                        remaining_cells <- tri_cells[!tri_cells %in% cell]
                        if (all(c(cell_1, cell_2) %in% remaining_cells)) {
                            table_test[1,1] <- table_test[1,1] + dictionary[[tri]]
                        }
                    }
                    
                    # Update other cells of the contingency table
                    if (cell_1 %in% tri_cells) {table_test[1,2] <- table_test[1,2] + dictionary[[tri]]}
                    if (cell_2 %in% tri_cells) {table_test[2,1] <- table_test[2,1] + dictionary[[tri]]}
                    table_test[2,2] <- table_test[2,2] + dictionary[[tri]]
                }
                
                # Add small constant to avoid zero cells
                table_test <- table_test + 0.5
                
                pvalue <- c(pvalue, fisher.test(table_test)$p.value)
                effect_size[i] <- rcompanion::cramerV(table_test)
            }
            
            list(pvalue = pvalue, effect_size = effect_size)
        }, mc.cores = num_cores)
        
        # Process results
        all_metap <- sapply(results, function(x) metap::meanz(x$pvalue)$p)
        all_effect_size <- lapply(results, function(x) x$effect_size)
        
        data$meta_p <- unlist(all_metap)
        data$effect_size_1 <- sapply(all_effect_size, function(x) x[[1]])
        data$effect_size_2 <- sapply(all_effect_size, function(x) x[[2]])
        data$effect_size_3 <- sapply(all_effect_size, function(x) x[[3]])
        
        # Apply p-value cutoff and effect size threshold
        data <- data[data$meta_p <= p_value_cutoff & 
                    (data$effect_size_1 >= effect_size_threshold | 
                     data$effect_size_2 >= effect_size_threshold | 
                     data$effect_size_3 >= effect_size_threshold), ]
        
        # Apply multiple testing correction if specified
        if (multiple_testing_correction != "none") {
            data$meta_p <- p.adjust(data$meta_p, method = multiple_testing_correction)
        }
        
        file_path <- file.path(output_dir, sprintf('AD_%d_control_%d_tri_ES&P.csv', m, r))
        write.csv(data, file = file_path, append = FALSE, row.names=FALSE)
        print(sprintf("Completed processing month %d, replicate %d", m, r))
        
    }, error = function(e) {
        warning(sprintf("Error processing files for month %d, replicate %d: %s", m, r, e$message))
    })
}

# Process all combinations in parallel
mclapply(months, function(m) {
    lapply(control_group, function(c) {
        lapply(replicates, function(r) {
            process_replicate(m, c, r)
        })
    })
}, mc.cores = num_cores)