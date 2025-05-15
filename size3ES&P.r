library(gtools)
library(metap)
library(rcompanion)
library(parallel)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
timepoints <- c(8, 13)
control_groups <- c('control_8', 'control_13')
sample_groups <- c('sample_8', 'sample_13')
replicates <- c(1, 2)
p_value_cutoff <- 0.05
effect_size_threshold <- 0.1
multiple_testing_correction <- "BH"
output_dir <- "output/size3es"

# Parse arguments
for (arg in args) {
    if (grepl("--size3es-timepoints=", arg)) {
        timepoints <- as.numeric(strsplit(sub("--size3es-timepoints=", "", arg), ",")[[1]])
    } else if (grepl("--size3es-control-groups=", arg)) {
        control_groups <- strsplit(sub("--size3es-control-groups=", "", arg), ",")[[1]]
    } else if (grepl("--size3es-sample-groups=", arg)) {
        sample_groups <- strsplit(sub("--size3es-sample-groups=", "", arg), ",")[[1]]
    } else if (grepl("--size3es-replicates=", arg)) {
        replicates <- as.numeric(strsplit(sub("--size3es-replicates=", "", arg), ",")[[1]])
    } else if (grepl("--size3es-pvalue-cutoff=", arg)) {
        p_value_cutoff <- as.numeric(sub("--size3es-pvalue-cutoff=", "", arg))
    } else if (grepl("--size3es-effect-size-threshold=", arg)) {
        effect_size_threshold <- as.numeric(sub("--size3es-effect-size-threshold=", "", arg))
    } else if (grepl("--size3es-multiple-testing-correction=", arg)) {
        multiple_testing_correction <- sub("--size3es-multiple-testing-correction=", "", arg)
    } else if (grepl("--size3es-output-dir=", arg)) {
        output_dir <- sub("--size3es-output-dir=", "", arg)
    }
}

# Validate input parameters
validate_parameters <- function() {
    # Check timepoints
    if (!all(timepoints %in% c(8, 13))) {
        stop("Invalid timepoints. Must be 8 and/or 13.")
    }
    
    # Check control groups
    expected_control_groups <- paste0("control_", timepoints)
    if (!all(control_groups %in% expected_control_groups)) {
        stop("Invalid control groups. Must match timepoints.")
    }
    
    # Check sample groups
    expected_sample_groups <- paste0("sample_", timepoints)
    if (!all(sample_groups %in% expected_sample_groups)) {
        stop("Invalid sample groups. Must match timepoints.")
    }
    
    # Check replicates
    if (!all(replicates %in% c(1, 2))) {
        stop("Invalid replicates. Must be 1 and/or 2.")
    }
    
    # Check p-value cutoff
    if (p_value_cutoff < 0 || p_value_cutoff > 1) {
        stop("Invalid p-value cutoff. Must be between 0 and 1.")
    }
    
    # Check effect size threshold
    if (effect_size_threshold < 0 || effect_size_threshold > 1) {
        stop("Invalid effect size threshold. Must be between 0 and 1.")
    }
    
    # Check multiple testing correction
    valid_corrections <- c("BH", "bonferroni", "holm", "hochberg", "none")
    if (!multiple_testing_correction %in% valid_corrections) {
        stop("Invalid multiple testing correction method.")
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

process_replicate <- function(timepoint, group, replicate) {
    print(sprintf("Processing timepoint %d, group %s, replicate %d", timepoint, group, replicate))
    
    # Check if input files exist
    node_rank_file <- sprintf('result_rank/%s_%d_control_node_rank.csv', group, replicate)
    tri_rank_file <- sprintf('result_rank/%s_%d_control_tri_rank.csv', group, replicate)
    
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
        
        # Create output filename with new naming convention
        output_file <- sprintf('%s_%d_control_%d_tri_ES&P.csv', 
                             ifelse(grepl("control", group), "control", "sample"),
                             timepoint, 
                             replicate)
        file_path <- file.path(output_dir, output_file)
        
        write.csv(data, file = file_path, append = FALSE, row.names=FALSE)
        print(sprintf("Completed processing timepoint %d, group %s, replicate %d", 
                     timepoint, group, replicate))
        
    }, error = function(e) {
        warning(sprintf("Error processing files for timepoint %d, group %s, replicate %d: %s", 
                       timepoint, group, replicate, e$message))
    })
}

# Validate parameters before processing
validate_parameters()

# Process all combinations in parallel
mclapply(timepoints, function(tp) {
    # Process control groups
    lapply(control_groups[grepl(paste0("_", tp), control_groups)], function(cg) {
        lapply(replicates, function(r) {
            process_replicate(tp, cg, r)
        })
    })
    
    # Process sample groups
    lapply(sample_groups[grepl(paste0("_", tp), sample_groups)], function(sg) {
        lapply(replicates, function(r) {
            process_replicate(tp, sg, r)
        })
    })
}, mc.cores = num_cores)