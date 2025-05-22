library(gtools)
library(metap)
library(rcompanion)

# Function to calculate effect size and p-value for a single dataset
calculate_effect_size_and_p <- function(data, cell_types, p_value_cutoff = 0.05, 
                                      effect_size_threshold = 0.1, 
                                      multiple_testing_correction = "BH") {
  # Create all possible triplets
  triplets <- combinations(length(cell_types), 3, cell_types)
  
  # Initialize results dataframe
  results <- data.frame(
    tri_celltype = character(),
    occurrence_num = integer(),
    rank = integer(),
    meta_p = numeric(),
    effect_size_1 = numeric(),
    effect_size_2 = numeric(),
    effect_size_3 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each triplet
  for (i in 1:nrow(triplets)) {
    triplet <- triplets[i,]
    
    # Create contingency tables for each cell type
    tables <- lapply(triplet, function(cell_type) {
      table(data$label == cell_type, data$top_lev)
    })
    
    # Calculate effect sizes using Cramer's V
    effect_sizes <- sapply(tables, function(tab) {
      if (all(dim(tab) >= 2)) {
        cramerV(tab)
      } else {
        NA
      }
    })
    
    # Calculate p-values using Fisher's exact test
    p_values <- sapply(tables, function(tab) {
      if (all(dim(tab) >= 2)) {
        fisher.test(tab)$p.value
      } else {
        NA
      }
    })
    
    # Combine p-values using meanz method
    if (all(!is.na(p_values))) {
      meta_p <- meanz(p_values)$p
    } else {
      meta_p <- NA
    }
    
    # Count occurrences
    occurrence_num <- sum(data$label %in% triplet)
    
    # Add to results if meets criteria
    if (!is.na(meta_p) && all(!is.na(effect_sizes)) && meta_p < p_value_cutoff && all(effect_sizes >= effect_size_threshold)) {
      results <- rbind(results, data.frame(
        tri_celltype = paste(triplet, collapse = "&"),
        occurrence_num = occurrence_num,
        rank = NA,  # Will be calculated after all triplets are processed
        meta_p = meta_p,
        effect_size_1 = effect_sizes[1],
        effect_size_2 = effect_sizes[2],
        effect_size_3 = effect_sizes[3],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply multiple testing correction if specified
  if (multiple_testing_correction != "none") {
    results$meta_p <- p.adjust(results$meta_p, method = multiple_testing_correction)
  }
  
  # Calculate ranks based on meta_p and effect sizes
  if (nrow(results) > 0) {
    results$rank <- rank(-(results$meta_p * (results$effect_size_1 + 
                                            results$effect_size_2 + 
                                            results$effect_size_3)))
  }
  
  return(results)
}

# Main function to process a dataset
process_dataset <- function(input_file, output_file, cell_types, 
                          p_value_cutoff = 0.05, 
                          effect_size_threshold = 0.1,
                          multiple_testing_correction = "BH") {
  # Read input data
  data <- read.csv(input_file)
  
  # Calculate effect sizes and p-values
  results <- calculate_effect_size_and_p(
    data = data,
    cell_types = cell_types,
    p_value_cutoff = p_value_cutoff,
    effect_size_threshold = effect_size_threshold,
    multiple_testing_correction = multiple_testing_correction
  )
  
  # Write results
  write.csv(results, output_file, row.names = FALSE)
  
  return(results)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
input_file <- NULL
output_file <- NULL
p_value_cutoff <- 0.05
effect_size_threshold <- 0.1
multiple_testing_correction <- "BH"

# Parse arguments
for (i in seq_along(args)) {
  if (args[i] == "--input-file") {
    input_file <- args[i+1]
  } else if (args[i] == "--output-file") {
    output_file <- args[i+1]
  } else if (args[i] == "--pvalue-cutoff") {
    p_value_cutoff <- as.numeric(args[i+1])
  } else if (args[i] == "--effect-size-threshold") {
    effect_size_threshold <- as.numeric(args[i+1])
  } else if (args[i] == "--multiple-testing-correction") {
    multiple_testing_correction <- args[i+1]
  }
}

if (is.null(input_file) || is.null(output_file)) {
  stop("Both --input-file and --output-file must be specified.")
}

# Define cell types (edit as needed for your dataset)
cell_types <- c("Astro", "CA1", "CTX-Ex", "Inh", "Micro", "OPC", "SMC")

# Run the analysis
process_dataset(
  input_file = input_file,
  output_file = output_file,
  cell_types = cell_types,
  p_value_cutoff = p_value_cutoff,
  effect_size_threshold = effect_size_threshold,
  multiple_testing_correction = multiple_testing_correction
)