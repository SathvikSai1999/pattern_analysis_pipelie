#!/usr/bin/env Rscript

# Load required libraries
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("--DEG-method", "-m"), type="character", default="both",
              help="DEG analysis method: 'deseq2', 'wilcox', or 'both' [default: %default]"),
  make_option(c("--DEG-p-adj-cutoff", "-p"), type="double", default=0.05,
              help="Adjusted p-value cutoff [default: %default]"),
  make_option(c("--log2fc-cutoff", "-l"), type="double", default=1.0,
              help="Log2 fold change cutoff [default: %default]"),
  make_option(c("--deg-input-dir", "-i"), type="character", default="data",
              help="Input directory containing data files [default: %default]"),
  make_option(c("--output_dir", "-o"), type="character", default="output/DEG",
              help="Output directory for results [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Print debug information
print("Debug: Command line arguments:")
print(opt)
print(paste("Debug: Working directory:", getwd()))
print(paste("Debug: Output directory:", opt$output_dir))

# Function to run DEG analysis
run_deg_analysis <- function(method) {
  tryCatch({
    if (method == "deseq2") {
      system2("Rscript", c("DEG_deseq2.r", opt$"deg-input-dir", opt$output_dir))
    } else if (method == "wilcox") {
      system2("Rscript", c("DEG_wilcox.r", opt$"deg-input-dir", opt$output_dir))
    }
  }, error = function(e) {
    message(sprintf("Error in %s analysis: %s", method, e$message))
    return(FALSE)
  })
  return(TRUE)
}

# Main execution
main <- function() {
  # Create output directory if it doesn't exist
  print(paste("Debug: Creating directory:", opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = TRUE)
  
  # Run selected method(s)
  if (opt$"DEG-method" == "both") {
    message("Running both DESeq2 and Wilcoxon analyses...")
    deseq2_success <- run_deg_analysis("deseq2")
    wilcox_success <- run_deg_analysis("wilcox")
    
    if (!deseq2_success || !wilcox_success) {
      message("Warning: One or more analyses failed. Check the logs for details.")
    }
  } else {
    message(sprintf("Running %s analysis...", opt$"DEG-method"))
    success <- run_deg_analysis(opt$"DEG-method")
    if (!success) {
      stop(sprintf("%s analysis failed. Check the logs for details.", opt$"DEG-method"))
    }
  }
  
  message("DEG analysis completed successfully!")
}

# Run main function
main() 