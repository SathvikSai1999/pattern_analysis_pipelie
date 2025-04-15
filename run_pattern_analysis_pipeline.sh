#!/bin/bash

# Set up logging
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to run a step with error handling
run_step() {
    local step_name=$1
    local command=$2
    local log_file="$LOG_DIR/${step_name}_$(date +%Y%m%d_%H%M%S).log"
    
    log "🚀 Starting $step_name..."
    log "Command: $command"
    
    if eval "$command" > "$log_file" 2>&1; then
        log "✅ $step_name completed successfully"
        return 0
    else
        log "❌ Error in $step_name. Check $log_file for details"
        return 1
    fi
}

# 🚀 Ensure correct environment is activated
log "🔁 Activating 'trimnn' environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /Users/sathviksai/miniconda/envs/trimnn

log "==================================="
log "🚀 Starting TrimNN Analysis Pipeline"
log "==================================="

# Set number of parallel jobs (adjust based on your system)
MAX_JOBS=4

# Run steps in parallel where possible
log "🔍 Running parallel steps..."
(
    run_step "extract_motifs" "python extract_motifs.py" &
    run_step "DEG" "Rscript DEG.r" &
    wait
) | tee -a "$LOG_FILE"

# Run dependent steps sequentially
log "🔍 Running sequential steps..."
run_step "extract_deg" "python extract_deg.py"
run_step "go_and_pathway" "Rscript go_and_pathway.r"
run_step "cellchat" "Rscript cellchat.r"
run_step "size3ES_and_P" "Rscript 'size3ES&P.r'"

log "==============================="
log "Pipeline execution completed!"
log "Check $LOG_DIR for detailed logs"
log "==============================="

