#!/bin/bash

# Ensure the script continues running even if terminal closes
# and set up proper environment handling
if [ -z "$NOHUP_ACTIVE" ]; then
    export NOHUP_ACTIVE=1
    LOG_DIR="logs"
    mkdir -p "$LOG_DIR"
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOG_FILE="$LOG_DIR/pipeline_${TIMESTAMP}.log"
    
    # Run the script in background but pipe output to both log file and terminal
    exec "$0" "$@" 2>&1 | tee "$LOG_FILE" &
    echo "Pipeline started in background. Progress will be shown in terminal and saved to $LOG_FILE"
    exit 0
fi

# Default parameter values
HOP_DISTANCE="3"                    # Default hop distance for motif extraction
DEG_P_ADJ_CUTOFF="0.05"            # Default adjusted p-value cutoff for DEG analysis
GO_PVALUE_CUTOFF="0.01"            # Default p-value cutoff for GO analysis
PATHWAY_DB="Reactome"              # Default pathway database
CELLCHAT_TYPE="truncatedMean"      # Default CellChat communication probability type
CELLCHAT_TRIM="0.1"               # Default CellChat trimming factor
CELLCHAT_SEARCH="Secreted Signaling" # Default CellChat search type
GO_QVALUE_CUTOFF="0.05"           # Default q-value cutoff for GO analysis
GO_ONTOLOGY_TYPES="BP,MF,CC"      # Default GO ontology types
GO_SHOW_CATEGORY="20"             # Default number of categories to show in GO plots
SIZE3ES_MONTHS="8,13"             # Default months for size analysis
SIZE3ES_REPLICATES="1,2"          # Default replicates for size analysis
SIZE3ES_CONTROL_GROUP="control"    # Default control group
SIZE3ES_PVALUE_CUTOFF="0.05"      # Default p-value cutoff for size analysis
SIZE3ES_EFFECT_SIZE_THRESHOLD="0.1" # Default effect size threshold
SIZE3ES_MULTIPLE_TESTING_CORRECTION="BH" # Default multiple testing correction method

# Function to validate numeric parameters
validate_numeric() {
    local param_name=$1
    local param_value=$2
    local min_value=$3
    local max_value=$4
    
    if ! [[ "$param_value" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "Error: $param_name must be a numeric value"
        exit 1
    fi
    
    if (( $(echo "$param_value < $min_value" | bc -l) )); then
        echo "Error: $param_name must be greater than or equal to $min_value"
        exit 1
    fi
    
    if (( $(echo "$param_value > $max_value" | bc -l) )); then
        echo "Error: $param_name must be less than or equal to $max_value"
        exit 1
    fi
}

# Function to validate list parameters
validate_list() {
    local param_name=$1
    local param_value=$2
    local valid_values=$3
    
    IFS=',' read -ra values <<< "$param_value"
    for value in "${values[@]}"; do
        if [[ ! " ${valid_values[@]} " =~ " ${value} " ]]; then
            echo "Error: Invalid value '$value' for $param_name. Valid values are: ${valid_values[*]}"
            exit 1
        fi
    done
}

# Function to display usage
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --hop-distance VALUE        Hop distance for motif extraction [default: 3]"
    echo "  --DEG-p-adj-cutoff VALUE   Adjusted p-value cutoff for DEG analysis [default: 0.05]"
    echo "  --go-pvalue-cutoff VALUE   P-value cutoff for GO analysis [default: 0.01]"
    echo "  --pathway-db VALUE         Pathway database selection [default: Reactome]"
    echo "  --cellchat-type TYPE       CellChat communication probability type [default: truncatedMean]"
    echo "  --cellchat-trim VALUE      CellChat trimming factor [default: 0.1]"
    echo "  --cellchat-search VALUE    CellChat search type [default: Secreted Signaling]"
    echo "  --go-qvalue-cutoff VALUE   Q-value cutoff for GO analysis [default: 0.05]"
    echo "  --go-ontology-types LIST   Comma-separated list of GO ontology types [default: BP,MF,CC]"
    echo "  --go-show-category VALUE   Number of categories to show in GO plots [default: 20]"
    echo "  --size3es-months LIST      Comma-separated list of months [default: 8,13]"
    echo "  --size3es-replicates LIST  Comma-separated list of replicates [default: 1,2]"
    echo "  --size3es-control-group VALUE Control group name [default: control]"
    echo "  --size3es-pvalue-cutoff VALUE P-value cutoff [default: 0.05]"
    echo "  --size3es-effect-size-threshold VALUE Effect size threshold [default: 0.1]"
    echo "  --size3es-multiple-testing-correction VALUE Multiple testing correction method [default: BH]"
    echo "  -h, --help                 Display this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --hop-distance)
            HOP_DISTANCE="$2"
            validate_numeric "hop-distance" "$HOP_DISTANCE" 1 10
            shift 2
            ;;
        --DEG-p-adj-cutoff)
            DEG_P_ADJ_CUTOFF="$2"
            validate_numeric "DEG-p-adj-cutoff" "$DEG_P_ADJ_CUTOFF" 0 1
            shift 2
            ;;
        --go-pvalue-cutoff)
            GO_PVALUE_CUTOFF="$2"
            validate_numeric "go-pvalue-cutoff" "$GO_PVALUE_CUTOFF" 0 1
            shift 2
            ;;
        --pathway-db)
            PATHWAY_DB="$2"
            validate_list "pathway-db" "$PATHWAY_DB" "Reactome KEGG GO"
            shift 2
            ;;
        --cellchat-type)
            CELLCHAT_TYPE="$2"
            validate_list "cellchat-type" "$CELLCHAT_TYPE" "truncatedMean mean median"
            shift 2
            ;;
        --cellchat-trim)
            CELLCHAT_TRIM="$2"
            validate_numeric "cellchat-trim" "$CELLCHAT_TRIM" 0 1
            shift 2
            ;;
        --cellchat-search)
            CELLCHAT_SEARCH="$2"
            validate_list "cellchat-search" "$CELLCHAT_SEARCH" "Secreted Signaling ECM-Receptor Cell-Cell Contact"
            shift 2
            ;;
        --go-qvalue-cutoff)
            GO_QVALUE_CUTOFF="$2"
            validate_numeric "go-qvalue-cutoff" "$GO_QVALUE_CUTOFF" 0 1
            shift 2
            ;;
        --go-ontology-types)
            GO_ONTOLOGY_TYPES="$2"
            validate_list "go-ontology-types" "$GO_ONTOLOGY_TYPES" "BP MF CC"
            shift 2
            ;;
        --go-show-category)
            GO_SHOW_CATEGORY="$2"
            validate_numeric "go-show-category" "$GO_SHOW_CATEGORY" 1 100
            shift 2
            ;;
        --size3es-months)
            SIZE3ES_MONTHS="$2"
            validate_list "size3es-months" "$SIZE3ES_MONTHS" "8 13"
            shift 2
            ;;
        --size3es-replicates)
            SIZE3ES_REPLICATES="$2"
            validate_list "size3es-replicates" "$SIZE3ES_REPLICATES" "1 2"
            shift 2
            ;;
        --size3es-control-group)
            SIZE3ES_CONTROL_GROUP="$2"
            validate_list "size3es-control-group" "$SIZE3ES_CONTROL_GROUP" "control"
            shift 2
            ;;
        --size3es-pvalue-cutoff)
            SIZE3ES_PVALUE_CUTOFF="$2"
            validate_numeric "size3es-pvalue-cutoff" "$SIZE3ES_PVALUE_CUTOFF" 0 1
            shift 2
            ;;
        --size3es-effect-size-threshold)
            SIZE3ES_EFFECT_SIZE_THRESHOLD="$2"
            validate_numeric "size3es-effect-size-threshold" "$SIZE3ES_EFFECT_SIZE_THRESHOLD" 0 1
            shift 2
            ;;
        --size3es-multiple-testing-correction)
            SIZE3ES_MULTIPLE_TESTING_CORRECTION="$2"
            validate_list "size3es-multiple-testing-correction" "$SIZE3ES_MULTIPLE_TESTING_CORRECTION" "BH bonferroni holm hochberg none"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Ensure conda environment activation works in detached state
source "$(conda info --base)/etc/profile.d/conda.sh"

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
    
    log "Starting $step_name..."
    log "Command: $command"
    
    if eval "$command" > "$log_file" 2>&1; then
        log "$step_name completed successfully"
        return 0
    else
        log "Error in $step_name. Check $log_file for details"
        return 1
    fi
}

# Ensure correct environment is activated
log "Activating 'trimnn' environment..."
conda activate trimnn

# Get the Python interpreter from the trimnn environment
PYTHON_INTERPRETER=$(conda run -n trimnn which python)
log "Using Python interpreter: $PYTHON_INTERPRETER"

log "==================================="
log "Starting TrimNN Analysis Pipeline"
log "==================================="

# Set number of parallel jobs (adjust based on your system)
MAX_JOBS=4

# Build command strings with parameters
MOTIF_CMD="$PYTHON_INTERPRETER extract_motifs.py"
if [ ! -z "$HOP_DISTANCE" ]; then
    MOTIF_CMD="$MOTIF_CMD --hop-distance=$HOP_DISTANCE"
fi

DEG_CMD="Rscript DEG.r"
if [ ! -z "$DEG_P_ADJ_CUTOFF" ]; then
    DEG_CMD="$DEG_CMD --p-adj-cutoff=$DEG_P_ADJ_CUTOFF"
fi

GO_CMD="Rscript go_and_pathway.r"
if [ ! -z "$GO_PVALUE_CUTOFF" ]; then
    GO_CMD="$GO_CMD --pvalue-cutoff=$GO_PVALUE_CUTOFF"
fi
if [ ! -z "$PATHWAY_DB" ]; then
    GO_CMD="$GO_CMD --pathway-db=$PATHWAY_DB"
fi
if [ ! -z "$GO_QVALUE_CUTOFF" ]; then
    GO_CMD="$GO_CMD --qvalue-cutoff=$GO_QVALUE_CUTOFF"
fi
if [ ! -z "$GO_ONTOLOGY_TYPES" ]; then
    GO_CMD="$GO_CMD --ontology-types=$GO_ONTOLOGY_TYPES"
fi
if [ ! -z "$GO_SHOW_CATEGORY" ]; then
    GO_CMD="$GO_CMD --show-category=$GO_SHOW_CATEGORY"
fi

CELLCHAT_CMD="Rscript cellchat.r"
if [ ! -z "$CELLCHAT_TYPE" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --type=$CELLCHAT_TYPE"
fi
if [ ! -z "$CELLCHAT_TRIM" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --trim=$CELLCHAT_TRIM"
fi
if [ ! -z "$CELLCHAT_SEARCH" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --search=$CELLCHAT_SEARCH"
fi

SIZE3ES_CMD="Rscript 'size3ES&P.r'"
if [ ! -z "$SIZE3ES_MONTHS" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --months=$SIZE3ES_MONTHS"
fi
if [ ! -z "$SIZE3ES_REPLICATES" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --replicates=$SIZE3ES_REPLICATES"
fi
if [ ! -z "$SIZE3ES_CONTROL_GROUP" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --control-group=$SIZE3ES_CONTROL_GROUP"
fi
if [ ! -z "$SIZE3ES_PVALUE_CUTOFF" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --p-value-cutoff=$SIZE3ES_PVALUE_CUTOFF"
fi
if [ ! -z "$SIZE3ES_EFFECT_SIZE_THRESHOLD" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --effect-size-threshold=$SIZE3ES_EFFECT_SIZE_THRESHOLD"
fi
if [ ! -z "$SIZE3ES_MULTIPLE_TESTING_CORRECTION" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --multiple-testing-correction=$SIZE3ES_MULTIPLE_TESTING_CORRECTION"
fi

# Run steps in parallel where possible
log "Running parallel steps..."
(
    run_step "extract_motifs" "$MOTIF_CMD" &
    run_step "DEG" "$DEG_CMD" &
    wait
) | tee -a "$LOG_FILE"

# Run dependent steps sequentially
log "Running sequential steps..."
run_step "extract_deg" "python extract_deg.py"
run_step "go_and_pathway" "$GO_CMD"
run_step "cellchat" "$CELLCHAT_CMD"

# Generate rank files after go_and_pathway and cellchat, but before size3ES&P
log "Generating rank files..."
run_step "generate_rank_files" "$PYTHON_INTERPRETER generate_rank_files.py"

run_step "size3ES_and_P" "$SIZE3ES_CMD"

log "==============================="
log "Pipeline execution completed!"
log "Check $LOG_DIR for detailed logs"
log "==============================="

