#!/bin/bash

# Default parameter values
RUN_ONLY="all"                    # Default to running all steps
VERSION="1.0.0"                   # Pipeline version

# Input directories
MOTIF_INPUT_DIR="data/motifs"     # Input directory for motif extraction
DEG_INPUT_DIR="data/DEG"          # Input directory for DEG analysis
GO_INPUT_DIR="data/GO"            # Input directory for GO analysis
CELLCHAT_INPUT_DIR="data/cellchat" # Input directory for CellChat
SIZE3ES_INPUT_DIR="data/size3es"  # Input directory for Size3ES

# Output directories
MOTIF_OUTPUT_DIR="output/motifs"     # Output directory for motifs
DEG_OUTPUT_DIR="output/DEG"          # Output directory for DEG
GO_OUTPUT_DIR="output/GO"            # Output directory for GO
CELLCHAT_OUTPUT_DIR="output/cellchat" # Output directory for CellChat
SIZE3ES_OUTPUT_DIR="output/size3es"   # Output directory for Size3ES

DEG_P_ADJ_CUTOFF="0.05"            # Default adjusted p-value cutoff for DEG analysis
GO_PVALUE_CUTOFF="0.01"            # Default p-value cutoff for GO analysis
PATHWAY_DB="Reactome"              # Default pathway database
CELLCHAT_TYPE="truncatedMean"      # Default CellChat communication probability type
CELLCHAT_TRIM="0.1"               # Default CellChat trimming factor
CELLCHAT_SEARCH="Secreted Signaling" # Default CellChat search type
GO_QVALUE_CUTOFF="0.05"           # Default q-value cutoff for GO analysis
GO_SHOW_CATEGORY="20"             # Default number of categories to show in GO plots
SIZE3ES_TIMEPOINTS="8,13"             # Default timepoints for size analysis
SIZE3ES_CONTROL_GROUPS="control_8,control_13" # Default control groups
SIZE3ES_SAMPLE_GROUPS="sample_8,sample_13"   # Default sample groups
SIZE3ES_REPLICATES="1,2"              # Default replicates for size analysis
SIZE3ES_PVALUE_CUTOFF="0.05"          # Default p-value cutoff for size analysis
SIZE3ES_EFFECT_SIZE_THRESHOLD="0.1"    # Default effect size threshold
SIZE3ES_MULTIPLE_TESTING_CORRECTION="BH" # Default multiple testing correction method
DEG_METHOD="both"                # Default DEG method

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

# Function to display usage with detailed examples
usage() {
    echo "Pattern Analysis Pipeline - Advanced Usage Guide"
    echo ""
    echo "Basic Usage:"
    echo "  $0 [options]"
    echo ""
    echo "Examples:"
    echo "1. Basic Run with Default Parameters:"
    echo "    $0"
    echo ""
    echo "2. Run Specific Analysis Steps:"
    echo "    $0 --run-only deg,go"
    echo ""
    echo "3. Full Pipeline with Custom Parameters:"
    echo "    $0 \\"
    echo "      --motif-input-dir /path/to/motif/input \\"
    echo "      --motif-output-dir /path/to/motif/output \\"
    echo "      --deg-input-dir /path/to/deg/input \\"
    echo "      --deg-output-dir /path/to/deg/output \\"
    echo "      --go-input-dir /path/to/go/input \\"
    echo "      --go-output-dir /path/to/go/output \\"
    echo "      --cellchat-input-dir /path/to/cellchat/input \\"
    echo "      --cellchat-output-dir /path/to/cellchat/output \\"
    echo "      --size3es-input-dir /path/to/size3es/input \\"
    echo "      --size3es-output-dir /path/to/size3es/output"
    echo ""
    echo "4. Advanced DEG Analysis:"
    echo "    $0 --run-only deg \\"
    echo "      --deg-input-dir /path/to/deg/input \\"
    echo "      --deg-output-dir /path/to/deg/output \\"
    echo "      --DEG-method deseq2 \\"
    echo "      --DEG-p-adj-cutoff 0.01 \\"
    echo "      --log2fc-cutoff 1.5"
    echo ""
    echo "5. Advanced GO Analysis:"
    echo "    $0 --run-only go \\"
    echo "      --go-input-dir /path/to/go/input \\"
    echo "      --go-output-dir /path/to/go/output \\"
    echo "      --go-pvalue-cutoff 0.005 \\"
    echo "      --pathway-db KEGG \\"
    echo "      --go-qvalue-cutoff 0.01 \\"
    echo "      --go-show-category 30"
    echo ""
    echo "6. Advanced CellChat Analysis:"
    echo "    $0 --run-only cellchat \\"
    echo "      --cellchat-input-dir /path/to/cellchat/input \\"
    echo "      --cellchat-output-dir /path/to/cellchat/output \\"
    echo "      --cellchat-type mean \\"
    echo "      --cellchat-trim 0.2 \\"
    echo "      --cellchat-search 'ECM-Receptor'"
    echo ""
    echo "7. Advanced Size3ES Analysis:"
    echo "    $0 --run-only size3es \\"
    echo "      --size3es-input-dir /path/to/size3es/input \\"
    echo "      --size3es-output-dir /path/to/size3es/output \\"
    echo "      --size3es-timepoints '8,13' \\"
    echo "      --size3es-control-groups 'control_8,control_13' \\"
    echo "      --size3es-sample-groups 'sample_8,sample_13' \\"
    echo "      --size3es-replicates '1,2' \\"
    echo "      --size3es-pvalue-cutoff 0.01 \\"
    echo "      --size3es-effect-size-threshold 0.15 \\"
    echo "      --size3es-multiple-testing-correction bonferroni"
    echo ""
    echo "Available Options:"
    echo ""
    echo "General Options:"
    echo "  --run-only VALUE           Specify which analyses to run [default: all]"
    echo "                             Options: all, deg, go, cellchat, size3es"
    echo "                             Multiple values allowed, comma-separated"
    echo "  --version                  Display pipeline version"
    echo ""
    echo "Directory Options:"
    echo "  --motif-input-dir PATH     Input directory for motif extraction"
    echo "  --motif-output-dir PATH    Output directory for motif extraction"
    echo "  --deg-input-dir PATH       Input directory for DEG analysis"
    echo "  --deg-output-dir PATH      Output directory for DEG analysis"
    echo "  --go-input-dir PATH        Input directory for GO analysis"
    echo "  --go-output-dir PATH       Output directory for GO analysis"
    echo "  --cellchat-input-dir PATH  Input directory for CellChat analysis"
    echo "  --cellchat-output-dir PATH Output directory for CellChat analysis"
    echo "  --size3es-input-dir PATH   Input directory for Size3ES analysis"
    echo "  --size3es-output-dir PATH  Output directory for Size3ES analysis"
    echo ""
    echo "DEG Analysis Options:"
    echo "  --DEG-method VALUE         Analysis method [default: both]"
    echo "                             Options: deseq2, wilcox, both"
    echo "  --DEG-p-adj-cutoff VALUE   Adjusted p-value cutoff [default: 0.05]"
    echo "  --log2fc-cutoff VALUE      Log2 fold change cutoff [default: 1.0]"
    echo ""
    echo "GO Analysis Options:"
    echo "  --go-pvalue-cutoff VALUE   P-value cutoff [default: 0.01]"
    echo "  --pathway-db VALUE         Pathway database [default: Reactome]"
    echo "                             Options: Reactome, KEGG, GO"
    echo "  --go-qvalue-cutoff VALUE   Q-value cutoff [default: 0.05]"
    echo "  --go-show-category VALUE   Number of categories to show [default: 20]"
    echo ""
    echo "CellChat Options:"
    echo "  --cellchat-type VALUE      Communication probability type [default: truncatedMean]"
    echo "                             Options: truncatedMean, mean, median"
    echo "  --cellchat-trim VALUE      Trimming factor [default: 0.1]"
    echo "  --cellchat-search VALUE    Search type [default: Secreted Signaling]"
    echo "                             Options: Secreted Signaling, ECM-Receptor, Cell-Cell Contact"
    echo ""
    echo "Size3ES Options:"
    echo "  --size3es-timepoints VALUE Comma-separated list of timepoints [default: 8,13]"
    echo "  --size3es-control-groups VALUE Comma-separated list of control groups [default: control_8,control_13]"
    echo "  --size3es-sample-groups VALUE Comma-separated list of sample groups [default: sample_8,sample_13]"
    echo "  --size3es-replicates VALUE Comma-separated list of replicates [default: 1,2]"
    echo "  --size3es-pvalue-cutoff VALUE P-value cutoff [default: 0.05]"
    echo "  --size3es-effect-size-threshold VALUE Effect size threshold [default: 0.1]"
    echo "  --size3es-multiple-testing-correction VALUE"
    echo "                             Multiple testing correction method [default: BH]"
    echo "                             Options: BH, bonferroni, holm, hochberg, none"
    echo ""
    echo "For detailed documentation, please refer to the README.md file."
    exit 1
}

# Function to display version
show_version() {
    echo "Pattern Analysis Pipeline version $VERSION"
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --run-only)
            RUN_ONLY="$2"
            validate_list "run-only" "$RUN_ONLY" "all deg go cellchat size3es"
            shift 2
            ;;
        --motif-input-dir)
            MOTIF_INPUT_DIR="$2"
            shift 2
            ;;
        --motif-output-dir)
            MOTIF_OUTPUT_DIR="$2"
            shift 2
            ;;
        --deg-input-dir)
            DEG_INPUT_DIR="$2"
            shift 2
            ;;
        --deg-output-dir)
            DEG_OUTPUT_DIR="$2"
            shift 2
            ;;
        --go-input-dir)
            GO_INPUT_DIR="$2"
            shift 2
            ;;
        --go-output-dir)
            GO_OUTPUT_DIR="$2"
            shift 2
            ;;
        --cellchat-input-dir)
            CELLCHAT_INPUT_DIR="$2"
            shift 2
            ;;
        --cellchat-output-dir)
            CELLCHAT_OUTPUT_DIR="$2"
            shift 2
            ;;
        --size3es-input-dir)
            SIZE3ES_INPUT_DIR="$2"
            shift 2
            ;;
        --size3es-output-dir)
            SIZE3ES_OUTPUT_DIR="$2"
            shift 2
            ;;
        --DEG-p-adj-cutoff)
            DEG_P_ADJ_CUTOFF="$2"
            validate_numeric "DEG-p-adj-cutoff" "$DEG_P_ADJ_CUTOFF" 0 1
            shift 2
            ;;
        --DEG-method)
            DEG_METHOD="$2"
            validate_list "DEG-method" "$DEG_METHOD" "deseq2 wilcox both"
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
        --go-show-category)
            GO_SHOW_CATEGORY="$2"
            validate_numeric "go-show-category" "$GO_SHOW_CATEGORY" 1 100
            shift 2
            ;;
        --size3es-timepoints)
            SIZE3ES_TIMEPOINTS="$2"
            shift 2
            ;;
        --size3es-control-groups)
            SIZE3ES_CONTROL_GROUPS="$2"
            shift 2
            ;;
        --size3es-sample-groups)
            SIZE3ES_SAMPLE_GROUPS="$2"
            shift 2
            ;;
        --size3es-replicates)
            SIZE3ES_REPLICATES="$2"
            shift 2
            ;;
        --size3es-pvalue-cutoff)
            SIZE3ES_PVALUE_CUTOFF="$2"
            shift 2
            ;;
        --size3es-effect-size-threshold)
            SIZE3ES_EFFECT_SIZE_THRESHOLD="$2"
            shift 2
            ;;
        --size3es-multiple-testing-correction)
            SIZE3ES_MULTIPLE_TESTING_CORRECTION="$2"
            shift 2
            ;;
        --version)
            show_version
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

# Ensure correct environment is activated
echo "Activating 'trimnn' environment..."
conda activate trimnn

# Get the Python interpreter from the trimnn environment
PYTHON_INTERPRETER=$(conda run -n trimnn which python)
echo "Using Python interpreter: $PYTHON_INTERPRETER"

echo "==================================="
echo "Starting TrimNN Analysis Pipeline"
echo "==================================="

# Build command strings with parameters
MOTIF_CMD="$PYTHON_INTERPRETER extract_motifs.py --input-dir $MOTIF_INPUT_DIR --output-dir $MOTIF_OUTPUT_DIR"

DEG_CMD="Rscript DEG_wrapper.r --DEG-method ${DEG_METHOD:-both} --DEG-p-adj-cutoff ${DEG_P_ADJ_CUTOFF:-0.01} --deg-input-dir ${DEG_INPUT_DIR} --output_dir ${DEG_OUTPUT_DIR}"

GO_CMD="Rscript go_and_pathway.r --input-dir $GO_INPUT_DIR --output-dir $GO_OUTPUT_DIR"
if [ ! -z "$GO_PVALUE_CUTOFF" ]; then
    GO_CMD="$GO_CMD --pvalue-cutoff=$GO_PVALUE_CUTOFF"
fi
if [ ! -z "$PATHWAY_DB" ]; then
    GO_CMD="$GO_CMD --pathway-db=$PATHWAY_DB"
fi
if [ ! -z "$GO_QVALUE_CUTOFF" ]; then
    GO_CMD="$GO_CMD --qvalue-cutoff=$GO_QVALUE_CUTOFF"
fi
if [ ! -z "$GO_SHOW_CATEGORY" ]; then
    GO_CMD="$GO_CMD --show-category=$GO_SHOW_CATEGORY"
fi

CELLCHAT_CMD="Rscript cellchat.r --input-dir $CELLCHAT_INPUT_DIR --output-dir $CELLCHAT_OUTPUT_DIR"
if [ ! -z "$CELLCHAT_TYPE" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --type=$CELLCHAT_TYPE"
fi
if [ ! -z "$CELLCHAT_TRIM" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --trim=$CELLCHAT_TRIM"
fi
if [ ! -z "$CELLCHAT_SEARCH" ]; then
    CELLCHAT_CMD="$CELLCHAT_CMD --search=$CELLCHAT_SEARCH"
fi

SIZE3ES_CMD="Rscript 'size3ES&P.r' --input-dir $SIZE3ES_INPUT_DIR --output-dir $SIZE3ES_OUTPUT_DIR"
if [ ! -z "$SIZE3ES_TIMEPOINTS" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-timepoints=$SIZE3ES_TIMEPOINTS"
fi
if [ ! -z "$SIZE3ES_CONTROL_GROUPS" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-control-groups=$SIZE3ES_CONTROL_GROUPS"
fi
if [ ! -z "$SIZE3ES_SAMPLE_GROUPS" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-sample-groups=$SIZE3ES_SAMPLE_GROUPS"
fi
if [ ! -z "$SIZE3ES_REPLICATES" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-replicates=$SIZE3ES_REPLICATES"
fi
if [ ! -z "$SIZE3ES_PVALUE_CUTOFF" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-pvalue-cutoff=$SIZE3ES_PVALUE_CUTOFF"
fi
if [ ! -z "$SIZE3ES_EFFECT_SIZE_THRESHOLD" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-effect-size-threshold=$SIZE3ES_EFFECT_SIZE_THRESHOLD"
fi
if [ ! -z "$SIZE3ES_MULTIPLE_TESTING_CORRECTION" ]; then
    SIZE3ES_CMD="$SIZE3ES_CMD --size3es-multiple-testing-correction=$SIZE3ES_MULTIPLE_TESTING_CORRECTION"
fi

# Main pipeline execution
echo "Starting pipeline execution..."

# Run DEG analysis if requested
if [[ "$RUN_ONLY" == "all" || "$RUN_ONLY" == "deg" ]]; then
    echo "Running DEG analysis..."
    if ! eval "$DEG_CMD"; then
        echo "DEG analysis failed."
        exit 1
    fi
fi

# Run GO analysis if requested
if [[ "$RUN_ONLY" == "all" || "$RUN_ONLY" == "go" ]]; then
    echo "Running GO analysis..."
    if ! eval "$GO_CMD"; then
        echo "GO analysis failed."
        exit 1
    fi
fi

# Run CellChat analysis if requested
if [[ "$RUN_ONLY" == "all" || "$RUN_ONLY" == "cellchat" ]]; then
    echo "Running CellChat analysis..."
    if ! eval "$CELLCHAT_CMD"; then
        echo "CellChat analysis failed."
        exit 1
    fi
fi

# Run Size3ES analysis if requested
if [[ "$RUN_ONLY" == "all" || "$RUN_ONLY" == "size3es" ]]; then
    echo "Running Size3ES analysis..."
    if ! eval "$SIZE3ES_CMD"; then
        echo "Size3ES analysis failed."
        exit 1
    fi
fi

echo "==============================="
echo "Pipeline execution completed!"
echo "==============================="

