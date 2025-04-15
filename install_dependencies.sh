#!/bin/bash

set -e  # Exit if any command fails

echo "Setting executable permissions for all script files..."
chmod +x *.py *.r *.sh
echo "Permissions set successfully!"

echo "🔍 Checking if Conda is installed..."

if ! command -v conda &> /dev/null; then
    echo "Conda not found. Preparing to install Miniconda..."

    INSTALL_DIR="$HOME/miniconda"

    if [ -d "$INSTALL_DIR" ]; then
        echo "Existing Miniconda directory found at $INSTALL_DIR"
        read -p "Do you want to [D]elete and reinstall or [U]pdate existing install? (D/U): " CHOICE

        if [[ "$CHOICE" =~ ^[Dd]$ ]]; then
            echo "🧹 Deleting existing Miniconda directory..."
            rm -rf "$INSTALL_DIR"
        elif [[ "$CHOICE" =~ ^[Uu]$ ]]; then
            echo "⚙️ Updating Miniconda..."
            INSTALL_FLAG="-u"
        else
            echo "Invalid option. Aborting."
            exit 1
        fi
    fi

    # Detect OS and architecture
    OS_TYPE="$(uname -s)"
    ARCH_TYPE="$(uname -m)"

    if [[ "$OS_TYPE" == "Darwin" && "$ARCH_TYPE" == "arm64" ]]; then
        INSTALLER_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    elif [[ "$OS_TYPE" == "Darwin" && "$ARCH_TYPE" == "x86_64" ]]; then
        INSTALLER_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    elif [[ "$OS_TYPE" == "Linux" ]]; then
        INSTALLER_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    else
        echo "Unsupported OS or architecture: $OS_TYPE $ARCH_TYPE"
        exit 1
    fi

    echo "📦 Downloading Miniconda from $INSTALLER_URL"
    curl -o miniconda.sh "$INSTALLER_URL"
    bash miniconda.sh -b $INSTALL_FLAG -p "$INSTALL_DIR"

    echo "⚙️ Initializing Conda..."
    eval "$($INSTALL_DIR/bin/conda shell.bash hook)"
    conda init

    echo "Miniconda installed. Please restart your shell and rerun the script."
    exit 0
else
    echo "Conda is already installed."
    eval "$(conda shell.bash hook)"
fi

# Get conda base directory
CONDA_BASE=$(conda info --base)
ENV_NAME="trimnn"
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"

# Check if we're already in the trimnn environment
if [[ "$CONDA_DEFAULT_ENV" == "$ENV_NAME" ]]; then
    echo "Already in '$ENV_NAME' environment. Skipping activation."
else
    # Environment setup
    echo "📦 Checking for existing '$ENV_NAME' environment..."
    if conda env list | grep -q "^$ENV_NAME "; then
        echo "Environment '$ENV_NAME' already exists. Skipping creation."
    else
        echo "Creating Conda environment '$ENV_NAME' with Python 3.9..."
        conda create -n "$ENV_NAME" python=3.9 -y
    fi

    echo "Activating '$ENV_NAME' environment..."
    conda activate "$ENV_NAME"
fi

# Python dependencies
echo "Installing Python dependencies..."
pip install numpy pandas scipy scikit-learn matplotlib seaborn

# R dependencies
echo "Installing R and Bioconductor dependencies..."
Rscript -e '
# List of required packages
packages <- c(
    "CellChat",
    "gtools",
    "metap",
    "rcompanion",
    "parallel",
    "ggplot2",
    "gridExtra",
    "knitr",
    "rmarkdown"
)

# Function to install packages if not already installed
install_if_missing <- function(package) {
    if (!require(package, character.only = TRUE)) {
        if (package == "CellChat") {
            # CellChat requires BiocManager
            if (!require("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install("CellChat")
        } else {
            install.packages(package)
        }
    }
}

# Install each package
for (package in packages) {
    cat(sprintf("\nInstalling package: %s\n", package))
    install_if_missing(package)
}'

echo "All dependencies installed successfully in the '$ENV_NAME' environment!"

