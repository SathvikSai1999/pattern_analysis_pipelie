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
echo "📦 Installing Python dependencies in '$ENV_NAME' environment..."
conda install -n "$ENV_NAME" -y \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    python-dateutil \
    pytz \
    tzdata

# Verify Python package installation
echo "Verifying Python package installation..."
conda run -n "$ENV_NAME" python -c "
import sys
print(f'Python version: {sys.version}')
try:
    import numpy
    print(f'Numpy version: {numpy.__version__}')
except ImportError:
    print('Numpy not installed')
try:
    import pandas
    print(f'Pandas version: {pandas.__version__}')
except ImportError:
    print('Pandas not installed')
try:
    import scipy
    print(f'Scipy version: {scipy.__version__}')
except ImportError:
    print('Scipy not installed')
try:
    import sklearn
    print(f'Scikit-learn version: {sklearn.__version__}')
except ImportError:
    print('Scikit-learn not installed')
try:
    import matplotlib
    print(f'Matplotlib version: {matplotlib.__version__}')
except ImportError:
    print('Matplotlib not installed')
try:
    import seaborn
    print(f'Seaborn version: {seaborn.__version__}')
except ImportError:
    print('Seaborn not installed')
"

# R dependencies
echo "📦 Installing R and Bioconductor dependencies..."
Rscript -e '
# Set CRAN mirror explicitly
repos <- c(CRAN = "https://cloud.r-project.org/")
options(repos = repos)

# Function to safely install packages
install_package <- function(package_name) {
    tryCatch({
        if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
            cat(sprintf("\nInstalling package: %s\n", package_name))
            install.packages(package_name, repos = "https://cloud.r-project.org/", dependencies = TRUE)
        } else {
            cat(sprintf("\nPackage already installed: %s\n", package_name))
        }
    }, error = function(e) {
        cat(sprintf("\nError installing %s: %s\n", package_name, e$message))
    })
}

# Install basic CRAN packages first
cran_packages <- c(
    "gtools",
    "metap",
    "rcompanion",
    "parallel",
    "ggplot2",
    "gridExtra",
    "knitr",
    "rmarkdown",
    "devtools",  # Required for GitHub installation
    "openxlsx",  # Added for go_and_pathway.r
    "ComplexUpset", # Added for go_and_pathway.r
    "ggvenn"     # Added for go_and_pathway.r
)

# Install CRAN packages
for (pkg in cran_packages) {
    install_package(pkg)
}

# Make sure devtools is available
if (!require("devtools", quietly = TRUE)) {
    cat("\nInstalling devtools package which is required for GitHub installations...\n")
    install.packages("devtools", repos = "https://cloud.r-project.org/", dependencies = TRUE)
}

# Install BioConductor dependencies
cat("\nInstalling Bioconductor dependencies...\n")
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
# Install required Bioconductor packages
BiocManager::install(c(
    "multtest",       # Required by metap
    "org.Mm.eg.db",   # Mouse genome annotations
    "org.Hs.eg.db",   # Human genome annotations 
    "ReactomePA",     # For pathway analysis
    "ComplexHeatmap"  # For heatmap visualization
), update = FALSE, ask = FALSE)

# Install CellChat directly from GitHub since it has compatibility issues with Bioconductor
cat("\nInstalling CellChat directly from GitHub...\n")
tryCatch({
    if (!require("CellChat", quietly = TRUE)) {
        # Load devtools
        library(devtools)
        
        # Install dependencies first
        if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        
        BiocManager::install(c("ComplexHeatmap", "circlize"), update = FALSE, ask = FALSE)
        
        # Install CellChat from GitHub
        devtools::install_github("sqjin/CellChat")
        
        # Check if installation succeeded
        if (require("CellChat", quietly = TRUE)) {
            cat("\nCellChat successfully installed from GitHub!\n")
        } else {
            cat("\nFailed to install CellChat from GitHub. You may need to install it manually.\n")
        }
    } else {
        cat("\nCellChat is already installed.\n")
    }
}, error = function(e) {
    cat(sprintf("\nError installing CellChat from GitHub: %s\n", e$message))
    cat("\nPlease install CellChat manually following instructions at: https://github.com/sqjin/CellChat\n")
})

# Verify installation of all packages
all_packages <- c(cran_packages, "CellChat")
installed_packages <- rownames(installed.packages())
missing_packages <- all_packages[!all_packages %in% installed_packages]

if (length(missing_packages) > 0) {
    cat("\n\nWarning: The following packages could not be installed:\n")
    cat(paste(missing_packages, collapse = ", "), "\n")
    cat("You may need to install them manually.\n")
} else {
    cat("\n\nAll required R packages have been successfully installed.\n")
}
'

echo "R dependencies installation attempt completed. Check the output for any warnings or errors."
echo "All dependencies installed successfully in the '$ENV_NAME' environment!"
