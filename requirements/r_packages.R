#!/usr/bin/env Rscript

###############################################################################
# Install R Packages for Diapause GWAS Analysis
#
# Usage: Rscript r_packages.R
#
# This script installs all required R packages for the analysis pipeline
###############################################################################

# Function to check and install packages
install_if_needed <- function(package_name, bioc = FALSE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    cat("Installing", package_name, "...\n")
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package_name)
    } else {
      install.packages(package_name)
    }
    cat(package_name, "installed successfully\n\n")
  } else {
    cat(package_name, "already installed\n")
  }
}

cat("================================================\n")
cat("Installing R Packages\n")
cat("R version:", R.version$version.string, "\n")
cat("================================================\n\n")

# CRAN packages
cat("Installing CRAN packages...\n\n")
cran_packages <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "gridExtra",
  "gplots",
  "igraph",
  "plotly"
)

for (pkg in cran_packages) {
  install_if_needed(pkg, bioc = FALSE)
}

# Bioconductor packages
cat("\nInstalling Bioconductor packages...\n\n")
bioc_packages <- c(
  "adegenet",
  "rrBLUP",
  "vcfR",
  "goseq",
  "GenomicRanges",
  "GenomicFeatures"
)

for (pkg in bioc_packages) {
  install_if_needed(pkg, bioc = TRUE)
}

cat("\n================================================\n")
cat("Package installation complete!\n")
cat("================================================\n")

# Print session info
cat("\nR Session Info:\n")
sessionInfo()

# Verify key packages
cat("\n\nVerifying installation...\n")
required_packages <- c(
  "dplyr", "tidyr", "stringr", "ggplot2",
  "adegenet", "rrBLUP", "vcfR", "goseq"
)

all_loaded <- TRUE
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("WARNING: Failed to install or load", pkg, "\n")
    all_loaded <- FALSE
  }
}

if (all_loaded) {
  cat("\n✓ All packages installed successfully!\n")
  cat("  Ready to run analysis pipeline\n")
} else {
  cat("\n✗ Some packages failed to install\n")
  cat("  Please check error messages above\n")
}
