################################################################################
#                CAMK2D Analysis Pipeline - Dependency Installer              #
#                           One-time setup script                             #
################################################################################

cat("================================================================================\n")
cat("                  CAMK2D Analysis Pipeline - Dependency Installer              \n")
cat("================================================================================\n\n")

# Function to install package if not already installed
install_if_missing <- function(pkg, type = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from", type, "...\n")
    if (type == "CRAN") {
      install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
    } else if (type == "Bioconductor") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    # Check if installation succeeded
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "installed successfully\n")
      return(TRUE)
    } else {
      cat("  ✗ Failed to install", pkg, "\n")
      return(FALSE)
    }
  } else {
    cat("  ✓", pkg, "already installed\n")
    return(TRUE)
  }
}

# Track failed installations
failed_packages <- character()

################################################################################
# STEP 1: Install BiocManager (required for Bioconductor packages)
################################################################################

cat("\n[1/5] Installing BiocManager...\n")
cat("--------------------------------\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Suppress BiocManager messages
suppressMessages(library(BiocManager))

################################################################################
# STEP 2: Core R Packages (CRAN)
################################################################################

cat("\n[2/5] Installing Core R Packages...\n")
cat("------------------------------------\n")

cran_packages <- c(
  # Data manipulation and visualization
  "tidyverse",      # Includes: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
  "lubridate",      # Date/time manipulation
  "glue",           # String interpolation
  "here",           # Path management
  
  # Tables and reporting
  "knitr",          # Report generation
  "rmarkdown",      # R Markdown documents
  "kableExtra",     # Enhanced tables
  "gt",             # Grammar of tables
  "DT",             # Interactive data tables
  "openxlsx",       # Excel file handling
  
  # Visualization
  "plotly",         # Interactive plots
  "corrplot",       # Correlation plots
  "ComplexHeatmap", # Advanced heatmaps
  "circlize",       # Circular visualizations
  "VennDiagram",    # Venn diagrams
  "UpSetR",         # UpSet plots
  "visNetwork",     # Network visualization
  "networkD3",      # D3.js network viz
  "ggraph",         # Network plots with ggplot2
  "igraph",         # Network analysis
  "tidygraph",      # Tidy network data
  
  # Dashboard and Shiny
  "shiny",          # Interactive web apps
  "shinydashboard", # Dashboard layouts
  "flexdashboard",  # Flexible dashboards
  
  # Statistical analysis
  "metafor",        # Meta-analysis
  "meta",           # Meta-analysis tools
  "bayesmeta",      # Bayesian meta-analysis
  "dmetar",         # Meta-analysis tools
  "boot",           # Bootstrap methods
  "broom",          # Tidy model outputs
  "emmeans",        # Estimated marginal means
  "lme4",           # Mixed effects models
  
  # Machine learning
  "caret",          # Classification and regression
  "randomForest",   # Random forests
  "e1071",          # Support vector machines
  "glmnet",         # Regularized regression
  
  # Text mining
  "tm",             # Text mining framework
  "tidytext",       # Tidy text mining
  "wordcloud",      # Word clouds
  "wordcloud2",     # HTML word clouds
  "topicmodels",    # Topic modeling
  "SnowballC",      # Text stemming
  
  # Parallel processing
  "parallel",       # Parallel computing
  "doParallel",     # Parallel backend
  
  # Bioinformatics (CRAN)
  "rentrez",        # NCBI Entrez access
  "seqinr",         # Sequence analysis
  "Peptides",       # Peptide properties
  "Biostrings"      # Biological sequences
)

cat("Installing", length(cran_packages), "CRAN packages...\n\n")

for (pkg in cran_packages) {
  success <- install_if_missing(pkg, "CRAN")
  if (!success) {
    failed_packages <- c(failed_packages, pkg)
  }
}

################################################################################
# STEP 3: Bioconductor Packages
################################################################################

cat("\n[3/5] Installing Bioconductor Packages...\n")
cat("------------------------------------------\n")

bioc_packages <- c(
  # Core Bioconductor
  "Biobase",          # Base functions
  
  # Database access
  "GEOquery",         # Gene Expression Omnibus
  "biomaRt",          # BioMart database access
  
  # Annotation
  "org.Hs.eg.db",     # Human gene annotation
  "clusterProfiler",  # Functional enrichment
  "enrichplot",       # Enrichment visualization
  
  # Protein databases
  "UniProt.ws",       # UniProt database
  "STRINGdb",         # STRING protein interactions
  
  # PubMed access
  "RISmed"            # PubMed/MEDLINE interface
)

cat("Installing", length(bioc_packages), "Bioconductor packages...\n\n")

for (pkg in bioc_packages) {
  success <- install_if_missing(pkg, "Bioconductor")
  if (!success) {
    failed_packages <- c(failed_packages, pkg)
  }
}

################################################################################
# STEP 4: Version Verification
################################################################################

cat("\n[4/5] Verifying Package Versions...\n")
cat("------------------------------------\n")

# Key packages to check versions
key_packages <- c(
  "tidyverse", "shiny", "rmarkdown", "metafor", 
  "GEOquery", "UniProt.ws", "STRINGdb"
)

cat("\nKey package versions:\n")
for (pkg in key_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("  %-15s : %s\n", pkg, version))
  }
}

################################################################################
# STEP 5: Final Status Report
################################################################################

cat("\n[5/5] Installation Summary\n")
cat("--------------------------\n")

# Count installed packages
all_packages <- c(cran_packages, bioc_packages)
installed_count <- sum(sapply(all_packages, requireNamespace, quietly = TRUE))

cat("\nTotal packages required:", length(all_packages), "\n")
cat("Successfully installed:", installed_count, "\n")

if (length(failed_packages) > 0) {
  cat("\n⚠️  WARNING: The following packages failed to install:\n")
  for (pkg in failed_packages) {
    cat("  -", pkg, "\n")
  }
  cat("\nTry installing these manually with:\n")
  cat("  install.packages('package_name') for CRAN packages\n")
  cat("  BiocManager::install('package_name') for Bioconductor packages\n")
} else {
  cat("\n✅ All packages installed successfully!\n")
}

################################################################################
# R Session Info
################################################################################

cat("\n================================================================================\n")
cat("                              System Information                               \n")
cat("================================================================================\n")

cat("\nR Version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("BiocManager Version:", as.character(packageVersion("BiocManager")), "\n")

# Check memory
mem_info <- gc()
cat("\nMemory Usage:\n")
cat("  Used:", round(sum(mem_info[,2]), 1), "MB\n")
cat("  Max used:", round(sum(mem_info[,6]), 1), "MB\n")

################################################################################
# Test Critical Functions
################################################################################

cat("\n================================================================================\n")
cat("                         Testing Critical Functions                            \n")
cat("================================================================================\n\n")

# Test rmarkdown
cat("Testing rmarkdown render function... ")
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("✓ OK\n")
} else {
  cat("✗ FAILED\n")
}

# Test Shiny
cat("Testing Shiny... ")
if (requireNamespace("shiny", quietly = TRUE)) {
  cat("✓ OK\n")
} else {
  cat("✗ FAILED\n")
}

# Test GEOquery
cat("Testing GEOquery... ")
if (requireNamespace("GEOquery", quietly = TRUE)) {
  cat("✓ OK\n")
} else {
  cat("✗ FAILED\n")
}

# Test UniProt.ws
cat("Testing UniProt.ws... ")
if (requireNamespace("UniProt.ws", quietly = TRUE)) {
  cat("✓ OK\n")
} else {
  cat("✗ FAILED\n")
}

################################################################################
# Final Instructions
################################################################################

cat("\n================================================================================\n")
cat("                            INSTALLATION COMPLETE                              \n")
cat("================================================================================\n\n")

if (length(failed_packages) == 0) {
  cat("✅ Your system is ready to run the CAMK2D analysis pipeline!\n\n")
  cat("Next steps:\n")
  cat("1. Open 'MASTER_ANALYSIS.Rmd' in RStudio\n")
  cat("2. Click the 'Knit' button to run all analyses\n")
  cat("   OR\n")
  cat("   Run: source('run_all_analyses.R')\n\n")
} else {
  cat("⚠️  Some packages failed to install. Please install them manually before running the pipeline.\n\n")
  cat("For help with specific packages:\n")
  cat("- CRAN packages: https://cran.r-project.org/\n")
  cat("- Bioconductor: https://bioconductor.org/packages/\n\n")
}

cat("For detailed instructions, see: HOW_TO_RUN.md\n")
cat("\n================================================================================\n")