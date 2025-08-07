# CAMK2D Comprehensive Analysis Pipeline

A production-ready R-based bioinformatics pipeline for analyzing CAMK2D (Calcium/Calmodulin Dependent Protein Kinase II Delta) in cardiac diseases.

## 🎯 Overview

This pipeline performs comprehensive analysis of CAMK2D through five integrated modules:

1. **Literature Mining** - PubMed analysis and research trends
2. **GEO Transcriptomics** - Gene expression meta-analysis
3. **Phosphoproteomic Analysis & Protein Interactions** - Protein analysis and network studies
4. **Integrated Statistical Analysis & Meta-Analysis** - Comprehensive statistical integration
5. **Automated Reporting & Dashboard Creation** - Interactive reporting and visualization

## Research Objectives

### Primary Goals
- Identify CAMK family members and their roles in aFIB and HF
- Find transcriptomic evidence in disease states  
- Identify CAMK2D phosphorylation targets and inhibitor responses
- Characterize phosphoproteomic signatures

### Key Features
- **Comprehensive Analysis**: End-to-end research pipeline from literature to results
- **Reproducible Research**: All analyses documented in R Markdown with explanations
- **Professional Reporting**: Generate publication-ready figures and tables
- **Interactive Dashboards**: Explore data with interactive visualizations
- **Modular Design**: Each component can be run independently or as part of the full pipeline

## Quick Start

### Prerequisites
- **R** (version 4.0.0 or higher)
- **RStudio** (recommended IDE)
- Internet connection for database access

### Installation

1. **Clone or download this repository**
   ```bash
   git clone [repository-url]
   cd CAMK2D
   ```

2. **Install required packages**
   ```r
   # In R or RStudio, run:
   source("install_packages.R")      # Install CRAN packages
   source("install_bioconductor.R")  # Install Bioconductor packages
   ```

3. **Start with literature mining**
   ```r
   # Open and run in RStudio:
   rmarkdown::render("01_literature_mining.Rmd")
   ```

## Usage Instructions

### Running Individual Components

Each R Markdown document can be run independently:

```r
# Literature mining and analysis
rmarkdown::render("01_literature_mining.Rmd")

# Transcriptomic analysis
rmarkdown::render("02_geo_transcriptomics.Rmd")

# Phosphoproteomic analysis
rmarkdown::render("03_phosphoproteomics.Rmd")

# Meta-analysis and integration
rmarkdown::render("04_meta_analysis.Rmd")

# Interactive dashboard
rmarkdown::run("05_dashboard.Rmd")
```

### Running the Complete Pipeline

```r
# Run all analyses in sequence
documents <- c(
  "01_literature_mining.Rmd",
  "02_geo_transcriptomics.Rmd",
  "03_phosphoproteomics.Rmd",
  "04_meta_analysis.Rmd"
)

for (doc in documents) {
  rmarkdown::render(doc)
}

# Launch interactive dashboard
rmarkdown::run("05_dashboard.Rmd")
```

## Directory Structure

```
CAMK2D/
├── README.md                     # This file
├── promts.md                    # Original project specifications
├── user_manual.md               # Detailed usage instructions
├── install_packages.R           # CRAN package installer
├── install_bioconductor.R       # Bioconductor package installer
├── 01_literature_mining.Rmd     # Literature analysis
├── 02_geo_transcriptomics.Rmd   # Gene expression analysis
├── 03_phosphoproteomics.Rmd     # Protein analysis
├── 04_meta_analysis.Rmd         # Statistical integration
├── 05_dashboard.Rmd             # Interactive dashboard
├── data/                        # Input and intermediate data
├── results/                     # Analysis results
└── figures/                     # Generated visualizations
```

## R Package Requirements

### Core Packages (CRAN)
- **Data manipulation**: `tidyverse`, `dplyr`, `ggplot2`
- **Reporting**: `rmarkdown`, `knitr`, `flexdashboard`, `shiny`
- **Statistics**: `metafor`, `meta`, `broom`
- **Text mining**: `RISmed`, `rentrez`, `tm`, `wordcloud`
- **Networks**: `igraph`, `visNetwork`

### Bioconductor Packages
- **Genomics**: `GEOquery`, `limma`, `edgeR`, `DESeq2`
- **Annotation**: `org.Hs.eg.db`, `GO.db`, `clusterProfiler`
- **Visualization**: `ComplexHeatmap`, `EnhancedVolcano`
- **Proteomics**: `UniProt.ws`, `STRINGdb`, `Biostrings`

## Expected Outputs

Each analysis generates:
- **HTML Reports**: Interactive documents with embedded plots and tables
- **PDF Documents**: Publication-ready formatted documents  
- **Data Files**: Excel, CSV, and RData exports
- **Figures**: High-resolution images for presentations and publications
- **Interactive Dashboards**: Shiny applications for data exploration

## Troubleshooting

### Common Issues

1. **Package installation errors**
   - Ensure you have the latest R version
   - Install packages individually if batch installation fails
   - Check system dependencies (especially for Bioconductor packages)

2. **Memory issues with large datasets**
   - Increase R memory limit: `memory.limit(size = 8000)` (Windows)
   - Use data subsampling for testing
   - Consider running analyses on server with more memory

3. **Network issues accessing databases**
   - Check internet connection
   - Some analyses may require VPN or institutional access
   - Cached results are stored to minimize repeated downloads

### Getting Help

1. Check the `user_manual.md` for detailed instructions
2. Review individual R Markdown documents for specific guidance
3. Examine error messages and traceback information
4. Test package installation with `sessionInfo()`

## Contributing

This is a research automation system designed for CAMK2D studies but can be adapted for:
- Other protein kinase families
- Different disease contexts
- Alternative data sources

To extend the system:
1. Follow the existing R Markdown structure
2. Document all methods and interpretations
3. Include data previews and diagnostic plots
4. Test with sample data before running full analyses

## Timeline Estimation

- **Initial setup**: 1-2 hours (package installation and testing)
- **Each component**: 2-4 hours runtime (depending on data size)
- **Complete pipeline**: 8-12 hours total runtime
- **Results review**: 2-3 hours per component

## License

This research automation system is provided for academic and research use. Please cite appropriately if used in publications.

## Contact

For questions about this research automation system, please refer to the documentation or create an issue in the repository.