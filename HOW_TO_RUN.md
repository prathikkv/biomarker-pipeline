# How to Run CAMK2D Analysis Pipeline in RStudio

## 🚀 Quick Start (Recommended)

### Option 1: Run Everything with Master Analysis Notebook
1. Open RStudio
2. Open file: `MASTER_ANALYSIS.Rmd`
3. Click **"Knit"** button OR **"Run All"** button
4. Wait ~10-15 minutes for complete analysis

### Option 2: Run Everything with R Script
1. Open RStudio
2. Open file: `run_all_analyses.R`
3. Click **"Source"** button (top-right of script panel)
4. Wait ~10-15 minutes for complete analysis

## 📋 Step-by-Step Instructions

### Method 1: Using RStudio Interface (Visual)

1. **Open RStudio**
   - Set working directory: `setwd("/Users/macbookpro/Desktop/CAMK2D")`

2. **Run analyses in order:**
   ```r
   # Step 1: Literature Mining
   rmarkdown::render("01_literature_mining.Rmd")
   
   # Step 2: GEO Transcriptomics
   rmarkdown::render("02_geo_transcriptomics_meta_analysis.Rmd")
   
   # Step 3: Phosphoproteomics
   rmarkdown::render("03_phosphoproteomics.Rmd")
   
   # Step 4: Meta-Analysis
   rmarkdown::render("04_meta_analysis.Rmd")
   
   # Step 5: Dashboard
   rmarkdown::render("05_dashboard.Rmd")
   ```

### Method 2: Using Master Script (Automated)

1. **In RStudio Console:**
   ```r
   # Set working directory
   setwd("/Users/macbookpro/Desktop/CAMK2D")
   
   # Run all analyses
   source("run_all_analyses.R")
   ```

### Method 3: Using RMarkdown Notebook (Interactive)

1. **Open `MASTER_ANALYSIS.Rmd` in RStudio**
2. **Choose one:**
   - Click **"Knit"** to run all and generate report
   - Click **"Run All"** (Ctrl/Cmd + Alt + R) to run all chunks
   - Run chunks individually with green arrow buttons

## 🖥️ View Interactive Dashboard

After analysis completes:

### Static Version:
- Open `dashboard.html` in any web browser

### Interactive Shiny Version:
```r
# In RStudio Console:
rmarkdown::run("05_dashboard.Rmd")
```
Then open browser to: http://127.0.0.1:XXXX (port shown in console)

## ⏱️ Expected Runtime

| Step | Approximate Time |
|------|-----------------|
| Literature Mining | 2-3 minutes |
| GEO Transcriptomics | 3-4 minutes |
| Phosphoproteomics | 2-3 minutes |
| Meta-Analysis | 2-3 minutes |
| Dashboard | 1-2 minutes |
| **Total** | **10-15 minutes** |

## 📁 Output Files

After completion, you'll have:

### HTML Reports:
- `01_literature_mining.html`
- `02_geo_transcriptomics_meta_analysis.html`
- `03_phosphoproteomics.html`
- `04_meta_analysis.html`
- `dashboard.html`

### Data Files:
- `data/dashboard_metrics.rds`
- `data/phospho_integration_data.rds`
- `results/meta_analysis_summary.rds`

### Excel Reports:
- `results/CAMK2D_Literature_Analysis_*.xlsx`
- `results/Meta_Analysis_Results_*.xlsx`
- `results/Phosphoproteomics_Analysis_*.xlsx`

## 🔧 Troubleshooting

### If you get package errors:
```r
# Install required packages
install.packages(c("tidyverse", "rmarkdown", "flexdashboard", 
                   "plotly", "DT", "visNetwork"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("UniProt.ws", "STRINGdb"))
```

### If a step fails:
- Check error message in console
- Run that individual step again
- Continue with next steps (they may still work)

### Memory issues:
```r
# Clear memory between steps
gc()

# Or restart R session
# Session > Restart R
```

## 🎯 Individual Step Re-runs

To re-run only specific analyses:

```r
# Just literature mining
rmarkdown::render("01_literature_mining.Rmd")

# Just dashboard
rmarkdown::render("05_dashboard.Rmd")
```

## 💡 Tips

1. **First time**: Run `MASTER_ANALYSIS.Rmd` - it shows progress nicely
2. **Subsequent runs**: Use `run_all_analyses.R` - it's faster
3. **Development**: Run steps individually to test changes
4. **Production**: Use the shell script for automated runs

## 📊 Verify Success

After running, check:
```r
# In R Console:
# Check if all HTML files exist
list.files(pattern = "*.html")

# Check if data files were created
file.exists("data/dashboard_metrics.rds")
file.exists("results/meta_analysis_summary.rds")

# Launch dashboard to verify
rmarkdown::run("05_dashboard.Rmd")
```

## ❓ Help

- Check `CLAUDE.md` for analysis guidelines
- Review individual `.Rmd` files for step details
- Dashboard documentation: `user_manual.html`