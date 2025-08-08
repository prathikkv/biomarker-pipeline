################################################################################
#                     CAMK2D COMPLETE ANALYSIS PIPELINE                       #
#                         Master Script to Run All Steps                      #
################################################################################

# This script runs all CAMK2D analyses in sequence
# Run this in RStudio to execute the complete pipeline

cat("=========================================\n")
cat("   CAMK2D ANALYSIS PIPELINE STARTING    \n")
cat("=========================================\n\n")

# Record start time
start_time <- Sys.time()

# Set working directory (adjust if needed)
if (!file.exists("01_literature_mining.Rmd")) {
  cat("Please set working directory to CAMK2D folder\n")
  cat("Use: setwd('/Users/macbookpro/Desktop/CAMK2D')\n")
  stop("Working directory not set correctly")
}

################################################################################
# STEP 1: LITERATURE MINING
################################################################################

cat("STEP 1: LITERATURE MINING\n")
cat("-------------------------\n")
tryCatch({
  cat("Running literature mining analysis...\n")
  rmarkdown::render(
    "01_literature_mining.Rmd",
    output_file = "01_literature_mining.html",
    quiet = TRUE
  )
  cat("✅ Literature mining completed\n\n")
}, error = function(e) {
  cat("❌ Error in literature mining:", e$message, "\n")
  cat("Continuing with next step...\n\n")
})

################################################################################
# STEP 2: GEO TRANSCRIPTOMICS META-ANALYSIS
################################################################################

cat("STEP 2: GEO TRANSCRIPTOMICS\n")
cat("---------------------------\n")
tryCatch({
  cat("Running GEO transcriptomics analysis...\n")
  rmarkdown::render(
    "02_geo_transcriptomics_meta_analysis.Rmd",
    output_file = "02_geo_transcriptomics_meta_analysis.html",
    quiet = TRUE
  )
  cat("✅ GEO transcriptomics completed\n\n")
}, error = function(e) {
  cat("❌ Error in GEO analysis:", e$message, "\n")
  cat("Continuing with next step...\n\n")
})

################################################################################
# STEP 3: PHOSPHOPROTEOMICS ANALYSIS
################################################################################

cat("STEP 3: PHOSPHOPROTEOMICS\n")
cat("-------------------------\n")
tryCatch({
  cat("Running phosphoproteomics analysis...\n")
  rmarkdown::render(
    "03_phosphoproteomics.Rmd",
    output_file = "03_phosphoproteomics.html",
    quiet = TRUE
  )
  cat("✅ Phosphoproteomics completed\n\n")
}, error = function(e) {
  cat("❌ Error in phosphoproteomics:", e$message, "\n")
  cat("Continuing with next step...\n\n")
})

################################################################################
# STEP 4: MODULAR META-ANALYSIS INTEGRATION
################################################################################

cat("STEP 4: META-ANALYSIS INTEGRATION (MODULAR)\n")
cat("--------------------------------------------\n")

# STEP 4a: Core Meta-Analysis (Essential - guaranteed to work)
cat("4a. Running core meta-analysis...\n")
tryCatch({
  rmarkdown::render(
    "04a_basic_meta_analysis.Rmd",
    output_file = "04a_basic_meta_analysis.html",
    quiet = TRUE
  )
  cat("✅ Core meta-analysis completed\n")
}, error = function(e) {
  cat("❌ Error in core meta-analysis:", e$message, "\n")
  cat("Continuing with next step...\n")
})

# STEP 4b: Publication Bias Analysis (Optional)
cat("4b. Running publication bias analysis...\n")
tryCatch({
  rmarkdown::render(
    "04b_publication_bias.Rmd",
    output_file = "04b_publication_bias.html",
    quiet = TRUE
  )
  cat("✅ Publication bias analysis completed\n")
}, error = function(e) {
  cat("⚠️ Publication bias analysis failed:", e$message, "\n")
  cat("Skipping advanced bias testing...\n")
})

# STEP 4c: Advanced Statistics (Optional)
cat("4c. Running advanced statistics...\n")
tryCatch({
  rmarkdown::render(
    "04c_advanced_statistics.Rmd", 
    output_file = "04c_advanced_statistics.html",
    quiet = TRUE
  )
  cat("✅ Advanced statistics completed\n")
}, error = function(e) {
  cat("⚠️ Advanced statistics failed:", e$message, "\n")
  cat("Skipping bootstrap/Bayesian analysis...\n")
})

# STEP 4d: Machine Learning Analysis (Optional)
cat("4d. Running machine learning analysis...\n")
tryCatch({
  rmarkdown::render(
    "04d_machine_learning.Rmd",
    output_file = "04d_machine_learning.html",
    quiet = TRUE
  )
  cat("✅ Machine learning analysis completed\n")
}, error = function(e) {
  cat("⚠️ Machine learning analysis failed:", e$message, "\n")
  cat("Skipping ML analysis...\n")
})

# STEP 4e: Advanced Visualizations (Optional)
cat("4e. Running advanced visualizations...\n")
tryCatch({
  rmarkdown::render(
    "04e_visualizations.Rmd",
    output_file = "04e_visualizations.html", 
    quiet = TRUE
  )
  cat("✅ Advanced visualizations completed\n")
}, error = function(e) {
  cat("⚠️ Visualizations failed:", e$message, "\n")
  cat("Skipping advanced plots...\n")
})

# STEP 4f: Integration Summary (Optional)
cat("4f. Running integration summary...\n")
tryCatch({
  rmarkdown::render(
    "04f_integration_summary.Rmd",
    output_file = "04f_integration_summary.html",
    quiet = TRUE
  )
  cat("✅ Integration summary completed\n")
}, error = function(e) {
  cat("⚠️ Integration summary failed:", e$message, "\n")
  cat("Skipping final integration...\n")
})

cat("Meta-analysis integration completed (modular system with graceful failures)\n\n")

################################################################################
# STEP 5: DASHBOARD GENERATION
################################################################################

cat("STEP 5: DASHBOARD GENERATION\n")
cat("---------------------------\n")
tryCatch({
  cat("Building interactive dashboard...\n")
  rmarkdown::render(
    "05_dashboard.Rmd",
    output_file = "dashboard.html",
    quiet = TRUE
  )
  cat("✅ Dashboard generated\n\n")
}, error = function(e) {
  cat("❌ Error in dashboard generation:", e$message, "\n\n")
})

################################################################################
# RESULTS SUMMARY
################################################################################

cat("=========================================\n")
cat("        ANALYSIS PIPELINE COMPLETE       \n")
cat("=========================================\n\n")

# Calculate runtime
end_time <- Sys.time()
if (exists("start_time")) {
  runtime <- difftime(end_time, start_time, units = "mins")
  cat("Total runtime:", round(runtime, 2), "minutes\n\n")
} else {
  cat("Runtime tracking not available\n\n")
}

# Check generated files
cat("GENERATED FILES:\n")
cat("---------------\n")

output_files <- c(
  "01_literature_mining.html",
  "02_geo_transcriptomics_meta_analysis.html", 
  "03_phosphoproteomics.html",
  "04a_basic_meta_analysis.html",
  "04b_publication_bias.html",
  "04c_advanced_statistics.html",
  "04d_machine_learning.html",
  "04e_visualizations.html",
  "04f_integration_summary.html",
  "dashboard.html"
)

for (file in output_files) {
  if (file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "(not found)\n")
  }
}

# Check data files
cat("\nDATA FILES:\n")
cat("-----------\n")

data_files <- c(
  "data/dashboard_metrics.rds",
  "data/phospho_integration_data.rds",
  "results/core_meta_analysis.rds",
  "results/publication_bias_analysis.rds",
  "results/advanced_statistics.rds",
  "results/machine_learning_analysis.rds",
  "results/advanced_visualizations.rds",
  "results/final_integrated_results.rds"
)

for (file in data_files) {
  if (file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "(not found)\n")
  }
}

# Check Excel reports
cat("\nEXCEL REPORTS:\n")
cat("--------------\n")

excel_files <- list.files("results", pattern = "*.xlsx", full.names = TRUE)
if (length(excel_files) > 0) {
  for (file in excel_files) {
    cat("✅", basename(file), "\n")
  }
} else {
  cat("No Excel files found\n")
}

cat("\n=========================================\n")
cat("          PIPELINE COMPLETE!             \n")
cat("=========================================\n\n")

cat("Next steps:\n")
cat("1. Review HTML reports for each analysis step\n")
cat("2. Open dashboard.html in browser for interactive exploration\n")
cat("3. Check results/ folder for Excel reports\n\n")

cat("To launch interactive dashboard:\n")
cat("rmarkdown::run('05_dashboard.Rmd')\n")