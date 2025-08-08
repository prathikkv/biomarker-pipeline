# Advanced Dataset Validation for CAMK2D Research Proposal
# Rigorous scientific validation to identify highest quality datasets

cat("=== ADVANCED CAMK2D RESEARCH VALIDATION ===\n")
cat("Implementing rigorous dataset validation for scientific excellence\n\n")

# Load required functions
source("config/datasets.R")
source("functions/advanced_dataset_validator.R")

# Get all current datasets
all_configs <- get_all_dataset_configs()
all_dataset_ids <- c()
dataset_categories <- c()
dataset_info <- data.frame()

cat("=== CURRENT DATASET INVENTORY ===\n")
for (config_name in names(all_configs)) {
  config <- all_configs[[config_name]]
  cat(config_name, ":", nrow(config), "datasets\n")
  
  all_dataset_ids <- c(all_dataset_ids, config$GSE_ID)
  dataset_categories <- c(dataset_categories, rep(config_name, nrow(config)))
  
  # Create comprehensive info dataframe
  config_info <- data.frame(
    GSE_ID = config$GSE_ID,
    Category = config_name,
    Title = config$Title,
    Sample_Count = config$Sample_Count,
    Platform = config$Platform,
    Disease_Focus = config$Disease_Focus,
    Quality_Score = config$Quality_Score,
    Verified = config$Verified,
    stringsAsFactors = FALSE
  )
  dataset_info <- rbind(dataset_info, config_info)
}

cat("\nTotal datasets to validate:", length(all_dataset_ids), "\n")
cat("Total potential samples:", sum(dataset_info$Sample_Count), "\n\n")

# Run advanced validation on all datasets
cat("=== RUNNING ADVANCED VALIDATION ===\n")
cat("This will perform rigorous quality assessment...\n\n")

# Advanced batch validation (without downloading data for now - can be enabled later)
advanced_results <- batch_advanced_validation(all_dataset_ids, download_data = FALSE)

# Extract validation results
validation_results <- advanced_results$validation_results
summary_stats <- advanced_results$summary_stats

# Create comprehensive results dataframe
results_df <- data.frame(
  GSE_ID = character(),
  Category = character(),
  Title = character(),
  Original_Sample_Count = numeric(),
  Platform = character(),
  Overall_Score = numeric(),
  CAMK_Score = numeric(),
  Relevance_Score = numeric(),
  Power_Score = numeric(),
  CAMK2D_Detected = logical(),
  Tier = character(),
  Ready_for_Analysis = logical(),
  Exclusion_Reason = character(),
  stringsAsFactors = FALSE
)

# Process validation results
for (gse_id in names(validation_results)) {
  result <- validation_results[[gse_id]]
  dataset_row <- dataset_info[dataset_info$GSE_ID == gse_id, ]
  
  if (!is.null(result$overall_score) && result$overall_score > 0) {
    # Determine tier
    tier <- if (result$overall_score >= 8.0) "Tier 1 (Excellent)" 
            else if (result$overall_score >= 6.5) "Tier 2 (Good)"
            else if (result$overall_score >= 5.0) "Tier 3 (Acceptable)"
            else "Excluded"
    
    results_df <- rbind(results_df, data.frame(
      GSE_ID = gse_id,
      Category = dataset_row$Category[1],
      Title = substr(dataset_row$Title[1], 1, 50),  # Truncate for display
      Original_Sample_Count = dataset_row$Sample_Count[1],
      Platform = dataset_row$Platform[1],
      Overall_Score = round(result$overall_score, 2),
      CAMK_Score = round(result$camk_detection$detection_score, 1),
      Relevance_Score = round(result$scientific_relevance$overall_relevance, 1),
      Power_Score = round(result$statistical_power$power_score, 1),
      CAMK2D_Detected = result$camk_detection$camk2d_status$detected,
      Tier = tier,
      Ready_for_Analysis = result$ready_for_analysis %||% FALSE,
      Exclusion_Reason = result$exclusion_reason %||% "",
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by overall score (descending)
results_df <- results_df[order(-results_df$Overall_Score), ]

cat("\n", rep("=", 100), "\n")
cat("=== VALIDATION RESULTS SUMMARY ===\n")
cat("Total datasets processed:", nrow(results_df), "\n")
cat("Tier 1 (Excellent, ≥8.0):", sum(results_df$Tier == "Tier 1 (Excellent)"), "\n")
cat("Tier 2 (Good, 6.5-7.9):", sum(results_df$Tier == "Tier 2 (Good)"), "\n")
cat("Tier 3 (Acceptable, 5.0-6.4):", sum(results_df$Tier == "Tier 3 (Acceptable)"), "\n")
cat("Excluded (<5.0):", sum(results_df$Tier == "Excluded"), "\n")
cat("Ready for analysis:", sum(results_df$Ready_for_Analysis), "\n")
cat("Average score:", round(mean(results_df$Overall_Score), 2), "\n\n")

# Display top performing datasets
cat("=== TOP 15 DATASETS FOR CAMK2D RESEARCH ===\n")
top_datasets <- head(results_df, 15)
for (i in 1:nrow(top_datasets)) {
  dataset <- top_datasets[i, ]
  cat(sprintf("%2d. %s (%s) - Score: %.1f - %s samples - %s\n",
              i, dataset$GSE_ID, dataset$Category, dataset$Overall_Score,
              dataset$Original_Sample_Count, dataset$Platform))
}

cat("\n=== DATASETS RECOMMENDED FOR EXCLUSION ===\n")
excluded_datasets <- results_df[results_df$Tier == "Excluded", ]
if (nrow(excluded_datasets) > 0) {
  for (i in 1:nrow(excluded_datasets)) {
    dataset <- excluded_datasets[i, ]
    cat(sprintf("- %s (Score: %.1f) - Reason: %s\n",
                dataset$GSE_ID, dataset$Overall_Score, 
                ifelse(dataset$Exclusion_Reason == "", "Low quality", dataset$Exclusion_Reason)))
  }
} else {
  cat("No datasets recommended for exclusion - all meet minimum standards!\n")
}

# Focus on largest datasets for each disease category
cat("\n=== LARGEST DATASETS BY DISEASE CATEGORY ===\n")

# Atrial Fibrillation datasets
af_datasets <- results_df[grepl("Atrial.Fibrillation|AF", results_df$Category) | 
                          grepl("atrial.fibrillation|AF", results_df$Title, ignore.case = TRUE), ]
af_datasets <- af_datasets[order(-af_datasets$Original_Sample_Count), ]

cat("ATRIAL FIBRILLATION - Top 5 largest datasets:\n")
if (nrow(af_datasets) > 0) {
  top_af <- head(af_datasets, 5)
  for (i in 1:nrow(top_af)) {
    dataset <- top_af[i, ]
    cat(sprintf("  %d. %s - %d samples (Score: %.1f) - %s\n",
                i, dataset$GSE_ID, dataset$Original_Sample_Count, 
                dataset$Overall_Score, dataset$Platform))
  }
} else {
  cat("  No AF-specific datasets found in current configuration\n")
}

# Heart Failure datasets
hf_datasets <- results_df[grepl("Heart.Failure|HF", results_df$Category) | 
                          grepl("heart.failure|HF", results_df$Title, ignore.case = TRUE), ]
hf_datasets <- hf_datasets[order(-hf_datasets$Original_Sample_Count), ]

cat("\nHEART FAILURE - Top 5 largest datasets:\n")
if (nrow(hf_datasets) > 0) {
  top_hf <- head(hf_datasets, 5)
  for (i in 1:nrow(top_hf)) {
    dataset <- top_hf[i, ]
    cat(sprintf("  %d. %s - %d samples (Score: %.1f) - %s\n",
                i, dataset$GSE_ID, dataset$Original_Sample_Count, 
                dataset$Overall_Score, dataset$Platform))
  }
} else {
  cat("  No HF-specific datasets found in current configuration\n")
}

# Identify datasets needing expression data validation
cat("\n=== DATASETS REQUIRING EXPRESSION DATA VALIDATION ===\n")
cat("The following high-scoring datasets should have their actual expression data downloaded and analyzed:\n\n")

priority_datasets <- results_df[results_df$Overall_Score >= 6.0 & results_df$Ready_for_Analysis, ]
priority_datasets <- head(priority_datasets[order(-priority_datasets$Overall_Score), ], 10)

for (i in 1:nrow(priority_datasets)) {
  dataset <- priority_datasets[i, ]
  cat(sprintf("Priority %d: %s (%s) - Score: %.1f\n",
              i, dataset$GSE_ID, dataset$Platform, dataset$Overall_Score))
  cat(sprintf("  Samples: %d, Category: %s\n",
              dataset$Original_Sample_Count, dataset$Category))
  cat(sprintf("  Action: Download expression data and validate CAMK2D presence\n\n"))
}

# Calculate optimized dataset portfolio statistics
ready_datasets <- results_df[results_df$Ready_for_Analysis, ]
total_optimized_samples <- sum(ready_datasets$Original_Sample_Count)
rna_seq_datasets <- sum(ready_datasets$Platform == "RNA-seq")
modern_datasets <- nrow(ready_datasets)  # All are considered modern after filtering

cat("\n", rep("=", 100), "\n")
cat("=== OPTIMIZED DATASET PORTFOLIO STATISTICS ===\n")
cat("Recommended datasets for analysis:", nrow(ready_datasets), "\n")
cat("Total samples in optimized portfolio:", total_optimized_samples, "\n")
cat("RNA-seq datasets (preferred):", rna_seq_datasets, "/", nrow(ready_datasets), 
    "(", round(rna_seq_datasets/nrow(ready_datasets)*100, 1), "%)\n")
cat("Average quality score:", round(mean(ready_datasets$Overall_Score), 2), "/10\n")
cat("Tier 1 datasets:", sum(ready_datasets$Tier == "Tier 1 (Excellent)"), "\n")
cat("Tier 2 datasets:", sum(ready_datasets$Tier == "Tier 2 (Good)"), "\n")

# Export results for further analysis
write.csv(results_df, "data/advanced_validation_results.csv", row.names = FALSE)
write.csv(priority_datasets, "data/priority_datasets_for_analysis.csv", row.names = FALSE)

cat("\nResults exported to:\n")
cat("- data/advanced_validation_results.csv\n")
cat("- data/priority_datasets_for_analysis.csv\n")

cat("\n=== NEXT STEPS RECOMMENDED ===\n")
cat("1. Download expression data for top 10 priority datasets\n")
cat("2. Validate actual CAMK2D expression levels\n")
cat("3. Perform comprehensive database search for additional large AF/HF datasets\n")
cat("4. Implement DGE analysis pipeline on validated datasets\n")
cat("5. Integrate results with phosphoproteomic analysis\n")

cat("\n=== ADVANCED VALIDATION COMPLETE ===\n")
cat("Pipeline ready for high-impact CAMK2D research!\n")