# Simplified Dataset Validation for CAMK2D Research
# Focus on identifying the best datasets for rigorous analysis

cat("=== SIMPLIFIED CAMK2D DATASET VALIDATION ===\n")

# Load configuration
source("config/datasets.R")

# Helper functions
calculate_relevance_score <- function(row, config_name) {
  score <- 5  # Base score
  
  # Disease relevance
  if (grepl("atrial.*fibrillation", tolower(row$Disease_Focus))) score <- score + 3
  if (grepl("heart.*failure", tolower(row$Disease_Focus))) score <- score + 3
  if (grepl("cardiomyopathy", tolower(row$Disease_Focus))) score <- score + 2.5
  
  # CAMK2D specific studies get bonus
  if (grepl("CAMK2D|camk2d", row$Title) || grepl("camk2d_specific", config_name)) score <- score + 2
  
  # Platform bonus
  if (row$Platform == "RNA-seq") score <- score + 1
  if (row$Platform == "Single-cell RNA-seq") score <- score + 1.5
  
  return(min(score, 10))
}

calculate_power_score <- function(sample_count) {
  if (sample_count >= 100) return(10)
  if (sample_count >= 75) return(9)
  if (sample_count >= 50) return(8)
  if (sample_count >= 30) return(6)
  if (sample_count >= 20) return(4)
  return(2)
}

assess_camk2d_suitability <- function(row, config_name) {
  score <- 6  # Base assumption that CAMK2D is present
  
  # Boost for CAMK2D-specific studies
  if (grepl("CAMK2D|camk2d", row$Title)) score <- score + 3
  if (grepl("camk2d_specific", config_name)) score <- score + 2
  
  # Boost for cardiac tissue
  if (grepl("cardiac|heart|atrial", tolower(row$Title))) score <- score + 1
  
  # Platform considerations
  if (row$Platform == "RNA-seq") score <- score + 0.5
  
  return(min(score, 10))
}

determine_tier <- function(overall_score) {
  if (overall_score >= 8.5) return("Tier 1 (Excellent)")
  if (overall_score >= 7.0) return("Tier 2 (Good)")
  if (overall_score >= 5.5) return("Tier 3 (Acceptable)")
  return("Exclude")
}

generate_recommendation <- function(overall_score, row) {
  if (overall_score >= 8.5) return("Priority for analysis")
  if (overall_score >= 7.0) return("Suitable for meta-analysis")
  if (overall_score >= 5.5) return("Use with caution")
  return("Exclude from analysis")
}

# Get all datasets
all_configs <- get_all_dataset_configs()

# Create comprehensive dataset info
dataset_analysis <- data.frame()

for (config_name in names(all_configs)) {
  config <- all_configs[[config_name]]
  
  for (i in 1:nrow(config)) {
    row <- config[i, ]
    
    # Calculate scientific relevance score
    relevance_score <- calculate_relevance_score(row, config_name)
    
    # Calculate statistical power score
    power_score <- calculate_power_score(row$Sample_Count)
    
    # Assess CAMK2D suitability
    camk2d_score <- assess_camk2d_suitability(row, config_name)
    
    # Overall assessment
    overall_score <- (relevance_score * 0.4 + power_score * 0.3 + camk2d_score * 0.3)
    
    # Determine tier and recommendation
    tier <- determine_tier(overall_score)
    recommendation <- generate_recommendation(overall_score, row)
    
    dataset_analysis <- rbind(dataset_analysis, data.frame(
      GSE_ID = row$GSE_ID,
      Category = config_name,
      Title = substr(row$Title, 1, 50),
      Sample_Count = row$Sample_Count,
      Platform = row$Platform,
      Disease_Focus = row$Disease_Focus,
      Quality_Score = row$Quality_Score,
      Relevance_Score = round(relevance_score, 1),
      Power_Score = round(power_score, 1),
      CAMK2D_Score = round(camk2d_score, 1),
      Overall_Score = round(overall_score, 1),
      Tier = tier,
      Recommendation = recommendation,
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by overall score
dataset_analysis <- dataset_analysis[order(-dataset_analysis$Overall_Score), ]

# Display results
cat("\n=== DATASET VALIDATION RESULTS ===\n")
cat("Total datasets analyzed:", nrow(dataset_analysis), "\n")

tier_summary <- table(dataset_analysis$Tier)
for (tier in names(tier_summary)) {
  cat(tier, ":", tier_summary[tier], "\n")
}

cat("\nAverage overall score:", round(mean(dataset_analysis$Overall_Score), 2), "\n")

# Top 15 datasets
cat("\n=== TOP 15 DATASETS FOR CAMK2D RESEARCH ===\n")
top_15 <- head(dataset_analysis, 15)
for (i in 1:nrow(top_15)) {
  row <- top_15[i, ]
  cat(sprintf("%2d. %-12s (%-20s) Score: %4.1f Samples: %3d %-15s %s\n",
              i, row$GSE_ID, substr(row$Category, 1, 20), row$Overall_Score, 
              row$Sample_Count, row$Platform, row$Tier))
}

# Largest datasets by disease
cat("\n=== LARGEST DATASETS BY DISEASE CATEGORY ===\n")

# AF datasets
af_datasets <- dataset_analysis[grepl("atrial_fib|AF", dataset_analysis$Category) | 
                                grepl("atrial fibrillation", dataset_analysis$Disease_Focus), ]
af_datasets <- af_datasets[order(-af_datasets$Sample_Count), ]

cat("\nATRIAL FIBRILLATION - Top 5:\n")
if (nrow(af_datasets) > 0) {
  top_af <- head(af_datasets, 5)
  for (i in 1:min(5, nrow(top_af))) {
    row <- top_af[i, ]
    cat(sprintf("  %d. %-12s - %3d samples (Score: %4.1f) %s\n",
                i, row$GSE_ID, row$Sample_Count, row$Overall_Score, row$Platform))
  }
}

# HF datasets  
hf_datasets <- dataset_analysis[grepl("heart_failure|HF", dataset_analysis$Category) | 
                                grepl("heart failure", dataset_analysis$Disease_Focus), ]
hf_datasets <- hf_datasets[order(-hf_datasets$Sample_Count), ]

cat("\nHEART FAILURE - Top 5:\n")
if (nrow(hf_datasets) > 0) {
  top_hf <- head(hf_datasets, 5)
  for (i in 1:min(5, nrow(top_hf))) {
    row <- top_hf[i, ]
    cat(sprintf("  %d. %-12s - %3d samples (Score: %4.1f) %s\n",
                i, row$GSE_ID, row$Sample_Count, row$Overall_Score, row$Platform))
  }
}

# Identify datasets for exclusion
excluded_datasets <- dataset_analysis[dataset_analysis$Tier == "Exclude", ]
cat("\n=== DATASETS RECOMMENDED FOR EXCLUSION ===\n")
if (nrow(excluded_datasets) > 0) {
  for (i in 1:nrow(excluded_datasets)) {
    row <- excluded_datasets[i, ]
    cat(sprintf("- %s (Score: %.1f) - %s\n", row$GSE_ID, row$Overall_Score, row$Recommendation))
  }
} else {
  cat("No datasets recommended for exclusion!\n")
}

# Calculate optimized portfolio stats
recommended_datasets <- dataset_analysis[dataset_analysis$Overall_Score >= 5.5, ]
total_samples <- sum(recommended_datasets$Sample_Count)
rna_seq_count <- sum(recommended_datasets$Platform == "RNA-seq")

cat("\n=== OPTIMIZED PORTFOLIO STATISTICS ===\n")
cat("Recommended datasets:", nrow(recommended_datasets), "\n")
cat("Total samples:", total_samples, "\n")
cat("RNA-seq datasets:", rna_seq_count, "/", nrow(recommended_datasets), 
    "(", round(rna_seq_count/nrow(recommended_datasets)*100, 1), "%)\n")
cat("Average score:", round(mean(recommended_datasets$Overall_Score), 2), "\n")

# Export results
write.csv(dataset_analysis, "data/dataset_validation_results.csv", row.names = FALSE)
write.csv(recommended_datasets, "data/recommended_datasets.csv", row.names = FALSE)

cat("\n=== RESULTS EXPORTED ===\n")
cat("- data/dataset_validation_results.csv (all datasets)\n")
cat("- data/recommended_datasets.csv (recommended for analysis)\n")

# Specific recommendations for Activity 2
cat("\n=== RECOMMENDATIONS FOR ACTIVITY 2 (LARGEST DATASETS) ===\n")

if (nrow(af_datasets) > 0) {
  largest_af <- af_datasets[1, ]
  cat("LARGEST AF DATASET:", largest_af$GSE_ID, "with", largest_af$Sample_Count, "samples\n")
  cat("  Platform:", largest_af$Platform, "Score:", largest_af$Overall_Score, "\n")
}

if (nrow(hf_datasets) > 0) {
  largest_hf <- hf_datasets[1, ]
  cat("LARGEST HF DATASET:", largest_hf$GSE_ID, "with", largest_hf$Sample_Count, "samples\n") 
  cat("  Platform:", largest_hf$Platform, "Score:", largest_hf$Overall_Score, "\n")
}

cat("\n=== NEXT STEPS ===\n")
cat("1. Download expression data for top-scoring datasets\n")
cat("2. Validate actual CAMK2D expression levels\n")
cat("3. Perform DGE analysis on largest AF and HF datasets\n")
cat("4. Search databases for additional large-scale studies\n")
cat("5. Implement phosphoproteomic analysis pipeline\n")

cat("\n=== VALIDATION COMPLETE ===\n")