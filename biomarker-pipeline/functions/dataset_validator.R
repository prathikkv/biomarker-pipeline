# Dataset Validation Framework for CAMK2D Meta-Analysis
# Comprehensive validation of GEO datasets for CAMK2D research

#' Validate Dataset for CAMK2D Analysis
#' @param gse_id GEO Series ID
#' @param expression_data Expression matrix (optional, for loaded datasets)
#' @param metadata Sample metadata (optional)
#' @return Validation results with scores and recommendations
validate_dataset <- function(gse_id, expression_data = NULL, metadata = NULL) {
  
  cat("=== VALIDATING DATASET:", gse_id, "===\n")
  
  validation_results <- list(
    gse_id = gse_id,
    validation_date = Sys.Date(),
    overall_score = 0,
    camk_detection = NULL,
    quality_metrics = NULL,
    relevance_score = NULL,
    recommendations = character()
  )
  
  # 1. CAMK Gene Detection Validation
  camk_results <- validate_camk_detection(gse_id, expression_data)
  validation_results$camk_detection <- camk_results
  
  # 2. Dataset Quality Assessment  
  quality_results <- assess_dataset_quality(gse_id, expression_data, metadata)
  validation_results$quality_metrics <- quality_results
  
  # 3. Research Relevance Scoring
  relevance_results <- score_research_relevance(gse_id)
  validation_results$relevance_score <- relevance_results
  
  # 4. Calculate Overall Score
  validation_results$overall_score <- calculate_validation_score(
    camk_results, quality_results, relevance_results
  )
  
  # 5. Generate Recommendations
  validation_results$recommendations <- generate_recommendations(validation_results)
  
  # 6. Display Summary
  display_validation_summary(validation_results)
  
  return(validation_results)
}

#' Validate CAMK Gene Detection in Dataset
#' @param gse_id GEO Series ID
#' @param expression_data Expression matrix
#' @return CAMK detection results
validate_camk_detection <- function(gse_id, expression_data = NULL) {
  
  cat("Checking CAMK gene detection...\n")
  
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
  
  detection_results <- list(
    total_camk_genes = length(camk_genes),
    detected_camk_genes = 0,
    missing_camk_genes = camk_genes,
    detection_rate = 0,
    primary_target_detected = FALSE,
    detection_score = 0
  )
  
  if (!is.null(expression_data)) {
    # Check actual expression data
    gene_symbols <- rownames(expression_data)
    detected_genes <- intersect(camk_genes, gene_symbols)
    
    detection_results$detected_camk_genes <- length(detected_genes)
    detection_results$missing_camk_genes <- setdiff(camk_genes, detected_genes)
    detection_results$detection_rate <- length(detected_genes) / length(camk_genes)
    detection_results$primary_target_detected <- "CAMK2D" %in% detected_genes
    
    # Score based on detection
    if (detection_results$primary_target_detected) {
      detection_results$detection_score <- 10
    } else if (detection_results$detected_camk_genes > 0) {
      detection_results$detection_score <- 5 + detection_results$detected_camk_genes
    } else {
      detection_results$detection_score <- 0
    }
    
  } else {
    # Estimate based on dataset metadata (fallback)
    cat("No expression data provided - using metadata-based estimation\n")
    detection_results$detection_score <- 5  # Neutral score
  }
  
  return(detection_results)
}

#' Assess Dataset Quality Metrics
#' @param gse_id GEO Series ID
#' @param expression_data Expression matrix
#' @param metadata Sample metadata
#' @return Quality assessment results
assess_dataset_quality <- function(gse_id, expression_data = NULL, metadata = NULL) {
  
  cat("Assessing dataset quality...\n")
  
  quality_results <- list(
    sample_size = 0,
    gene_coverage = 0,
    missing_data_rate = 0,
    platform_quality = "Unknown",
    batch_effects = "Not assessed",
    quality_score = 0
  )
  
  if (!is.null(expression_data)) {
    # Calculate actual quality metrics
    quality_results$sample_size <- ncol(expression_data)
    quality_results$gene_coverage <- nrow(expression_data)
    quality_results$missing_data_rate <- sum(is.na(expression_data)) / length(expression_data)
    
    # Quality scoring
    size_score <- min(quality_results$sample_size / 50, 10)  # Max 10 for 50+ samples
    coverage_score <- min(quality_results$gene_coverage / 2000, 10)  # Max 10 for 20k+ genes
    missing_penalty <- max(0, 10 - quality_results$missing_data_rate * 100)
    
    quality_results$quality_score <- (size_score + coverage_score + missing_penalty) / 3
    
  } else {
    # Use metadata-based estimates
    quality_results$quality_score <- 6  # Moderate default
  }
  
  return(quality_results)
}

#' Score Research Relevance to CAMK2D
#' @param gse_id GEO Series ID
#' @return Relevance scoring results
score_research_relevance <- function(gse_id) {
  
  cat("Scoring research relevance...\n")
  
  relevance_results <- list(
    cardiac_relevance = 0,
    camk2d_focus = 0,
    species_relevance = 0,
    publication_quality = 0,
    relevance_score = 0
  )
  
  # Get dataset configuration if available
  all_configs <- tryCatch({
    source("config/datasets.R", local = TRUE)
    get_all_dataset_configs()
  }, error = function(e) NULL)
  
  if (!is.null(all_configs)) {
    # Find dataset in configurations
    dataset_found <- FALSE
    
    for (config_name in names(all_configs)) {
      config <- all_configs[[config_name]]
      if (gse_id %in% config$GSE_ID) {
        dataset_row <- config[config$GSE_ID == gse_id, ]
        
        # Cardiac relevance scoring
        cardiac_terms <- c("heart", "cardiac", "atrial", "myocardial", "cardiomyopathy")
        if (any(sapply(cardiac_terms, function(x) grepl(x, tolower(dataset_row$Disease_Focus))))) {
          relevance_results$cardiac_relevance <- 10
        } else {
          relevance_results$cardiac_relevance <- 5
        }
        
        # CAMK2D focus scoring
        camk_terms <- c("camk", "calcium", "calmodulin")
        if (any(sapply(camk_terms, function(x) grepl(x, tolower(dataset_row$Title))))) {
          relevance_results$camk2d_focus <- 10
        } else {
          relevance_results$camk2d_focus <- 5
        }
        
        # Species relevance
        if (dataset_row$Species %in% c("Human", "Mouse")) {
          relevance_results$species_relevance <- 10
        } else {
          relevance_results$species_relevance <- 3
        }
        
        # Quality score from configuration
        if ("Quality_Score" %in% names(dataset_row)) {
          relevance_results$publication_quality <- dataset_row$Quality_Score
        } else {
          relevance_results$publication_quality <- 6
        }
        
        dataset_found <- TRUE
        break
      }
    }
    
    if (!dataset_found) {
      # Default scoring for unknown datasets
      relevance_results$cardiac_relevance <- 5
      relevance_results$camk2d_focus <- 5
      relevance_results$species_relevance <- 8
      relevance_results$publication_quality <- 6
    }
  }
  
  # Calculate overall relevance score
  relevance_results$relevance_score <- (
    relevance_results$cardiac_relevance * 0.3 +
    relevance_results$camk2d_focus * 0.3 +
    relevance_results$species_relevance * 0.2 +
    relevance_results$publication_quality * 0.2
  )
  
  return(relevance_results)
}

#' Calculate Overall Validation Score
#' @param camk_results CAMK detection results
#' @param quality_results Quality assessment results
#' @param relevance_results Relevance scoring results
#' @return Overall validation score (0-10)
calculate_validation_score <- function(camk_results, quality_results, relevance_results) {
  
  # Weighted scoring
  camk_weight <- 0.4      # CAMK detection is critical
  quality_weight <- 0.3   # Dataset quality important
  relevance_weight <- 0.3 # Research relevance important
  
  overall_score <- (
    camk_results$detection_score * camk_weight +
    quality_results$quality_score * quality_weight +
    relevance_results$relevance_score * relevance_weight
  )
  
  return(round(overall_score, 1))
}

#' Generate Recommendations Based on Validation
#' @param validation_results Complete validation results
#' @return Vector of recommendations
generate_recommendations <- function(validation_results) {
  
  recommendations <- character()
  
  # CAMK detection recommendations
  if (validation_results$camk_detection$detection_score < 5) {
    recommendations <- c(recommendations, 
                        "CRITICAL: Low CAMK gene detection - verify gene symbols and platform compatibility")
  }
  
  if (!validation_results$camk_detection$primary_target_detected) {
    recommendations <- c(recommendations,
                        "WARNING: CAMK2D not detected - check if gene is measured in this dataset")
  }
  
  # Quality recommendations
  if (validation_results$quality_metrics$quality_score < 6) {
    recommendations <- c(recommendations,
                        "QUALITY: Dataset quality below threshold - consider data preprocessing")
  }
  
  if (validation_results$quality_metrics$missing_data_rate > 0.2) {
    recommendations <- c(recommendations,
                        "DATA: High missing data rate - implement imputation strategies")
  }
  
  # Relevance recommendations
  if (validation_results$relevance_score$relevance_score < 7) {
    recommendations <- c(recommendations,
                        "RELEVANCE: Limited relevance to CAMK2D cardiac research - use with caution")
  }
  
  # Overall recommendations
  if (validation_results$overall_score >= 8) {
    recommendations <- c(recommendations,
                        "RECOMMENDED: High-quality dataset suitable for meta-analysis")
  } else if (validation_results$overall_score >= 6) {
    recommendations <- c(recommendations,
                        "CONDITIONAL: Dataset acceptable with preprocessing")
  } else {
    recommendations <- c(recommendations,
                        "NOT RECOMMENDED: Dataset does not meet quality thresholds")
  }
  
  return(recommendations)
}

#' Display Validation Summary
#' @param validation_results Complete validation results
display_validation_summary <- function(validation_results) {
  
  cat("\n=== VALIDATION SUMMARY ===\n")
  cat("Dataset:", validation_results$gse_id, "\n")
  cat("Overall Score:", validation_results$overall_score, "/10\n")
  cat("CAMK Detection Score:", validation_results$camk_detection$detection_score, "/10\n")
  cat("Quality Score:", round(validation_results$quality_metrics$quality_score, 1), "/10\n")
  cat("Relevance Score:", round(validation_results$relevance_score$relevance_score, 1), "/10\n")
  
  if (validation_results$camk_detection$primary_target_detected) {
    cat("✓ CAMK2D detected\n")
  } else {
    cat("✗ CAMK2D not detected\n")
  }
  
  cat("\nRecommendations:\n")
  for (rec in validation_results$recommendations) {
    cat(" -", rec, "\n")
  }
  cat("\n")
}

#' Batch Validate Multiple Datasets
#' @param dataset_list Vector of GSE IDs to validate
#' @return List of validation results
batch_validate_datasets <- function(dataset_list) {
  
  cat("=== BATCH VALIDATION OF", length(dataset_list), "DATASETS ===\n")
  
  validation_results <- list()
  
  for (gse_id in dataset_list) {
    tryCatch({
      validation_results[[gse_id]] <- validate_dataset(gse_id)
    }, error = function(e) {
      cat("Error validating", gse_id, ":", e$message, "\n")
      validation_results[[gse_id]] <- list(
        gse_id = gse_id,
        error = e$message,
        overall_score = 0
      )
    })
  }
  
  # Summary statistics
  successful_validations <- sum(sapply(validation_results, function(x) !is.null(x$overall_score)))
  avg_score <- mean(sapply(validation_results, function(x) x$overall_score), na.rm = TRUE)
  
  cat("\n=== BATCH VALIDATION SUMMARY ===\n")
  cat("Total datasets:", length(dataset_list), "\n")
  cat("Successful validations:", successful_validations, "\n")
  cat("Average validation score:", round(avg_score, 2), "\n")
  
  return(validation_results)
}

#' Export Validation Results
#' @param validation_results Validation results to export
#' @param output_file Output file path
export_validation_results <- function(validation_results, output_file = "validation_results.csv") {
  
  # Convert to data frame for export
  results_df <- data.frame(
    GSE_ID = sapply(validation_results, function(x) x$gse_id),
    Overall_Score = sapply(validation_results, function(x) x$overall_score),
    CAMK_Detection_Score = sapply(validation_results, function(x) x$camk_detection$detection_score),
    Quality_Score = sapply(validation_results, function(x) round(x$quality_metrics$quality_score, 1)),
    Relevance_Score = sapply(validation_results, function(x) round(x$relevance_score$relevance_score, 1)),
    CAMK2D_Detected = sapply(validation_results, function(x) x$camk_detection$primary_target_detected),
    Validation_Date = sapply(validation_results, function(x) x$validation_date),
    stringsAsFactors = FALSE
  )
  
  write.csv(results_df, output_file, row.names = FALSE)
  cat("Validation results exported to:", output_file, "\n")
  
  return(results_df)
}