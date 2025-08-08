# Enhanced Differential Gene Expression Analysis Pipeline
# Scientific rigor focused implementation for CAMK2D biomarker study

#' Enhanced Differential Gene Expression Analysis
#' Implements comprehensive DGE workflow with statistical validation
#' @param expression_matrix Expression data matrix (genes x samples)
#' @param group_vector Group assignment vector (disease vs control)
#' @param dataset_id Dataset identifier for tracking
#' @param platform Platform type (RNA-seq vs microarray)
#' @param method Statistical method ("limma", "deseq2", "edger")
#' @return Comprehensive DGE results with quality metrics
perform_rigorous_dge <- function(expression_matrix, group_vector, dataset_id, 
                                platform = "RNA-seq", method = "limma") {
  
  cat("🧬 Performing rigorous DGE analysis for", dataset_id, "\n")
  
  # Input validation and quality control
  qc_results <- validate_dge_inputs(expression_matrix, group_vector)
  if (!qc_results$valid) {
    warning("QC failed for ", dataset_id, ": ", qc_results$issues)
    return(NULL)
  }
  
  # Data preprocessing based on platform
  processed_data <- preprocess_expression_data(expression_matrix, platform)
  
  # Statistical analysis
  if (method == "limma") {
    dge_results <- perform_limma_analysis(processed_data$expression, group_vector)
  } else if (method == "deseq2" && platform == "RNA-seq") {
    dge_results <- perform_deseq2_analysis(processed_data$counts, group_vector)
  } else {
    dge_results <- perform_limma_analysis(processed_data$expression, group_vector)
  }
  
  # Extract CAMK family results with enhanced annotation
  camk_results <- extract_camk_family_results(dge_results, dataset_id)
  
  # Quality assessment and statistical validation
  quality_metrics <- assess_dge_quality(dge_results, processed_data, group_vector)
  
  # Return comprehensive results
  list(
    dataset_id = dataset_id,
    platform = platform,
    method = method,
    sample_size = list(
      total = length(group_vector),
      disease = sum(group_vector == "disease"),
      control = sum(group_vector == "control")
    ),
    dge_results = dge_results,
    camk_results = camk_results,
    quality_metrics = quality_metrics,
    preprocessing_info = processed_data$info
  )
}

#' Comprehensive Input Validation for DGE Analysis
#' @param expression_matrix Expression data matrix
#' @param group_vector Group assignment vector
#' @return List with validation status and issues
validate_dge_inputs <- function(expression_matrix, group_vector) {
  issues <- character(0)
  
  # Check dimensions
  if (ncol(expression_matrix) != length(group_vector)) {
    issues <- c(issues, "Sample count mismatch between expression and groups")
  }
  
  # Check group balance
  group_table <- table(group_vector)
  if (any(group_table < 3)) {
    issues <- c(issues, "Insufficient samples in one or more groups (min 3 required)")
  }
  
  # Check for missing data
  na_prop <- sum(is.na(expression_matrix)) / prod(dim(expression_matrix))
  if (na_prop > 0.1) {
    issues <- c(issues, sprintf("High missing data proportion: %.1f%%", na_prop * 100))
  }
  
  # Check expression range (detect log-transformed data)
  expr_range <- range(expression_matrix, na.rm = TRUE)
  if (expr_range[2] - expr_range[1] > 50) {
    issues <- c(issues, "Possible non-log transformed data detected")
  }
  
  list(
    valid = length(issues) == 0,
    issues = paste(issues, collapse = "; ")
  )
}

#' Platform-Specific Data Preprocessing
#' @param expression_matrix Raw expression matrix
#' @param platform Platform type
#' @return Preprocessed data with metadata
preprocess_expression_data <- function(expression_matrix, platform) {
  
  preprocessing_info <- list()
  
  if (platform == "RNA-seq") {
    # RNA-seq specific preprocessing
    
    # Filter low-expressed genes (CPM > 1 in at least 3 samples)
    keep_genes <- rowSums(expression_matrix > 1) >= 3
    filtered_matrix <- expression_matrix[keep_genes, ]
    preprocessing_info$genes_filtered <- sum(!keep_genes)
    
    # Log2 transformation if data appears to be raw counts
    if (max(filtered_matrix, na.rm = TRUE) > 100) {
      processed_matrix <- log2(filtered_matrix + 1)
      preprocessing_info$log_transformed <- TRUE
    } else {
      processed_matrix <- filtered_matrix
      preprocessing_info$log_transformed <- FALSE
    }
    
  } else {
    # Microarray preprocessing
    
    # Remove genes with consistently low expression
    gene_means <- rowMeans(expression_matrix, na.rm = TRUE)
    keep_genes <- gene_means > quantile(gene_means, 0.1, na.rm = TRUE)
    processed_matrix <- expression_matrix[keep_genes, ]
    preprocessing_info$genes_filtered <- sum(!keep_genes)
    
    # Quantile normalization if needed
    if (requireNamespace("limma", quietly = TRUE)) {
      processed_matrix <- limma::normalizeQuantiles(processed_matrix)
      preprocessing_info$quantile_normalized <- TRUE
    }
  }
  
  preprocessing_info$final_genes <- nrow(processed_matrix)
  preprocessing_info$final_samples <- ncol(processed_matrix)
  
  list(
    expression = processed_matrix,
    counts = if (platform == "RNA-seq") round(2^processed_matrix - 1) else NULL,
    info = preprocessing_info
  )
}

#' LIMMA-based Differential Expression Analysis
#' @param expression_matrix Preprocessed expression matrix
#' @param group_vector Group assignments
#' @return LIMMA results with comprehensive statistics
perform_limma_analysis <- function(expression_matrix, group_vector) {
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package required for differential expression analysis")
  }
  
  # Create design matrix
  design_matrix <- model.matrix(~ group_vector)
  colnames(design_matrix) <- c("Intercept", "Disease_vs_Control")
  
  # Fit linear model
  fit <- limma::lmFit(expression_matrix, design_matrix)
  fit <- limma::eBayes(fit, trend = TRUE)
  
  # Extract results
  results <- limma::topTable(fit, coef = "Disease_vs_Control", 
                            number = Inf, sort.by = "P")
  
  # Add additional statistics
  results$Gene_Symbol <- rownames(results)
  results$Significant <- results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5
  results$Direction <- ifelse(results$logFC > 0, "Up", "Down")
  results$Effect_Size <- abs(results$logFC)
  
  # Calculate confidence intervals
  results$CI_Lower <- results$logFC - 1.96 * sqrt(fit$s2.post) * fit$stdev.unscaled[, "Disease_vs_Control"]
  results$CI_Upper <- results$logFC + 1.96 * sqrt(fit$s2.post) * fit$stdev.unscaled[, "Disease_vs_Control"]
  
  return(results)
}

#' Extract CAMK Family Results with Enhanced Annotation
#' @param dge_results Full DGE results
#' @param dataset_id Dataset identifier
#' @return CAMK-specific results with annotations
extract_camk_family_results <- function(dge_results, dataset_id) {
  
  # Define comprehensive CAMK gene list
  camk_genes <- c(
    # Primary CAMK2 family
    "CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G",
    # Related CAMK family
    "CAMK1", "CAMK4", "CAMK1D", "CAMK1G",
    # Upstream regulators
    "CAMKK1", "CAMKK2",
    # Alternative names and variants
    "CAMK2", "CAMKII", "CAMKIID"
  )
  
  # Extract CAMK results
  camk_mask <- dge_results$Gene_Symbol %in% camk_genes | 
               grepl("CAMK", dge_results$Gene_Symbol, ignore.case = TRUE)
  
  camk_results <- dge_results[camk_mask, ]
  
  if (nrow(camk_results) > 0) {
    # Add annotation and prioritization
    camk_results$Dataset_ID <- dataset_id
    camk_results$CAMK_Type <- classify_camk_gene(camk_results$Gene_Symbol)
    camk_results$Clinical_Priority <- assign_clinical_priority(camk_results)
    camk_results$Evidence_Strength <- calculate_evidence_strength(camk_results)
    
    # Sort by significance and effect size
    camk_results <- camk_results[order(camk_results$adj.P.Val, -abs(camk_results$logFC)), ]
  }
  
  return(camk_results)
}

#' Classify CAMK Gene Types
#' @param gene_symbols Vector of gene symbols
#' @return Classification vector
classify_camk_gene <- function(gene_symbols) {
  sapply(gene_symbols, function(gene) {
    if (gene %in% c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G")) {
      return("CAMK2_Primary")
    } else if (gene %in% c("CAMK1", "CAMK4", "CAMK1D", "CAMK1G")) {
      return("CAMK_Related")
    } else if (gene %in% c("CAMKK1", "CAMKK2")) {
      return("CAMK_Upstream")
    } else {
      return("CAMK_Other")
    }
  })
}

#' Assign Clinical Priority Scores
#' @param camk_results CAMK-specific DGE results
#' @return Priority scores
assign_clinical_priority <- function(camk_results) {
  
  priority_scores <- numeric(nrow(camk_results))
  
  for (i in seq_len(nrow(camk_results))) {
    gene <- camk_results$Gene_Symbol[i]
    pval <- camk_results$adj.P.Val[i]
    effect_size <- abs(camk_results$logFC[i])
    
    # Base score from statistical significance
    if (pval < 0.001) score <- 10
    else if (pval < 0.01) score <- 8
    else if (pval < 0.05) score <- 6
    else score <- 3
    
    # Boost for large effect sizes
    if (effect_size > 1.5) score <- score + 3
    else if (effect_size > 1.0) score <- score + 2
    else if (effect_size > 0.5) score <- score + 1
    
    # Gene-specific clinical relevance
    if (gene == "CAMK2D") score <- score + 5  # Primary target
    else if (gene %in% c("CAMK2A", "CAMK2B")) score <- score + 3
    else if (gene == "CAMKK2") score <- score + 2  # Druggable upstream target
    
    priority_scores[i] <- min(score, 20)  # Cap at 20
  }
  
  return(priority_scores)
}

#' Calculate Evidence Strength
#' @param camk_results CAMK-specific results
#' @return Evidence strength scores
calculate_evidence_strength <- function(camk_results) {
  
  # Combine p-value, effect size, and statistical power
  strength_scores <- with(camk_results, {
    
    # Statistical significance component (0-5)
    sig_component <- pmax(0, 5 * (1 - adj.P.Val))
    
    # Effect size component (0-3)
    effect_component <- pmin(3, abs(logFC))
    
    # Statistical power component based on t-statistic (0-2)
    power_component <- pmin(2, abs(t) / 5)
    
    # Combine components
    total_score <- sig_component + effect_component + power_component
    
    # Normalize to 0-10 scale
    10 * total_score / max(total_score, na.rm = TRUE)
  })
  
  return(round(strength_scores, 2))
}

#' Assess DGE Analysis Quality
#' @param dge_results Full DGE results
#' @param processed_data Preprocessed data
#' @param group_vector Group assignments
#' @return Quality assessment metrics
assess_dge_quality <- function(dge_results, processed_data, group_vector) {
  
  quality_metrics <- list()
  
  # Statistical power assessment
  quality_metrics$sample_size <- length(group_vector)
  quality_metrics$group_balance <- min(table(group_vector)) / max(table(group_vector))
  
  # Multiple testing burden
  quality_metrics$genes_tested <- nrow(dge_results)
  quality_metrics$significant_genes <- sum(dge_results$adj.P.Val < 0.05, na.rm = TRUE)
  quality_metrics$discovery_rate <- quality_metrics$significant_genes / quality_metrics$genes_tested
  
  # Effect size distribution
  quality_metrics$median_effect_size <- median(abs(dge_results$logFC), na.rm = TRUE)
  quality_metrics$large_effects <- sum(abs(dge_results$logFC) > 1, na.rm = TRUE)
  
  # P-value distribution (should be uniform under null hypothesis)
  pval_hist <- hist(dge_results$P.Value, breaks = 20, plot = FALSE)
  quality_metrics$pvalue_inflation <- max(pval_hist$counts) / min(pval_hist$counts)
  
  # Data quality indicators
  quality_metrics$preprocessing_quality <- processed_data$info
  quality_metrics$overall_score <- calculate_overall_quality_score(quality_metrics)
  
  return(quality_metrics)
}

#' Calculate Overall Quality Score
#' @param metrics Quality metrics list
#' @return Overall quality score (0-10)
calculate_overall_quality_score <- function(metrics) {
  
  # Sample size score (0-3)
  size_score <- pmin(3, metrics$sample_size / 20)
  
  # Balance score (0-2)
  balance_score <- 2 * metrics$group_balance
  
  # Discovery score (0-3)
  discovery_score <- pmin(3, 10 * metrics$discovery_rate)
  
  # Effect size score (0-2)
  effect_score <- pmin(2, 2 * metrics$median_effect_size)
  
  total_score <- size_score + balance_score + discovery_score + effect_score
  return(round(total_score, 2))
}

#' Generate DGE Summary Report
#' @param dge_analysis_list List of DGE analysis results
#' @return Formatted summary table
generate_dge_summary <- function(dge_analysis_list) {
  
  if (length(dge_analysis_list) == 0) {
    return(data.frame(Message = "No successful DGE analyses to summarize"))
  }
  
  summary_data <- do.call(rbind, lapply(dge_analysis_list, function(analysis) {
    if (is.null(analysis)) return(NULL)
    
    camk_sig <- sum(analysis$camk_results$Significant, na.rm = TRUE)
    
    data.frame(
      Dataset_ID = analysis$dataset_id,
      Platform = analysis$platform,
      Method = analysis$method,
      Total_Samples = analysis$sample_size$total,
      Disease_Samples = analysis$sample_size$disease,
      Control_Samples = analysis$sample_size$control,
      Genes_Tested = nrow(analysis$dge_results),
      Significant_Genes = sum(analysis$dge_results$Significant, na.rm = TRUE),
      CAMK_Genes_Found = nrow(analysis$camk_results),
      CAMK_Significant = camk_sig,
      Quality_Score = analysis$quality_metrics$overall_score,
      stringsAsFactors = FALSE
    )
  }))
  
  # Add summary statistics
  if (nrow(summary_data) > 0) {
    summary_data$Success_Rate <- ifelse(summary_data$CAMK_Genes_Found > 0, "Success", "No CAMK")
    summary_data$Power_Estimate <- with(summary_data, 
      pmin(100, 100 * Total_Samples / 50))  # Rough power estimate
  }
  
  return(summary_data)
}