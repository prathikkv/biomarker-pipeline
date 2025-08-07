# Meta-Analysis Functions
# Advanced statistical meta-analysis functions for multi-dataset integration

#' Perform Meta-Analysis of CAMK Gene Expression
#' @param datasets_list List of datasets with expression data
#' @param target_gene Target gene for meta-analysis
#' @param method Meta-analysis method ("random" or "fixed")
#' @return Meta-analysis results object
perform_camk_meta_analysis <- function(datasets_list, target_gene = "CAMK2D", method = "random") {
  cat("Performing meta-analysis for", target_gene, "using", method, "effects model...\n")
  
  # Check if metafor is available
  if (!requireNamespace("metafor", quietly = TRUE)) {
    warning("metafor package not available. Using simplified meta-analysis.")
    return(perform_simple_meta_analysis(datasets_list, target_gene))
  }
  
  library(metafor)
  
  # Extract effect sizes from each dataset
  effect_sizes <- extract_effect_sizes(datasets_list, target_gene)
  
  if (nrow(effect_sizes) < 3) {
    warning("Insufficient datasets for robust meta-analysis: ", nrow(effect_sizes))
    return(NULL)
  }
  
  # Perform meta-analysis
  tryCatch({
    meta_result <- rma(yi = effect_sizes$effect_size,
                       vi = effect_sizes$variance,
                       method = ifelse(method == "random", "REML", "FE"),
                       data = effect_sizes)
    
    # Calculate additional statistics
    heterogeneity_stats <- calculate_heterogeneity_stats(meta_result)
    
    # Perform publication bias tests
    bias_tests <- test_publication_bias(effect_sizes)
    
    # Create comprehensive results object
    results <- list(
      target_gene = target_gene,
      method = method,
      meta_result = meta_result,
      effect_sizes = effect_sizes,
      heterogeneity = heterogeneity_stats,
      publication_bias = bias_tests,
      n_studies = nrow(effect_sizes),
      total_samples = sum(effect_sizes$sample_size, na.rm = TRUE)
    )
    
    cat("✓ Meta-analysis completed for", target_gene, "\n")
    cat("  - Studies included:", nrow(effect_sizes), "\n")
    cat("  - Total samples:", sum(effect_sizes$sample_size, na.rm = TRUE), "\n")
    cat("  - Overall effect:", round(meta_result$beta, 3), "\n")
    cat("  - Heterogeneity I²:", round(heterogeneity_stats$I2, 1), "%\n")
    
    return(results)
    
  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(perform_simple_meta_analysis(datasets_list, target_gene))
  })
}

#' Extract Effect Sizes from Datasets
#' @param datasets_list List of datasets
#' @param target_gene Target gene
#' @return Data frame with effect sizes and variances
extract_effect_sizes <- function(datasets_list, target_gene) {
  effect_data <- data.frame(
    study_id = character(),
    effect_size = numeric(),
    variance = numeric(),
    sample_size = numeric(),
    species = character(),
    platform = character(),
    disease_focus = character(),
    stringsAsFactors = FALSE
  )
  
  for (study_id in names(datasets_list)) {
    dataset <- datasets_list[[study_id]]
    
    # Check if target gene was detected
    if (is.null(dataset$camk_detection_results) || 
        !target_gene %in% names(dataset$camk_detection_results) ||
        !dataset$camk_detection_results[[target_gene]]$found) {
      next
    }
    
    # Calculate effect size (mean expression as proxy)
    gene_indices <- dataset$camk_detection_results[[target_gene]]$indices
    if (length(gene_indices) == 0) next
    
    gene_expression <- dataset$expression[gene_indices, , drop = FALSE]
    if (nrow(gene_expression) == 1) {
      expr_values <- as.numeric(gene_expression[1, ])
    } else {
      expr_values <- colMeans(gene_expression, na.rm = TRUE)
    }
    
    # Calculate effect size metrics
    mean_expr <- mean(expr_values, na.rm = TRUE)
    sd_expr <- sd(expr_values, na.rm = TRUE)
    n_samples <- length(expr_values)
    
    # Use standardized mean as effect size (log-transformed)
    effect_size <- log2(mean_expr + 1)  # Log2 transformation
    variance <- (sd_expr^2) / n_samples  # Approximate variance
    
    # Add to effect data
    effect_data <- rbind(effect_data, data.frame(
      study_id = study_id,
      effect_size = effect_size,
      variance = variance,
      sample_size = n_samples,
      species = if (!is.null(dataset$species)) dataset$species else "Unknown",
      platform = if (!is.null(dataset$platform)) dataset$platform else "Unknown",
      disease_focus = if (!is.null(dataset$disease_focus)) dataset$disease_focus else "Unknown",
      stringsAsFactors = FALSE
    ))
  }
  
  # Remove invalid entries
  effect_data <- effect_data[!is.na(effect_data$effect_size) & 
                           !is.na(effect_data$variance) & 
                           effect_data$variance > 0, ]
  
  cat("Extracted effect sizes for", nrow(effect_data), "studies\n")
  
  return(effect_data)
}

#' Calculate Heterogeneity Statistics
#' @param meta_result Meta-analysis result from rma()
#' @return List with heterogeneity statistics
calculate_heterogeneity_stats <- function(meta_result) {
  # Extract heterogeneity statistics
  I2 <- max(0, (meta_result$QE - meta_result$k + 1) / meta_result$QE * 100)
  tau2 <- meta_result$tau2
  Q_stat <- meta_result$QE
  Q_pval <- meta_result$QEp
  
  # Interpret heterogeneity
  heterogeneity_interpretation <- ifelse(I2 < 25, "Low",
                                       ifelse(I2 < 50, "Moderate",
                                             ifelse(I2 < 75, "Substantial", "Considerable")))
  
  return(list(
    I2 = I2,
    tau2 = tau2,
    Q_statistic = Q_stat,
    Q_pvalue = Q_pval,
    interpretation = heterogeneity_interpretation
  ))
}

#' Test for Publication Bias
#' @param effect_sizes Data frame with effect sizes
#' @return List with publication bias test results
test_publication_bias <- function(effect_sizes) {
  if (nrow(effect_sizes) < 10) {
    return(list(
      egger_test = "Not performed (< 10 studies)",
      rank_correlation = "Not performed (< 10 studies)",
      funnel_plot_data = effect_sizes
    ))
  }
  
  tryCatch({
    # Egger's regression test
    egger_result <- lm(effect_sizes$effect_size / sqrt(effect_sizes$variance) ~ 
                       I(1/sqrt(effect_sizes$variance)))
    egger_pvalue <- summary(egger_result)$coefficients[1, 4]
    
    # Begg's rank correlation test
    rank_test <- cor.test(effect_sizes$effect_size, effect_sizes$variance, method = "kendall")
    
    return(list(
      egger_test = list(
        coefficient = coef(egger_result)[1],
        pvalue = egger_pvalue,
        significant = egger_pvalue < 0.05
      ),
      rank_correlation = list(
        tau = rank_test$estimate,
        pvalue = rank_test$p.value,
        significant = rank_test$p.value < 0.05
      ),
      funnel_plot_data = effect_sizes
    ))
    
  }, error = function(e) {
    return(list(
      egger_test = paste("Failed:", e$message),
      rank_correlation = paste("Failed:", e$message),
      funnel_plot_data = effect_sizes
    ))
  })
}

#' Perform Simplified Meta-Analysis (Fallback)
#' @param datasets_list List of datasets
#' @param target_gene Target gene
#' @return Simplified meta-analysis results
perform_simple_meta_analysis <- function(datasets_list, target_gene) {
  cat("Performing simplified meta-analysis (metafor not available)...\n")
  
  effect_sizes <- extract_effect_sizes(datasets_list, target_gene)
  
  if (nrow(effect_sizes) < 2) {
    return(NULL)
  }
  
  # Simple weighted average
  weights <- 1 / effect_sizes$variance
  weighted_mean <- sum(effect_sizes$effect_size * weights) / sum(weights)
  
  # Simple heterogeneity assessment
  Q_stat <- sum(weights * (effect_sizes$effect_size - weighted_mean)^2)
  I2 <- max(0, (Q_stat - nrow(effect_sizes) + 1) / Q_stat * 100)
  
  return(list(
    target_gene = target_gene,
    method = "simple_weighted",
    effect_sizes = effect_sizes,
    overall_effect = weighted_mean,
    heterogeneity = list(I2 = I2, Q_statistic = Q_stat),
    n_studies = nrow(effect_sizes),
    total_samples = sum(effect_sizes$sample_size, na.rm = TRUE)
  ))
}

#' Perform Cross-Species Meta-Analysis
#' @param datasets_list List of datasets (human and mouse)
#' @param target_gene Target gene
#' @param mapping_info Cross-species mapping information
#' @return Cross-species meta-analysis results
perform_cross_species_meta_analysis <- function(datasets_list, target_gene, mapping_info) {
  cat("Performing cross-species meta-analysis for", target_gene, "...\n")
  
  # Separate human and mouse datasets
  human_indices <- sapply(datasets_list, function(x) {
    !is.null(x$species) && x$species == "Human"
  })
  
  mouse_indices <- sapply(datasets_list, function(x) {
    !is.null(x$species) && x$species == "Mouse"
  })
  
  human_datasets <- datasets_list[human_indices]
  mouse_datasets <- datasets_list[mouse_indices]
  
  # Perform species-specific meta-analyses
  human_meta <- perform_camk_meta_analysis(human_datasets, target_gene)
  mouse_meta <- perform_camk_meta_analysis(mouse_datasets, target_gene)
  
  # Calculate cross-species concordance
  concordance <- calculate_cross_species_concordance(human_meta, mouse_meta)
  
  return(list(
    target_gene = target_gene,
    human_meta_analysis = human_meta,
    mouse_meta_analysis = mouse_meta,
    cross_species_concordance = concordance,
    combined_sample_size = ifelse(!is.null(human_meta), human_meta$total_samples, 0) + 
                          ifelse(!is.null(mouse_meta), mouse_meta$total_samples, 0)
  ))
}

#' Calculate Cross-Species Concordance
#' @param human_meta Human meta-analysis results
#' @param mouse_meta Mouse meta-analysis results  
#' @return Concordance assessment
calculate_cross_species_concordance <- function(human_meta, mouse_meta) {
  if (is.null(human_meta) || is.null(mouse_meta)) {
    return(list(
      concordance = "Cannot assess - missing species data",
      correlation = NA,
      direction_agreement = NA
    ))
  }
  
  # Extract effect sizes
  human_effect <- if (!is.null(human_meta$meta_result)) human_meta$meta_result$beta else human_meta$overall_effect
  mouse_effect <- if (!is.null(mouse_meta$meta_result)) mouse_meta$meta_result$beta else mouse_meta$overall_effect
  
  # Check direction agreement
  direction_agreement <- (human_effect > 0 && mouse_effect > 0) || (human_effect < 0 && mouse_effect < 0)
  
  # Simple correlation assessment (would need individual study data for true correlation)
  effect_correlation <- NA  # Placeholder
  
  concordance_level <- ifelse(direction_agreement, 
                             ifelse(abs(human_effect - mouse_effect) < 1, "High", "Moderate"),
                             "Low")
  
  return(list(
    concordance = concordance_level,
    direction_agreement = direction_agreement,
    human_effect_size = human_effect,
    mouse_effect_size = mouse_effect,
    effect_difference = abs(human_effect - mouse_effect)
  ))
}

#' Create Forest Plot Data
#' @param meta_results Meta-analysis results
#' @return Data frame formatted for forest plot creation
create_forest_plot_data <- function(meta_results) {
  if (is.null(meta_results) || is.null(meta_results$effect_sizes)) {
    return(NULL)
  }
  
  effect_data <- meta_results$effect_sizes
  
  # Add confidence intervals
  effect_data$ci_lower <- effect_data$effect_size - 1.96 * sqrt(effect_data$variance)
  effect_data$ci_upper <- effect_data$effect_size + 1.96 * sqrt(effect_data$variance)
  
  # Add weights
  if (!is.null(meta_results$meta_result)) {
    effect_data$weight <- weights(meta_results$meta_result)
  } else {
    effect_data$weight <- 1 / effect_data$variance
    effect_data$weight <- effect_data$weight / sum(effect_data$weight) * 100
  }
  
  # Add overall summary row
  if (!is.null(meta_results$meta_result)) {
    summary_row <- data.frame(
      study_id = "SUMMARY",
      effect_size = meta_results$meta_result$beta,
      variance = meta_results$meta_result$se^2,
      sample_size = sum(effect_data$sample_size),
      species = "Combined",
      platform = "Meta-analysis",
      disease_focus = "Combined",
      ci_lower = meta_results$meta_result$ci.lb,
      ci_upper = meta_results$meta_result$ci.ub,
      weight = 100,
      stringsAsFactors = FALSE
    )
    
    effect_data <- rbind(effect_data, summary_row)
  }
  
  return(effect_data)
}

#' Generate Meta-Analysis Summary Statistics
#' @param meta_results Meta-analysis results
#' @return Summary statistics data frame
generate_meta_summary_stats <- function(meta_results) {
  if (is.null(meta_results)) {
    return(data.frame(Metric = "No results", Value = "N/A"))
  }
  
  summary_stats <- data.frame(
    Metric = c(
      "Target Gene",
      "Number of Studies", 
      "Total Sample Size",
      "Overall Effect Size",
      "95% CI Lower",
      "95% CI Upper",
      "Heterogeneity I²",
      "Q-statistic P-value",
      "Method Used"
    ),
    Value = c(
      meta_results$target_gene,
      meta_results$n_studies,
      meta_results$total_samples,
      if (!is.null(meta_results$meta_result)) round(meta_results$meta_result$beta, 3) else round(meta_results$overall_effect, 3),
      if (!is.null(meta_results$meta_result)) round(meta_results$meta_result$ci.lb, 3) else "N/A",
      if (!is.null(meta_results$meta_result)) round(meta_results$meta_result$ci.ub, 3) else "N/A",
      paste0(round(meta_results$heterogeneity$I2, 1), "%"),
      if (!is.null(meta_results$heterogeneity$Q_pvalue)) round(meta_results$heterogeneity$Q_pvalue, 3) else "N/A",
      meta_results$method
    ),
    stringsAsFactors = FALSE
  )
  
  return(summary_stats)
}