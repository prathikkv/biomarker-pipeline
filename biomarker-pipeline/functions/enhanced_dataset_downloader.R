# Enhanced Dataset Downloader and Validator
# Scientifically rigorous data acquisition and validation

#' Enhanced Dataset Download with Comprehensive Validation
#' Downloads and validates datasets with quality control
#' @param geo_id GEO dataset identifier
#' @param force_download Force re-download even if cached
#' @param validate_immediately Perform immediate validation
#' @return Enhanced dataset object with validation results
download_and_validate_dataset <- function(geo_id, force_download = FALSE, 
                                        validate_immediately = TRUE) {
  
  cat("📥 Downloading and validating", geo_id, "\n")
  
  tryCatch({
    
    # Download dataset with enhanced error handling
    gse_data <- download_geo_dataset(geo_id, force_download)
    
    if (is.null(gse_data)) {
      return(create_failed_dataset_result(geo_id, "Download failed"))
    }
    
    # Immediate validation if requested
    if (validate_immediately) {
      validation_results <- validate_dataset_comprehensively(gse_data, geo_id)
      gse_data$validation <- validation_results
    }
    
    # Extract and enhance metadata
    enhanced_metadata <- extract_enhanced_metadata(gse_data, geo_id)
    gse_data$enhanced_metadata <- enhanced_metadata
    
    cat("✓ Successfully processed", geo_id, "\n")
    return(gse_data)
    
  }, error = function(e) {
    cat("✗ Error processing", geo_id, ":", e$message, "\n")
    return(create_failed_dataset_result(geo_id, e$message))
  })
}

#' Download GEO Dataset with Enhanced Error Handling
#' @param geo_id GEO identifier
#' @param force_download Force re-download
#' @return GEO dataset object or NULL
download_geo_dataset <- function(geo_id, force_download = FALSE) {
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEOquery package required for data download")
  }
  
  # Set up caching directory
  cache_dir <- file.path(getwd(), "data", "geo_cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Check cache first
  cache_file <- file.path(cache_dir, paste0(geo_id, ".rds"))
  
  if (!force_download && file.exists(cache_file)) {
    cat("🔄 Loading", geo_id, "from cache\n")
    return(readRDS(cache_file))
  }
  
  # Download from GEO
  cat("🌐 Downloading", geo_id, "from GEO database\n")
  
  # Multiple attempts with different strategies
  for (attempt in 1:3) {
    
    tryCatch({
      
      # Primary download method
      if (attempt == 1) {
        gse_data <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE, 
                                    AnnotGPL = TRUE, destdir = cache_dir)
      }
      # Fallback method without annotation
      else if (attempt == 2) {
        gse_data <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE, 
                                    AnnotGPL = FALSE, destdir = cache_dir)
      }
      # Last resort - minimal download
      else {
        gse_data <- GEOquery::getGEO(geo_id, GSEMatrix = FALSE, destdir = cache_dir)
      }
      
      # Validate download
      if (is.null(gse_data) || length(gse_data) == 0) {
        stop("Empty dataset returned")
      }
      
      # Cache successful download
      saveRDS(gse_data, cache_file)
      cat("💾 Cached", geo_id, "for future use\n")
      
      return(gse_data)
      
    }, error = function(e) {
      cat("⚠ Attempt", attempt, "failed:", e$message, "\n")
      if (attempt == 3) {
        return(NULL)
      }
      Sys.sleep(2^attempt)  # Exponential backoff
    })
  }
  
  return(NULL)
}

#' Comprehensive Dataset Validation
#' @param gse_data Downloaded GEO dataset
#' @param geo_id Dataset identifier
#' @return Comprehensive validation results
validate_dataset_comprehensively <- function(gse_data, geo_id) {
  
  validation_results <- list(
    dataset_id = geo_id,
    timestamp = Sys.time(),
    overall_valid = TRUE,
    issues = character(0),
    warnings = character(0)
  )
  
  tryCatch({
    
    # Extract first GSE matrix (most common case)
    if (is.list(gse_data) && length(gse_data) > 0) {
      gse_matrix <- gse_data[[1]]
    } else {
      gse_matrix <- gse_data
    }
    
    # Basic structure validation
    structure_results <- validate_dataset_structure(gse_matrix, geo_id)
    validation_results <- merge_validation_results(validation_results, structure_results)
    
    # Expression data validation
    if (validation_results$overall_valid) {
      expression_results <- validate_expression_data(gse_matrix, geo_id)
      validation_results <- merge_validation_results(validation_results, expression_results)
    }
    
    # Phenotype data validation
    if (validation_results$overall_valid) {
      pheno_results <- validate_phenotype_data(gse_matrix, geo_id)
      validation_results <- merge_validation_results(validation_results, pheno_results)
    }
    
    # CAMK gene detection
    if (validation_results$overall_valid) {
      camk_results <- detect_camk_genes_in_dataset(gse_matrix, geo_id)
      validation_results$camk_detection <- camk_results
    }
    
    # Calculate overall score
    validation_results$quality_score <- calculate_dataset_quality_score(validation_results)
    
  }, error = function(e) {
    validation_results$overall_valid <- FALSE
    validation_results$issues <- c(validation_results$issues, 
                                 paste("Validation error:", e$message))
  })
  
  return(validation_results)
}

#' Validate Dataset Structure
#' @param gse_matrix GEO matrix object
#' @param geo_id Dataset identifier
#' @return Structure validation results
validate_dataset_structure <- function(gse_matrix, geo_id) {
  
  results <- list(overall_valid = TRUE, issues = character(0), warnings = character(0))
  
  # Check if it's an ExpressionSet
  if (!methods::is(gse_matrix, "ExpressionSet")) {
    results$issues <- c(results$issues, "Not a valid ExpressionSet object")
    results$overall_valid <- FALSE
    return(results)
  }
  
  # Check dimensions
  expr_matrix <- Biobase::exprs(gse_matrix)
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    results$issues <- c(results$issues, "Empty expression matrix")
    results$overall_valid <- FALSE
    return(results)
  }
  
  # Check sample size
  n_samples <- ncol(expr_matrix)
  if (n_samples < 6) {
    results$issues <- c(results$issues, paste("Too few samples:", n_samples, "< 6"))
    results$overall_valid <- FALSE
  } else if (n_samples < 20) {
    results$warnings <- c(results$warnings, paste("Small sample size:", n_samples))
  }
  
  # Check gene count
  n_genes <- nrow(expr_matrix)
  if (n_genes < 1000) {
    results$warnings <- c(results$warnings, paste("Few genes detected:", n_genes))
  }
  
  results$dimensions <- list(samples = n_samples, genes = n_genes)
  return(results)
}

#' Validate Expression Data Quality
#' @param gse_matrix GEO matrix object
#' @param geo_id Dataset identifier
#' @return Expression validation results
validate_expression_data <- function(gse_matrix, geo_id) {
  
  results <- list(overall_valid = TRUE, issues = character(0), warnings = character(0))
  
  expr_matrix <- Biobase::exprs(gse_matrix)
  
  # Check for missing data
  na_proportion <- sum(is.na(expr_matrix)) / length(expr_matrix)
  if (na_proportion > 0.5) {
    results$issues <- c(results$issues, paste("High missing data:", round(na_proportion * 100, 1), "%"))
    results$overall_valid <- FALSE
  } else if (na_proportion > 0.1) {
    results$warnings <- c(results$warnings, paste("Moderate missing data:", round(na_proportion * 100, 1), "%"))
  }
  
  # Check expression range and distribution
  expr_range <- range(expr_matrix, na.rm = TRUE)
  expr_median <- median(expr_matrix, na.rm = TRUE)
  
  if (diff(expr_range) < 0.1) {
    results$issues <- c(results$issues, "Very low expression variability")
    results$overall_valid <- FALSE
  }
  
  # Detect platform type based on expression values
  if (expr_range[2] > 50 && expr_median > 5) {
    platform_type <- "Likely raw counts (RNA-seq)"
  } else if (expr_range[1] >= 0 && expr_range[2] < 20) {
    platform_type <- "Likely log-transformed (Microarray/RNA-seq)"
  } else {
    platform_type <- "Unknown platform type"
    results$warnings <- c(results$warnings, "Cannot determine platform type from expression values")
  }
  
  results$expression_stats <- list(
    range = expr_range,
    median = expr_median,
    na_proportion = na_proportion,
    platform_type = platform_type
  )
  
  return(results)
}

#' Validate Phenotype Data
#' @param gse_matrix GEO matrix object
#' @param geo_id Dataset identifier
#' @return Phenotype validation results
validate_phenotype_data <- function(gse_matrix, geo_id) {
  
  results <- list(overall_valid = TRUE, issues = character(0), warnings = character(0))
  
  pheno_data <- Biobase::pData(gse_matrix)
  
  if (nrow(pheno_data) == 0) {
    results$warnings <- c(results$warnings, "No phenotype data available")
    return(results)
  }
  
  # Check for disease/condition information
  pheno_cols <- colnames(pheno_data)
  disease_indicators <- grepl("disease|condition|group|treatment|tissue|diagnosis", 
                             pheno_cols, ignore.case = TRUE)
  
  if (!any(disease_indicators)) {
    results$warnings <- c(results$warnings, "No obvious disease/condition columns found")
  }
  
  # Check for demographic information
  demo_indicators <- grepl("age|sex|gender|ethnicity", pheno_cols, ignore.case = TRUE)
  if (any(demo_indicators)) {
    results$demographic_info <- TRUE
  } else {
    results$warnings <- c(results$warnings, "Limited demographic information")
  }
  
  results$phenotype_stats <- list(
    n_variables = ncol(pheno_data),
    disease_cols = sum(disease_indicators),
    demo_cols = sum(demo_indicators)
  )
  
  return(results)
}

#' Detect CAMK Genes in Dataset
#' @param gse_matrix GEO matrix object
#' @param geo_id Dataset identifier
#' @return CAMK detection results
detect_camk_genes_in_dataset <- function(gse_matrix, geo_id) {
  
  # Get feature data for gene mapping
  feature_data <- Biobase::fData(gse_matrix)
  expr_matrix <- Biobase::exprs(gse_matrix)
  
  # Define CAMK genes to search for
  target_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", 
                   "CAMK1", "CAMK4", "CAMKK1", "CAMKK2")
  
  # Search strategies
  found_genes <- list()
  
  # Strategy 1: Direct rowname matching
  rowname_matches <- intersect(rownames(expr_matrix), target_genes)
  if (length(rowname_matches) > 0) {
    found_genes$rowname_matches <- rowname_matches
  }
  
  # Strategy 2: Search in feature annotation
  if (nrow(feature_data) > 0) {
    for (col in colnames(feature_data)) {
      if (is.character(feature_data[[col]])) {
        for (gene in target_genes) {
          matches <- grep(gene, feature_data[[col]], ignore.case = TRUE)
          if (length(matches) > 0) {
            if (is.null(found_genes$annotation_matches)) {
              found_genes$annotation_matches <- list()
            }
            found_genes$annotation_matches[[gene]] <- rownames(feature_data)[matches]
          }
        }
      }
    }
  }
  
  # Strategy 3: Fuzzy matching for gene symbols
  fuzzy_matches <- list()
  all_symbols <- c(rownames(expr_matrix), unlist(feature_data))
  for (gene in target_genes) {
    fuzzy_hits <- grep(gene, all_symbols, ignore.case = TRUE, value = TRUE)
    if (length(fuzzy_hits) > 0) {
      fuzzy_matches[[gene]] <- unique(fuzzy_hits)
    }
  }
  if (length(fuzzy_matches) > 0) {
    found_genes$fuzzy_matches <- fuzzy_matches
  }
  
  # Summary
  total_found <- length(unique(c(
    found_genes$rowname_matches,
    unlist(found_genes$annotation_matches),
    unlist(found_genes$fuzzy_matches)
  )))
  
  list(
    dataset_id = geo_id,
    target_genes = target_genes,
    found_genes = found_genes,
    total_found = total_found,
    detection_rate = total_found / length(target_genes)
  )
}

#' Calculate Dataset Quality Score
#' @param validation_results Validation results object
#' @return Quality score (0-10)
calculate_dataset_quality_score <- function(validation_results) {
  
  if (!validation_results$overall_valid) {
    return(0)
  }
  
  score <- 10  # Start with perfect score
  
  # Penalize for issues and warnings
  score <- score - length(validation_results$issues) * 2
  score <- score - length(validation_results$warnings) * 0.5
  
  # Sample size bonus/penalty
  if (!is.null(validation_results$dimensions)) {
    n_samples <- validation_results$dimensions$samples
    if (n_samples >= 50) score <- score + 1
    else if (n_samples < 20) score <- score - 1
  }
  
  # CAMK detection bonus
  if (!is.null(validation_results$camk_detection)) {
    detection_rate <- validation_results$camk_detection$detection_rate
    score <- score + detection_rate * 2
  }
  
  return(max(0, min(10, score)))
}

#' Create Failed Dataset Result
#' @param geo_id Dataset identifier
#' @param error_message Error description
#' @return Failed result object
create_failed_dataset_result <- function(geo_id, error_message) {
  
  list(
    dataset_id = geo_id,
    success = FALSE,
    error = error_message,
    validation = list(
      overall_valid = FALSE,
      issues = error_message,
      quality_score = 0
    )
  )
}

#' Merge Validation Results
#' @param main_results Main validation results
#' @param new_results New validation component
#' @return Merged results
merge_validation_results <- function(main_results, new_results) {
  
  # Merge validity status
  main_results$overall_valid <- main_results$overall_valid && new_results$overall_valid
  
  # Merge issues and warnings
  main_results$issues <- c(main_results$issues, new_results$issues)
  main_results$warnings <- c(main_results$warnings, new_results$warnings)
  
  # Merge additional components
  for (name in names(new_results)) {
    if (!name %in% c("overall_valid", "issues", "warnings")) {
      main_results[[name]] <- new_results[[name]]
    }
  }
  
  return(main_results)
}

#' Extract Enhanced Metadata
#' @param gse_data GEO dataset object
#' @param geo_id Dataset identifier
#' @return Enhanced metadata object
extract_enhanced_metadata <- function(gse_data, geo_id) {
  
  if (is.list(gse_data)) {
    gse_matrix <- gse_data[[1]]
  } else {
    gse_matrix <- gse_data
  }
  
  metadata <- list(
    dataset_id = geo_id,
    title = Biobase::experimentData(gse_matrix)@title,
    summary = Biobase::experimentData(gse_matrix)@abstract,
    organism = Biobase::experimentData(gse_matrix)@other$organism,
    platform = Biobase::annotation(gse_matrix),
    sample_count = ncol(Biobase::exprs(gse_matrix)),
    gene_count = nrow(Biobase::exprs(gse_matrix)),
    extraction_date = Sys.Date()
  )
  
  return(metadata)
}