# GEO Data Loading Functions  
# Clean, modular functions for loading and preprocessing GEO datasets

#' Load Single GEO Dataset
#' @param gse_id GEO series ID (e.g. "GSE57338")
#' @param detect_camk Whether to detect CAMK genes (default TRUE)
#' @return List with dataset components
load_single_geo_dataset <- function(gse_id, detect_camk = TRUE) {
  cat("Loading dataset:", gse_id, "\n")
  
  tryCatch({
    # Download from GEO
    gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
    
    if (length(gse_list) == 0) {
      stop("No data retrieved for ", gse_id)
    }
    
    gse <- gse_list[[1]]
    
    # Extract core data
    expression_matrix <- exprs(gse)
    phenotype_data <- pData(gse)
    feature_data <- fData(gse)
    
    cat("Data extracted: ", nrow(expression_matrix), " features x ", 
        ncol(expression_matrix), " samples\n")
    
    # Extract gene symbols
    gene_symbols <- extract_gene_symbols(feature_data, expression_matrix)
    
    # Process phenotype data
    phenotype_data <- process_phenotype_data(phenotype_data)
    
    # Detect CAMK genes if requested
    camk_results <- NULL
    if (detect_camk && exists("detect_camk_genes")) {
      camk_results <- detect_camk_genes(gene_symbols, feature_data, expression_matrix)
    }
    
    # Create clean dataset object
    dataset <- list(
      gse_id = gse_id,
      expression = expression_matrix,
      phenotypes = phenotype_data,
      feature_data = feature_data,
      gene_symbols = gene_symbols,
      camk_detection_results = camk_results
    )
    
    cat("✓ Successfully loaded:", gse_id, "\n")
    return(dataset)
    
  }, error = function(e) {
    cat("✗ Failed to load", gse_id, ":", e$message, "\n")
    return(NULL)
  })
}

#' Process Phenotype Data
#' @param phenotype_data Raw phenotype data from GEO
#' @return Processed phenotype data with standardized disease status
process_phenotype_data <- function(phenotype_data) {
  # Look for disease status columns
  status_cols <- grep("status|disease|condition|diagnosis|characteristics", 
                     names(phenotype_data), ignore.case = TRUE, value = TRUE)
  
  if (length(status_cols) > 0) {
    disease_status <- as.character(phenotype_data[[status_cols[1]]])
    
    # Standardize disease status
    control_patterns <- c("control", "normal", "healthy", "donor")
    disease_patterns <- c("heart.*failure", "cardiomyopathy", "diseased", "failing", "patient")
    
    is_control <- grepl(paste(control_patterns, collapse = "|"), disease_status, ignore.case = TRUE)
    is_disease <- grepl(paste(disease_patterns, collapse = "|"), disease_status, ignore.case = TRUE)
    
    phenotype_data$Disease_Status <- ifelse(is_control, "Control", 
                                          ifelse(is_disease, "Disease", "Unknown"))
  } else {
    phenotype_data$Disease_Status <- "Unknown"
  }
  
  return(phenotype_data)
}

#' Load Multiple GEO Datasets
#' @param dataset_config Data frame with GSE_ID and metadata
#' @param max_datasets Maximum number of datasets to load
#' @param species Species filter (optional)
#' @return List of successfully loaded datasets
load_multiple_geo_datasets <- function(dataset_config, max_datasets = 3, species = NULL) {
  cat("Loading multiple datasets...\n")
  
  # Filter by species if specified
  if (!is.null(species)) {
    dataset_config <- dataset_config[dataset_config$Species == species, ]
  }
  
  successfully_loaded <- list()
  failed_datasets <- character()
  
  n_to_load <- min(nrow(dataset_config), max_datasets)
  
  for (i in 1:n_to_load) {
    gse_id <- dataset_config$GSE_ID[i]
    
    dataset <- load_single_geo_dataset(gse_id)
    
    if (!is.null(dataset)) {
      # Add metadata
      dataset$platform <- dataset_config$Platform[i]
      dataset$disease_focus <- dataset_config$Disease_Focus[i]
      dataset$sample_count <- ncol(dataset$expression)
      
      successfully_loaded[[gse_id]] <- dataset
      cat("✓ Added:", gse_id, "with", ncol(dataset$expression), "samples\n")
    } else {
      failed_datasets <- c(failed_datasets, gse_id)
    }
    
    # Polite delay
    if (i < n_to_load) Sys.sleep(2)
  }
  
  cat("Successfully loaded", length(successfully_loaded), "datasets\n")
  if (length(failed_datasets) > 0) {
    cat("Failed:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  return(successfully_loaded)
}

#' Quick Dataset Summary
#' @param dataset Dataset object
#' @return Summary statistics
summarize_dataset <- function(dataset) {
  if (is.null(dataset)) return(NULL)
  
  n_samples <- ncol(dataset$expression)
  n_genes <- nrow(dataset$expression)
  n_controls <- sum(dataset$phenotypes$Disease_Status == "Control", na.rm = TRUE)
  n_cases <- sum(dataset$phenotypes$Disease_Status == "Disease", na.rm = TRUE)
  
  camk_found <- 0
  if (!is.null(dataset$camk_detection_results)) {
    camk_found <- safe_count_camk_found(dataset$camk_detection_results)
  }
  
  list(
    gse_id = dataset$gse_id,
    samples = n_samples,
    genes = n_genes,
    controls = n_controls,
    cases = n_cases,
    camk_genes = camk_found
  )
}