# Dataset Configuration
# Clean, verified GEO dataset definitions for CAMK2D analysis

#' Get Human Heart Failure Dataset Configuration (Enhanced Meta-Analysis)
#' @return Data frame with expanded human heart failure datasets for meta-analysis
get_human_heart_failure_datasets <- function() {
  data.frame(
    GSE_ID = c("GSE57338", "GSE116250", "GSE5406", "GSE1145", "GSE79962", "GSE21610", "GSE26887", "GSE42955"),
    Title = c(
      "Human Heart Failure RNA-seq (Primary Discovery)",
      "Heart Failure RNA-seq (Etiology-specific)", 
      "Human Heart Failure Expression Profiles (Large Cohort)",
      "Heart Failure Etiology Comparison Study",
      "Multi-etiology Heart Failure Analysis",
      "Ischemic vs Non-ischemic Cardiomyopathy",
      "Dilated Cardiomyopathy Expression Analysis",
      "Heart Failure Progression Study"
    ),
    Sample_Count = c(313, 64, 210, 37, 48, 30, 25, 42),
    Platform = c("RNA-seq", "RNA-seq", "Microarray", "Microarray", "Microarray", "Microarray", "Microarray", "Microarray"),
    Disease_Focus = c("Heart Failure", "Heart Failure", "Heart Failure", "Heart Failure", "Heart Failure", "Heart Failure", "Heart Failure", "Heart Failure"),
    Quality_Score = c(9, 8, 9, 7, 7, 6, 6, 7),
    Species = c("Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human"),
    Verified = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Core datasets verified, others need validation
    Meta_Analysis_Priority = c(1, 1, 1, 2, 2, 3, 3, 3),  # Priority for meta-analysis inclusion
    stringsAsFactors = FALSE
  )
}

#' Get Human Atrial Fibrillation Dataset Configuration (Enhanced Meta-Analysis)
#' @return Data frame with expanded human atrial fibrillation datasets for meta-analysis
get_human_atrial_fib_datasets <- function() {
  data.frame(
    GSE_ID = c("GSE203367", "GSE226283", "GSE226282", "GSE41177", "GSE79768", "GSE14975", 
               "GSE244414", "GSE208494", "GSE237003", "GSE163754", "GSE205741", "GSE179132"),
    Title = c(
      "Multi-regional Atrial Fibrillation RNA-seq",
      "Atrial Fibrillation Cardiac Tissue Analysis",
      "CATCH-ME Consortium AF RNA-seq Series",
      "Established Atrial Fibrillation Dataset",
      "AF Validation Microarray Study", 
      "Atrial Tissue-Specific Analysis",
      "Sex-specific Atrial Fibrillation Transcriptomics (2024)",
      "Human AF with Anti-arrhythmic Drug Evaluation (2024)",
      "Comprehensive 75-sample AF Bulk RNA-seq (2023)",
      "Human Left Atrial AF Samples RNA-seq (2021)",
      "Human Atrial AF Gene Expression Profiling (2022)",
      "Multi-omics Human Left Atrial Appendage (2021)"
    ),
    Sample_Count = c(21, 21, 35, 28, 24, 18, 40, 32, 75, 45, 50, 66),
    Platform = c("RNA-seq", "RNA-seq", "RNA-seq", "Microarray", "Microarray", "Microarray",
                 "RNA-seq", "RNA-seq", "RNA-seq", "RNA-seq", "RNA-seq", "Multi-omics"),
    Disease_Focus = c("Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation",
                      "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation", "Atrial Fibrillation"),
    Quality_Score = c(7, 7, 8, 8, 7, 6, 9, 8, 9, 8, 8, 9),
    Species = c("Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Human"),
    Verified = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),  # Recent high-quality studies verified
    Meta_Analysis_Priority = c(2, 2, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1),  # Recent studies are high priority
    stringsAsFactors = FALSE
  )
}

#' Get Mouse Cardiac Dataset Configuration (Enhanced Cross-Species Meta-Analysis)
#' @return Data frame with expanded mouse cardiac datasets for cross-species validation
get_mouse_cardiac_datasets <- function() {
  data.frame(
    GSE_ID = c("GSE48760", "GSE1869", "GSE36074", "GSE18224", "GSE23294", "GSE71733", 
               "GSE201018", "GSE266652", "GSE151156", "GSE160043"),
    Title = c(
      "Isoproterenol-induced Cardiac Hypertrophy Model",
      "TAC Mouse Model (Pressure Overload)", 
      "Mouse Pressure Overload-induced Heart Failure",
      "CAMK2D Knockout Mouse Cardiac Study",
      "CAMK2D Transgenic Overexpression Model",
      "Aging-induced AF Mouse Model",
      "RBM20 I536T variant affects CAMK2D splicing in cardiac tissue",
      "Heart-specific CAMK2D gene silencing via AAV-base editing",
      "Mouse Cardiac Hypertrophy with Calcium Signaling Focus (2020)",
      "Single-cell Transcriptome of Atrial Appendage AF vs non-AF (2020)"
    ),
    Sample_Count = c(24, 18, 36, 16, 20, 22, 12, 30, 35, 42),
    Platform = c("Microarray", "Microarray", "Microarray", "Microarray", "Microarray", "RNA-seq", 
                 "RNA-seq", "Amplicon-seq", "RNA-seq", "Single-cell RNA-seq"),
    Disease_Focus = c("Cardiac Hypertrophy", "Heart Failure", "Heart Failure", "CAMK2D Knockout", "CAMK2D Overexpression", 
                      "Atrial Fibrillation", "CAMK2D Splicing", "CAMK2D Gene Editing", 
                      "Cardiac Hypertrophy", "Atrial Fibrillation"),
    Quality_Score = c(7, 8, 8, 9, 8, 7, 9, 10, 8, 8),
    Species = c("Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse"),
    Verified = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),  # Updated based on validation results
    Meta_Analysis_Priority = c(2, 1, 1, 1, 1, 3, 1, 1, 1, 1),  # CAMK2D-specific studies high priority
    Model_Type = c("Pharmacological", "Surgical", "Surgical", "Genetic", "Genetic", "Aging", 
                   "Genetic", "Gene Editing", "Pharmacological", "Genetic"),
    stringsAsFactors = FALSE
  )
}

#' Get CAMK2D-Specific Dataset Configuration (High Priority)
#' @return Data frame with CAMK2D-focused datasets from recent studies
get_camk2d_specific_datasets <- function() {
  data.frame(
    GSE_ID = c("GSE225336", "GSE228762", "GSE161032", "GSE154030", "GSE244117"),
    Title = c(
      "Myocardial Fibrosis Genetics - UK Biobank 41,505 participants",
      "CAMK2D CRISPR Knockout in mpkCCD cells",
      "CAMK2D in Acute Myeloid Leukemia",
      "CAMK2D Splice Variants in Human Hippocampus",
      "Single-nucleus RNA-seq of Human Atrial Tissue Biopsies (2023)"
    ),
    Sample_Count = c(100, 12, 48, 24, 38),  # Estimated for GSE225336 tissue samples
    Platform = c("RNA-seq", "RNA-seq", "RNA-seq", "RNA-seq", "Single-nucleus RNA-seq"),
    Disease_Focus = c("Myocardial Fibrosis", "Kidney/Vasopressin", "Leukemia", "Neural", "Atrial Fibrillation"),
    Quality_Score = c(10, 8, 6, 7, 9),
    Species = c("Human", "Mouse", "Human", "Human", "Human"),
    Verified = c(TRUE, FALSE, FALSE, FALSE, TRUE),  # UK Biobank and single-cell studies verified
    Meta_Analysis_Priority = c(1, 2, 3, 3, 1),  # Fibrosis study is highest priority
    Study_Type = c("Population Genetics", "CRISPR Knockout", "Cancer", "Splice Variants", "Single-cell"),
    CAMK2D_Focus = c("Genetic Loci", "Direct Knockout", "Expression", "Splice Variants", "Cellular Resolution"),
    stringsAsFactors = FALSE
  )
}

#' Get All Dataset Configurations Combined
#' @return List with all dataset configurations
get_all_dataset_configs <- function() {
  list(
    human_heart_failure = get_human_heart_failure_datasets(),
    human_atrial_fib = get_human_atrial_fib_datasets(),
    mouse_cardiac = get_mouse_cardiac_datasets(),
    camk2d_specific = get_camk2d_specific_datasets()
  )
}

#' Get High-Priority Datasets Only  
#' @return Data frame with only the most reliable datasets
get_priority_datasets <- function() {
  # Get all configs
  all_configs <- get_all_dataset_configs()
  
  # Get common columns across all configs
  all_column_names <- lapply(all_configs, names)
  common_columns <- Reduce(intersect, all_column_names)
  
  # Filter for verified datasets only
  priority_datasets <- data.frame()
  
  for (dataset_type in names(all_configs)) {
    config <- all_configs[[dataset_type]]
    verified <- config[config$Verified == TRUE, ]
    
    if (nrow(verified) > 0) {
      # Use only common columns
      verified_subset <- verified[, common_columns]
      priority_datasets <- rbind(priority_datasets, verified_subset)
    }
  }
  
  # Sort by quality score (descending) if available
  if ("Quality_Score" %in% names(priority_datasets)) {
    priority_datasets <- priority_datasets[order(-priority_datasets$Quality_Score), ]
  }
  
  return(priority_datasets)
}

#' Get Enhanced Meta-Analysis Configuration
#' @return List with enhanced analysis parameters for multi-dataset meta-analysis
get_analysis_config <- function() {
  list(
    # Enhanced dataset loading settings for meta-analysis
    max_human_hf_datasets = 8,      # All available HF datasets
    max_human_af_datasets = 6,      # All available AF datasets  
    max_mouse_datasets = 8,         # All available mouse datasets
    
    # Meta-analysis settings
    meta_analysis_enabled = TRUE,
    cross_species_validation = TRUE,
    min_datasets_for_meta = 3,      # Minimum datasets required for meta-analysis
    meta_analysis_method = "random", # "random" or "fixed" effects
    
    # CAMK gene targets (expanded)
    target_genes = c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4"),
    primary_target = "CAMK2D",      # Primary focus for cross-species validation
    
    # Enhanced expression validation thresholds
    min_expression_threshold = 0.5,
    min_detection_percentage = 10,
    cross_species_correlation_threshold = 0.6,  # Minimum correlation between species
    
    # Advanced integration settings
    min_common_genes = 100,
    preferred_integration_method = "gene_symbols",
    batch_correction_method = "ComBat",
    cross_species_mapping_method = "biomaRt",
    
    # Meta-analysis statistical settings
    heterogeneity_threshold = 50,   # I² threshold for acceptable heterogeneity
    publication_bias_alpha = 0.05,  # Significance level for bias testing
    effect_size_method = "SMD",     # Standardized Mean Difference
    confidence_level = 0.95,        # For confidence intervals
    
    # Quality control settings
    min_sample_size = 10,           # Minimum samples per dataset
    max_missing_data = 0.3,         # Maximum proportion missing data
    outlier_detection_method = "IQR", # Outlier detection method
    
    # Enhanced output settings
    cache_downloads = TRUE,
    generate_detailed_reports = TRUE,
    create_visualizations = TRUE,
    export_meta_analysis = TRUE,
    create_forest_plots = TRUE,
    cross_species_networks = TRUE,
    
    # Therapeutic target prioritization
    therapeutic_scoring = TRUE,
    druggability_assessment = TRUE,
    conservation_weighting = 0.4,   # Weight for evolutionary conservation
    effect_size_weighting = 0.3,    # Weight for effect size magnitude  
    replication_weighting = 0.3     # Weight for cross-dataset replication
  )
}

#' Add New Dataset to Configuration
#' @param gse_id GEO Series ID
#' @param title Dataset title
#' @param sample_count Number of samples
#' @param platform Technology platform
#' @param disease_focus Primary disease/condition
#' @param species Organism
#' @param dataset_type Type of dataset (human_heart_failure, human_atrial_fib, mouse_cardiac, camk2d_specific)
#' @return Message indicating success or failure
add_dataset <- function(gse_id, title, sample_count, platform, disease_focus, species, dataset_type) {
  valid_types <- c("human_heart_failure", "human_atrial_fib", "mouse_cardiac", "camk2d_specific")
  
  if (!dataset_type %in% valid_types) {
    stop("Invalid dataset_type. Must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  cat("Adding dataset", gse_id, "to", dataset_type, "configuration\n")
  cat("Title:", title, "\n")
  cat("Samples:", sample_count, "\n")
  cat("Platform:", platform, "\n")
  cat("Disease:", disease_focus, "\n")
  cat("Species:", species, "\n")
  
  return(paste("Dataset", gse_id, "prepared for addition to", dataset_type))
}

#' Search for Datasets by Criteria
#' @param keyword Search term for title or disease focus
#' @param min_samples Minimum number of samples
#' @param species Filter by species
#' @param verified_only Only return verified datasets
#' @return Data frame of matching datasets
search_datasets <- function(keyword = NULL, min_samples = 0, species = NULL, verified_only = FALSE) {
  all_configs <- get_all_dataset_configs()
  
  # Get common columns across all configs
  all_column_names <- lapply(all_configs, names)
  common_columns <- Reduce(intersect, all_column_names)
  
  # Combine all datasets using only common columns
  all_datasets <- data.frame()
  for (config_name in names(all_configs)) {
    config <- all_configs[[config_name]]
    config$Dataset_Type <- config_name
    
    # Select only common columns plus Dataset_Type
    config_subset <- config[, c(common_columns, "Dataset_Type")]
    all_datasets <- rbind(all_datasets, config_subset)
  }
  
  # Apply filters
  if (!is.null(keyword)) {
    keyword_lower <- tolower(keyword)
    matching <- grepl(keyword_lower, tolower(all_datasets$Title)) | 
                grepl(keyword_lower, tolower(all_datasets$Disease_Focus))
    all_datasets <- all_datasets[matching, ]
  }
  
  if (min_samples > 0) {
    all_datasets <- all_datasets[all_datasets$Sample_Count >= min_samples, ]
  }
  
  if (!is.null(species)) {
    all_datasets <- all_datasets[all_datasets$Species == species, ]
  }
  
  if (verified_only) {
    all_datasets <- all_datasets[all_datasets$Verified == TRUE, ]
  }
  
  return(all_datasets)
}

#' Get Dataset Statistics
#' @return Summary statistics of all datasets
get_dataset_stats <- function() {
  all_configs <- get_all_dataset_configs()
  
  stats <- list()
  total_datasets <- 0
  total_samples <- 0
  total_verified <- 0
  
  for (config_name in names(all_configs)) {
    config <- all_configs[[config_name]]
    n_datasets <- nrow(config)
    n_samples <- sum(config$Sample_Count)
    n_verified <- sum(config$Verified)
    
    stats[[config_name]] <- list(
      datasets = n_datasets,
      samples = n_samples,
      verified = n_verified,
      verification_rate = round(n_verified / n_datasets * 100, 1)
    )
    
    total_datasets <- total_datasets + n_datasets
    total_samples <- total_samples + n_samples
    total_verified <- total_verified + n_verified
  }
  
  stats$total <- list(
    datasets = total_datasets,
    samples = total_samples,
    verified = total_verified,
    verification_rate = round(total_verified / total_datasets * 100, 1)
  )
  
  return(stats)
}