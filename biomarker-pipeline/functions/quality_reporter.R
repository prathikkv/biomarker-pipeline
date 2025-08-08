# Quality Reporting Functions
# Clean functions for generating data quality and CAMK gene presence reports

#' Create Data Quality Report
#' @param datasets_list List of successfully loaded datasets
#' @param failed_datasets Vector of failed dataset IDs
#' @param dataset_config Original dataset configuration
#' @return Data frame with quality metrics
create_quality_report <- function(datasets_list, failed_datasets = NULL, dataset_config = NULL) {
  quality_report <- data.frame(
    Dataset = character(),
    Status = character(),
    Samples = numeric(),
    Genes = numeric(),
    CAMK_Found = numeric(),
    Platform = character(),
    stringsAsFactors = FALSE
  )
  
  # Add successful datasets
  for (gse_id in names(datasets_list)) {
    dataset <- datasets_list[[gse_id]]
    
    camk_count <- 0
    if (!is.null(dataset$camk_detection_results)) {
      # Safely count CAMK genes found with error handling
      tryCatch({
        camk_count <- sum(sapply(dataset$camk_detection_results, function(x) {
          if (!is.null(x) && !is.null(x$found)) x$found else FALSE
        }))
      }, error = function(e) {
        camk_count <- 0
      })
    }
    
    platform <- if (!is.null(dataset$platform)) dataset$platform else "Unknown"
    
    quality_report <- rbind(quality_report, data.frame(
      Dataset = gse_id,
      Status = "Loaded",
      Samples = ncol(dataset$expression),
      Genes = nrow(dataset$expression),
      CAMK_Found = camk_count,
      Platform = platform,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add failed datasets
  if (!is.null(failed_datasets) && !is.null(dataset_config)) {
    for (gse_id in failed_datasets) {
      platform <- "Unknown"
      # Safely access Platform column
      if (gse_id %in% dataset_config$GSE_ID && "Platform" %in% names(dataset_config)) {
        platform_matches <- dataset_config$Platform[dataset_config$GSE_ID == gse_id]
        if (length(platform_matches) > 0 && !is.na(platform_matches[1])) {
          platform <- platform_matches[1]
        }
      }
      
      quality_report <- rbind(quality_report, data.frame(
        Dataset = gse_id,
        Status = "Failed",
        Samples = 0,
        Genes = 0,
        CAMK_Found = 0,
        Platform = platform,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(quality_report)
}

#' Display Quality Report Summary
#' @param quality_report Data quality report
display_quality_summary <- function(quality_report) {
  cat("=== DATA QUALITY SUMMARY ===\n")
  
  if (is.null(quality_report) || nrow(quality_report) == 0) {
    cat("No quality data available\n")
    return()
  }
  
  print(quality_report)
  
  cat("\nSummary Statistics:\n")
  cat("- Total datasets attempted:", nrow(quality_report), "\n")
  cat("- Successfully loaded:", sum(quality_report$Status == "Loaded"), "\n")
  cat("- Failed to load:", sum(quality_report$Status == "Failed"), "\n")
  cat("- Total samples available:", sum(quality_report$Samples), "\n")
  cat("- Datasets with CAMK genes:", sum(quality_report$CAMK_Found > 0), "\n")
}

#' Create CAMK Gene Presence Matrix
#' @param datasets_list List of successfully loaded datasets
#' @return Data frame showing CAMK gene presence across datasets
create_camk_presence_matrix <- function(datasets_list) {
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
  
  presence_matrix <- data.frame(
    Dataset = character(),
    stringsAsFactors = FALSE
  )
  
  # Add columns for each CAMK gene
  for (gene in camk_genes) {
    presence_matrix[[gene]] <- character()
  }
  
  presence_matrix$Total_Found <- numeric()
  presence_matrix$Platform <- character()
  
  # Process each dataset
  for (gse_id in names(datasets_list)) {
    dataset <- datasets_list[[gse_id]]
    
    if (is.null(dataset$camk_detection_results)) {
      next
    }
    
    results <- dataset$camk_detection_results
    
    # Create row for this dataset
    row_data <- list(Dataset = gse_id)
    
    total_found <- 0
    for (gene in camk_genes) {
      if (gene %in% names(results) && 
          !is.null(results[[gene]]) && 
          !is.null(results[[gene]]$found) && 
          results[[gene]]$found) {
        count_val <- if (!is.null(results[[gene]]$count)) results[[gene]]$count else "?"
        row_data[[gene]] <- paste0("✓ (", count_val, ")")
        total_found <- total_found + 1
      } else {
        row_data[[gene]] <- "✗"
      }
    }
    
    row_data$Total_Found <- total_found
    row_data$Platform <- if (!is.null(dataset$platform)) dataset$platform else "Unknown"
    
    # Add row to matrix
    presence_matrix <- rbind(presence_matrix, data.frame(row_data, stringsAsFactors = FALSE))
  }
  
  return(presence_matrix)
}

#' Display CAMK Presence Report
#' @param presence_matrix CAMK gene presence matrix
display_camk_presence_report <- function(presence_matrix) {
  cat("=== CAMK GENE PRESENCE REPORT ===\n")
  
  if (is.null(presence_matrix) || nrow(presence_matrix) == 0) {
    cat("No CAMK detection results available\n")
    return()
  }
  
  cat("Detection Results (✓ = Found, ✗ = Not Found, numbers = probe count):\n\n")
  print(presence_matrix)
  
  # Calculate summary statistics
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
  total_datasets <- nrow(presence_matrix)
  
  cat("\n=== CAMK GENE DETECTION SUMMARY ===\n")
  for (gene in camk_genes) {
    if (gene %in% names(presence_matrix)) {
      found_count <- sum(grepl("✓", presence_matrix[[gene]]))
      percentage <- round(found_count / total_datasets * 100, 1)
      cat(sprintf("%-7s: %d/%d datasets (%s%%)\n", gene, found_count, total_datasets, percentage))
    }
  }
  
  # Highlight cardiac-specific genes
  if ("CAMK2D" %in% names(presence_matrix) && "CAMK2G" %in% names(presence_matrix)) {
    camk2d_found <- sum(grepl("✓", presence_matrix$CAMK2D))
    camk2g_found <- sum(grepl("✓", presence_matrix$CAMK2G))
    
    cat("\n** Primary cardiac isoforms **\n")
    cat("CAMK2D (delta): ", camk2d_found, "/", total_datasets, " datasets\n")
    cat("CAMK2G (gamma): ", camk2g_found, "/", total_datasets, " datasets\n")
  }
}

#' Generate Complete Quality Report
#' @param datasets_list List of successfully loaded datasets
#' @param failed_datasets Vector of failed dataset IDs (optional)
#' @param dataset_config Original dataset configuration (optional)
generate_complete_report <- function(datasets_list, failed_datasets = NULL, dataset_config = NULL) {
  # Data quality report
  quality_report <- create_quality_report(datasets_list, failed_datasets, dataset_config)
  display_quality_summary(quality_report)
  
  cat("\n")
  
  # CAMK presence report
  presence_matrix <- create_camk_presence_matrix(datasets_list)
  display_camk_presence_report(presence_matrix)
  
  return(list(
    quality_report = quality_report,
    presence_matrix = presence_matrix
  ))
}