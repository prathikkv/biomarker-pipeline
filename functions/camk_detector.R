# CAMK Gene Detection Functions
# Clean, modular functions for detecting CAMK genes in GEO datasets

#' CAMK Gene Aliases Configuration
#' @return Named list of CAMK gene aliases
get_camk_aliases <- function() {
  list(
    "CAMK2D" = c("CAMK2D", "CAMKIID", "CaMKIIdelta"),
    "CAMK2A" = c("CAMK2A", "CAMKIIA", "CaMKIIalpha"), 
    "CAMK2B" = c("CAMK2B", "CAMKIIB", "CaMKIIbeta"),
    "CAMK2G" = c("CAMK2G", "CAMKIIG", "CaMKIIgamma"),
    "CAMK1" = c("CAMK1", "CAMKI", "CaMKI"),
    "CAMK4" = c("CAMK4", "CAMKIV", "CaMKIV")
  )
}

#' Extract Gene Symbols from Feature Data
#' @param feature_data Feature data from GEO dataset
#' @param expression_matrix Expression matrix for fallback row names
#' @return Vector of gene symbols
extract_gene_symbols <- function(feature_data, expression_matrix) {
  # Try standard gene symbol columns
  if ("Gene Symbol" %in% names(feature_data)) {
    return(feature_data$`Gene Symbol`)
  }
  
  if ("GENE_SYMBOL" %in% names(feature_data)) {
    return(feature_data$GENE_SYMBOL)
  }
  
  # Try gene assignment column
  if ("gene_assignment" %in% names(feature_data)) {
    gene_symbols <- gsub(".*gene_symbol:", "", feature_data$gene_assignment)
    gene_symbols <- gsub("/.*", "", gene_symbols)
    return(gene_symbols)
  }
  
  # Try any column with "gene" or "symbol" in name
  gene_cols <- grep("gene|symbol", names(feature_data), ignore.case = TRUE, value = TRUE)
  if (length(gene_cols) > 0) {
    return(feature_data[[gene_cols[1]]])
  }
  
  # Fallback to row names
  return(rownames(expression_matrix))
}

#' Detect Single CAMK Gene
#' @param target_gene Name of target CAMK gene
#' @param aliases Vector of aliases for the gene  
#' @param gene_symbols Vector of gene symbols from dataset
#' @param feature_data Feature data for extended search
#' @return List with detection results
detect_single_camk_gene <- function(target_gene, aliases, gene_symbols, feature_data = NULL) {
  found_indices <- integer()
  found_descriptions <- character()
  
  # Search in gene symbols
  for (alias in aliases) {
    if (is.character(alias) && nchar(alias) > 0) {
      # Exact matching for main gene names
      matches <- which(toupper(gene_symbols) %in% toupper(alias))
      
      if (length(matches) > 0) {
        found_indices <- unique(c(found_indices, matches))
        found_descriptions <- c(found_descriptions, gene_symbols[matches[1]])
      }
    }
  }
  
  # Extended search in feature data if no matches found
  if (length(found_indices) == 0 && !is.null(feature_data)) {
    for (col in names(feature_data)) {
      if (is.character(feature_data[[col]])) {
        for (alias in aliases) {
          matches <- grep(alias, feature_data[[col]], ignore.case = TRUE, fixed = TRUE)
          if (length(matches) > 0) {
            found_indices <- unique(c(found_indices, matches))
            found_descriptions <- c(found_descriptions, feature_data[[col]][matches[1]])
            break
          }
        }
        if (length(found_indices) > 0) break
      }
    }
  }
  
  return(list(
    found = length(found_indices) > 0,
    indices = found_indices,
    count = length(found_indices),
    descriptions = unique(found_descriptions)
  ))
}

#' Validate Gene Expression Levels
#' @param gene_indices Indices of detected genes
#' @param expression_matrix Expression matrix
#' @return List with validation results
validate_gene_expression <- function(gene_indices, expression_matrix) {
  if (length(gene_indices) == 0) {
    return(list(validated = FALSE, mean = 0, non_zero_pct = 0))
  }
  
  gene_expressions <- expression_matrix[gene_indices, , drop = FALSE]
  
  if (nrow(gene_expressions) == 1) {
    expr_values <- as.numeric(gene_expressions[1, ])
  } else {
    expr_values <- colMeans(gene_expressions, na.rm = TRUE)
  }
  
  mean_expr <- mean(expr_values, na.rm = TRUE)
  non_zero_pct <- sum(expr_values > 0.1, na.rm = TRUE) / length(expr_values) * 100
  
  return(list(
    validated = mean_expr > 0.5 && non_zero_pct > 10,
    mean = mean_expr,
    non_zero_pct = non_zero_pct
  ))
}

#' Detect All CAMK Genes in Dataset
#' @param gene_symbols Vector of gene symbols
#' @param feature_data Feature data from GEO dataset
#' @param expression_matrix Expression matrix for validation
#' @return List with detection results for all CAMK genes
detect_camk_genes <- function(gene_symbols, feature_data = NULL, expression_matrix = NULL) {
  camk_aliases <- get_camk_aliases()
  camk_detection_results <- list()
  
  cat("=== CAMK GENE DETECTION ===\n")
  
  for (target_gene in names(camk_aliases)) {
    aliases <- camk_aliases[[target_gene]]
    
    # Detect gene
    detection_result <- detect_single_camk_gene(target_gene, aliases, gene_symbols, feature_data)
    
    # Validate expression if matrix provided
    if (!is.null(expression_matrix) && detection_result$found) {
      validation <- validate_gene_expression(detection_result$indices, expression_matrix)
      detection_result <- c(detection_result, validation)
    }
    
    camk_detection_results[[target_gene]] <- detection_result
    
    # Report results
    if (detection_result$found) {
      status <- if (!is.null(detection_result$validated)) {
        ifelse(detection_result$validated, "✓ FOUND & VALIDATED", "⚠ FOUND (low signal)")
      } else {
        "✓ FOUND"
      }
      cat(target_gene, ":", status, "(", detection_result$count, "probes)\n")
    } else {
      cat("✗", target_gene, ": NOT FOUND\n")
    }
  }
  
  # Summary
  genes_found <- safe_count_camk_found(camk_detection_results)
  cat("Summary:", genes_found, "/", length(camk_aliases), "CAMK genes detected\n")
  
  return(camk_detection_results)
}

#' Safe Count of CAMK Genes Found
#' @param detection_results CAMK detection results object
#' @return Integer count of genes found (safely handles any data type issues)
safe_count_camk_found <- function(detection_results) {
  if (is.null(detection_results) || length(detection_results) == 0) {
    return(0)
  }
  
  total_found <- 0
  tryCatch({
    for (gene_result in detection_results) {
      if (!is.null(gene_result) && 
          !is.null(gene_result$found) && 
          is.logical(gene_result$found) && 
          length(gene_result$found) == 1 &&
          gene_result$found == TRUE) {
        total_found <- total_found + 1
      }
    }
  }, error = function(e) {
    # If any error occurs, fall back to safer manual counting
    for (gene_result in detection_results) {
      if (is.list(gene_result) && "found" %in% names(gene_result)) {
        if (isTRUE(gene_result$found)) {
          total_found <- total_found + 1
        }
      }
    }
  })
  
  return(total_found)
}