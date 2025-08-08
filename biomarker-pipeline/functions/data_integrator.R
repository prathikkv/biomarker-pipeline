# Data Integration Functions
# Functions for integrating multiple GEO datasets with batch correction

#' Find Common Genes Across Datasets
#' @param datasets_list List of dataset objects
#' @return List with common genes and mapping info
find_common_genes <- function(datasets_list) {
  if (length(datasets_list) == 0) return(NULL)
  
  # Try gene symbol mapping first
  gene_mappings <- list()
  
  for (gse_id in names(datasets_list)) {
    dataset <- datasets_list[[gse_id]]
    if (!is.null(dataset$gene_symbols)) {
      mapping <- data.frame(
        probe_id = rownames(dataset$expression),
        gene_symbol = dataset$gene_symbols,
        stringsAsFactors = FALSE
      )
      # Remove empty symbols
      mapping <- mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol != "", ]
      gene_mappings[[gse_id]] <- mapping
    }
  }
  
  # Find common gene symbols
  if (length(gene_mappings) > 0) {
    common_symbols <- Reduce(intersect, lapply(gene_mappings, function(x) unique(x$gene_symbol)))
    cat("Found", length(common_symbols), "common gene symbols\n")
    
    if (length(common_symbols) > 100) {
      return(list(
        method = "gene_symbols",
        common_genes = common_symbols,
        mappings = gene_mappings
      ))
    }
  }
  
  # Fallback to probe ID matching
  expression_list <- lapply(datasets_list, function(x) x$expression)
  common_probes <- Reduce(intersect, lapply(expression_list, rownames))
  
  cat("Found", length(common_probes), "common probe IDs\n")
  
  return(list(
    method = "probe_ids",
    common_genes = common_probes,
    mappings = NULL
  ))
}

#' Create Integrated Expression Matrix
#' @param datasets_list List of dataset objects
#' @param common_info Common gene information from find_common_genes
#' @return Integrated expression matrix
create_integrated_expression <- function(datasets_list, common_info) {
  if (common_info$method == "gene_symbols") {
    return(create_symbol_based_matrix(datasets_list, common_info))
  } else {
    return(create_probe_based_matrix(datasets_list, common_info))
  }
}

#' Create Expression Matrix Using Gene Symbols
#' @param datasets_list List of dataset objects
#' @param common_info Common gene information
#' @return Integrated expression matrix with gene symbols as row names
create_symbol_based_matrix <- function(datasets_list, common_info) {
  common_symbols <- common_info$common_genes
  gene_mappings <- common_info$mappings
  
  integrated_matrices <- list()
  
  for (gse_id in names(gene_mappings)) {
    mapping <- gene_mappings[[gse_id]]
    expr_data <- datasets_list[[gse_id]]$expression
    
    # Create gene symbol matrix
    expr_symbols <- matrix(NA, nrow = length(common_symbols), ncol = ncol(expr_data))
    rownames(expr_symbols) <- common_symbols
    colnames(expr_symbols) <- colnames(expr_data)
    
    for (i in 1:length(common_symbols)) {
      symbol <- common_symbols[i]
      probes <- mapping$probe_id[mapping$gene_symbol == symbol]
      
      if (length(probes) == 1) {
        expr_symbols[i, ] <- expr_data[probes, ]
      } else if (length(probes) > 1) {
        # Average multiple probes
        expr_symbols[i, ] <- colMeans(expr_data[probes, , drop = FALSE], na.rm = TRUE)
      }
    }
    
    integrated_matrices[[gse_id]] <- expr_symbols
  }
  
  # Combine all datasets
  combined_matrix <- do.call(cbind, integrated_matrices)
  return(combined_matrix)
}

#' Create Expression Matrix Using Probe IDs  
#' @param datasets_list List of dataset objects
#' @param common_info Common gene information
#' @return Integrated expression matrix with probe IDs as row names
create_probe_based_matrix <- function(datasets_list, common_info) {
  common_probes <- common_info$common_genes
  
  expression_list <- lapply(datasets_list, function(x) x$expression[common_probes, ])
  combined_matrix <- do.call(cbind, expression_list)
  
  return(combined_matrix)
}

#' Integrate Multiple Datasets
#' @param datasets_list List of successfully loaded datasets
#' @return Integrated dataset object or single best dataset
integrate_datasets <- function(datasets_list) {
  if (length(datasets_list) == 0) {
    cat("No datasets to integrate\n")
    return(NULL)
  }
  
  if (length(datasets_list) == 1) {
    cat("Only one dataset available, using without integration\n")
    return(datasets_list[[1]])
  }
  
  cat("Integrating", length(datasets_list), "datasets...\n")
  
  # Find common genes
  common_info <- find_common_genes(datasets_list)
  
  if (length(common_info$common_genes) == 0) {
    cat("No common genes found, using largest dataset\n")
    dataset_sizes <- sapply(datasets_list, function(x) ncol(x$expression))
    largest_idx <- which.max(dataset_sizes)
    return(datasets_list[[largest_idx]])
  }
  
  # Create integrated matrix
  tryCatch({
    integrated_expression <- create_integrated_expression(datasets_list, common_info)
    
    # Create integrated phenotype data
    integrated_phenotypes <- create_integrated_phenotypes(datasets_list)
    
    # Create batch information
    batch_info <- create_batch_info(datasets_list)
    
    cat("✓ Successfully integrated", length(datasets_list), "datasets\n")
    cat("  Final dimensions:", nrow(integrated_expression), "x", ncol(integrated_expression), "\n")
    
    return(list(
      expression = integrated_expression,
      phenotypes = integrated_phenotypes,
      batch_info = batch_info,
      common_genes = common_info$common_genes,
      integration_method = common_info$method,
      original_datasets = datasets_list
    ))
    
  }, error = function(e) {
    cat("Integration failed:", e$message, "\n")
    cat("Using largest dataset instead\n")
    
    dataset_sizes <- sapply(datasets_list, function(x) ncol(x$expression))
    largest_idx <- which.max(dataset_sizes)
    return(datasets_list[[largest_idx]])
  })
}

#' Create Integrated Phenotype Data
#' @param datasets_list List of dataset objects
#' @return Combined phenotype data frame
create_integrated_phenotypes <- function(datasets_list) {
  phenotype_list <- list()
  
  for (gse_id in names(datasets_list)) {
    dataset <- datasets_list[[gse_id]]
    pheno <- dataset$phenotypes
    
    if (!is.null(pheno)) {
      pheno$dataset_id <- gse_id
      pheno$sample_id <- rownames(pheno)
      phenotype_list[[gse_id]] <- pheno
    }
  }
  
  if (length(phenotype_list) > 0) {
    return(do.call(rbind, phenotype_list))
  } else {
    return(NULL)
  }
}

#' Create Batch Information
#' @param datasets_list List of dataset objects  
#' @return Vector of batch assignments
create_batch_info <- function(datasets_list) {
  batch_vector <- character()
  
  for (gse_id in names(datasets_list)) {
    n_samples <- ncol(datasets_list[[gse_id]]$expression)
    batch_vector <- c(batch_vector, rep(gse_id, n_samples))
  }
  
  return(factor(batch_vector))
}