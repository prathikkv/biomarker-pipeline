# Cross-Species Mapping Functions
# Functions for human-mouse ortholog mapping and cross-species validation

#' Setup Cross-Species Gene Mapping
#' @param species_pair Vector of species (default: c("human", "mouse"))
#' @return List with mapping configurations and marts
setup_cross_species_mapping <- function(species_pair = c("human", "mouse")) {
  cat("Setting up cross-species gene mapping...\n")
  
  tryCatch({
    # Check if biomaRt is available
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      warning("biomaRt package not available. Using static CAMK gene mappings.")
      return(get_static_camk_mappings())
    }
    
    library(biomaRt)
    
    # Human gene mapping
    human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Mouse gene mapping  
    mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # Get CAMK family mappings between human and mouse
    camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
    
    human_to_mouse <- getBM(
      attributes = c("hgnc_symbol", "mmusculus_homolog_associated_gene_name"),
      filters = "hgnc_symbol",
      values = camk_genes,
      mart = human_mart
    )
    
    cat("✓ Cross-species gene mapping initialized\n")
    cat("  - Human mart: hsapiens_gene_ensembl\n")
    cat("  - Mouse mart: mmusculus_gene_ensembl\n")
    cat("  - CAMK mappings found:", nrow(human_to_mouse), "\n")
    
    return(list(
      human_mart = human_mart,
      mouse_mart = mouse_mart,
      gene_mapping = human_to_mouse,
      species_pair = species_pair,
      mapping_method = "biomaRt"
    ))
    
  }, error = function(e) {
    warning("biomaRt setup failed: ", e$message)
    cat("Using static CAMK gene mappings as fallback\n")
    return(get_static_camk_mappings())
  })
}

#' Get Static CAMK Gene Mappings (Fallback)
#' @return List with static human-mouse CAMK gene mappings
get_static_camk_mappings <- function() {
  # Static mappings for CAMK family genes
  gene_mapping <- data.frame(
    hgnc_symbol = c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4"),
    mmusculus_homolog_associated_gene_name = c("Camk2d", "Camk2a", "Camk2b", "Camk2g", "Camk1", "Camk4"),
    conservation_score = c(0.95, 0.98, 0.96, 0.94, 0.92, 0.89),
    stringsAsFactors = FALSE
  )
  
  return(list(
    human_mart = NULL,
    mouse_mart = NULL, 
    gene_mapping = gene_mapping,
    species_pair = c("human", "mouse"),
    mapping_method = "static"
  ))
}

#' Map Gene Symbols Between Species
#' @param gene_symbols Vector of gene symbols to map
#' @param from_species Source species ("human" or "mouse")
#' @param to_species Target species ("mouse" or "human")  
#' @param mapping_info Cross-species mapping information
#' @return Data frame with mapped gene symbols
map_gene_symbols_cross_species <- function(gene_symbols, from_species = "human", to_species = "mouse", mapping_info) {
  cat("Mapping", length(gene_symbols), "genes from", from_species, "to", to_species, "\n")
  
  if (mapping_info$mapping_method == "biomaRt") {
    return(map_with_biomart(gene_symbols, from_species, to_species, mapping_info))
  } else {
    return(map_with_static(gene_symbols, from_species, to_species, mapping_info))
  }
}

#' Map Genes Using biomaRt
#' @param gene_symbols Vector of gene symbols
#' @param from_species Source species
#' @param to_species Target species
#' @param mapping_info Mapping configuration
#' @return Data frame with mapped genes
map_with_biomart <- function(gene_symbols, from_species, to_species, mapping_info) {
  tryCatch({
    if (from_species == "human" && to_species == "mouse") {
      mapped_genes <- getBM(
        attributes = c("hgnc_symbol", "mmusculus_homolog_associated_gene_name"),
        filters = "hgnc_symbol",
        values = gene_symbols,
        mart = mapping_info$human_mart
      )
      colnames(mapped_genes) <- c("from_gene", "to_gene")
    } else if (from_species == "mouse" && to_species == "human") {
      mapped_genes <- getBM(
        attributes = c("mgi_symbol", "hsapiens_homolog_associated_gene_name"),
        filters = "mgi_symbol", 
        values = gene_symbols,
        mart = mapping_info$mouse_mart
      )
      colnames(mapped_genes) <- c("from_gene", "to_gene")
    } else {
      stop("Unsupported species mapping: ", from_species, " to ", to_species)
    }
    
    # Remove empty mappings
    mapped_genes <- mapped_genes[mapped_genes$to_gene != "", ]
    
    cat("✓ Mapped", nrow(mapped_genes), "genes using biomaRt\n")
    return(mapped_genes)
    
  }, error = function(e) {
    warning("biomaRt mapping failed: ", e$message)
    return(map_with_static(gene_symbols, from_species, to_species, mapping_info))
  })
}

#' Map Genes Using Static Mappings
#' @param gene_symbols Vector of gene symbols
#' @param from_species Source species
#' @param to_species Target species  
#' @param mapping_info Mapping configuration
#' @return Data frame with mapped genes
map_with_static <- function(gene_symbols, from_species, to_species, mapping_info) {
  if (from_species == "human" && to_species == "mouse") {
    matched_genes <- mapping_info$gene_mapping[
      mapping_info$gene_mapping$hgnc_symbol %in% gene_symbols, 
      c("hgnc_symbol", "mmusculus_homolog_associated_gene_name")
    ]
    colnames(matched_genes) <- c("from_gene", "to_gene")
  } else if (from_species == "mouse" && to_species == "human") {
    matched_genes <- mapping_info$gene_mapping[
      mapping_info$gene_mapping$mmusculus_homolog_associated_gene_name %in% gene_symbols,
      c("mmusculus_homolog_associated_gene_name", "hgnc_symbol")
    ]
    colnames(matched_genes) <- c("from_gene", "to_gene")
  } else {
    stop("Unsupported species mapping")
  }
  
  cat("✓ Mapped", nrow(matched_genes), "genes using static mappings\n")
  return(matched_genes)
}

#' Calculate Cross-Species Gene Expression Correlation
#' @param human_expression Human gene expression matrix
#' @param mouse_expression Mouse gene expression matrix
#' @param mapping_info Cross-species mapping information
#' @return List with correlation results
calculate_cross_species_correlation <- function(human_expression, mouse_expression, mapping_info) {
  cat("Calculating cross-species expression correlations...\n")
  
  # Get common genes between species
  human_genes <- rownames(human_expression)
  mouse_genes <- rownames(mouse_expression)
  
  # Map human genes to mouse
  h2m_mapping <- map_gene_symbols_cross_species(human_genes, "human", "mouse", mapping_info)
  
  # Find genes present in both datasets
  common_human <- h2m_mapping$from_gene[h2m_mapping$to_gene %in% mouse_genes]
  common_mouse <- h2m_mapping$to_gene[h2m_mapping$from_gene %in% human_genes]
  
  if (length(common_human) < 10) {
    warning("Very few common genes found between species: ", length(common_human))
    return(NULL)
  }
  
  cat("Found", length(common_human), "common genes for correlation analysis\n")
  
  # Calculate mean expression for each gene in each species
  human_means <- rowMeans(human_expression[common_human, ], na.rm = TRUE)
  mouse_means <- rowMeans(mouse_expression[common_mouse, ], na.rm = TRUE)
  
  # Calculate correlation
  correlation <- cor(human_means, mouse_means, use = "complete.obs")
  
  # Create results object
  correlation_results <- list(
    correlation = correlation,
    common_genes_count = length(common_human),
    human_genes = common_human,
    mouse_genes = common_mouse,
    human_means = human_means,
    mouse_means = mouse_means,
    mapping_used = h2m_mapping
  )
  
  cat("✓ Cross-species correlation calculated: r =", round(correlation, 3), "\n")
  
  return(correlation_results)
}

#' Assess Evolutionary Conservation of CAMK Genes
#' @param datasets_list List of datasets (human and mouse)
#' @param mapping_info Cross-species mapping information
#' @return Conservation assessment results
assess_evolutionary_conservation <- function(datasets_list, mapping_info) {
  cat("Assessing evolutionary conservation of CAMK genes...\n")
  
  # Separate human and mouse datasets
  human_indices <- sapply(datasets_list, function(x) !is.null(x$species) && x$species == "Human")
  mouse_indices <- sapply(datasets_list, function(x) !is.null(x$species) && x$species == "Mouse")
  
  human_datasets <- datasets_list[human_indices]
  mouse_datasets <- datasets_list[mouse_indices]
  
  if (length(human_datasets) == 0 || length(mouse_datasets) == 0) {
    warning("Need both human and mouse datasets for conservation analysis")
    return(NULL)
  }
  
  conservation_results <- list()
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
  
  for (gene in camk_genes) {
    # Check detection in human datasets
    human_detection <- sapply(human_datasets, function(x) {
      if (!is.null(x$camk_detection_results) && gene %in% names(x$camk_detection_results)) {
        x$camk_detection_results[[gene]]$found
      } else {
        FALSE
      }
    })
    
    # Check detection in mouse datasets  
    mouse_detection <- sapply(mouse_datasets, function(x) {
      if (!is.null(x$camk_detection_results) && gene %in% names(x$camk_detection_results)) {
        x$camk_detection_results[[gene]]$found
      } else {
        FALSE
      }
    })
    
    # Calculate conservation metrics
    human_detection_rate <- mean(human_detection, na.rm = TRUE)
    mouse_detection_rate <- mean(mouse_detection, na.rm = TRUE)
    conservation_score <- min(human_detection_rate, mouse_detection_rate)
    
    conservation_results[[gene]] <- list(
      human_detection_rate = human_detection_rate,
      mouse_detection_rate = mouse_detection_rate,
      conservation_score = conservation_score,
      human_datasets_detected = sum(human_detection),
      mouse_datasets_detected = sum(mouse_detection),
      conservation_category = ifelse(conservation_score > 0.7, "High", 
                                    ifelse(conservation_score > 0.4, "Medium", "Low"))
    )
    
    cat(gene, ": Conservation =", round(conservation_score, 2), 
        "(", conservation_results[[gene]]$conservation_category, ")\n")
  }
  
  return(conservation_results)
}

#' Create Cross-Species Expression Comparison Matrix
#' @param datasets_list List of all datasets
#' @param mapping_info Cross-species mapping information  
#' @return Cross-species comparison matrix
create_cross_species_comparison_matrix <- function(datasets_list, mapping_info) {
  cat("Creating cross-species expression comparison matrix...\n")
  
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMK1", "CAMK4")
  
  # Initialize comparison matrix
  comparison_matrix <- data.frame(
    Gene = camk_genes,
    stringsAsFactors = FALSE
  )
  
  # Add columns for each dataset
  for (dataset_id in names(datasets_list)) {
    dataset <- datasets_list[[dataset_id]]
    species <- if (!is.null(dataset$species)) dataset$species else "Unknown"
    
    col_name <- paste0(dataset_id, "_", species)
    comparison_matrix[[col_name]] <- character(length(camk_genes))
    
    for (i in 1:length(camk_genes)) {
      gene <- camk_genes[i]
      
      if (!is.null(dataset$camk_detection_results) && gene %in% names(dataset$camk_detection_results)) {
        result <- dataset$camk_detection_results[[gene]]
        if (result$found) {
          comparison_matrix[[col_name]][i] <- paste0("✓ (", result$count, ")")
        } else {
          comparison_matrix[[col_name]][i] <- "✗"
        }
      } else {
        comparison_matrix[[col_name]][i] <- "−"  # No data
      }
    }
  }
  
  cat("✓ Cross-species comparison matrix created with", ncol(comparison_matrix) - 1, "datasets\n")
  
  return(comparison_matrix)
}