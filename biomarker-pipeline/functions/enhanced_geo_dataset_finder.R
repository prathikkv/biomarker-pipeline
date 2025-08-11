# Enhanced GEO Dataset Finder Module
# Ultra-detailed dataset discovery for high-quality CAMK2D transcriptomics analysis
# Author: Claude Code Assistant

#' Enhanced GEO Dataset Search with Dual Database Strategy
#' @param search_terms Vector of search terms
#' @param organism Filter by organism
#' @param min_samples Minimum sample requirement
#' @param max_results Maximum results per database
#' @param recent_only Focus on recent datasets (2020+)
#' @return Enhanced dataset information data frame
enhanced_geo_search <- function(search_terms, 
                               organism = NULL,
                               min_samples = 10,
                               max_results = 100,
                               recent_only = TRUE) {
  
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("rentrez package required. Install with: install.packages('rentrez')")
  }
  
  cat("🔍 Enhanced GEO Dataset Discovery\n")
  cat("  Strategy: Dual database search (GDS + GSE)\n")
  cat("  Focus: High-quality CAMK2D transcriptomics datasets\n\n")
  
  all_datasets <- list()
  
  # Search Strategy 1: GDS Database (Curated datasets)
  cat("📊 Searching GDS (curated datasets)...\n")
  gds_results <- search_gds_database(search_terms, organism, max_results, recent_only)
  if (nrow(gds_results) > 0) {
    gds_results$database_source <- "GDS"
    all_datasets$gds <- gds_results
    cat("  ✓ Found", nrow(gds_results), "GDS datasets\n")
  }
  
  # Search Strategy 2: GSE Database (Series)
  cat("📊 Searching GSE (series datasets)...\n")
  gse_results <- search_gse_database(search_terms, organism, max_results, recent_only)
  if (nrow(gse_results) > 0) {
    gse_results$database_source <- "GSE"
    all_datasets$gse <- gse_results
    cat("  ✓ Found", nrow(gse_results), "GSE datasets\n")
  }
  
  # Search Strategy 3: PubMed-linked datasets
  cat("📊 Searching PubMed-linked datasets...\n")
  pmid_results <- search_pmid_linked_datasets(search_terms, organism, recent_only)
  if (nrow(pmid_results) > 0) {
    pmid_results$database_source <- "PMID"
    all_datasets$pmid <- pmid_results
    cat("  ✓ Found", nrow(pmid_results), "PMID-linked datasets\n")
  }
  
  # Combine all results
  if (length(all_datasets) > 0) {
    combined_datasets <- do.call(rbind, all_datasets)
    rownames(combined_datasets) <- NULL
    
    # Remove duplicates based on GEO ID
    combined_datasets <- combined_datasets[!duplicated(combined_datasets$geo_id), ]
    
    # Filter by minimum samples
    combined_datasets <- combined_datasets[combined_datasets$sample_count >= min_samples, ]
    
    cat("\n📈 Search Summary:\n")
    cat("  Total unique datasets:", nrow(combined_datasets), "\n")
    cat("  Meeting sample criteria:", nrow(combined_datasets), "\n")
    
    # Enhance with detailed metadata
    enhanced_datasets <- enhance_dataset_metadata(combined_datasets)
    
    return(enhanced_datasets)
  } else {
    return(data.frame())
  }
}

#' Search GDS Database with Enhanced Query
#' @param search_terms Search terms
#' @param organism Organism filter
#' @param max_results Maximum results
#' @param recent_only Focus on recent datasets
#' @return GDS dataset results
search_gds_database <- function(search_terms, organism, max_results, recent_only) {
  
  # Construct enhanced GDS query
  query_parts <- c()
  
  # Primary search terms with field restrictions
  if (length(search_terms) > 0) {
    camk_terms <- grep("CAMK|calcium.*calmodulin|CaMK", search_terms, value = TRUE, ignore.case = TRUE)
    disease_terms <- setdiff(search_terms, camk_terms)
    
    if (length(camk_terms) > 0) {
      camk_query <- paste0("(", paste(camk_terms, collapse = " OR "), ")")
      query_parts <- c(query_parts, camk_query)
    }
    
    if (length(disease_terms) > 0) {
      disease_query <- paste0("(", paste(disease_terms, collapse = " OR "), ")")
      query_parts <- c(query_parts, disease_query)
    }
  }
  
  # Organism filter
  if (!is.null(organism)) {
    org_query <- paste0('"', organism, '"[Organism]')
    query_parts <- c(query_parts, org_query)
  }
  
  # Recent datasets filter
  if (recent_only) {
    query_parts <- c(query_parts, '2020:2024[Publication Date]')
  }
  
  # Database type
  query_parts <- c(query_parts, 'gds[Filter]')
  
  full_query <- paste(query_parts, collapse = " AND ")
  
  tryCatch({
    search_results <- rentrez::entrez_search(
      db = "gds",
      term = full_query,
      retmax = max_results
    )
    
    if (length(search_results$ids) > 0) {
      return(fetch_enhanced_gds_summaries(search_results$ids))
    } else {
      return(data.frame())
    }
  }, error = function(e) {
    cat("  ⚠ GDS search error:", e$message, "\n")
    return(data.frame())
  })
}

#' Search GSE Database with Series Focus
#' @param search_terms Search terms
#' @param organism Organism filter  
#' @param max_results Maximum results
#' @param recent_only Focus on recent datasets
#' @return GSE dataset results
search_gse_database <- function(search_terms, organism, max_results, recent_only) {
  
  # Enhanced GSE query construction
  query_parts <- c()
  
  # Multi-level search strategy
  primary_terms <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "calcium calmodulin kinase")
  secondary_terms <- c("cardiac", "heart", "myocardial", "atrial", "ventricular")
  
  # Primary search (high specificity)
  if (any(primary_terms %in% search_terms)) {
    primary_match <- intersect(primary_terms, search_terms)
    primary_query <- paste0("(", paste(primary_match, collapse = " OR "), ")")
    query_parts <- c(query_parts, primary_query)
  }
  
  # Secondary search (broader cardiac focus)
  secondary_match <- intersect(secondary_terms, search_terms)
  if (length(secondary_match) > 0) {
    secondary_query <- paste0("(", paste(secondary_match, collapse = " OR "), ")")
    query_parts <- c(query_parts, secondary_query)
  }
  
  # Technology filters for transcriptomics
  tech_query <- "(RNA-seq OR microarray OR transcriptome OR gene expression OR single cell)"
  query_parts <- c(query_parts, tech_query)
  
  # Organism and date filters
  if (!is.null(organism)) {
    query_parts <- c(query_parts, paste0('"', organism, '"[Organism]'))
  }
  
  if (recent_only) {
    query_parts <- c(query_parts, '2020:2024[Publication Date]')
  }
  
  # GSE series only
  query_parts <- c(query_parts, 'gse[Filter]')
  
  full_query <- paste(query_parts, collapse = " AND ")
  
  tryCatch({
    search_results <- rentrez::entrez_search(
      db = "gds",
      term = full_query,
      retmax = max_results
    )
    
    if (length(search_results$ids) > 0) {
      return(fetch_enhanced_gse_summaries(search_results$ids))
    } else {
      return(data.frame())
    }
  }, error = function(e) {
    cat("  ⚠ GSE search error:", e$message, "\n")
    return(data.frame())
  })
}

#' Search PubMed-linked Datasets
#' @param search_terms Search terms
#' @param organism Organism filter
#' @param recent_only Focus on recent datasets
#' @return PubMed-linked dataset results
search_pmid_linked_datasets <- function(search_terms, organism, recent_only) {
  
  # Search PubMed for papers with GEO datasets
  pmid_query <- paste(c(
    search_terms,
    "(GEO OR Gene Expression Omnibus OR GSE OR transcriptome OR RNA-seq)",
    if (recent_only) "2020:2024[Publication Date]" else NULL
  ), collapse = " AND ")
  
  tryCatch({
    pmid_search <- rentrez::entrez_search(
      db = "pubmed",
      term = pmid_query,
      retmax = 50
    )
    
    if (length(pmid_search$ids) > 0) {
      return(extract_geo_from_pmids(pmid_search$ids))
    } else {
      return(data.frame())
    }
  }, error = function(e) {
    cat("  ⚠ PubMed search error:", e$message, "\n")
    return(data.frame())
  })
}

#' Extract GEO Datasets from PubMed IDs
#' @param pmids Vector of PubMed IDs
#' @return Data frame with GEO datasets found in publications
extract_geo_from_pmids <- function(pmids) {
  
  dataset_list <- list()
  batch_size <- 10
  
  for (i in seq(1, length(pmids), batch_size)) {
    batch_end <- min(i + batch_size - 1, length(pmids))
    batch_pmids <- pmids[i:batch_end]
    
    # Fetch abstracts
    abstracts <- rentrez::entrez_fetch(
      db = "pubmed",
      id = batch_pmids,
      rettype = "abstract",
      retmode = "text"
    )
    
    # Extract GEO accessions from abstracts
    geo_matches <- regmatches(abstracts, gregexpr("G(SE|DS)\\d{4,}", abstracts, ignore.case = TRUE))[[1]]
    
    if (length(geo_matches) > 0) {
      # Create dataset entries for found GEO IDs
      for (geo_id in unique(geo_matches)) {
        dataset_info <- list(
          geo_id = toupper(geo_id),
          accession = toupper(geo_id),
          title = paste("Dataset from PMID", batch_pmids[1]),
          summary_text = "Identified from PubMed literature search",
          organism = "Unknown",
          platform = "",
          platform_title = "",
          sample_count = 0,
          publication_date = "",
          pmid = paste(batch_pmids, collapse = ";"),
          source_database = "PMID"
        )
        dataset_list[[geo_id]] <- dataset_info
      }
    }
    
    Sys.sleep(0.5)  # Rate limiting
  }
  
  if (length(dataset_list) > 0) {
    return(do.call(rbind, lapply(dataset_list, as.data.frame, stringsAsFactors = FALSE)))
  } else {
    return(data.frame())
  }
}

#' Fetch Enhanced GDS Summaries with Detailed Metadata
#' @param geo_ids Vector of GEO IDs
#' @return Enhanced dataset information
fetch_enhanced_gds_summaries <- function(geo_ids) {
  
  dataset_list <- list()
  batch_size <- 20
  
  for (i in seq(1, length(geo_ids), batch_size)) {
    batch_end <- min(i + batch_size - 1, length(geo_ids))
    batch_ids <- geo_ids[i:batch_end]
    
    # Fetch summaries
    summaries <- rentrez::entrez_summary(db = "gds", id = batch_ids)
    
    if (length(batch_ids) == 1) {
      summaries <- list(summaries)
      names(summaries) <- batch_ids
    }
    
    # Enhanced parsing
    for (id in names(summaries)) {
      summary <- summaries[[id]]
      parsed_info <- parse_enhanced_geo_summary(summary, id, "GDS")
      
      if (!is.null(parsed_info)) {
        dataset_list[[id]] <- parsed_info
      }
    }
    
    Sys.sleep(0.5)  # Rate limiting
  }
  
  if (length(dataset_list) > 0) {
    return(do.call(rbind, lapply(dataset_list, as.data.frame, stringsAsFactors = FALSE)))
  } else {
    return(data.frame())
  }
}

#' Fetch Enhanced GSE Summaries
#' @param geo_ids Vector of GEO IDs  
#' @return Enhanced GSE dataset information
fetch_enhanced_gse_summaries <- function(geo_ids) {
  
  # Similar structure to GDS but with GSE-specific parsing
  dataset_list <- list()
  batch_size <- 20
  
  for (i in seq(1, length(geo_ids), batch_size)) {
    batch_end <- min(i + batch_size - 1, length(geo_ids))
    batch_ids <- geo_ids[i:batch_end]
    
    summaries <- rentrez::entrez_summary(db = "gds", id = batch_ids)
    
    if (length(batch_ids) == 1) {
      summaries <- list(summaries)
      names(summaries) <- batch_ids
    }
    
    for (id in names(summaries)) {
      summary <- summaries[[id]]
      parsed_info <- parse_enhanced_geo_summary(summary, id, "GSE")
      
      if (!is.null(parsed_info)) {
        dataset_list[[id]] <- parsed_info
      }
    }
    
    Sys.sleep(0.5)
  }
  
  if (length(dataset_list) > 0) {
    return(do.call(rbind, lapply(dataset_list, as.data.frame, stringsAsFactors = FALSE)))
  } else {
    return(data.frame())
  }
}

#' Parse Enhanced GEO Summary with Ultra-Detailed Metadata
#' @param summary GEO summary object
#' @param geo_id Dataset ID
#' @param source_db Source database (GDS or GSE)
#' @return Enhanced dataset information list
parse_enhanced_geo_summary <- function(summary, geo_id, source_db) {
  
  # Core information
  info <- list(
    geo_id = geo_id,
    accession = ifelse(!is.null(summary$accession), summary$accession, geo_id),
    title = ifelse(!is.null(summary$title), summary$title, ""),
    summary_text = ifelse(!is.null(summary$summary), summary$summary, ""),
    organism = ifelse(!is.null(summary$taxon), summary$taxon, "Unknown"),
    platform = ifelse(!is.null(summary$gpl), summary$gpl, ""),
    platform_title = ifelse(!is.null(summary$platformtitle), summary$platformtitle, ""),
    sample_count = ifelse(!is.null(summary$n_samples), as.numeric(summary$n_samples), 0),
    publication_date = ifelse(!is.null(summary$pdat), summary$pdat, ""),
    pmid = ifelse(!is.null(summary$pubmedids) && length(summary$pubmedids) > 0,
                 paste(summary$pubmedids, collapse = ";"), ""),
    source_database = source_db
  )
  
  # Enhanced technology classification
  info$modality <- classify_enhanced_technology(info$platform_title, info$title, info$summary_text)
  
  # Enhanced disease classification
  disease_info <- classify_enhanced_disease(info$title, info$summary_text)
  info <- c(info, disease_info)
  
  # Sample group analysis
  sample_info <- analyze_enhanced_sample_groups(summary, info$title, info$summary_text)
  info <- c(info, sample_info)
  
  # CAMK2D relevance scoring
  info$camk2d_relevance <- score_camk2d_relevance(info$title, info$summary_text, info$platform_title)
  
  # DGE suitability assessment
  info$dge_suitability <- assess_dge_suitability(info$modality, info$sample_count, sample_info)
  
  # Quality metrics
  info$data_completeness <- calculate_enhanced_completeness(summary)
  info$publication_impact <- assess_publication_impact(info$pmid, info$publication_date)
  
  return(info)
}

#' Classify Enhanced Technology Types
#' @param platform_title Platform title
#' @param title Dataset title
#' @param summary_text Summary text
#' @return Enhanced technology classification
classify_enhanced_technology <- function(platform_title, title, summary_text) {
  
  combined_text <- tolower(paste(platform_title, title, summary_text))
  
  # RNA-seq variants
  if (grepl("single.?cell|sc.?rna|10x|chromium|drop.?seq", combined_text)) {
    return("Single-cell RNA-seq")
  } else if (grepl("spatial|visium|slide.?seq|merfish", combined_text)) {
    return("Spatial transcriptomics")
  } else if (grepl("long.?read|pacbio|nanopore|isoform", combined_text)) {
    return("Long-read RNA-seq")
  } else if (grepl("rna.?seq|rnaseq|illumina.*seq|nextseq|hiseq|novaseq", combined_text)) {
    return("RNA-seq")
  }
  
  # Epigenomics
  else if (grepl("atac.?seq|assay.*transposase", combined_text)) {
    return("ATAC-seq")
  } else if (grepl("chip.?seq|chromatin.*immunoprecip", combined_text)) {
    return("ChIP-seq")
  } else if (grepl("cut.?tag|cut.?run|cleavage.*tagmentation", combined_text)) {
    return("CUT&TAG")
  } else if (grepl("bisulfite|methylation|rrbs|wgbs", combined_text)) {
    return("Methylation")
  }
  
  # Specialized sequencing
  else if (grepl("amplicon.?seq|targeted.*seq|panel", combined_text)) {
    return("Amplicon-Seq")
  } else if (grepl("ribosome.*profiling|ribo.?seq|rpf", combined_text)) {
    return("Ribosome profiling")
  } else if (grepl("cage.?seq|cap.*analysis", combined_text)) {
    return("CAGE-seq")
  }
  
  # Multi-omics
  else if (grepl("proteom|mass.*spec|lc.?ms", combined_text)) {
    return("Proteomics")
  } else if (grepl("metabolom|gc.?ms|nmr", combined_text)) {
    return("Metabolomics")
  }
  
  # Arrays
  else if (grepl("affymetrix|agilent|illumina.*array|microarray|beadchip", combined_text)) {
    return("Microarray")
  }
  
  # Default
  else {
    return("Other")
  }
}

#' Classify Enhanced Disease Types
#' @param title Dataset title
#' @param summary_text Summary text
#' @return Disease classification information
classify_enhanced_disease <- function(title, summary_text) {
  
  combined_text <- tolower(paste(title, summary_text))
  
  disease_info <- list(
    disease = "Unknown",
    disease_category = "Unknown",
    disease_specificity = 0
  )
  
  # Cardiac disease classification
  cardiac_diseases <- list(
    "Heart failure" = c("heart.failure", "hf", "failing.heart", "cardiac.failure"),
    "Atrial fibrillation" = c("atrial.fibrillation", "af", "atrial.arrhythmia", "fibrillation"),
    "Cardiomyopathy" = c("cardiomyopathy", "dcm", "icm", "hcm", "dilated.cardiomyopathy"),
    "Myocardial infarction" = c("myocardial.infarction", "mi", "heart.attack", "acute.mi"),
    "Cardiac hypertrophy" = c("cardiac.hypertrophy", "ventricular.hypertrophy", "lvh"),
    "Arrhythmia" = c("arrhythmia", "arrhythmic", "cardiac.rhythm", "bradycardia", "tachycardia"),
    "Ischemic heart disease" = c("ischemic", "ischemia", "coronary.artery", "cad"),
    "Heart disease" = c("heart.disease", "cardiac.disease", "cardiovascular.disease")
  )
  
  # Check for specific diseases
  for (disease_name in names(cardiac_diseases)) {
    patterns <- cardiac_diseases[[disease_name]]
    matches <- sum(sapply(patterns, function(p) grepl(p, combined_text)))
    
    if (matches > 0) {
      disease_info$disease <- disease_name
      disease_info$disease_category <- "Cardiac"
      disease_info$disease_specificity <- matches
      break
    }
  }
  
  # Check for non-cardiac diseases
  if (disease_info$disease == "Unknown") {
    other_diseases <- list(
      "Diabetes" = c("diabetic", "diabetes", "hyperglycemia", "insulin"),
      "Kidney disease" = c("kidney", "renal", "nephropathy"),
      "Neurological" = c("brain", "neuron", "hippocampus", "cortex"),
      "Cancer" = c("cancer", "tumor", "leukemia", "carcinoma", "aml")
    )
    
    for (disease_name in names(other_diseases)) {
      patterns <- other_diseases[[disease_name]]
      matches <- sum(sapply(patterns, function(p) grepl(p, combined_text)))
      
      if (matches > 0) {
        disease_info$disease <- disease_name
        disease_info$disease_category <- "Non-cardiac"
        disease_info$disease_specificity <- matches
        break
      }
    }
  }
  
  # If still unknown, check for general health/disease terms
  if (disease_info$disease == "Unknown") {
    if (grepl("disease|pathology|disorder|syndrome", combined_text)) {
      disease_info$disease <- "Disease (unspecified)"
      disease_info$disease_category <- "General"
    } else if (grepl("healthy|normal|control|wild.type", combined_text)) {
      disease_info$disease <- "Healthy/Control"
      disease_info$disease_category <- "Control"
    }
  }
  
  return(disease_info)
}

#' Analyze Enhanced Sample Groups
#' @param summary GEO summary object
#' @param title Dataset title  
#' @param summary_text Summary text
#' @return Enhanced sample group information
analyze_enhanced_sample_groups <- function(summary, title, summary_text) {
  
  combined_text <- tolower(paste(title, summary_text))
  
  sample_info <- list(
    n_disease = NA,
    n_control = NA,
    total_samples = ifelse(!is.null(summary$n_samples), as.numeric(summary$n_samples), 0),
    experimental_design = "Unknown",
    has_replicates = FALSE,
    sample_balance = "Unknown"
  )
  
  # Extract sample numbers from text
  sample_numbers <- extract_sample_numbers(combined_text)
  sample_info <- c(sample_info, sample_numbers)
  
  # Determine experimental design
  if (grepl("case.control|patient.*control|disease.*control", combined_text)) {
    sample_info$experimental_design <- "Case-control"
  } else if (grepl("time.course|temporal|longitudinal", combined_text)) {
    sample_info$experimental_design <- "Time-course"
  } else if (grepl("dose.response|treatment.*dose", combined_text)) {
    sample_info$experimental_design <- "Dose-response"
  } else if (grepl("knockout|transgenic|overexpression", combined_text)) {
    sample_info$experimental_design <- "Genetic model"
  } else if (grepl("single.cell", combined_text)) {
    sample_info$experimental_design <- "Single-cell"
  }
  
  # Check for replicates
  if (grepl("replicate|biological.*replicate|n\\s*=\\s*[3-9]", combined_text)) {
    sample_info$has_replicates <- TRUE
  }
  
  # Assess sample balance
  if (!is.na(sample_info$n_disease) && !is.na(sample_info$n_control)) {
    ratio <- min(sample_info$n_disease, sample_info$n_control) / max(sample_info$n_disease, sample_info$n_control)
    if (ratio >= 0.8) {
      sample_info$sample_balance <- "Well-balanced"
    } else if (ratio >= 0.5) {
      sample_info$sample_balance <- "Moderately balanced"
    } else {
      sample_info$sample_balance <- "Imbalanced"
    }
  }
  
  return(sample_info)
}

#' Extract Sample Numbers from Text
#' @param text Combined text to search
#' @return List with extracted sample numbers
extract_sample_numbers <- function(text) {
  
  sample_numbers <- list(
    n_disease = NA,
    n_control = NA,
    n_treated = NA,
    n_untreated = NA
  )
  
  # Pattern matching for sample numbers
  patterns <- list(
    disease = c("(\\d+)\\s*(patient|case|disease|affected)",
                "n\\s*=\\s*(\\d+)\\s*(patient|case|disease)"),
    control = c("(\\d+)\\s*(control|normal|healthy)",
                "n\\s*=\\s*(\\d+)\\s*(control|normal)"),
    treated = c("(\\d+)\\s*(treated|treatment)",
                "n\\s*=\\s*(\\d+)\\s*(treated)"),
    untreated = c("(\\d+)\\s*(untreated|vehicle|sham)",
                  "n\\s*=\\s*(\\d+)\\s*(untreated|vehicle)")
  )
  
  for (group in names(patterns)) {
    for (pattern in patterns[[group]]) {
      matches <- regmatches(text, regexec(pattern, text, ignore.case = TRUE))[[1]]
      if (length(matches) > 1) {
        number <- as.numeric(matches[2])
        if (!is.na(number)) {
          if (group == "disease") sample_numbers$n_disease <- number
          else if (group == "control") sample_numbers$n_control <- number
          else if (group == "treated") sample_numbers$n_treated <- number
          else if (group == "untreated") sample_numbers$n_untreated <- number
          break
        }
      }
    }
  }
  
  return(sample_numbers)
}

#' Score CAMK2D Relevance
#' @param title Dataset title
#' @param summary_text Summary text
#' @param platform_title Platform title
#' @return CAMK2D relevance score (0-10)
score_camk2d_relevance <- function(title, summary_text, platform_title) {
  
  combined_text <- tolower(paste(title, summary_text, platform_title))
  
  score <- 0
  
  # Direct CAMK mentions (high score)
  camk_direct <- c("camk2d", "camk2a", "camk2b", "camk2g", "camkii", "camk2", "cam.kinase.ii")
  direct_matches <- sum(sapply(camk_direct, function(term) grepl(term, combined_text)))
  score <- score + direct_matches * 2
  
  # Calcium/calmodulin pathway (medium score)
  calcium_terms <- c("calcium", "calmodulin", "ca2\\+", "calcium.signaling", "calcium.channel")
  calcium_matches <- sum(sapply(calcium_terms, function(term) grepl(term, combined_text)))
  score <- score + calcium_matches * 1
  
  # Cardiac context (medium score)
  cardiac_terms <- c("cardiac", "heart", "myocardial", "atrial", "ventricular", "cardiomyocyte")
  cardiac_matches <- sum(sapply(cardiac_terms, function(term) grepl(term, combined_text)))
  score <- score + cardiac_matches * 1
  
  # CAMK2D substrates (high score)
  substrates <- c("phospholamban", "pln", "ryanodine", "ryr2", "troponin", "tnni3")
  substrate_matches <- sum(sapply(substrates, function(term) grepl(term, combined_text)))
  score <- score + substrate_matches * 2
  
  # Relevant diseases (medium score)
  diseases <- c("heart.failure", "atrial.fibrillation", "arrhythmia", "cardiomyopathy")
  disease_matches <- sum(sapply(diseases, function(term) grepl(term, combined_text)))
  score <- score + disease_matches * 1
  
  # Cap at 10
  return(min(10, score))
}

#' Assess DGE Suitability
#' @param modality Technology type
#' @param sample_count Number of samples
#' @param sample_info Sample group information
#' @return DGE suitability score (0-10)
assess_dge_suitability <- function(modality, sample_count, sample_info) {
  
  score <- 0
  
  # Technology scoring
  tech_scores <- list(
    "RNA-seq" = 10,
    "Single-cell RNA-seq" = 8,
    "Spatial transcriptomics" = 7,
    "Long-read RNA-seq" = 9,
    "Microarray" = 6,
    "Other" = 2
  )
  
  score <- score + ifelse(modality %in% names(tech_scores), tech_scores[[modality]], 2)
  
  # Sample size scoring (out of 3)
  if (sample_count >= 50) score <- score + 3
  else if (sample_count >= 20) score <- score + 2
  else if (sample_count >= 10) score <- score + 1
  
  # Experimental design scoring
  if (!is.null(sample_info$experimental_design)) {
    if (sample_info$experimental_design %in% c("Case-control", "Genetic model")) {
      score <- score + 2
    } else if (sample_info$experimental_design != "Unknown") {
      score <- score + 1
    }
  }
  
  # Sample balance scoring
  if (!is.null(sample_info$sample_balance)) {
    if (sample_info$sample_balance == "Well-balanced") score <- score + 2
    else if (sample_info$sample_balance == "Moderately balanced") score <- score + 1
  }
  
  # Replicate bonus
  if (!is.null(sample_info$has_replicates) && sample_info$has_replicates) {
    score <- score + 1
  }
  
  return(min(10, score))
}

#' Calculate Enhanced Data Completeness
#' @param summary GEO summary object
#' @return Completeness score (0-100)
calculate_enhanced_completeness <- function(summary) {
  
  score <- 0
  
  # Basic metadata
  if (!is.null(summary$title) && nchar(summary$title) > 20) score <- score + 15
  if (!is.null(summary$summary) && nchar(summary$summary) > 100) score <- score + 20
  if (!is.null(summary$n_samples) && summary$n_samples > 0) score <- score + 15
  if (!is.null(summary$taxon)) score <- score + 10
  if (!is.null(summary$gpl)) score <- score + 10
  
  # Publication information
  if (!is.null(summary$pubmedids) && length(summary$pubmedids) > 0) score <- score + 20
  
  # Data availability
  if (!is.null(summary$suppfile)) score <- score + 5
  if (!is.null(summary$ftplink)) score <- score + 5
  
  return(score)
}

#' Assess Publication Impact
#' @param pmid Publication ID
#' @param publication_date Publication date
#' @return Publication impact score (0-5)
assess_publication_impact <- function(pmid, publication_date) {
  
  score <- 0
  
  # Has publication
  if (!is.null(pmid) && pmid != "") {
    score <- score + 2
    
    # Recent publication bonus
    if (!is.null(publication_date) && nchar(publication_date) >= 4) {
      year <- as.numeric(substr(publication_date, 1, 4))
      if (!is.na(year)) {
        current_year <- as.numeric(format(Sys.Date(), "%Y"))
        if (year >= (current_year - 2)) score <- score + 2  # Last 2 years
        else if (year >= (current_year - 5)) score <- score + 1  # Last 5 years
      }
    }
  }
  
  return(min(5, score))
}

#' Enhance Dataset Metadata with Advanced Analytics
#' @param dataset_df Combined dataset results
#' @return Enhanced dataset information
enhance_dataset_metadata <- function(dataset_df) {
  
  if (nrow(dataset_df) == 0) {
    return(dataset_df)
  }
  
  cat("🚀 Enhancing dataset metadata with advanced analytics\n")
  
  # Add comprehensive quality metrics
  for (i in 1:nrow(dataset_df)) {
    # Network analysis - check for related datasets
    if (!is.null(dataset_df$pmid[i]) && dataset_df$pmid[i] != "") {
      dataset_df$related_studies[i] <- count_related_studies(dataset_df$pmid[i])
    } else {
      dataset_df$related_studies[i] <- 0
    }
    
    # Temporal analysis - dataset age
    if (!is.null(dataset_df$publication_date[i]) && dataset_df$publication_date[i] != "") {
      dataset_df$years_since_publication[i] <- calculate_dataset_age(dataset_df$publication_date[i])
    } else {
      dataset_df$years_since_publication[i] <- NA
    }
    
    # Complexity scoring
    dataset_df$study_complexity[i] <- assess_study_complexity(
      dataset_df$title[i], 
      dataset_df$summary_text[i],
      dataset_df$sample_count[i],
      dataset_df$modality[i]
    )
  }
  
  # Calculate overall enhanced score
  dataset_df$enhanced_quality_score <- calculate_enhanced_quality_score(dataset_df)
  
  # Final ranking
  dataset_df <- dataset_df[order(-dataset_df$enhanced_quality_score), ]
  
  cat("✓ Enhanced metadata for", nrow(dataset_df), "datasets\n")
  
  return(dataset_df)
}

#' Count Related Studies for Network Analysis
#' @param pmid Publication ID
#' @return Number of related studies
count_related_studies <- function(pmid) {
  
  if (is.null(pmid) || pmid == "") {
    return(0)
  }
  
  tryCatch({
    # Simple heuristic - split multiple PMIDs
    pmids <- strsplit(pmid, "[;,]")[[1]]
    return(length(pmids))
  }, error = function(e) {
    return(0)
  })
}

#' Calculate Dataset Age in Years
#' @param publication_date Publication date string
#' @return Years since publication
calculate_dataset_age <- function(publication_date) {
  
  if (is.null(publication_date) || publication_date == "") {
    return(NA)
  }
  
  tryCatch({
    # Extract year from various date formats
    year_match <- regmatches(publication_date, regexpr("\\d{4}", publication_date))
    if (length(year_match) > 0) {
      pub_year <- as.numeric(year_match[1])
      current_year <- as.numeric(format(Sys.Date(), "%Y"))
      return(current_year - pub_year)
    } else {
      return(NA)
    }
  }, error = function(e) {
    return(NA)
  })
}

#' Assess Study Complexity Score
#' @param title Dataset title
#' @param summary_text Summary text
#' @param sample_count Number of samples
#' @param modality Technology type
#' @return Complexity score (0-10)
assess_study_complexity <- function(title, summary_text, sample_count, modality) {
  
  complexity_score <- 0
  combined_text <- tolower(paste(title, summary_text))
  
  # Multi-condition studies
  if (grepl("time.course|temporal|longitudinal", combined_text)) complexity_score <- complexity_score + 2
  if (grepl("dose.response|concentration", combined_text)) complexity_score <- complexity_score + 2
  if (grepl("multi.tissue|cross.tissue", combined_text)) complexity_score <- complexity_score + 2
  
  # Technical complexity
  if (!is.null(modality)) {
    if (modality %in% c("Single-cell RNA-seq", "Spatial transcriptomics", "Long-read RNA-seq")) {
      complexity_score <- complexity_score + 2
    } else if (modality == "RNA-seq") {
      complexity_score <- complexity_score + 1
    }
  }
  
  # Sample size complexity
  if (!is.null(sample_count) && sample_count > 100) complexity_score <- complexity_score + 1
  if (!is.null(sample_count) && sample_count > 200) complexity_score <- complexity_score + 1
  
  # Multi-omics indicators
  if (grepl("proteom|metabolom|multi.omics", combined_text)) complexity_score <- complexity_score + 2
  
  return(min(10, complexity_score))
}

#' Calculate Enhanced Quality Score
#' @param dataset_df Dataset data frame
#' @return Enhanced quality scores
calculate_enhanced_quality_score <- function(dataset_df) {
  
  # Base DGE suitability (40%)
  dge_component <- dataset_df$dge_suitability * 0.4
  
  # CAMK2D relevance (30%)
  camk_component <- dataset_df$camk2d_relevance * 0.3
  
  # Publication impact (15%)
  pub_component <- dataset_df$publication_impact * 0.15
  
  # Data completeness (10%)
  complete_component <- (dataset_df$data_completeness / 100) * 10 * 0.1
  
  # Study complexity bonus (5%)
  complexity_component <- ifelse(!is.null(dataset_df$study_complexity), 
                                dataset_df$study_complexity * 0.05, 0)
  
  # Calculate final score
  enhanced_score <- dge_component + camk_component + pub_component + 
                   complete_component + complexity_component
  
  return(round(enhanced_score, 2))
}

#' Enhanced Export to Excel with Specialized Worksheets
#' @param dataset_df Enhanced dataset data frame
#' @param output_file Excel output file path
#' @return Success status
export_enhanced_dataset_results <- function(dataset_df, output_file) {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx package required. Install with: install.packages('openxlsx')")
  }
  
  if (nrow(dataset_df) == 0) {
    cat("⚠ No datasets to export\n")
    return(FALSE)
  }
  
  cat("📊 Creating enhanced Excel export with specialized worksheets\n")
  
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  # Sheet 1: Executive Summary
  openxlsx::addWorksheet(wb, "Executive Summary")
  
  summary_data <- data.frame(
    Metric = c(
      "Total Datasets Discovered", "Ultra High-Quality (9-10)", "High-Quality (7-8)", 
      "DGE-Ready Datasets", "Human Studies", "Mouse Studies", "Rat Studies",
      "RNA-seq Datasets", "Single-cell Datasets", "Recent Studies (2020+)",
      "With Publications", "Multi-condition Studies", "Average CAMK2D Relevance",
      "Average DGE Suitability", "Recommended for Analysis"
    ),
    Value = c(
      nrow(dataset_df),
      sum(dataset_df$enhanced_quality_score >= 9, na.rm = TRUE),
      sum(dataset_df$enhanced_quality_score >= 7 & dataset_df$enhanced_quality_score < 9, na.rm = TRUE),
      sum(dataset_df$dge_suitability >= 8, na.rm = TRUE),
      sum(dataset_df$organism == "Homo sapiens", na.rm = TRUE),
      sum(dataset_df$organism == "Mus musculus", na.rm = TRUE),
      sum(dataset_df$organism == "Rattus norvegicus", na.rm = TRUE),
      sum(dataset_df$modality == "RNA-seq", na.rm = TRUE),
      sum(dataset_df$modality == "Single-cell RNA-seq", na.rm = TRUE),
      sum(dataset_df$years_since_publication <= 4, na.rm = TRUE),
      sum(dataset_df$publication_impact >= 2, na.rm = TRUE),
      sum(dataset_df$study_complexity >= 5, na.rm = TRUE),
      round(mean(dataset_df$camk2d_relevance, na.rm = TRUE), 2),
      round(mean(dataset_df$dge_suitability, na.rm = TRUE), 2),
      sum(dataset_df$enhanced_quality_score >= 7, na.rm = TRUE)
    )
  )
  
  openxlsx::writeData(wb, "Executive Summary", summary_data)
  
  # Sheet 2: Top Recommendations (Ultra-filtered)
  openxlsx::addWorksheet(wb, "Top Recommendations")
  
  top_datasets <- head(dataset_df[dataset_df$enhanced_quality_score >= 7, ], 25)
  top_columns <- c(
    "geo_id", "title", "organism", "sample_count", "modality",
    "disease", "enhanced_quality_score", "dge_suitability", 
    "camk2d_relevance", "publication_impact", "pmid"
  )
  
  if (nrow(top_datasets) > 0) {
    openxlsx::writeData(wb, "Top Recommendations", 
                       top_datasets[, top_columns[top_columns %in% names(top_datasets)]])
  } else {
    # Write placeholder message if no top datasets
    openxlsx::writeData(wb, "Top Recommendations", 
                       data.frame(Message = "No datasets found with quality score >= 7. Lower threshold or check data sources."))
  }
  
  # Sheet 3: DGE Analysis Ready
  openxlsx::addWorksheet(wb, "DGE Ready")
  
  dge_ready <- dataset_df[dataset_df$dge_suitability >= 8, ]
  openxlsx::writeData(wb, "DGE Ready", dge_ready)
  
  # Sheet 4: CAMK2D Specific
  openxlsx::addWorksheet(wb, "CAMK2D Specific")
  
  camk_specific <- dataset_df[dataset_df$camk2d_relevance >= 5, ]
  openxlsx::writeData(wb, "CAMK2D Specific", camk_specific)
  
  # Sheet 5: Technology Breakdown
  openxlsx::addWorksheet(wb, "Technology Analysis")
  
  # Technology analysis without dplyr dependency
  tech_analysis <- aggregate(
    cbind(enhanced_quality_score, sample_count, dge_ready = as.numeric(dataset_df$dge_suitability >= 8)) ~ modality + organism,
    data = dataset_df,
    FUN = function(x) c(count = length(x), avg = round(mean(x, na.rm = TRUE), 2))
  )
  
  # Flatten the results
  tech_summary <- data.frame(
    modality = tech_analysis$modality,
    organism = tech_analysis$organism,
    count = tech_analysis$enhanced_quality_score[, "count"],
    avg_quality = tech_analysis$enhanced_quality_score[, "avg"],
    avg_samples = tech_analysis$sample_count[, "avg"],
    dge_ready = tech_analysis$dge_ready[, "count"]
  )
  
  openxlsx::writeData(wb, "Technology Analysis", tech_summary)
  
  # Sheet 6: Complete Dataset Catalog
  openxlsx::addWorksheet(wb, "Complete Catalog")
  openxlsx::writeData(wb, "Complete Catalog", dataset_df)
  
  # Sheet 7: Dataset URLs and Download Instructions
  openxlsx::addWorksheet(wb, "Download Instructions")
  
  # Check if we have datasets for download instructions
  if (nrow(top_datasets) == 0) {
    # Create placeholder data if no high-quality datasets found
    download_data <- data.frame(
      GEO_ID = "No datasets found",
      GEO_URL = "Please lower quality threshold",
      R_Command = "# No datasets meet criteria",
      Quality_Score = 0,
      DGE_Score = 0,
      stringsAsFactors = FALSE
    )
  } else {
    # Use actual top datasets
    top_10_ids <- head(top_datasets$geo_id, 10)
    download_data <- data.frame(
      GEO_ID = top_10_ids,
      GEO_URL = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", top_10_ids),
      R_Command = paste0('gse <- getGEO("', top_10_ids, '", GSEMatrix=TRUE)'),
      Quality_Score = head(top_datasets$enhanced_quality_score, length(top_10_ids)),
      DGE_Score = head(top_datasets$dge_suitability, length(top_10_ids)),
      stringsAsFactors = FALSE
    )
  }
  
  openxlsx::writeData(wb, "Download Instructions", download_data)
  
  # Apply formatting
  apply_enhanced_formatting(wb, dataset_df)
  
  # Save workbook
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("✅ Enhanced results exported to:", output_file, "\n")
  cat("  - Total datasets:", nrow(dataset_df), "\n")
  cat("  - High-quality datasets:", sum(dataset_df$enhanced_quality_score >= 7), "\n")
  cat("  - DGE-ready datasets:", sum(dataset_df$dge_suitability >= 8), "\n")
  cat("  - Worksheets created:", length(openxlsx::worksheetNames(wb)), "\n")
  
  return(TRUE)
}

#' Apply Enhanced Formatting to Excel Workbook
#' @param wb Workbook object
#' @param dataset_df Dataset data frame
apply_enhanced_formatting <- function(wb, dataset_df) {
  
  # Header style
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#4F81BD",
    fontColour = "white"
  )
  
  # High quality style
  high_quality_style <- openxlsx::createStyle(fgFill = "#C6EFCE")
  
  # Apply to each worksheet
  for (sheet in openxlsx::worksheetNames(wb)) {
    if (sheet != "Executive Summary" && sheet != "Download Instructions") {
      # Apply header style
      openxlsx::addStyle(wb, sheet, header_style, rows = 1, 
                        cols = 1:min(20, ncol(dataset_df)), gridExpand = TRUE)
      
      # Apply conditional formatting based on quality scores if present
      if ("enhanced_quality_score" %in% names(dataset_df)) {
        high_rows <- which(dataset_df$enhanced_quality_score >= 8) + 1
        if (length(high_rows) > 0) {
          openxlsx::addStyle(wb, sheet, high_quality_style, 
                           rows = high_rows, cols = 1:ncol(dataset_df), gridExpand = TRUE)
        }
      }
    }
  }
}

cat("✓ Enhanced GEO Dataset Finder with ultra-detailed metadata loaded\n")