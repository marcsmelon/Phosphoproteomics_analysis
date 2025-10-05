

################
# VERSION 8
###############
################################################################################
# PHOSPHOPROTEOMICS ANALYSIS: HEALTHY vs DISEASE
# Custom script for your specific data
################################################################################

# STEP 1: INSTALL AND LOAD PACKAGES ============================================
# Run once to install (uncomment if needed):
# install.packages(c("tidyverse", "reshape2", "pheatmap", "htmltools", "DT", "base64enc", "scales", "viridis", "ggrepel", "VennDiagram"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(htmltools)
  library(DT)
  library(base64enc)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(scales)
  library(viridis)
  library(VennDiagram)
})

# STEP 2: YULIA'S PREPROCESSING FUNCTION ========================================
# (Same as before - processes multiply phosphorylated sites)

process_phosphosites <- function(data, log_file = "phosphosite_processing.log") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
  })
  
  write_log <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_message <- sprintf("%s: %s\n", timestamp, message)
    cat(log_message, file = log_file, append = TRUE)
    cat(message, "\n")
  }
  
  cat("", file = log_file)
  write_log("Starting phosphosite processing...")
  
  data <- data %>% 
    dplyr::select(-any_of(c("PTM.CollapseKey", "rowname")))
  
  # Find intensity columns - modified to detect your column pattern
  intensity_cols <- grep("(Disease_|Healthy_)", names(data), value = TRUE)
  
  if(length(intensity_cols) == 0) {
    stop("Could not identify intensity columns. Looking for patterns containing 'Disease_' or 'Healthy_'")
  }
  
  meta_cols <- setdiff(names(data), intensity_cols)
  
  extract_base_name <- function(genephos) {
    cleaned <- gsub("[*]", "@STAR@", genephos)
    parts <- str_split(cleaned, "_[STY]\\d+_M\\d+$")[[1]]
    if(length(parts) > 0) {
      base <- gsub("@STAR@", "*", parts[1])
      return(base)
    }
    return(NA_character_)
  }
  
  m1_sites <- data %>% 
    filter(str_detect(GenePhos, "_M1$"))
  write_log(sprintf("Found %d M1 sites", nrow(m1_sites)))
  
  higher_mult_sites <- data %>% 
    filter(!str_detect(GenePhos, "_M1$")) %>%
    mutate(
      base_name = sapply(GenePhos, extract_base_name),
      site_aa = str_extract(GenePhos, "[STY](?=\\d+_M\\d+$)"),
      position = as.numeric(str_extract(GenePhos, "\\d+(?=_M\\d+$)")),
      multiplicity = as.numeric(str_extract(GenePhos, "(?<=M)\\d+$")),
      across(!all_of(c(intensity_cols, "position", "multiplicity")), as.character),
      across(all_of(intensity_cols), as.numeric)
    ) %>%
    filter(!is.na(position), !is.na(multiplicity), !is.na(base_name))
  
  write_log(sprintf("Found %d higher multiplicity sites", nrow(higher_mult_sites)))
  
  paired_results <- list()
  gene_groups <- split(higher_mult_sites, higher_mult_sites$base_name)
  total_genes <- length(gene_groups)
  write_log(sprintf("Processing %d protein groups", total_genes))
  
  for(gene_idx in seq_along(gene_groups)) {
    if(gene_idx %% 10 == 0) {
      cat(sprintf("\rProcessing protein group %d of %d (%.1f%%)", 
                  gene_idx, total_genes, 100*gene_idx/total_genes))
      flush.console()
    }
    
    sites <- gene_groups[[gene_idx]]
    if(nrow(sites) < 2) next
    
    sites <- sites[order(sites$position),]
    used_indices <- logical(nrow(sites))
    
    i <- 1
    while(i <= nrow(sites)) {
      if(used_indices[i]) {
        i <- i + 1
        next
      }
      
      site1 <- sites[i,]
      remaining_sites <- sites[!used_indices & 
                                 seq_len(nrow(sites)) > i & 
                                 sites$multiplicity == site1$multiplicity,, drop=FALSE]
      
      if(site1$multiplicity == 2) {
        if(nrow(remaining_sites) == 0) {
          i <- i + 1
          next
        }
        
        distances <- abs(remaining_sites$position - site1$position)
        potential_sites <- remaining_sites[distances <= 10,, drop=FALSE]
        
        if(nrow(potential_sites) > 0) {
          best_pair <- NULL
          best_shared_count <- -1
          
          for(j in 1:nrow(potential_sites)) {
            site2 <- potential_sites[j,]
            shared_count <- 0
            
            for(col in intensity_cols) {
              int1 <- site1[[col]]
              int2 <- site2[[col]]
              if(!is.na(int1) && !is.na(int2) && int1 != 0 && int2 != 0) {
                shared_count <- shared_count + 1
              }
            }
            
            if(shared_count > best_shared_count) {
              best_shared_count <- shared_count
              best_pair <- list(site = site2, index = j)
            }
          }
          
          if(!is.null(best_pair)) {
            site2 <- best_pair$site
            write_log(sprintf("Found M2 pair: %s and %s", site1$GenePhos, site2$GenePhos))
            
            result <- site1
            positions <- sort(c(site1$position, site2$position))
            site_aas <- c(site1$site_aa, site2$site_aa)[order(c(site1$position, site2$position))]
            
            result$GenePhos <- paste0(
              site1$base_name, "_",
              paste0(site_aas, positions, collapse = "_")
            )
            
            if("PTM.SiteAA" %in% names(result)) {
              result$PTM.SiteAA <- paste(site_aas, collapse = ";")
            }
            if("PTM.SiteLocation" %in% names(result)) {
              result$PTM.SiteLocation <- paste(positions, collapse = ";")
            }
            if("PTM.FlankingRegion" %in% names(result)) {
              result$PTM.FlankingRegion <- paste(site1$PTM.FlankingRegion, 
                                                 site2$PTM.FlankingRegion, sep = ";")
            }
            
            for(col in intensity_cols) {
              int1 <- site1[[col]]
              int2 <- site2[[col]]
              if(!is.na(int1) && !is.na(int2) && int1 != 0 && int2 != 0) {
                result[[col]] <- min(int1, int2)
              } else {
                result[[col]] <- NA_real_
              }
            }
            
            paired_results[[length(paired_results) + 1]] <- result
            used_indices[i] <- TRUE
            used_indices[which(!used_indices & seq_len(nrow(sites)) > i)[best_pair$index]] <- TRUE
          }
        }
      } else if(site1$multiplicity == 3) {
        if(nrow(remaining_sites) < 2) {
          i <- i + 1
          next
        }
        
        distances <- abs(remaining_sites$position - site1$position)
        potential_sites <- remaining_sites[distances <= 10,, drop=FALSE]
        
        if(nrow(potential_sites) >= 2) {
          best_triplet <- NULL
          best_shared_count <- -1
          
          for(j in 1:(nrow(potential_sites)-1)) {
            site2 <- potential_sites[j,]
            for(k in (j+1):nrow(potential_sites)) {
              site3 <- potential_sites[k,]
              if(abs(site3$position - site2$position) <= 10) {
                shared_count <- 0
                for(col in intensity_cols) {
                  int1 <- site1[[col]]
                  int2 <- site2[[col]]
                  int3 <- site3[[col]]
                  if(!is.na(int1) && !is.na(int2) && !is.na(int3) && 
                     int1 != 0 && int2 != 0 && int3 != 0) {
                    shared_count <- shared_count + 1
                  }
                }
                if(shared_count > best_shared_count) {
                  best_shared_count <- shared_count
                  best_triplet <- list(sites = list(site2, site3), indices = c(j, k))
                }
              }
            }
          }
          
          if(!is.null(best_triplet)) {
            write_log(sprintf("Found M3 triplet: %s and %s and %s",
                              site1$GenePhos,
                              best_triplet$sites[[1]]$GenePhos,
                              best_triplet$sites[[2]]$GenePhos))
            
            result <- site1
            positions <- sort(c(site1$position,
                                best_triplet$sites[[1]]$position,
                                best_triplet$sites[[2]]$position))
            site_aas <- c(site1$site_aa,
                          best_triplet$sites[[1]]$site_aa,
                          best_triplet$sites[[2]]$site_aa)[order(c(site1$position,
                                                                   best_triplet$sites[[1]]$position,
                                                                   best_triplet$sites[[2]]$position))]
            
            result$GenePhos <- paste0(site1$base_name, "_",
                                      paste0(site_aas, positions, collapse = "_"))
            
            if("PTM.SiteAA" %in% names(result)) {
              result$PTM.SiteAA <- paste(site_aas, collapse = ";")
            }
            if("PTM.SiteLocation" %in% names(result)) {
              result$PTM.SiteLocation <- paste(positions, collapse = ";")
            }
            if("PTM.FlankingRegion" %in% names(result)) {
              result$PTM.FlankingRegion <- paste(
                site1$PTM.FlankingRegion,
                best_triplet$sites[[1]]$PTM.FlankingRegion,
                best_triplet$sites[[2]]$PTM.FlankingRegion, sep = ";")
            }
            
            for(col in intensity_cols) {
              int1 <- site1[[col]]
              int2 <- best_triplet$sites[[1]][[col]]
              int3 <- best_triplet$sites[[2]][[col]]
              if(!is.na(int1) && !is.na(int2) && !is.na(int3) &&
                 int1 != 0 && int2 != 0 && int3 != 0) {
                result[[col]] <- min(int1, int2, int3)
              } else {
                result[[col]] <- NA_real_
              }
            }
            
            paired_results[[length(paired_results) + 1]] <- result
            used_indices[i] <- TRUE
            orig_indices <- which(!used_indices & seq_len(nrow(sites)) > i)[best_triplet$indices]
            used_indices[orig_indices] <- TRUE
          }
        }
      }
      i <- i + 1
    }
  }
  
  write_log("Processing results...")
  cat("\nProcessing results...\n")
  
  if(length(paired_results) > 0) {
    write_log(sprintf("Found %d pairs/triplets to collapse", length(paired_results)))
    collapsed_results <- do.call(rbind, paired_results) %>%
      dplyr::select(-c(base_name, site_aa, position, multiplicity)) %>%
      mutate(
        across(all_of(intensity_cols), as.numeric),
        across(!all_of(intensity_cols), as.character)
      )
    
    m1_sites_processed <- m1_sites %>%
      mutate(
        GenePhos = str_replace(GenePhos, "_M1$", ""),
        across(all_of(intensity_cols), as.numeric),
        across(!all_of(intensity_cols), as.character)
      )
    
    final_results <- bind_rows(m1_sites_processed, collapsed_results)
  } else {
    write_log("No pairs/triplets found, processing only M1 sites")
    final_results <- m1_sites %>%
      mutate(
        GenePhos = str_replace(GenePhos, "_M1$", ""),
        across(all_of(intensity_cols), as.numeric),
        across(!all_of(intensity_cols), as.character)
      )
  }
  
  write_log(sprintf("Final dataset contains %d rows", nrow(final_results)))
  write_log("Processing completed")
  
  return(final_results)
}

# STEP 3: MAIN ANALYSIS FUNCTION (MODIFIED FOR YOUR DATA) ======================

phospho_analysis_healthy_disease <- function(phospho_data,
                                             use_advanced_preprocessing = TRUE,
                                             filtering_threshold = 3,
                                             imputation_method = "min",
                                             normalization_method = "median_center",
                                             p_value_threshold = 0.1,
                                             create_plots = TRUE) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("        HEALTHY vs DISEASE PHOSPHOPROTEOMICS ANALYSIS\n")
  cat("================================================================================\n\n")
  
  message("Step 1/10: Loading data...")
  
  if(is.character(phospho_data)) {
    SN <- read.csv(phospho_data, header = TRUE, sep = "\t")
    if(nrow(SN) > 0 && any(grepl("^X", colnames(SN)))) {
      colnames(SN) <- gsub("^X", "", colnames(SN))
    }
    if(nrow(SN) > 1 && all(is.na(SN[1,]))) {
      SN <- SN[-1,]
    }
  } else {
    SN <- phospho_data
  }
  
  unfilt <- nrow(SN)
  message(paste("  ✓ Original number of entries:", unfilt))
  
  # Step 2: Preprocessing
  if(use_advanced_preprocessing) {
    message("\nStep 2/10: Advanced phosphosite preprocessing...")
    message("  → Collapsing multiply phosphorylated sites")
    message("  → Grouping sites by proximity (≤10 amino acids)")
    message("  → Taking minimum intensities for grouped sites")
    SN_processed <- process_phosphosites(SN, log_file = "preprocessing.log")
    SN_processed <- SN_processed %>% tibble::column_to_rownames("GenePhos")
    message(paste("  ✓ After preprocessing:", nrow(SN_processed), "phosphosites"))
    message(paste("  ✓ Reduction:", round((1 - nrow(SN_processed)/unfilt)*100, 1), "%"))
  } else {
    message("\nStep 2/10: Basic preprocessing...")
    if(!("GenePhos" %in% colnames(SN))) {
      if(all(c("PG_Genes", "PTM_CollapseKey") %in% colnames(SN))) {
        SN$phos <- gsub("^[^_]*", "", SN$PTM_CollapseKey)
        SN$GenePhos <- paste0(SN$PG_Genes, SN$phos)
      } else {
        warning("Unable to create GenePhos column. Using row numbers as IDs.")
        SN$GenePhos <- paste0("ID_", 1:nrow(SN))
      }
    }
    SN_processed <- SN %>%
      dplyr::distinct(GenePhos, .keep_all = TRUE) %>%
      tibble::column_to_rownames("GenePhos")
    message(paste("  ✓ After basic preprocessing:", nrow(SN_processed), "phosphosites"))
  }
  
  message("\nStep 3/10: Extracting intensity data...")
  
  # Find intensity columns (Disease_ or Healthy_)
  intensity_cols <- grep("(Disease_|Healthy_)", colnames(SN_processed))
  
  if(length(intensity_cols) == 0) {
    stop("No intensity columns found. Looking for patterns containing 'Disease_' or 'Healthy_'")
  }
  
  intensity_matrix <- as.matrix(SN_processed[, intensity_cols])
  intensity_matrix[intensity_matrix == "NaN"] <- NA
  storage.mode(intensity_matrix) <- "numeric"
  message(paste("  ✓ Found", ncol(intensity_matrix), "samples"))
  
  # Create condition labels from column names
  sample_names <- colnames(intensity_matrix)
  condition <- ifelse(grepl("Healthy_", sample_names), "Healthy", "Disease")
  condition <- factor(condition, levels = c("Healthy", "Disease"))
  
  message(paste("  ✓ Healthy samples:", sum(condition == "Healthy")))
  message(paste("  ✓ Disease samples:", sum(condition == "Disease")))
  
  # Create shortened sample names for better visualization
  sample_names_short <- paste0(condition, "_", ave(as.numeric(condition), condition, FUN = seq_along))
  names(sample_names_short) <- sample_names
  
  # Create sample mapping for reference
  sample_mapping <- data.frame(
    Short_Name = sample_names_short,
    Full_Name = names(sample_names_short),
    Condition = condition,
    stringsAsFactors = FALSE
  )
  
  message("\nStep 4/10: Filtering data...")
  
  if(is.numeric(filtering_threshold)) {
    is_valid <- rowSums(!is.na(intensity_matrix) & intensity_matrix > 0) >= filtering_threshold
    message(paste("  → Filtering rows present in at least", filtering_threshold, "samples"))
  } else if(filtering_threshold == "by_condition") {
    healthy_cols <- which(condition == "Healthy")
    disease_cols <- which(condition == "Disease")
    is_valid1 <- rowSums(!is.na(intensity_matrix[, healthy_cols]) & intensity_matrix[, healthy_cols] > 0) >= 2
    is_valid2 <- rowSums(!is.na(intensity_matrix[, disease_cols]) & intensity_matrix[, disease_cols] > 0) >= 2
    is_valid <- is_valid1 & is_valid2
    message("  → Filtering to phosphosites present in at least 2 samples per condition")
  } else {
    is_valid <- rep(TRUE, nrow(intensity_matrix))
  }
  
  intensity_matrix <- intensity_matrix[is_valid, ]
  message(paste("  ✓ Remaining phosphosites:", nrow(intensity_matrix)))
  
  message("\nStep 5/10: Imputing missing values...")
  intensity_matrix[intensity_matrix == 0] <- NA
  
  if(imputation_method == "min") {
    for(i in 1:nrow(intensity_matrix)) {
      if(all(is.na(intensity_matrix[i, ]))) next
      row_min <- min(intensity_matrix[i, ], na.rm = TRUE) / 2
      intensity_matrix[i, is.na(intensity_matrix[i, ])] <- row_min
    }
    message("  ✓ Minimum value imputation performed")
  } else if(imputation_method == "mean") {
    for(i in 1:nrow(intensity_matrix)) {
      if(all(is.na(intensity_matrix[i, ]))) next
      row_mean <- mean(intensity_matrix[i, ], na.rm = TRUE)
      intensity_matrix[i, is.na(intensity_matrix[i, ])] <- row_mean
    }
    message("  ✓ Mean imputation performed")
  }
  
  message("\nStep 6/10: Log2 transforming data...")
  min_nonzero <- min(intensity_matrix[intensity_matrix > 0], na.rm = TRUE)
  offset <- ifelse(min_nonzero < 1, 1, 0)
  log2_intensity_matrix <- log2(intensity_matrix + offset)
  message("  ✓ Log2 transformation complete")
  
  # Create distinct colors for Healthy vs Disease
  # Using very different colors: Blue for Healthy, Red/Orange for Disease
  condition_colors <- c("Healthy" = "#2E86AB", "Disease" = "#E63946")
  # Alternative bright options:
  # condition_colors <- c("Healthy" = "#0077BB", "Disease" = "#EE7733")  # Blue-Orange
  # condition_colors <- c("Healthy" = "#33A02C", "Disease" = "#E31A1C")  # Green-Red
  
  # Initialize plots list
  plots <- list()
  
  # QC GRAPHS ===================================================================
  message("\nGenerating QC graphs...")
  
  # Calculate PTM site stats (always, for results object)
  tmp <- gsub("_M.$", "", rownames(intensity_matrix))
  tmp <- gsub("^.*_", "", tmp)
  tmp <- gsub("[0-9]{1,9}$", "", tmp)
  ptm_site_counts <- table(tmp)
  ptm_site_percentages <- (ptm_site_counts / nrow(intensity_matrix)) * 100
  
  ptm_site_stats <- data.frame(
    SiteAA = names(ptm_site_counts),
    Count = as.integer(ptm_site_counts),
    Percent = as.numeric(ptm_site_percentages)
  )
  
  # Save PTM site stats
  write.csv(ptm_site_stats, "ptm_site_statistics.csv", row.names = FALSE)
  message("  ✓ PTM site statistics saved")
  
  if(create_plots) {
    # QC 1: PTM Site Composition Pie Chart
    message("  → Creating PTM composition plot...")
    png("ptm_site_composition.png", width = 800, height = 600)
    par(mar = c(5, 5, 4, 2) + 0.1)
    pie_colors <- c(S = "#4477AA", T = "#AA4488", Y = "#44AA77")
    pie(ptm_site_stats$Count, 
        labels = paste0(ptm_site_stats$SiteAA, "\n", 
                        ptm_site_stats$Count, " (", 
                        round(ptm_site_stats$Percent, 1), "%)"),
        col = pie_colors[ptm_site_stats$SiteAA],
        main = "PTM Site Composition")
    dev.off()
    plots$ptm_composition <- "ptm_site_composition.png"
    
    # QC 2: Total Phosphosites Detected per Sample
    message("  → Phosphosites detected per sample...")
    tmp_phospho <- data.frame(
      Sample = colnames(intensity_matrix),
      Sample_Short = sample_names_short[colnames(intensity_matrix)],
      Total_Phosphosites = colSums(!is.na(intensity_matrix) & intensity_matrix > 0),
      Condition = condition
    )
    tmp_phospho$Sample_Short <- factor(tmp_phospho$Sample_Short, 
                                       levels = sample_names_short[colnames(intensity_matrix)])
    
    p_phospho <- ggplot2::ggplot(tmp_phospho, ggplot2::aes(x = Sample_Short, y = Total_Phosphosites, fill = Condition)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = condition_colors) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Phosphosites with Intensity > 0") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90, size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        legend.position = "right"
      ) +
      ggplot2::labs(y = "Number of Phosphosites", x = "Sample") +
      ggplot2::geom_text(ggplot2::aes(label = Total_Phosphosites), 
                         vjust = -0.5, size = 3)
    
    ggplot2::ggsave("phosphosites_detected.png", p_phospho, width = 10, height = 6, dpi = 300)
    plots$phosphosites_detected <- "phosphosites_detected.png"
    
    # QC 3: Library Size (Total Intensities per Sample)
    message("  → Library size analysis...")
    tmp_lib <- data.frame(
      Sample = colnames(intensity_matrix),
      Sample_Short = sample_names_short[colnames(intensity_matrix)],
      Total_Intensities = colSums(intensity_matrix, na.rm = TRUE),
      Condition = condition
    )
    tmp_lib$Sample_Short <- factor(tmp_lib$Sample_Short, 
                                   levels = sample_names_short[colnames(intensity_matrix)])
    
    p_lib <- ggplot2::ggplot(tmp_lib, ggplot2::aes(x = Sample_Short, y = Total_Intensities, fill = Condition)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = condition_colors) +
      ggplot2::theme_classic() +
      ggplot2::ggtitle("Library Size (Total Raw Intensities)") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90, size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        legend.position = "right"
      ) +
      ggplot2::labs(y = "Total Intensity", x = "Sample") +
      ggplot2::scale_y_continuous(labels = scales::comma)
    
    ggplot2::ggsave("library_size.png", p_lib, width = 10, height = 6, dpi = 300)
    plots$library_size <- "library_size.png"
    
    # QC 4: Density Plot of Log2 Intensities (Before Normalization)
    message("  → Density plots...")
    png("density_plot_raw.png", width = 800, height = 600)
    par(mar = c(5, 5, 4, 5) + 0.1, bty = "L")
    limma::plotDensities(log2_intensity_matrix, 
                         group = condition,
                         col = condition_colors[as.character(condition)],
                         legend = "topright",
                         main = "Density Plot - Raw Log2 Intensities")
    dev.off()
    plots$density_raw <- "density_plot_raw.png"
    
    message("  ✓ QC graphs complete")
  }
  
  plots <- list()
  
  # Create shortened sample names for better visualization
  sample_names_short <- paste0(condition, "_", seq_along(condition))
  names(sample_names_short) <- colnames(log2_intensity_matrix)
  
  if(create_plots) {
    tmp_before <- reshape2::melt(t(log2_intensity_matrix))
    colnames(tmp_before) <- c("Sample", "Protein", "Intensity")
    tmp_before$Condition <- condition[match(tmp_before$Sample, colnames(log2_intensity_matrix))]
    tmp_before$Sample_Short <- sample_names_short[tmp_before$Sample]
    tmp_before <- tmp_before[order(tmp_before$Condition), ]
    tmp_before$Sample_Short <- factor(tmp_before$Sample_Short, levels = unique(tmp_before$Sample_Short))
    
    boxplot_before <- ggplot2::ggplot(tmp_before, ggplot2::aes(x = Sample_Short, y = Intensity, fill = Condition)) + 
      ggplot2::geom_boxplot(outlier.size = 0.5) + 
      ggplot2::scale_fill_manual(values = condition_colors) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90, size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10)
      ) + 
      ggplot2::ggtitle("Log2 Intensities Before Normalization") +
      ggplot2::labs(x = "Sample", y = "Log2 Intensity") +
      ggplot2::ylim(range(tmp_before$Intensity, na.rm = TRUE))
    
    ggplot2::ggsave("boxplot_before_normalization.png", boxplot_before, width = 10, height = 6, dpi = 300)
    plots$boxplot_before <- "boxplot_before_normalization.png"
  }
  
  message("\nStep 7/10: Normalizing data...")
  if(!is.null(normalization_method) && normalization_method == "quantile") {
    temp_elist <- new("EList", list(E = intensity_matrix))
    temp_elist <- limma::normalizeBetweenArrays(temp_elist, method = "quantile")
    min_nonzero <- min(temp_elist$E[temp_elist$E > 0], na.rm = TRUE)
    offset <- ifelse(min_nonzero < 1, 1, 0)
    normalized_matrix <- log2(temp_elist$E + offset)
    message("  ✓ Quantile normalization complete")
  } else if(!is.null(normalization_method) && normalization_method == "median_center") {
    normalized_matrix <- apply(log2_intensity_matrix, 2, function(y) {
      med <- median(y, na.rm = TRUE)
      y - med
    })
    message("  ✓ Median centering normalization complete")
  } else {
    normalized_matrix <- log2_intensity_matrix
    message("  ✓ No normalization performed")
  }
  
  if(create_plots) {
    tmp_after <- reshape2::melt(t(normalized_matrix))
    colnames(tmp_after) <- c("Sample", "Protein", "Intensity")
    tmp_after$Condition <- condition[match(tmp_after$Sample, colnames(normalized_matrix))]
    tmp_after$Sample_Short <- sample_names_short[tmp_after$Sample]
    tmp_after <- tmp_after[order(tmp_after$Condition), ]
    tmp_after$Sample_Short <- factor(tmp_after$Sample_Short, levels = unique(tmp_after$Sample_Short))
    
    boxplot_after <- ggplot2::ggplot(tmp_after, ggplot2::aes(x = Sample_Short, y = Intensity, fill = Condition)) + 
      ggplot2::geom_boxplot(outlier.size = 0.5) + 
      ggplot2::scale_fill_manual(values = condition_colors) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90, size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10)
      ) + 
      ggplot2::ggtitle("Normalized Intensities (Median Centered)") +
      ggplot2::labs(x = "Sample", y = "Normalized Intensity") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)
    
    ggplot2::ggsave("boxplot_after_normalization.png", boxplot_after, width = 10, height = 6, dpi = 300)
    plots$boxplot_after <- "boxplot_after_normalization.png"
    
    # Density plot after normalization
    png("density_plot_normalized.png", width = 800, height = 600)
    par(mar = c(5, 5, 4, 5) + 0.1, bty = "L")
    limma::plotDensities(normalized_matrix, 
                         group = condition,
                         col = condition_colors[as.character(condition)],
                         legend = "topright",
                         main = "Density Plot - Normalized Intensities")
    abline(v = 0, lty = 2, col = "gray50")
    dev.off()
    plots$density_normalized <- "density_plot_normalized.png"
    
    # QC 5: Pearson Correlation Heatmap
    message("  → Pearson correlation analysis...")
    correlation_matrix <- cor(normalized_matrix, method = "pearson", use = "pairwise.complete.obs")
    
    # Calculate average correlation (upper triangle only, excluding diagonal)
    upper_tri_values <- correlation_matrix[upper.tri(correlation_matrix)]
    average_R <- mean(upper_tri_values, na.rm = TRUE)
    message(paste("    Average Pearson correlation:", round(average_R, 3)))
    
    # Use short names for correlation matrix
    rownames(correlation_matrix) <- sample_names_short[rownames(correlation_matrix)]
    colnames(correlation_matrix) <- sample_names_short[colnames(correlation_matrix)]
    
    # Save correlation matrix
    write.csv(correlation_matrix, "pearson_correlation_matrix.csv", row.names = TRUE)
    
    # Melt correlation matrix for ggplot
    melted_cor <- reshape2::melt(correlation_matrix)
    
    # Create heatmap
    p_cor <- ggplot2::ggplot(data = melted_cor, ggplot2::aes(Var1, Var2, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_viridis_c(option = "C", 
                                    name = "Pearson\nCorrelation", 
                                    limits = c(0.5, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = "Pearson Correlation Heatmap (Normalized Data)", 
                    x = "", y = "") +
      ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), 
                         size = 2.5, color = "white")
    
    ggplot2::ggsave("pearson_correlation_heatmap.png", p_cor, width = 10, height = 9, dpi = 300)
    plots$correlation_heatmap <- "pearson_correlation_heatmap.png"
    
    # QC 6: Venn Diagram - Phosphosite Detection Overlap
    message("  → Creating Venn diagram...")
    
    # Identify which phosphosites are detected in each condition
    healthy_cols <- which(condition == "Healthy")
    disease_cols <- which(condition == "Disease")
    
    # Use more stringent criteria: detected in at least 50% of samples in that condition
    min_samples_healthy <- ceiling(length(healthy_cols) * 0.5)
    min_samples_disease <- ceiling(length(disease_cols) * 0.5)
    
    # Phosphosites consistently detected in Healthy (≥50% of samples)
    healthy_detected <- rownames(intensity_matrix)[
      rowSums(!is.na(intensity_matrix[, healthy_cols]) & 
                intensity_matrix[, healthy_cols] > 0) >= min_samples_healthy
    ]
    
    # Phosphosites consistently detected in Disease (≥50% of samples)
    disease_detected <- rownames(intensity_matrix)[
      rowSums(!is.na(intensity_matrix[, disease_cols]) & 
                intensity_matrix[, disease_cols] > 0) >= min_samples_disease
    ]
    
    # Calculate overlaps
    healthy_only <- setdiff(healthy_detected, disease_detected)
    disease_only <- setdiff(disease_detected, healthy_detected)
    both <- intersect(healthy_detected, disease_detected)
    
    message(paste("    ✓ Healthy-specific (≥50% samples):", length(healthy_only)))
    message(paste("    ✓ Both conditions (≥50% samples):", length(both)))
    message(paste("    ✓ Disease-specific (≥50% samples):", length(disease_only)))
    
    # Save Venn diagram data
    venn_stats <- data.frame(
      Category = c("Healthy Only", "Both Conditions", "Disease Only", 
                   "Total Detected in Healthy", "Total Detected in Disease"),
      Count = c(
        length(healthy_only),
        length(both),
        length(disease_only),
        length(healthy_detected),
        length(disease_detected)
      )
    )
    write.csv(venn_stats, "venn_diagram_statistics.csv", row.names = FALSE)
    
    # Save condition-specific phosphosites
    if(length(healthy_only) > 0) {
      write.csv(data.frame(GenePhos = healthy_only), 
                "phosphosites_healthy_specific.csv", row.names = FALSE)
    }
    if(length(disease_only) > 0) {
      write.csv(data.frame(GenePhos = disease_only), 
                "phosphosites_disease_specific.csv", row.names = FALSE)
    }
    if(length(both) > 0) {
      write.csv(data.frame(GenePhos = both), 
                "phosphosites_both_conditions.csv", row.names = FALSE)
    }
    
    # Create Venn diagram with proper positioning
    venn.plot <- VennDiagram::venn.diagram(
      x = list(
        Healthy = healthy_detected,
        Disease = disease_detected
      ),
      filename = NULL,
      category.names = c("Healthy", "Disease"),
      
      # Output
      output = TRUE,
      imagetype = "png",
      height = 2000,
      width = 2000,
      resolution = 300,
      compression = "lzw",
      
      # Circles - proper positioning
      lwd = 3,
      lty = 'solid',
      fill = c(condition_colors["Healthy"], condition_colors["Disease"]),
      alpha = 0.6,
      
      # Numbers
      cex = 2,
      fontface = "bold",
      fontfamily = "sans",
      
      # Category names
      cat.cex = 2,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = c(condition_colors["Healthy"], condition_colors["Disease"]),
      
      # Positioning to ensure circles overlap properly but not completely
      euler.d = TRUE,
      scaled = TRUE
    )
    
    # Save Venn diagram
    png("venn_diagram_phosphosites.png", width = 2000, height = 2000, res = 300, bg = "white")
    grid::grid.draw(venn.plot)
    dev.off()
    plots$venn_diagram <- "venn_diagram_phosphosites.png"
  }
  
  message("\nStep 8/10: Setting up experimental design...")
  # Simple two-group comparison (unpaired)
  design <- model.matrix(~0 + condition)
  colnames(design) <- levels(condition)
  rownames(design) <- colnames(normalized_matrix)
  contrast <- limma::makeContrasts(Disease - Healthy, levels = design)
  message("  ✓ Design matrix created (unpaired comparison)")
  
  message("\nStep 9/10: Fitting statistical model...")
  fit <- limma::lmFit(normalized_matrix, design = design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  message("  ✓ Model fitting complete")
  
  message("\nStep 10/10: Collecting results...")
  tt <- limma::topTable(fit, n = Inf, coef = 1)
  tt$GenePhos <- rownames(tt)
  write.csv(tt, "complete_phosphosite_statistics.csv", row.names = FALSE)
  
  sig_results <- tt[tt$adj.P.Val < p_value_threshold, ]
  DT <- limma::decideTests(fit, p.value = p_value_threshold)
  colnames(DT) <- "Disease vs Healthy"
  result_summary <- summary(DT)
  
  # Export normalized data matrix
  normalized_export <- as.data.frame(normalized_matrix)
  normalized_export$GenePhos <- rownames(normalized_export)
  normalized_export <- normalized_export[, c("GenePhos", colnames(normalized_matrix))]
  write.table(normalized_export, "normalized_phospho_data.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Export sample name mapping (already created earlier)
  write.csv(sample_mapping, "sample_name_mapping.csv", row.names = FALSE)
  
  message(paste("  ✓ Complete statistics saved to 'complete_phosphosite_statistics.csv'"))
  message(paste("  ✓ Normalized data saved to 'normalized_phospho_data.tsv'"))
  message(paste("  ✓ Sample mapping saved to 'sample_name_mapping.csv'"))
  
  if(create_plots) {
    message("\nCreating visualizations...")
    
    # Enhanced Volcano plot with top gene labels
    message("  → Creating enhanced volcano plot...")
    
    # Prepare data for volcano plot
    volcano_data <- tt %>%
      mutate(
        Expression = case_when(
          logFC >= log2(2) & adj.P.Val <= p_value_threshold ~ "Up-regulated in Disease",
          logFC <= -log2(2) & adj.P.Val <= p_value_threshold ~ "Up-regulated in Healthy",
          TRUE ~ "Unchanged"
        )
      )
    
    # Split data for layering
    data_unchanged <- volcano_data %>% filter(Expression == "Unchanged")
    data_updown <- volcano_data %>% filter(Expression != "Unchanged")
    
    # Identify top genes to label (top 10 each direction)
    top_n_genes <- 10
    top_genes <- bind_rows(
      data_updown %>%
        filter(Expression == "Up-regulated in Disease") %>%
        arrange(adj.P.Val, desc(abs(logFC))) %>%
        head(top_n_genes),
      data_updown %>%
        filter(Expression == "Up-regulated in Healthy") %>%
        arrange(adj.P.Val, desc(abs(logFC))) %>%
        head(top_n_genes)
    )
    
    # Create volcano plot
    p_volcano <- ggplot2::ggplot() +
      # Layer 1: Unchanged points (20% opacity, no stroke)
      ggplot2::geom_point(
        data = data_unchanged,
        ggplot2::aes(x = logFC, y = -log10(P.Value)),
        shape = 21,
        size = 1.5,
        fill = "gray60",
        color = "gray60",
        alpha = 0.2,
        stroke = 0
      ) +
      # Layer 2: Significant points (full opacity, black stroke)
      ggplot2::geom_point(
        data = data_updown,
        ggplot2::aes(x = logFC, y = -log10(P.Value), fill = Expression),
        shape = 21,
        size = 1.5,
        color = "black",
        stroke = 0.3,
        alpha = 0.8
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Up-regulated in Disease" = "#E63946",
          "Up-regulated in Healthy" = "#2E86AB"
        )
      ) +
      # Add threshold lines
      ggplot2::geom_hline(yintercept = -log10(p_value_threshold), 
                          linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)), 
                          linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      # Labels and theme
      ggplot2::xlab(expression("log"[2]*" Fold Change")) +
      ggplot2::ylab(expression("-log"[10]*" P-value")) +
      ggplot2::ggtitle("Volcano Plot: Disease vs Healthy") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10),
        # Transparent backgrounds
        panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
      ) +
      ggplot2::labs(fill = "Regulation")
    
    # Add labels for top genes if ggrepel is available
    if(requireNamespace("ggrepel", quietly = TRUE)) {
      p_volcano <- p_volcano +
        ggrepel::geom_label_repel(
          data = top_genes,
          ggplot2::aes(x = logFC, y = -log10(P.Value), label = GenePhos),
          size = 2.5,
          box.padding = 0.5,
          point.padding = 0.3,
          max.overlaps = 20,
          min.segment.length = 0,
          seed = 42
        )
      message("    ✓ Top genes labeled")
    } else {
      message("    Note: Install 'ggrepel' package for gene labels")
    }
    
    ggplot2::ggsave("volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 300, bg = "transparent")
    plots$volcano_plot <- "volcano_plot.png"
    
    # Also save as PDF (white background is fine for PDF)
    ggplot2::ggsave("volcano_plot.pdf", p_volcano, width = 10, height = 8)
    
    # Export top genes table
    write.csv(top_genes, "top_genes_volcano_plot.csv", row.names = FALSE)
    message("    ✓ Top genes saved to 'top_genes_volcano_plot.csv'")
    
    # MDS plot by condition
    png("mds_plot_condition.png", width = 800, height = 600)
    par(mar = c(5, 5, 4, 5) + 0.1)
    plot_colors <- condition_colors[as.character(condition)]
    limma::plotMDS(normalized_matrix, 
                   col = plot_colors,
                   main = "MDS Plot - Healthy vs Disease",
                   pch = 16, cex = 2)
    legend("topright", legend = levels(condition), 
           col = condition_colors,
           pch = 16, bg = "white", cex = 1.2)
    dev.off()
    plots$mds_condition <- "mds_plot_condition.png"
    
    message("  ✓ Visualizations saved")
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("                           ANALYSIS COMPLETE!\n")
  cat("================================================================================\n\n")
  cat("SUMMARY:\n")
  cat(sprintf("  • Total phosphosites analyzed: %d\n", nrow(normalized_matrix)))
  cat(sprintf("  • Significant phosphosites (adj.P.Val < %.2f): %d\n", p_value_threshold, nrow(sig_results)))
  cat(sprintf("  • Up-regulated in Disease: %d\n", sum(sig_results$logFC > 0)))
  cat(sprintf("  • Up-regulated in Healthy: %d\n", sum(sig_results$logFC < 0)))
  if(create_plots) {
    cat(sprintf("  • Average sample correlation: %.3f\n", average_R))
    both_count <- venn_stats$Count[venn_stats$Category == "Both Conditions"]
    healthy_only_count <- venn_stats$Count[venn_stats$Category == "Healthy Only"]
    disease_only_count <- venn_stats$Count[venn_stats$Category == "Disease Only"]
    total_detected <- both_count + healthy_only_count + disease_only_count
    if(total_detected > 0) {
      cat(sprintf("  • Consistently detected in both: %d (%.1f%%)\n", 
                  both_count, (both_count/total_detected)*100))
    }
  }
  cat("\nOUTPUT FILES:\n")
  cat("  • complete_phosphosite_statistics.csv\n")
  cat("  • normalized_phospho_data.tsv\n")
  cat("  • sample_name_mapping.csv\n")
  cat("  • ptm_site_statistics.csv\n")
  if(use_advanced_preprocessing) cat("  • preprocessing.log\n")
  if(create_plots) {
    cat("  • ptm_site_composition.png\n")
    cat("  • phosphosites_detected.png\n")
    cat("  • library_size.png\n")
    cat("  • density_plot_raw.png\n")
    cat("  • boxplot_before_normalization.png\n")
    cat("  • boxplot_after_normalization.png\n")
    cat("  • density_plot_normalized.png\n")
    cat("  • pearson_correlation_matrix.csv\n")
    cat("  • pearson_correlation_heatmap.png\n")
    cat("  • venn_diagram_statistics.csv\n")
    cat("  • venn_diagram_phosphosites.png\n")
    cat("  • phosphosites_healthy_specific.csv (if any)\n")
    cat("  • phosphosites_disease_specific.csv (if any)\n")
    cat("  • phosphosites_both_conditions.csv\n")
    cat("  • volcano_plot.png\n")
    cat("  • volcano_plot.pdf\n")
    cat("  • top_genes_volcano_plot.csv\n")
    cat("  • mds_plot_condition.png\n")
  }
  cat("\n")
  
  return(list(
    raw_data = intensity_matrix,
    log2_data = log2_intensity_matrix,
    normalized_data = normalized_matrix,
    correlation_matrix = if(create_plots) correlation_matrix else NULL,
    average_correlation = if(create_plots) average_R else NULL,
    venn_stats = if(create_plots) venn_stats else NULL,
    design = design,
    fit = fit,
    top_table = tt,
    significant_results = sig_results,
    top_genes_labeled = if(create_plots) top_genes else NULL,
    results_summary = result_summary,
    ptm_site_stats = ptm_site_stats,
    sample_mapping = sample_mapping,
    condition = condition,
    p_value_threshold = p_value_threshold,
    preprocessing_used = use_advanced_preprocessing,
    plots = plots
  ))
}

# STEP 4: RUN THE ANALYSIS ======================================================

# Set your working directory
setwd("~/Documents/Results/Phosphoproteomics/Timecourse_Experiment/Analysis_timecourse_phosphoproteomics/Phosphoproteomics_2D_hCO/Phosphoproteomics_2D_hCO/TNNT2/TNNT2_D15")

# Run the analysis
results <- phospho_analysis_healthy_disease(
  phospho_data = "Phospho(STY).tsv",  # <-- CHANGE THIS!
  use_advanced_preprocessing = TRUE,     # Use sophisticated preprocessing
  filtering_threshold = 3,               # Present in at least 3 samples
  imputation_method = "min",             # Half minimum imputation
  normalization_method = "median_center", # Median centering
  p_value_threshold = 0.1,               # FDR < 0.1
  create_plots = TRUE                    # Generate plots
)

# STEP 5: EXPLORE RESULTS =======================================================

# View top 20 most significant
head(results$top_table, 20)

# Export significant phosphosites
write.csv(results$significant_results, 
          "significant_phosphosites_disease_vs_healthy.csv", 
          row.names = FALSE)

# Look at specific proteins
proteins_of_interest <- c("PKP2", "GYS1", "ALPK3", "SRC", "PLN")
poi_results <- results$top_table %>%
  filter(grepl(paste(proteins_of_interest, collapse = "|"), GenePhos)) %>%
  arrange(adj.P.Val)
print(poi_results)

# Summary statistics
cat("\n=== RESULT SUMMARY ===\n")
cat("Up in Disease:", sum(results$significant_results$logFC > 0), "\n")
cat("Up in Healthy:", sum(results$significant_results$logFC < 0), "\n")

# View PTM site composition
cat("\n=== PTM SITE COMPOSITION ===\n")
print(results$ptm_site_stats)

# View sample name mapping
cat("\n=== SAMPLE NAME MAPPING ===\n")
cat("(Short names used in plots for clarity)\n\n")
print(results$sample_mapping)

# View correlation results
if(!is.null(results$average_correlation)) {
  cat("\n=== PEARSON CORRELATION ANALYSIS ===\n")
  cat(sprintf("Average correlation between samples: %.3f\n", results$average_correlation))
  cat("\nCorrelation Matrix:\n")
  print(round(results$correlation_matrix, 3))
  cat("\nInterpretation:\n")
  cat("  • High correlation (>0.9): Samples are very similar\n")
  cat("  • Medium correlation (0.7-0.9): Samples are moderately similar\n")
  cat("  • Low correlation (<0.7): Samples are different or may have quality issues\n")
}

# View top labeled genes from volcano plot
if(!is.null(results$top_genes_labeled)) {
  cat("\n=== TOP GENES LABELED IN VOLCANO PLOT ===\n")
  cat(sprintf("(Top %d most significant in each direction)\n\n", nrow(results$top_genes_labeled)/2))
  top_genes_display <- results$top_genes_labeled %>%
    select(GeneTPhos, logFC, P.Value, adj.P.Val, Expression) %>%
    arrange(Expression, adj.P.Val)
  print(as.data.frame(top_genes_display))
}

# View Venn diagram statistics
if(!is.null(results$venn_stats)) {
  cat("\n=== PHOSPHOSITE DETECTION OVERLAP ===\n")
  cat("(Detection = present in ≥50% of samples in that condition)\n\n")
  print(results$venn_stats)
  cat("\nInterpretation:\n")
  cat("  • Healthy Only: Consistently detected in healthy but not disease\n")
  cat("  • Both Conditions: Consistently detected in both conditions\n")
  cat("  • Disease Only: Consistently detected in disease but not healthy\n")
  
  # Calculate percentages based on total unique detected sites
  healthy_only <- results$venn_stats$Count[results$venn_stats$Category == "Healthy Only"]
  both <- results$venn_stats$Count[results$venn_stats$Category == "Both Conditions"]
  disease_only <- results$venn_stats$Count[results$venn_stats$Category == "Disease Only"]
  total_unique <- healthy_only + both + disease_only
  
  if(total_unique > 0) {
    cat(sprintf("\nPercentages (of %d total detected sites):\n", total_unique))
    cat(sprintf("  • Healthy Only: %d (%.1f%%)\n", healthy_only, (healthy_only/total_unique)*100))
    cat(sprintf("  • Both Conditions: %d (%.1f%%)\n", both, (both/total_unique)*100))
    cat(sprintf("  • Disease Only: %d (%.1f%%)\n", disease_only, (disease_only/total_unique)*100))
  }
}

# View QC plots
cat("\n=== QC PLOTS GENERATED ===\n")
if(length(results$plots) > 0) {
  for(plot_name in names(results$plots)) {
    cat(sprintf("  • %s: %s\n", plot_name, results$plots[[plot_name]]))
  }
}

cat("\n✓ Analysis complete! Check your working directory for outputs.\n")

