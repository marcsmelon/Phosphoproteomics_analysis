################################################################################
# PHOSPHOPROTEOMICS ANALYSIS: HEALTHY vs DISEASE
# Custom script for your specific data
################################################################################

# STEP 1: INSTALL AND LOAD PACKAGES ============================================
# Run once to install (uncomment if needed):
# install.packages(c("tidyverse", "reshape2", "pheatmap", "htmltools", "DT", "base64enc", "scales", "viridis", "ggrepel", "VennDiagram"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"))

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

# STEP 3: MAIN ANALYSIS FUNCTION ================================================

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
  
  intensity_cols <- grep("(Disease_|Healthy_)", colnames(SN_processed))
  
  if(length(intensity_cols) == 0) {
    stop("No intensity columns found. Looking for patterns containing 'Disease_' or 'Healthy_'")
  }
  
  intensity_matrix <- as.matrix(SN_processed[, intensity_cols])
  intensity_matrix[intensity_matrix == "NaN"] <- NA
  storage.mode(intensity_matrix) <- "numeric"
  message(paste("  ✓ Found", ncol(intensity_matrix), "samples"))
  
  sample_names <- colnames(intensity_matrix)
  condition <- ifelse(grepl("Healthy_", sample_names), "Healthy", "Disease")
  condition <- factor(condition, levels = c("Healthy", "Disease"))
  
  message(paste("  ✓ Healthy samples:", sum(condition == "Healthy")))
  message(paste("  ✓ Disease samples:", sum(condition == "Disease")))
  
  sample_names_short <- paste0(condition, "_", ave(as.numeric(condition), condition, FUN = seq_along))
  names(sample_names_short) <- sample_names
  
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
  
  condition_colors <- c("Healthy" = "#2E86AB", "Disease" = "#E63946")
  plots <- list()
  
  message("\nGenerating QC graphs...")
  
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
  
  write.csv(ptm_site_stats, "ptm_site_statistics.csv", row.names = FALSE)
  message("  ✓ PTM site statistics saved")
  
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
    
    message("  → Pearson correlation analysis...")
    correlation_matrix <- cor(normalized_matrix, method = "pearson", use = "pairwise.complete.obs")
    
    upper_tri_values <- correlation_matrix[upper.tri(correlation_matrix)]
    average_R <- mean(upper_tri_values, na.rm = TRUE)
    message(paste("    Average Pearson correlation:", round(average_R, 3)))
    
    rownames(correlation_matrix) <- sample_names_short[rownames(correlation_matrix)]
    colnames(correlation_matrix) <- sample_names_short[colnames(correlation_matrix)]
    
    write.csv(correlation_matrix, "pearson_correlation_matrix.csv", row.names = TRUE)
    
    melted_cor <- reshape2::melt(correlation_matrix)
    
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
    
    message("  → Creating Venn diagram...")
    
    healthy_cols <- which(condition == "Healthy")
    disease_cols <- which(condition == "Disease")
    
    min_samples_healthy <- ceiling(length(healthy_cols) * 0.5)
    min_samples_disease <- ceiling(length(disease_cols) * 0.5)
    
    healthy_detected <- rownames(intensity_matrix)[
      rowSums(!is.na(intensity_matrix[, healthy_cols]) & 
                intensity_matrix[, healthy_cols] > 0) >= min_samples_healthy
    ]
    
    disease_detected <- rownames(intensity_matrix)[
      rowSums(!is.na(intensity_matrix[, disease_cols]) & 
                intensity_matrix[, disease_cols] > 0) >= min_samples_disease
    ]
    
    healthy_only <- setdiff(healthy_detected, disease_detected)
    disease_only <- setdiff(disease_detected, healthy_detected)
    both <- intersect(healthy_detected, disease_detected)
    
    message(paste("    ✓ Healthy-specific (≥50% samples):", length(healthy_only)))
    message(paste("    ✓ Both conditions (≥50% samples):", length(both)))
    message(paste("    ✓ Disease-specific (≥50% samples):", length(disease_only)))
    
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
    
    venn.plot <- VennDiagram::venn.diagram(
      x = list(
        Healthy = healthy_detected,
        Disease = disease_detected
      ),
      filename = NULL,
      category.names = c("Healthy", "Disease"),
      output = TRUE,
      imagetype = "png",
      height = 2000,
      width = 2000,
      resolution = 300,
      compression = "lzw",
      lwd = 3,
      lty = 'solid',
      fill = c(condition_colors["Healthy"], condition_colors["Disease"]),
      alpha = 0.6,
      cex = 2,
      fontface = "bold",
      fontfamily = "sans",
      cat.cex = 2,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = c(condition_colors["Healthy"], condition_colors["Disease"]),
      euler.d = TRUE,
      scaled = TRUE
    )
    
    png("venn_diagram_phosphosites.png", width = 2000, height = 2000, res = 300, bg = "white")
    grid::grid.draw(venn.plot)
    dev.off()
    plots$venn_diagram <- "venn_diagram_phosphosites.png"
  }
  
  message("\nStep 8/10: Setting up experimental design...")
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
  
  normalized_export <- as.data.frame(normalized_matrix)
  normalized_export$GenePhos <- rownames(normalized_export)
  normalized_export <- normalized_export[, c("GenePhos", colnames(normalized_matrix))]
  write.table(normalized_export, "normalized_phospho_data.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.csv(sample_mapping, "sample_name_mapping.csv", row.names = FALSE)
  
  message(paste("  ✓ Complete statistics saved to 'complete_phosphosite_statistics.csv'"))
  message(paste("  ✓ Normalized data saved to 'normalized_phospho_data.tsv'"))
  message(paste("  ✓ Sample mapping saved to 'sample_name_mapping.csv'"))
  
  if(create_plots) {
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
    
    # ========================================================================
    # ENHANCED MDS PLOTS WITH SAMPLE LABELS FOR OUTLIER IDENTIFICATION
    # ========================================================================
    message("  → Creating MDS plots with sample labels...")
    
    plot_colors <- condition_colors[as.character(condition)]
    
    # Calculate MDS coordinates
    mds_result <- limma::plotMDS(normalized_matrix, plot = FALSE)
    
    # VERSION 1: MDS Plot with Short Sample Names (Cleaner)
    png("mds_plot_condition.png", width = 1000, height = 800)
    par(mar = c(5, 5, 4, 5) + 0.1)
    
    plot(mds_result$x, mds_result$y,
         col = plot_colors,
         pch = 16, cex = 2,
         xlab = paste("Leading logFC dim 1 (", round(mds_result$var.explained[1] * 100, 1), "%)", sep = ""),
         ylab = paste("Leading logFC dim 2 (", round(mds_result$var.explained[2] * 100, 1), "%)", sep = ""),
         main = "MDS Plot - Healthy vs Disease")
    
    # Add sample labels using short names
    text(mds_result$x, mds_result$y,
         labels = sample_names_short[colnames(normalized_matrix)],
         pos = 3, cex = 0.7, col = "black")
    
    legend("topright", legend = levels(condition), 
           col = condition_colors,
           pch = 16, bg = "white", cex = 1.2)
    
    grid(col = "gray80", lty = "dotted")
    dev.off()
    plots$mds_condition <- "mds_plot_condition.png"
    
    # VERSION 2: MDS Plot with Full Sample Names
    png("mds_plot_full_labels.png", width = 1200, height = 900)
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    plot(mds_result$x, mds_result$y,
         col = plot_colors,
         pch = 16, cex = 2,
         xlab = paste("Leading logFC dim 1 (", round(mds_result$var.explained[1] * 100, 1), "%)", sep = ""),
         ylab = paste("Leading logFC dim 2 (", round(mds_result$var.explained[2] * 100, 1), "%)", sep = ""),
         main = "MDS Plot - Full Sample Names")
    
    text(mds_result$x, mds_result$y,
         labels = colnames(normalized_matrix),
         pos = 3, cex = 0.6, col = "black")
    
    legend("topright", legend = levels(condition), 
           col = condition_colors,
           pch = 16, bg = "white", cex = 1.2)
    
    grid(col = "gray80", lty = "dotted")
    dev.off()
    plots$mds_full_labels <- "mds_plot_full_labels.png"
    
    # Create MDS coordinates data frame for outlier analysis
    mds_coordinates <- data.frame(
      Sample_Full = colnames(normalized_matrix),
      Sample_Short = sample_names_short[colnames(normalized_matrix)],
      Condition = condition,
      Dim1 = mds_result$x,
      Dim2 = mds_result$y,
      stringsAsFactors = FALSE
    )
    
    # Calculate group centroids
    healthy_centroid <- c(
      mean(mds_result$x[condition == "Healthy"]),
      mean(mds_result$y[condition == "Healthy"])
    )
    
    disease_centroid <- c(
      mean(mds_result$x[condition == "Disease"]),
      mean(mds_result$y[condition == "Disease"])
    )
    
    # Calculate distance from each sample to its group centroid
    mds_coordinates$Distance_to_Centroid <- sapply(1:nrow(mds_coordinates), function(i) {
      sample_coords <- c(mds_coordinates$Dim1[i], mds_coordinates$Dim2[i])
      if(mds_coordinates$Condition[i] == "Healthy") {
        sqrt(sum((sample_coords - healthy_centroid)^2))
      } else {
        sqrt(sum((sample_coords - disease_centroid)^2))
      }
    })
    
    # Identify potential outliers (samples > 2 SD from centroid)
    mean_dist <- mean(mds_coordinates$Distance_to_Centroid)
    sd_dist <- sd(mds_coordinates$Distance_to_Centroid)
    mds_coordinates$Potential_Outlier <- mds_coordinates$Distance_to_Centroid > (mean_dist + 2*sd_dist)
    
    # Save coordinates
    write.csv(mds_coordinates, "mds_sample_coordinates.csv", row.names = FALSE)
    
    # Print outlier information
    cat("\n")
    cat("================================================================================\n")
    cat("MDS OUTLIER ANALYSIS\n")
    cat("================================================================================\n\n")
    cat("MDS Sample Coordinates saved to: mds_sample_coordinates.csv\n\n")
    
    if(any(mds_coordinates$Potential_Outlier)) {
      cat("Potential outliers detected (>2 SD from group centroid):\n")
      outlier_info <- mds_coordinates[mds_coordinates$Potential_Outlier == TRUE, 
                                      c("Sample_Short", "Sample_Full", "Condition", "Distance_to_Centroid")]
      print(outlier_info)
      cat("\n")
    } else {
      cat("No obvious outliers detected (all samples within 2 SD of group centroid)\n\n")
    }
    
    message("  ✓ MDS plots with labels created")
  }
  
  if(create_plots) {
    message("\nCreating visualizations...")
    
    message("  → Creating enhanced volcano plot...")
    
    volcano_data <- tt %>%
      mutate(
        Expression = case_when(
          logFC >= log2(2) & adj.P.Val <= p_value_threshold ~ "Up-regulated in Disease",
          logFC <= -log2(2) & adj.P.Val <= p_value_threshold ~ "Up-regulated in Healthy",
          TRUE ~ "Unchanged"
        )
      )
    
    data_unchanged <- volcano_data %>% filter(Expression == "Unchanged")
    data_updown <- volcano_data %>% filter(Expression != "Unchanged")
    
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
    
    p_volcano <- ggplot2::ggplot() +
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
      ggplot2::geom_hline(yintercept = -log10(p_value_threshold), 
                          linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)), 
                          linetype = "dashed", color = "gray30", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
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
        panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
      ) +
      ggplot2::labs(fill = "Regulation")
    
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
    
    ggplot2::ggsave("volcano_plot.pdf", p_volcano, width = 10, height = 8)
    
    write.csv(top_genes, "top_genes_volcano_plot.csv", row.names = FALSE)
    message("    ✓ Top genes saved to 'top_genes_volcano_plot.csv'")
    
    message("  ✓ Visualizations saved")
  } else {
    top_genes <- NULL
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
    cat("  • mds_plot_full_labels.png\n")
    cat("  • mds_sample_coordinates.csv\n")
  }
  cat("\n")
  
  return(list(
    raw_data = intensity_matrix,
    log2_data = log2_intensity_matrix,
    normalized_data = normalized_matrix,
    correlation_matrix = if(create_plots) correlation_matrix else NULL,
    average_correlation = if(create_plots) average_R else NULL,
    venn_stats = if(create_plots) venn_stats else NULL,
    mds_coordinates = if(create_plots) mds_coordinates else NULL,
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

# STEP 3B: ENRICHMENT ANALYSIS FUNCTION =========================================

run_enrichment_analysis <- function(results, 
                                    p_value_threshold = 0.1,
                                    fc_threshold = 1,
                                    organism = "org.Hs.eg.db",
                                    run_go = TRUE,
                                    run_kegg = TRUE,
                                    run_reactome = TRUE) {
  
  required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if(length(missing_packages) > 0) {
    message("ERROR: Required packages not installed: ", paste(missing_packages, collapse = ", "))
    message("Install with: BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))")
    return(NULL)
  }
  
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
  })
  
  if(run_reactome) {
    if(!requireNamespace("ReactomePA", quietly = TRUE)) {
      message("Note: ReactomePA not installed. Install with: BiocManager::install('ReactomePA')")
      run_reactome <- FALSE
    } else {
      library(ReactomePA)
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    PATHWAY ENRICHMENT ANALYSIS\n")
  cat("================================================================================\n\n")
  
  message("Step 1/6: Extracting genes from phosphosites...")
  
  extract_gene <- function(genephos) {
    gene <- gsub("_[STY].*$", "", genephos)
    return(gene)
  }
  
  all_genes <- unique(sapply(rownames(results$normalized_data), extract_gene))
  sig_genes <- unique(sapply(results$significant_results$GenePhos, extract_gene))
  
  up_genes <- unique(sapply(
    results$significant_results$GenePhos[results$significant_results$logFC > fc_threshold],
    extract_gene
  ))
  
  down_genes <- unique(sapply(
    results$significant_results$GenePhos[results$significant_results$logFC < -fc_threshold],
    extract_gene
  ))
  
  message(paste("  ✓ Total genes detected:", length(all_genes)))
  message(paste("  ✓ Significant genes:", length(sig_genes)))
  message(paste("  ✓ Up-regulated genes:", length(up_genes)))
  message(paste("  ✓ Down-regulated genes:", length(down_genes)))
  
  message("\nStep 2/6: Converting gene symbols to Entrez IDs...")
  
  gene2entrez <- function(genes) {
    if(length(genes) == 0) {
      return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
    }
    
    entrez <- tryCatch({
      clusterProfiler::bitr(genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = organism)
    }, error = function(e) {
      message("    Warning: Error converting gene symbols: ", e$message)
      data.frame(SYMBOL = character(0), ENTREZID = character(0))
    })
    
    if(is.null(entrez) || nrow(entrez) == 0) {
      return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
    }
    
    return(entrez)
  }
  
  all_entrez <- gene2entrez(all_genes)
  sig_entrez <- gene2entrez(sig_genes)
  
  if(length(up_genes) > 0) {
    up_entrez <- gene2entrez(up_genes)
  } else {
    up_entrez <- data.frame(SYMBOL = character(0), ENTREZID = character(0))
  }
  
  if(length(down_genes) > 0) {
    down_entrez <- gene2entrez(down_genes)
  } else {
    down_entrez <- data.frame(SYMBOL = character(0), ENTREZID = character(0))
  }
  
  message(paste("  ✓ Converted", nrow(all_entrez), "/", length(all_genes), "background genes"))
  message(paste("  ✓ Converted", nrow(sig_entrez), "/", length(sig_genes), "significant genes"))
  message(paste("  ✓ Converted", nrow(up_entrez), "/", length(up_genes), "up-regulated genes"))
  message(paste("  ✓ Converted", nrow(down_entrez), "/", length(down_genes), "down-regulated genes"))
  
  if(nrow(all_entrez) < 100) {
    message("\nWARNING: Very few background genes converted. Enrichment may not be reliable.")
  }
  
  if(nrow(sig_entrez) < 5) {
    message("\nERROR: Too few significant genes converted to Entrez IDs.")
    message("Cannot proceed with enrichment analysis.")
    return(NULL)
  }
  
  enrichment_results <- list()
  
  if(run_go && nrow(sig_entrez) >= 5) {
    message("\nStep 3/6: Running GO enrichment...")
    
    message("  → GO Biological Process (all significant)...")
    go_bp <- tryCatch({
      clusterProfiler::enrichGO(
        gene = sig_entrez$ENTREZID,
        universe = all_entrez$ENTREZID,
        OrgDb = organism,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) {
      message("    Error in GO enrichment: ", e$message)
      NULL
    })
    
    if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
      write.csv(go_bp@result, "enrichment_GO_BP_all.csv", row.names = FALSE)
      enrichment_results$GO_BP_all <- go_bp
      
      if(nrow(go_bp@result) >= 10) {
        p <- enrichplot::dotplot(go_bp, showCategory = 20) + 
          ggtitle("GO Biological Process - All Significant Genes")
        ggsave("enrichment_GO_BP_all_dotplot.png", p, width = 10, height = 8, dpi = 300)
      }
      message(paste("    ✓ Found", nrow(go_bp@result), "enriched BP terms"))
    } else {
      message("    No significant GO BP terms found")
    }
    
    if(nrow(up_entrez) >= 5) {
      message("  → GO Biological Process (up-regulated)...")
      go_bp_up <- tryCatch({
        clusterProfiler::enrichGO(
          gene = up_entrez$ENTREZID,
          universe = all_entrez$ENTREZID,
          OrgDb = organism,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
      }, error = function(e) {
        message("    Error in GO enrichment: ", e$message)
        NULL
      })
      
      if(!is.null(go_bp_up) && nrow(go_bp_up@result) > 0) {
        write.csv(go_bp_up@result, "enrichment_GO_BP_upregulated.csv", row.names = FALSE)
        enrichment_results$GO_BP_up <- go_bp_up
        
        if(nrow(go_bp_up@result) >= 5) {
          p <- enrichplot::dotplot(go_bp_up, showCategory = 15) + 
            ggtitle("GO BP - Up-regulated in Disease")
          ggsave("enrichment_GO_BP_upregulated_dotplot.png", p, width = 10, height = 8, dpi = 300)
        }
        message(paste("    ✓ Found", nrow(go_bp_up@result), "enriched terms"))
      } else {
        message("    No significant GO BP terms found for up-regulated genes")
      }
    }
    
    if(nrow(down_entrez) >= 5) {
      message("  → GO Biological Process (down-regulated)...")
      go_bp_down <- tryCatch({
        clusterProfiler::enrichGO(
          gene = down_entrez$ENTREZID,
          universe = all_entrez$ENTREZID,
          OrgDb = organism,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
      }, error = function(e) {
        message("    Error in GO enrichment: ", e$message)
        NULL
      })
      
      if(!is.null(go_bp_down) && nrow(go_bp_down@result) > 0) {
        write.csv(go_bp_down@result, "enrichment_GO_BP_downregulated.csv", row.names = FALSE)
        enrichment_results$GO_BP_down <- go_bp_down
        
        if(nrow(go_bp_down@result) >= 5) {
          p <- enrichplot::dotplot(go_bp_down, showCategory = 15) + 
            ggtitle("GO BP - Up-regulated in Healthy")
          ggsave("enrichment_GO_BP_downregulated_dotplot.png", p, width = 10, height = 8, dpi = 300)
        }
        message(paste("    ✓ Found", nrow(go_bp_down@result), "enriched terms"))
      } else {
        message("    No significant GO BP terms found for down-regulated genes")
      }
    }
  }
  
  if(run_kegg && nrow(sig_entrez) >= 5) {
    message("\nStep 4/6: Running KEGG pathway enrichment...")
    
    kegg_all <- tryCatch({
      clusterProfiler::enrichKEGG(
        gene = sig_entrez$ENTREZID,
        universe = all_entrez$ENTREZID,
        organism = 'hsa',
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("    Error in KEGG enrichment: ", e$message)
      NULL
    })
    
    if(!is.null(kegg_all) && nrow(kegg_all@result) > 0) {
      write.csv(kegg_all@result, "enrichment_KEGG_all.csv", row.names = FALSE)
      enrichment_results$KEGG_all <- kegg_all
      
      if(nrow(kegg_all@result) >= 5) {
        p <- enrichplot::dotplot(kegg_all, showCategory = 15) + 
          ggtitle("KEGG Pathways - All Significant Genes")
        ggsave("enrichment_KEGG_all_dotplot.png", p, width = 10, height = 8, dpi = 300)
      }
      message(paste("    ✓ Found", nrow(kegg_all@result), "enriched KEGG pathways"))
    } else {
      message("    No significant KEGG pathways found")
    }
  }
  
  if(run_reactome && nrow(sig_entrez) >= 5) {
    message("\nStep 5/6: Running Reactome pathway enrichment...")
    
    reactome_all <- tryCatch({
      ReactomePA::enrichPathway(
        gene = sig_entrez$ENTREZID,
        universe = all_entrez$ENTREZID,
        organism = "human",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) {
      message("    Error in Reactome enrichment: ", e$message)
      NULL
    })
    
    if(!is.null(reactome_all) && nrow(reactome_all@result) > 0) {
      write.csv(reactome_all@result, "enrichment_Reactome_all.csv", row.names = FALSE)
      enrichment_results$Reactome_all <- reactome_all
      
      if(nrow(reactome_all@result) >= 5) {
        p <- enrichplot::dotplot(reactome_all, showCategory = 15) + 
          ggtitle("Reactome Pathways - All Significant Genes")
        ggsave("enrichment_Reactome_all_dotplot.png", p, width = 10, height = 8, dpi = 300)
      }
      message(paste("    ✓ Found", nrow(reactome_all@result), "enriched Reactome pathways"))
    } else {
      message("    No significant Reactome pathways found")
    }
  }
  
  message("\nStep 6/6: Enrichment analysis complete!")
  
  cat("\n")
  cat("================================================================================\n")
  cat("                   ENRICHMENT ANALYSIS COMPLETE!\n")
  cat("================================================================================\n\n")
  
  if(length(enrichment_results) > 0) {
    cat("OUTPUT FILES:\n")
    for(name in names(enrichment_results)) {
      result_name <- gsub("_", " ", name)
      cat(sprintf("  • enrichment_%s.csv\n", name))
      if(nrow(enrichment_results[[name]]@result) >= 5) {
        cat(sprintf("  • enrichment_%s_dotplot.png\n", name))
      }
    }
  } else {
    cat("No enrichment results generated.\n")
  }
  
  cat("\n")
  
  return(enrichment_results)
}

# STEP 4: RUN THE ANALYSIS ======================================================
# Set your working directory
setwd("/path/to/your/data")  # <-- CHANGE THIS to your data folder!

# Run the analysis
results <- phospho_analysis_healthy_disease(
  phospho_data = "Phospho(STY).tsv",     # <-- CHANGE THIS to your file name!
  use_advanced_preprocessing = TRUE,      # Use sophisticated preprocessing
  filtering_threshold = 3,                # Present in at least 3 samples
  imputation_method = "min",              # Half minimum imputation
  normalization_method = "median_center", # Median centering
  p_value_threshold = 0.1,                # FDR < 0.1
  create_plots = TRUE                     # Generate plots
)

# Optional: Run enrichment analysis (requires Bioconductor packages)
enrichment <- run_enrichment_analysis(
   results = results,
   p_value_threshold = 0.1,
   fc_threshold = 1,
   run_go = TRUE,
   run_kegg = TRUE,
   run_reactome = TRUE
 )

