################################################################################
# HUMAN DCM PROTEOMICS BENCHMARKING — v3
#
# Based on the processing approach from Proteomics_pipeline.R, adapted for
# protein-level data with the same 8-comparison benchmark design as the
# phosphoproteomics benchmark (Human_DCM_Phosphoproteomics_Benchmark_v3).
#
# Per-comparison outputs:
#   • QC plots: boxplots (before/after norm), density, correlation heatmap
#   • PCA/MDS with sample labels + outlier detection
#   • limma DE (Disease – Healthy contrast)
#   • Volcano plot
#   • clusterProfiler enrichment: GO BP (all/up/down), KEGG, Reactome + dotplots
#   • fgsea GSEA: Hallmarks, KEGG, Reactome (ranked by t-statistic)
#
# Input:  fullCohortProt_PivotTable.tsv
#         (columns: rowname, PG.GroupLabel [UniProt], 66 intensity columns)
#
# Contrast direction: Disease – Healthy
#   positive logFC = higher protein abundance in DCM
#   negative logFC = lower protein abundance in DCM
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PROTEOMICS_DATA_PATH <- "fullCohortProt_PivotTable.tsv"
OUTPUT_BASE_DIR      <- "Proteomics_Benchmark"

# Preprocessing parameters (matching Proteomics_pipeline.R defaults)
FILTERING_MIN_SAMPLES    <- 3       # min non-NA samples per protein
IMPUTATION_METHOD        <- "min"   # "min" = row min / 2
NORMALIZATION_METHOD     <- "median_center"

# Enrichment parameters
ENRICHMENT_P_THRESHOLD   <- 0.05    # FDR cutoff for selecting DE genes
ENRICHMENT_FC_THRESHOLD  <- 0       # logFC cutoff (0 = use all significant)

# ==============================================================================
# COLOUR SCHEME — applied to all plots
# ==============================================================================

COL_DCM     <- "#F5CD6A"   # Up-regulated in DCM (gold)
COL_HEALTHY <- "#3A2044"   # Down-regulated / Healthy (dark purple)
COL_NS      <- "grey75"    # Non-significant
COL_FDR     <- "#D4A03C"   # FDR < 0.05 but below FC cutoff (muted gold)

# Named vectors for ggplot scale_*_manual
CONDITION_COLORS <- c("Healthy" = COL_HEALTHY, "Disease" = COL_DCM)
DIRECTION_COLORS <- c("Up in DCM" = COL_DCM, "Down in DCM" = COL_HEALTHY)

# Benchmark barplot palette
BENCH_COLORS <- c("FDR < 0.01" = COL_HEALTHY,
                  "FDR < 0.05" = COL_DCM,
                  "FDR < 0.10" = "#E8D5A0")  # light gold

# Transparent background theme (appended to all ggplots)
THEME_TRANSPARENT <- theme(
  panel.background  = element_rect(fill = "transparent", color = NA),
  plot.background   = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA)
)

# ==============================================================================
# BLOOD MARKER GENES — excluded from heatmaps to avoid contamination artifacts
# ==============================================================================

BLOOD_MARKER_GENES <- c(
  # Hemoglobin
  "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBZ",
  # RBC enzymes
  "CA1", "CA2", "BPGM", "PRDX2", "BLVRB", "EPB41", "EPB42",
  "SLC4A1", "ANK1", "SPTA1", "SPTB",
  # Albumin + acute phase
  "ALB", "HP", "HPR", "ORM1", "ORM2", "CRP", "SAA1", "SAA2",
  "SERPINA1", "SERPINA3", "A2M",
  # Transport
  "TF", "TTR", "CP", "GC", "APOD", "RBP4", "SHBG",
  # Coagulation
  "FGA", "FGB", "FGG", "F2", "F5", "F7", "F9", "F10", "F12",
  "PLG", "SERPINC1", "SERPIND1", "PROC", "PROS1",
  # Apolipoproteins
  "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2",
  "APOC3", "APOE", "APOH", "APOL1", "APOM",
  # Complement
  "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C3", "C4A",
  "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP",
  # Immunoglobulins
  "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGHA1",
  "IGHA2", "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN",
  # Platelets
  "PF4", "PPBP", "GP1BA", "GP1BB", "GP5", "GP9", "ITGA2B",
  "ITGB3", "SELP", "THBS1", "VWF"
)

# ==============================================================================
# PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(reshape2)
  library(scales)
  library(VennDiagram)
})

# Bioconductor — loaded on demand in functions
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# library(clusterProfiler)
# library(enrichplot)
# library(ReactomePA)
# library(fgsea)
# library(msigdbr)

# ==============================================================================
# SAMPLE DEFINITIONS — identical to phosphoproteomics benchmark
# ==============================================================================

# ---- Healthy donors (14 columns) ----
healthy_cols <- c(
  "166-D", "246-D", "263-D", "270-D", "277-D", "285-D", "345-D",
  "382-D", "391-D", "450-D", "460-D", "480-D", "493-D", "508-D"
)

# ---- DCM_filtered from metadata (13 columns — curated set) ----
dcm_filtered_cols <- c(
  "120-V",  "0091-V", "63-T",   "272-V",  "309-V",  "328-V",
  "0062-V", "0021-V", "288-V",  "0019-V", "468-V",  "499-V",
  "545-T"
)

# ---- Myocarditis columns (excluded from DCM-only comparisons) ----
myocarditis_cols <- c(
  "0020-V", "0039-V", "0053-T", "0054-V", "0077-V", "72-V",
  "79-V",   "395-V"
)

# ---- Unknown / no metadata match ----
unknown_cols <- c("0058-T", "0096-V")

# ---- All disease columns (28 VAD + 24 Tx = 52) ----
all_vad_cols <- c(
  "0019-V", "0020-V", "0021-V", "0034-V", "0039-V", "0050-V",
  "0054-V", "0062-V", "0077-V", "0091-V", "0096-V", "104-V",
  "120-V",  "148-V",  "154-V",  "172-V",  "180-V",  "272-V",
  "288-V",  "309-V",  "328-V",  "395-V",  "426-V",  "468-V",
  "499-V",  "53-V",   "72-V",   "79-V"
)

all_tx_cols <- c(
  "0025-T", "0029-T", "0031-T", "0050-T", "0053-T", "0058-T",
  "154-T",  "172-T",  "223-T",  "262-T",  "269-T",  "272-T",
  "288-T",  "344-T",  "367-T",  "426-T",  "449-T",  "468-T",
  "479-T",  "490-T",  "499-T",  "53-T",   "545-T",  "63-T"
)

all_disease_cols <- c(all_vad_cols, all_tx_cols)

# ---- One-per-patient de-duplication ----
build_one_per_patient <- function(prefer = "V") {
  all_cols <- c(all_vad_cols, all_tx_cols)
  get_id <- function(col) sub("-(V|T|D)$", "", col)
  
  col_df <- data.frame(col = all_cols,
                       id  = sapply(all_cols, get_id),
                       type = ifelse(grepl("-V$", all_cols), "V", "T"),
                       stringsAsFactors = FALSE)
  
  selected <- col_df %>%
    group_by(id) %>%
    summarise(col = {
      if (n() == 1) col
      else if (prefer == "V" && any(type == "V")) col[type == "V"][1]
      else if (prefer == "T" && any(type == "T")) col[type == "T"][1]
      else col[1]
    }, .groups = "drop") %>%
    pull(col)
  
  return(selected)
}

one_per_patient_preferV <- build_one_per_patient("V")
one_per_patient_preferT <- build_one_per_patient("T")

dcm_one_per_patient_V <- setdiff(one_per_patient_preferV,
                                 c(myocarditis_cols, unknown_cols))
dcm_one_per_patient_T <- setdiff(one_per_patient_preferT,
                                 c(myocarditis_cols, unknown_cols))

all_dcm_naive <- setdiff(all_disease_cols, c(myocarditis_cols, unknown_cols))

# ---- Genetic subgroups ----
col_group <- c(
  "120-V"  = "Sarcomeric", "309-V"  = "Sarcomeric", "63-T"   = "Sarcomeric",
  "154-V"  = "Sarcomeric", "154-T"  = "Sarcomeric",
  "367-T"  = "Sarcomeric", "449-T"  = "Sarcomeric", "490-T"  = "Sarcomeric",
  "0019-V" = "Cytoskeletal", "0021-V" = "Cytoskeletal", "0034-V" = "Cytoskeletal",
  "0050-V" = "Cytoskeletal", "0050-T" = "Cytoskeletal",
  "172-V"  = "Cytoskeletal", "172-T"  = "Cytoskeletal",
  "223-T"  = "Cytoskeletal", "269-T"  = "Cytoskeletal", "328-V"  = "Cytoskeletal",
  "344-T"  = "Cytoskeletal", "479-T"  = "Cytoskeletal", "499-V"  = "Cytoskeletal",
  "499-T"  = "Cytoskeletal", "545-T"  = "Cytoskeletal",
  "0020-V" = "Idiopathic", "0062-V" = "Idiopathic", "0091-V" = "Idiopathic",
  "148-V"  = "Idiopathic", "180-V"  = "Idiopathic", "262-T"  = "Idiopathic",
  "272-V"  = "Idiopathic", "272-T"  = "Idiopathic", "288-V"  = "Idiopathic",
  "288-T"  = "Idiopathic", "426-V"  = "Idiopathic", "426-T"  = "Idiopathic",
  "468-V"  = "Idiopathic", "468-T"  = "Idiopathic",
  "0025-T" = "Idiopathic", "0029-T" = "Idiopathic", "0031-T" = "Idiopathic",
  "104-V"  = "Idiopathic", "53-V"   = "Idiopathic", "53-T"   = "Idiopathic"
)

sarcomeric_1pp   <- intersect(dcm_one_per_patient_V,
                              names(col_group[col_group == "Sarcomeric"]))
cytoskeletal_1pp <- intersect(dcm_one_per_patient_V,
                              names(col_group[col_group == "Cytoskeletal"]))
idiopathic_1pp   <- intersect(dcm_one_per_patient_V,
                              names(col_group[col_group == "Idiopathic"]))

# ==============================================================================
# DEFINE 8 COMPARISONS
# ==============================================================================

comparisons <- list(
  C1 = list(d = dcm_filtered_cols, h = healthy_cols,
            label = "C1: DCM_filtered curated (13 vs 14)",
            dir   = "C1_DCM_filtered"),
  C2 = list(d = dcm_one_per_patient_V, h = healthy_cols,
            label = "C2: All DCM 1-per-patient prefer-V (35 vs 14)",
            dir   = "C2_AllDCM_1pp_preferV"),
  C3 = list(d = dcm_one_per_patient_T, h = healthy_cols,
            label = "C3: All DCM 1-per-patient prefer-T (35 vs 14)",
            dir   = "C3_AllDCM_1pp_preferT"),
  C4 = list(d = all_dcm_naive, h = healthy_cols,
            label = "C4: All DCM naive incl. duplicates (44 vs 14) [REFERENCE]",
            dir   = "C4_AllDCM_naive"),
  C5 = list(d = sarcomeric_1pp, h = healthy_cols,
            label = "C5: Sarcomeric mutations (vs 14 Healthy)",
            dir   = "C5_Sarcomeric"),
  C6 = list(d = cytoskeletal_1pp, h = healthy_cols,
            label = "C6: Cytoskeletal/structural (vs 14 Healthy)",
            dir   = "C6_Cytoskeletal"),
  C7 = list(d = idiopathic_1pp, h = healthy_cols,
            label = "C7: Idiopathic DCM (vs 14 Healthy)",
            dir   = "C7_Idiopathic"),
  C8 = list(d = setdiff(one_per_patient_preferV, unknown_cols), h = healthy_cols,
            label = "C8: All disease incl. myocarditis 1pp (43 vs 14)",
            dir   = "C8_AllDisease_1pp")
)

# ==============================================================================
# UNIPROT → GENE SYMBOL MAPPING
# ==============================================================================

map_uniprot_to_gene <- function(uniprot_ids) {
  tryCatch({
    suppressPackageStartupMessages({
      library(org.Hs.eg.db)
      library(AnnotationDbi)
    })
    
    primary_ids <- sapply(strsplit(uniprot_ids, ";"), function(x) trimws(x[1]))
    
    gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys    = unique(primary_ids),
                                      keytype = "UNIPROT",
                                      columns = "SYMBOL")
    gene_map <- gene_map[!duplicated(gene_map$UNIPROT), ]
    
    symbol_vec <- setNames(gene_map$SYMBOL, gene_map$UNIPROT)
    result     <- symbol_vec[primary_ids]
    names(result) <- uniprot_ids
    
    n_mapped <- sum(!is.na(result))
    message(sprintf("  ✓ Mapped %d / %d UniProt IDs to gene symbols (%.1f%%)",
                    n_mapped, length(uniprot_ids),
                    100 * n_mapped / length(uniprot_ids)))
    return(result)
    
  }, error = function(e) {
    message("  ⚠ org.Hs.eg.db mapping failed: ", e$message)
    message("  → Returning UniProt IDs as-is (install org.Hs.eg.db for gene symbols)")
    result <- setNames(uniprot_ids, uniprot_ids)
    return(result)
  })
}

# ==============================================================================
# QC PLOTS FUNCTION (per comparison)
# Following Proteomics_pipeline.R: boxplots, density, correlation, Venn
# ==============================================================================

generate_qc_plots <- function(log2_mat, norm_mat, condition, sample_labels,
                              out_path, label) {
  
  condition_colors <- CONDITION_COLORS
  
  # ---- Boxplot BEFORE normalization ----
  tryCatch({
    tmp_before <- reshape2::melt(t(log2_mat))
    colnames(tmp_before) <- c("Sample", "Protein", "Intensity")
    tmp_before$Condition <- condition[match(tmp_before$Sample, colnames(log2_mat))]
    tmp_before$Label <- sample_labels[as.character(tmp_before$Sample)]
    tmp_before <- tmp_before[order(tmp_before$Condition), ]
    tmp_before$Label <- factor(tmp_before$Label, levels = unique(tmp_before$Label))
    
    p <- ggplot(tmp_before, aes(x = Label, y = Intensity, fill = Condition)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = condition_colors) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
      labs(title = paste0("Log2 Intensities Before Normalization\n", label),
           x = "Sample", y = "Log2 Intensity")
    
    ggsave(file.path(out_path, "boxplot_before_normalization.png"), p,
           width = max(8, length(sample_labels) * 0.25), height = 6, dpi = 300, bg = "transparent")
    message("    ✓ Boxplot (before norm)")
  }, error = function(e) message("    ⚠ Boxplot before: ", e$message))
  
  # ---- Boxplot AFTER normalization ----
  tryCatch({
    tmp_after <- reshape2::melt(t(norm_mat))
    colnames(tmp_after) <- c("Sample", "Protein", "Intensity")
    tmp_after$Condition <- condition[match(tmp_after$Sample, colnames(norm_mat))]
    tmp_after$Label <- sample_labels[as.character(tmp_after$Sample)]
    tmp_after <- tmp_after[order(tmp_after$Condition), ]
    tmp_after$Label <- factor(tmp_after$Label, levels = unique(tmp_after$Label))
    
    p <- ggplot(tmp_after, aes(x = Label, y = Intensity, fill = Condition)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = condition_colors) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
      labs(title = paste0("Normalized Intensities (Median Centered)\n", label),
           x = "Sample", y = "Normalized Intensity")
    
    ggsave(file.path(out_path, "boxplot_after_normalization.png"), p,
           width = max(8, length(sample_labels) * 0.25), height = 6, dpi = 300, bg = "transparent")
    message("    ✓ Boxplot (after norm)")
  }, error = function(e) message("    ⚠ Boxplot after: ", e$message))
  
  # ---- Density plot ----
  tryCatch({
    png(file.path(out_path, "density_plot_normalized.png"), width = 800, height = 600)
    par(mar = c(5, 5, 4, 5) + 0.1, bty = "L")
    limma::plotDensities(norm_mat,
                         group = condition,
                         col = condition_colors[as.character(condition)],
                         legend = "topright",
                         main = paste0("Density Plot — ", label))
    abline(v = 0, lty = 2, col = "gray50")
    dev.off()
    message("    ✓ Density plot")
  }, error = function(e) message("    ⚠ Density plot: ", e$message))
  
  # ---- Pearson correlation heatmap ----
  tryCatch({
    cor_mat <- cor(norm_mat, method = "pearson", use = "pairwise.complete.obs")
    avg_r   <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
    message(sprintf("    Average Pearson r: %.3f", avg_r))
    
    rownames(cor_mat) <- sample_labels[rownames(cor_mat)]
    colnames(cor_mat) <- sample_labels[colnames(cor_mat)]
    
    write.csv(cor_mat, file.path(out_path, "pearson_correlation_matrix.csv"), row.names = TRUE)
    
    melted <- reshape2::melt(cor_mat)
    p <- ggplot(melted, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "C", name = "Pearson\nCorrelation",
                           limits = c(min(0.5, min(melted$value, na.rm = TRUE)), 1)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
      labs(title = paste0("Pearson Correlation — ", label), x = "", y = "") +
      geom_text(aes(label = round(value, 2)), size = 1.8, color = "white")
    
    ggsave(file.path(out_path, "pearson_correlation_heatmap.png"), p,
           width = max(8, ncol(norm_mat) * 0.35), height = max(7, ncol(norm_mat) * 0.33),
           dpi = 300, bg = "transparent")
    message("    ✓ Correlation heatmap")
    
    return(avg_r)
  }, error = function(e) {
    message("    ⚠ Correlation: ", e$message)
    return(NA)
  })
}

# ==============================================================================
# PCA / MDS FUNCTION (per comparison)
# Following Proteomics_pipeline.R: MDS + outlier detection
# ==============================================================================

generate_pca_mds <- function(norm_mat, condition, sample_labels,
                             out_path, label) {
  
  condition_colors <- CONDITION_COLORS
  plot_colors <- condition_colors[as.character(condition)]
  
  # ---- MDS plot (limma) ----
  tryCatch({
    mds <- limma::plotMDS(norm_mat, plot = FALSE)
    
    # MDS with sample labels
    png(file.path(out_path, "mds_plot.png"), width = 1000, height = 800)
    par(mar = c(5, 5, 4, 5) + 0.1)
    
    plot(mds$x, mds$y,
         col = plot_colors, pch = 16, cex = 2,
         xlab = paste0("Leading logFC dim 1 (",
                       round(mds$var.explained[1] * 100, 1), "%)"),
         ylab = paste0("Leading logFC dim 2 (",
                       round(mds$var.explained[2] * 100, 1), "%)"),
         main = paste0("MDS Plot — ", label))
    
    text(mds$x, mds$y,
         labels = sample_labels[colnames(norm_mat)],
         pos = 3, cex = 0.6, col = "black")
    
    legend("topright", legend = levels(condition),
           col = condition_colors, pch = 16, cex = 1.2)
    grid(col = "gray80", lty = "dotted")
    dev.off()
    
    # Save MDS coordinates + outlier detection
    mds_df <- data.frame(
      Sample      = colnames(norm_mat),
      Label       = sample_labels[colnames(norm_mat)],
      Condition   = as.character(condition),
      Dim1        = mds$x,
      Dim2        = mds$y,
      stringsAsFactors = FALSE
    )
    
    # Centroid distances
    for (cond in c("Healthy", "Disease")) {
      idx <- mds_df$Condition == cond
      if (sum(idx) > 1) {
        cx <- mean(mds_df$Dim1[idx])
        cy <- mean(mds_df$Dim2[idx])
        mds_df$Dist_to_Centroid[idx] <- sqrt((mds_df$Dim1[idx] - cx)^2 +
                                               (mds_df$Dim2[idx] - cy)^2)
      }
    }
    
    mean_dist <- mean(mds_df$Dist_to_Centroid, na.rm = TRUE)
    sd_dist   <- sd(mds_df$Dist_to_Centroid, na.rm = TRUE)
    mds_df$Potential_Outlier <- mds_df$Dist_to_Centroid > (mean_dist + 2 * sd_dist)
    
    write.csv(mds_df, file.path(out_path, "mds_coordinates.csv"), row.names = FALSE)
    
    if (any(mds_df$Potential_Outlier, na.rm = TRUE)) {
      outliers <- mds_df[mds_df$Potential_Outlier == TRUE, ]
      cat(sprintf("    ⚠ Potential outliers (>2 SD from centroid): %s\n",
                  paste(outliers$Sample, collapse = ", ")))
    } else {
      cat("    ✓ No outliers detected (all within 2 SD of group centroid)\n")
    }
    
    message("    ✓ MDS plot + coordinates saved")
    
  }, error = function(e) message("    ⚠ MDS plot: ", e$message))
  
  # ---- PCA (prcomp) for additional perspective ----
  tryCatch({
    # Remove rows with any NA for PCA
    complete_rows <- complete.cases(norm_mat)
    pca_mat <- norm_mat[complete_rows, ]
    
    if (nrow(pca_mat) > 50) {
      pca <- prcomp(t(pca_mat), center = TRUE, scale. = FALSE)
      
      pca_df <- data.frame(
        Sample    = colnames(pca_mat),
        Label     = sample_labels[colnames(pca_mat)],
        Condition = as.character(condition),
        PC1       = pca$x[, 1],
        PC2       = pca$x[, 2],
        stringsAsFactors = FALSE
      )
      
      var_exp <- summary(pca)$importance[2, 1:2] * 100
      
      p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Label)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text_repel(size = 2.5, max.overlaps = 20) +
        scale_color_manual(values = CONDITION_COLORS) +
        labs(title    = paste0("PCA — ", label),
             x = sprintf("PC1 (%.1f%%)", var_exp[1]),
             y = sprintf("PC2 (%.1f%%)", var_exp[2]),
             color = NULL) +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 12),
              legend.position = "bottom") +
        THEME_TRANSPARENT
      
      ggsave(file.path(out_path, "pca_plot.png"), p_pca,
             width = 10, height = 8, dpi = 300, bg = "transparent")
      
      write.csv(pca_df, file.path(out_path, "pca_coordinates.csv"), row.names = FALSE)
      message("    ✓ PCA plot + coordinates saved")
    }
  }, error = function(e) message("    ⚠ PCA: ", e$message))
}

# ==============================================================================
# ENRICHMENT: clusterProfiler (GO BP, KEGG, Reactome) — per comparison
# Mirrors run_enrichment_analysis() from Proteomics_pipeline.R
# ==============================================================================

run_clusterprofiler_enrichment <- function(top_table, out_path, label,
                                           p_thresh = ENRICHMENT_P_THRESHOLD,
                                           fc_thresh = ENRICHMENT_FC_THRESHOLD) {
  tryCatch({
    suppressPackageStartupMessages({
      library(clusterProfiler)
      library(org.Hs.eg.db)
      library(enrichplot)
    })
    
    message("\n  ── clusterProfiler enrichment: ", label, " ──")
    
    # Background: all tested genes
    all_genes <- unique(top_table$Gene[!is.na(top_table$Gene) & top_table$Gene != ""])
    # Significant genes
    sig_tt  <- top_table[top_table$adj.P.Val < p_thresh & !is.na(top_table$Gene) &
                           top_table$Gene != "", ]
    sig_genes <- unique(sig_tt$Gene)
    up_genes  <- unique(sig_tt$Gene[sig_tt$logFC > fc_thresh])
    down_genes <- unique(sig_tt$Gene[sig_tt$logFC < -fc_thresh])
    
    message(sprintf("    Background: %d genes | Significant: %d | Up: %d | Down: %d",
                    length(all_genes), length(sig_genes), length(up_genes), length(down_genes)))
    
    # Convert to Entrez IDs
    gene2entrez <- function(genes) {
      if (length(genes) < 2) return(data.frame(SYMBOL = character(), ENTREZID = character()))
      tryCatch(
        clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db),
        error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
      )
    }
    
    all_entrez  <- gene2entrez(all_genes)
    sig_entrez  <- gene2entrez(sig_genes)
    up_entrez   <- gene2entrez(up_genes)
    down_entrez <- gene2entrez(down_genes)
    
    message(sprintf("    Entrez mapping: %d/%d background, %d/%d significant",
                    nrow(all_entrez), length(all_genes),
                    nrow(sig_entrez), length(sig_genes)))
    
    if (nrow(sig_entrez) < 5) {
      message("    ⚠ Too few genes for enrichment (<5). Skipping.")
      return(NULL)
    }
    
    results_list <- list()
    
    # ---- GO Biological Process: all significant ----
    tryCatch({
      go_bp <- clusterProfiler::enrichGO(
        gene     = sig_entrez$ENTREZID,
        universe = all_entrez$ENTREZID,
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable = TRUE
      )
      if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
        write.csv(go_bp@result, file.path(out_path, "enrichment_GO_BP_all.csv"),
                  row.names = FALSE)
        results_list$GO_BP_all <- go_bp
        if (nrow(go_bp@result) >= 5) {
          p <- enrichplot::dotplot(go_bp, showCategory = 20) +
            ggtitle(paste0("GO BP — All Significant\n", label))
          ggsave(file.path(out_path, "enrichment_GO_BP_all_dotplot.png"), p,
                 width = 10, height = 8, dpi = 300, bg = "transparent")
        }
        message(sprintf("    GO BP (all): %d terms", nrow(go_bp@result)))
      }
    }, error = function(e) message("    ⚠ GO BP all: ", e$message))
    
    # ---- GO BP: up-regulated ----
    if (nrow(up_entrez) >= 5) {
      tryCatch({
        go_up <- clusterProfiler::enrichGO(
          gene = up_entrez$ENTREZID, universe = all_entrez$ENTREZID,
          OrgDb = org.Hs.eg.db, ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
          readable = TRUE)
        if (!is.null(go_up) && nrow(go_up@result) > 0) {
          write.csv(go_up@result, file.path(out_path, "enrichment_GO_BP_up_in_DCM.csv"),
                    row.names = FALSE)
          results_list$GO_BP_up <- go_up
          if (nrow(go_up@result) >= 5) {
            p <- enrichplot::dotplot(go_up, showCategory = 15) +
              ggtitle(paste0("GO BP — Up in DCM\n", label))
            ggsave(file.path(out_path, "enrichment_GO_BP_up_dotplot.png"), p,
                   width = 10, height = 8, dpi = 300, bg = "transparent")
          }
          message(sprintf("    GO BP (up in DCM): %d terms", nrow(go_up@result)))
        }
      }, error = function(e) message("    ⚠ GO BP up: ", e$message))
    }
    
    # ---- GO BP: down-regulated ----
    if (nrow(down_entrez) >= 5) {
      tryCatch({
        go_down <- clusterProfiler::enrichGO(
          gene = down_entrez$ENTREZID, universe = all_entrez$ENTREZID,
          OrgDb = org.Hs.eg.db, ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
          readable = TRUE)
        if (!is.null(go_down) && nrow(go_down@result) > 0) {
          write.csv(go_down@result, file.path(out_path, "enrichment_GO_BP_down_in_DCM.csv"),
                    row.names = FALSE)
          results_list$GO_BP_down <- go_down
          if (nrow(go_down@result) >= 5) {
            p <- enrichplot::dotplot(go_down, showCategory = 15) +
              ggtitle(paste0("GO BP — Down in DCM (Up in Healthy)\n", label))
            ggsave(file.path(out_path, "enrichment_GO_BP_down_dotplot.png"), p,
                   width = 10, height = 8, dpi = 300, bg = "transparent")
          }
          message(sprintf("    GO BP (down in DCM): %d terms", nrow(go_down@result)))
        }
      }, error = function(e) message("    ⚠ GO BP down: ", e$message))
    }
    
    # ---- KEGG ----
    tryCatch({
      kegg <- clusterProfiler::enrichKEGG(
        gene = sig_entrez$ENTREZID, universe = all_entrez$ENTREZID,
        organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(kegg) && nrow(kegg@result) > 0) {
        write.csv(kegg@result, file.path(out_path, "enrichment_KEGG.csv"),
                  row.names = FALSE)
        results_list$KEGG <- kegg
        if (nrow(kegg@result) >= 5) {
          p <- enrichplot::dotplot(kegg, showCategory = 15) +
            ggtitle(paste0("KEGG Pathways\n", label))
          ggsave(file.path(out_path, "enrichment_KEGG_dotplot.png"), p,
                 width = 10, height = 8, dpi = 300, bg = "transparent")
        }
        message(sprintf("    KEGG: %d pathways", nrow(kegg@result)))
      }
    }, error = function(e) message("    ⚠ KEGG: ", e$message))
    
    # ---- Reactome ----
    tryCatch({
      if (requireNamespace("ReactomePA", quietly = TRUE)) {
        library(ReactomePA)
        reactome <- ReactomePA::enrichPathway(
          gene = sig_entrez$ENTREZID, universe = all_entrez$ENTREZID,
          organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
          readable = TRUE)
        if (!is.null(reactome) && nrow(reactome@result) > 0) {
          write.csv(reactome@result, file.path(out_path, "enrichment_Reactome.csv"),
                    row.names = FALSE)
          results_list$Reactome <- reactome
          if (nrow(reactome@result) >= 5) {
            p <- enrichplot::dotplot(reactome, showCategory = 15) +
              ggtitle(paste0("Reactome Pathways\n", label))
            ggsave(file.path(out_path, "enrichment_Reactome_dotplot.png"), p,
                   width = 10, height = 8, dpi = 300, bg = "transparent")
          }
          message(sprintf("    Reactome: %d pathways", nrow(reactome@result)))
        }
      } else {
        message("    ℹ ReactomePA not installed, skipping Reactome")
      }
    }, error = function(e) message("    ⚠ Reactome: ", e$message))
    
    return(results_list)
    
  }, error = function(e) {
    message("  ⚠ clusterProfiler enrichment failed: ", e$message)
    return(NULL)
  })
}

# ==============================================================================
# ENRICHMENT: fgsea GSEA (Hallmarks, KEGG, Reactome) — per comparison
# Complements clusterProfiler with ranked gene-set enrichment
# ==============================================================================

run_fgsea_enrichment <- function(top_table, out_path, label) {
  tryCatch({
    suppressPackageStartupMessages({
      library(fgsea)
      library(msigdbr)
    })
    
    cat(sprintf("\n  ── fgsea GSEA: %s ──\n", label))
    
    ranked <- top_table %>%
      filter(!is.na(Gene) & Gene != "") %>%
      arrange(P.Value) %>%
      filter(!duplicated(Gene)) %>%
      arrange(desc(t))
    
    gene_ranks <- setNames(ranked$t, ranked$Gene)
    cat(sprintf("    Ranked gene list: %d genes\n", length(gene_ranks)))
    
    # ---- Build gene sets — handle msigdbr version differences ----
    # msigdbr >= 7.5.1 renamed CP:KEGG → CP:KEGG_LEGACY / CP:KEGG_MEDICUS
    # and may have changed other subcategory names.
    # Strategy: pull all C2, then filter by gs_subcat using grep.
    
    # Hallmarks (stable across versions)
    hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
      split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
    cat(sprintf("    Hallmarks: %d gene sets\n", length(hallmarks)))
    
    # KEGG — try CP:KEGG first, fall back to grep for KEGG in subcategory
    kegg <- tryCatch({
      msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
        split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
    }, error = function(e) {
      cat("    ℹ CP:KEGG not found, trying KEGG_LEGACY/KEGG_MEDICUS...\n")
      all_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
      kegg_sets <- all_c2[grepl("KEGG", all_c2$gs_subcat, ignore.case = TRUE), ]
      if (nrow(kegg_sets) > 0) {
        kegg_sets %>% split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
      } else {
        cat("    ⚠ No KEGG gene sets found in msigdbr\n")
        list()
      }
    })
    cat(sprintf("    KEGG: %d gene sets\n", length(kegg)))
    
    # Reactome — same strategy
    reactome <- tryCatch({
      msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
        split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
    }, error = function(e) {
      cat("    ℹ CP:REACTOME not found, trying grep fallback...\n")
      all_c2 <- if (exists("all_c2")) all_c2 else msigdbr(species = "Homo sapiens", category = "C2")
      reac_sets <- all_c2[grepl("REACTOME", all_c2$gs_subcat, ignore.case = TRUE), ]
      if (nrow(reac_sets) > 0) {
        reac_sets %>% split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
      } else {
        cat("    ⚠ No Reactome gene sets found in msigdbr\n")
        list()
      }
    })
    cat(sprintf("    Reactome: %d gene sets\n", length(reactome)))
    
    # ---- Run fgsea on each gene set collection ----
    run_one <- function(gsets, name) {
      if (length(gsets) == 0) {
        cat(sprintf("    Skipping %s (no gene sets)\n", name))
        return(NULL)
      }
      res <- fgsea(pathways = gsets, stats = gene_ranks,
                   minSize = 10, maxSize = 500, nPermSimple = 10000)
      res <- res[order(res$pval), ]
      res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
      write.csv(as.data.frame(res), file.path(out_path, paste0("fgsea_", name, ".csv")),
                row.names = FALSE)
      n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
      cat(sprintf("    fgsea %s: %d tested, %d FDR<0.05\n", name, nrow(res), n_sig))
      return(res)
    }
    
    res_hall <- run_one(hallmarks, "Hallmarks")
    res_kegg <- run_one(kegg, "KEGG")
    res_reac <- run_one(reactome, "Reactome")
    
    # ---- Hallmarks barplot ----
    tryCatch({
      if (!is.null(res_hall) && nrow(res_hall) > 0) {
        top_hall <- as.data.frame(res_hall) %>%
          filter(padj < 0.25) %>%
          arrange(NES) %>%
          mutate(pathway = gsub("^HALLMARK_", "", pathway),
                 pathway = gsub("_", " ", pathway),
                 direction = ifelse(NES > 0, "Up in DCM", "Down in DCM"))
        
        if (nrow(top_hall) > 0) {
          p <- ggplot(top_hall, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
            geom_col(alpha = 0.85) +
            scale_fill_manual(values = DIRECTION_COLORS) +
            labs(title = paste0("Hallmark Pathways (fgsea) — ", label),
                 subtitle = "FDR < 0.25",
                 x = "Normalized Enrichment Score", y = NULL, fill = NULL) +
            theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  legend.position = "bottom")
          
          ggsave(file.path(out_path, "fgsea_hallmarks_barplot.png"), p,
                 width = 10, height = max(4, nrow(top_hall) * 0.3 + 2), dpi = 300, bg = "transparent")
          cat("    ✓ Hallmarks barplot saved\n")
        }
      }
    }, error = function(e) cat(sprintf("    ⚠ Hallmarks barplot: %s\n", e$message)))
    
    # ---- KEGG barplot ----
    tryCatch({
      if (!is.null(res_kegg) && nrow(res_kegg) > 0) {
        top_kegg <- as.data.frame(res_kegg) %>%
          filter(padj < 0.25) %>%
          arrange(NES) %>%
          mutate(pathway = gsub("^KEGG_", "", pathway),
                 pathway = gsub("_", " ", pathway),
                 direction = ifelse(NES > 0, "Up in DCM", "Down in DCM"))
        
        if (nrow(top_kegg) > 0) {
          p <- ggplot(top_kegg, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
            geom_col(alpha = 0.85) +
            scale_fill_manual(values = DIRECTION_COLORS) +
            labs(title = paste0("KEGG Pathways (fgsea) — ", label),
                 subtitle = "FDR < 0.25",
                 x = "Normalized Enrichment Score", y = NULL, fill = NULL) +
            theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  legend.position = "bottom")
          
          ggsave(file.path(out_path, "fgsea_kegg_barplot.png"), p,
                 width = 10, height = max(4, nrow(top_kegg) * 0.3 + 2), dpi = 300, bg = "transparent")
          cat("    ✓ KEGG barplot saved\n")
        }
      }
    }, error = function(e) cat(sprintf("    ⚠ KEGG barplot: %s\n", e$message)))
    
    # ---- Reactome barplot ----
    tryCatch({
      if (!is.null(res_reac) && nrow(res_reac) > 0) {
        top_reac <- as.data.frame(res_reac) %>%
          filter(padj < 0.25) %>%
          slice_head(n = 30) %>%  # Reactome can have many — cap at 30
          arrange(NES) %>%
          mutate(pathway = gsub("^REACTOME_", "", pathway),
                 pathway = gsub("_", " ", pathway),
                 direction = ifelse(NES > 0, "Up in DCM", "Down in DCM"))
        
        if (nrow(top_reac) > 0) {
          p <- ggplot(top_reac, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
            geom_col(alpha = 0.85) +
            scale_fill_manual(values = DIRECTION_COLORS) +
            labs(title = paste0("Reactome Pathways (fgsea) — ", label),
                 subtitle = "FDR < 0.25 | top 30",
                 x = "Normalized Enrichment Score", y = NULL, fill = NULL) +
            theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  legend.position = "bottom")
          
          ggsave(file.path(out_path, "fgsea_reactome_barplot.png"), p,
                 width = 11, height = max(4, nrow(top_reac) * 0.3 + 2), dpi = 300, bg = "transparent")
          cat("    ✓ Reactome barplot saved\n")
        }
      }
    }, error = function(e) cat(sprintf("    ⚠ Reactome barplot: %s\n", e$message)))
    
    return(list(hallmarks = res_hall, kegg = res_kegg, reactome = res_reac))
    
  }, error = function(e) {
    cat(sprintf("\n  ⚠⚠⚠ fgsea FAILED: %s\n", e$message))
    cat("    Check: library(fgsea); library(msigdbr)\n")
    cat("    Try:   msigdbr_collections() to see available subcategories\n\n")
    return(NULL)
  })
}

# ==============================================================================
# RUN_COMPARISON — main function per comparison
# ==============================================================================

run_comparison <- function(data_mat, gene_symbols, disease_cols, healthy_cols,
                           label, output_dir,
                           filtering_min = FILTERING_MIN_SAMPLES,
                           imputation    = IMPUTATION_METHOD,
                           normalization = NORMALIZATION_METHOD) {
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("  %s\n", label))
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  
  out_path <- file.path(OUTPUT_BASE_DIR, output_dir)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Subset columns ----
  d_cols <- intersect(disease_cols, colnames(data_mat))
  h_cols <- intersect(healthy_cols, colnames(data_mat))
  nD <- length(d_cols)
  nH <- length(h_cols)
  
  message(sprintf("  Samples: %d Disease vs %d Healthy", nD, nH))
  
  if (nD < 2 || nH < 2) {
    message("  ✗ Skipping: too few samples")
    return(NULL)
  }
  
  sub_mat <- data_mat[, c(d_cols, h_cols), drop = FALSE]
  condition <- factor(c(rep("Disease", nD), rep("Healthy", nH)),
                      levels = c("Healthy", "Disease"))
  names(condition) <- c(d_cols, h_cols)
  
  # Sample labels for plots
  sample_labels <- setNames(
    paste0(condition, "_", ave(as.numeric(condition), condition, FUN = seq_along)),
    names(condition)
  )
  
  # ---- Filter ----
  n_valid <- rowSums(!is.na(sub_mat) & sub_mat > 0)
  keep    <- n_valid >= filtering_min
  sub_mat <- sub_mat[keep, , drop = FALSE]
  message(sprintf("  Filtering (>= %d samples): %d → %d proteins",
                  filtering_min, nrow(data_mat), nrow(sub_mat)))
  
  # ---- Impute ----
  sub_mat[sub_mat == 0] <- NA
  if (imputation == "min") {
    for (i in seq_len(nrow(sub_mat))) {
      nas <- is.na(sub_mat[i, ])
      if (any(nas) && !all(nas)) {
        sub_mat[i, nas] <- min(sub_mat[i, !nas]) / 2
      }
    }
  }
  
  # ---- Log2 transform ----
  min_val <- min(sub_mat[sub_mat > 0], na.rm = TRUE)
  offset  <- ifelse(min_val < 1, 1, 0)
  log2_mat <- log2(sub_mat + offset)
  log2_mat[is.infinite(log2_mat)] <- NA
  
  # ---- Normalize ----
  if (normalization == "median_center") {
    norm_mat <- apply(log2_mat, 2, function(y) {
      med <- median(y, na.rm = TRUE)
      y - med
    })
    rownames(norm_mat) <- rownames(log2_mat)
  } else if (normalization == "quantile") {
    temp <- new("EList", list(E = 2^log2_mat))
    temp <- limma::normalizeBetweenArrays(temp, method = "quantile")
    norm_mat <- log2(temp$E)
    rownames(norm_mat) <- rownames(log2_mat)
  } else {
    norm_mat <- log2_mat
  }
  
  # =============== QC PLOTS ===============
  message("  QC plots...")
  avg_r <- generate_qc_plots(log2_mat, norm_mat, condition, sample_labels,
                             out_path, label)
  
  # =============== PCA / MDS ===============
  message("  PCA / MDS...")
  generate_pca_mds(norm_mat, condition, sample_labels, out_path, label)
  
  # =============== LIMMA ===============
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("Healthy", "Disease")
  contrast_mat <- limma::makeContrasts(Disease - Healthy, levels = design)
  
  fit  <- limma::lmFit(norm_mat, design)
  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2, robust = TRUE, trend = TRUE)
  
  top_table <- limma::topTable(fit2, number = Inf, sort.by = "none")
  top_table$Protein <- rownames(top_table)
  top_table$Gene <- gene_symbols[top_table$Protein]
  top_table$Gene[is.na(top_table$Gene)] <- top_table$Protein[is.na(top_table$Gene)]
  
  # ---- Summary stats ----
  n_total  <- nrow(top_table)
  n_fdr01  <- sum(top_table$adj.P.Val < 0.01, na.rm = TRUE)
  n_fdr05  <- sum(top_table$adj.P.Val < 0.05, na.rm = TRUE)
  n_fdr10  <- sum(top_table$adj.P.Val < 0.10, na.rm = TRUE)
  n_fdr20  <- sum(top_table$adj.P.Val < 0.20, na.rm = TRUE)
  n_pv01   <- sum(top_table$P.Value   < 0.01, na.rm = TRUE)
  n_up_fdr05   <- sum(top_table$adj.P.Val < 0.05 & top_table$logFC > 0, na.rm = TRUE)
  n_down_fdr05 <- sum(top_table$adj.P.Val < 0.05 & top_table$logFC < 0, na.rm = TRUE)
  pct_sig  <- round(100 * n_fdr05 / n_total, 1)
  
  cat(sprintf("  Total proteins tested : %d\n",  n_total))
  cat(sprintf("  FDR < 0.01           : %d\n",  n_fdr01))
  cat(sprintf("  FDR < 0.05           : %d (%d up / %d down)\n",
              n_fdr05, n_up_fdr05, n_down_fdr05))
  cat(sprintf("  FDR < 0.10           : %d\n",  n_fdr10))
  cat(sprintf("  FDR < 0.20           : %d\n",  n_fdr20))
  cat(sprintf("  raw p < 0.01         : %d\n",  n_pv01))
  cat(sprintf("  %% significant (FDR<0.05): %.1f%%\n", pct_sig))
  
  # ---- Save results ----
  write.csv(top_table, file.path(out_path, "top_table_all.csv"), row.names = FALSE)
  write.csv(top_table[top_table$adj.P.Val < 0.05, ],
            file.path(out_path, "DE_FDR005.csv"), row.names = FALSE)
  write.csv(top_table[top_table$adj.P.Val < 0.10, ],
            file.path(out_path, "DE_FDR010.csv"), row.names = FALSE)
  write.csv(top_table[top_table$adj.P.Val < 0.20, ],
            file.path(out_path, "DE_FDR020.csv"), row.names = FALSE)
  write.csv(top_table[top_table$P.Value < 0.01, ],
            file.path(out_path, "DE_rawP001.csv"), row.names = FALSE)
  
  # Export normalized matrix
  norm_export <- as.data.frame(norm_mat)
  norm_export$Protein <- rownames(norm_mat)
  norm_export$Gene    <- gene_symbols[norm_export$Protein]
  norm_export <- norm_export[, c("Protein", "Gene", colnames(norm_mat))]
  write.table(norm_export, file.path(out_path, "normalized_protein_data.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("  ✓ DE tables + normalized data saved")
  
  # ---- Volcano plots: raw p-value AND FDR ----
  # Square format | new colour scheme | transparent background
  tryCatch({
    plot_df <- top_table %>%
      mutate(
        neg_log10_p   = -log10(pmax(P.Value, 1e-50)),
        neg_log10_fdr = -log10(pmax(adj.P.Val, 1e-50)),
        sig = case_when(
          logFC >= log2(2)  & adj.P.Val < 0.05 ~ "Up in DCM",
          logFC <= -log2(2) & adj.P.Val < 0.05 ~ "Down in DCM",
          adj.P.Val < 0.05                     ~ "FDR < 0.05",
          TRUE                                  ~ "NS"
        ),
        label = ifelse(
          adj.P.Val < 0.01 &
            abs(logFC) > quantile(abs(logFC[adj.P.Val < 0.01]), 0.9, na.rm = TRUE),
          Gene, "")
      )
    
    volcano_colors <- c("Up in DCM"   = COL_DCM,
                        "Down in DCM" = COL_HEALTHY,
                        "FDR < 0.05"  = COL_FDR,
                        "NS"          = COL_NS)
    
    volcano_theme <- theme_minimal(base_size = 11) +
      theme(plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom",
            aspect.ratio = 1) +
      THEME_TRANSPARENT
    
    # ---- (A) Volcano: raw p-value on y-axis ----
    p_volcano_pval <- ggplot(plot_df, aes(x = logFC, y = neg_log10_p, color = sig)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 25,
                      fontface = "bold", color = "black") +
      scale_color_manual(values = volcano_colors) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = COL_HEALTHY,
                 linewidth = 0.6) +
      annotate("text", x = max(plot_df$logFC, na.rm = TRUE), y = -log10(0.05),
               label = "p = 0.05", hjust = 1, vjust = -0.5, size = 3,
               fontface = "italic", color = COL_HEALTHY) +
      geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      labs(title    = label,
           subtitle = sprintf("Proteomics | %d proteins | %d FDR<0.05 (%d up / %d down)",
                              n_total, n_fdr05, n_up_fdr05, n_down_fdr05),
           x = "log2 Fold Change (Disease – Healthy)",
           y = expression("-log"[10]*"(p-value)"), color = NULL) +
      volcano_theme
    
    ggsave(file.path(out_path, "volcano_pvalue.png"), p_volcano_pval,
           width = 8, height = 8, dpi = 300, bg = "transparent")
    ggsave(file.path(out_path, "volcano_pvalue.pdf"), p_volcano_pval,
           width = 8, height = 8, bg = "transparent")
    message("  ✓ Volcano (p-value) saved")
    
    # ---- (B) Volcano: FDR/adj.P.Val on y-axis ----
    p_volcano_fdr <- ggplot(plot_df, aes(x = logFC, y = neg_log10_fdr, color = sig)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 25,
                      fontface = "bold", color = "black") +
      scale_color_manual(values = volcano_colors) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = COL_HEALTHY,
                 linewidth = 0.6) +
      annotate("text", x = max(plot_df$logFC, na.rm = TRUE), y = -log10(0.05),
               label = "FDR = 0.05", hjust = 1, vjust = -0.5, size = 3,
               fontface = "italic", color = COL_HEALTHY) +
      geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      labs(title    = label,
           subtitle = sprintf("Proteomics | %d proteins | %d FDR<0.05 (%d up / %d down)",
                              n_total, n_fdr05, n_up_fdr05, n_down_fdr05),
           x = "log2 Fold Change (Disease – Healthy)",
           y = expression("-log"[10]*"(FDR / adj.P.Val)"), color = NULL) +
      volcano_theme
    
    ggsave(file.path(out_path, "volcano_FDR.png"), p_volcano_fdr,
           width = 8, height = 8, dpi = 300, bg = "transparent")
    ggsave(file.path(out_path, "volcano_FDR.pdf"), p_volcano_fdr,
           width = 8, height = 8, bg = "transparent")
    message("  ✓ Volcano (FDR) saved")
  }, error = function(e) message("  ⚠ Volcano: ", e$message))
  
  # ---- Top DE heatmap (excluding blood markers, square cells) ----
  tryCatch({
    message("  Top DE heatmap (blood-excluded)...")
    
    # Select top DE proteins: FDR < 0.05, exclude blood markers
    de_for_heatmap <- top_table %>%
      filter(adj.P.Val < 0.05,
             !Gene %in% BLOOD_MARKER_GENES,
             !is.na(Gene), Gene != "") %>%
      arrange(adj.P.Val) %>%
      filter(!duplicated(Gene)) %>%
      slice_head(n = 50)  # top 50 by significance
    
    if (nrow(de_for_heatmap) >= 5) {
      # Extract normalized intensities for these proteins
      heat_proteins <- de_for_heatmap$Protein
      heat_mat <- norm_mat[heat_proteins, , drop = FALSE]
      rownames(heat_mat) <- de_for_heatmap$Gene
      
      # Z-score per row (protein)
      heat_z <- t(scale(t(heat_mat)))
      heat_z[is.na(heat_z)] <- 0
      
      # Melt for ggplot
      melted <- reshape2::melt(as.matrix(heat_z))
      colnames(melted) <- c("Gene", "Sample", "Z_score")
      melted$Condition <- as.character(condition[as.character(melted$Sample)])
      
      # Order: genes by logFC, samples by condition then name
      gene_order <- de_for_heatmap %>% arrange(logFC) %>% pull(Gene)
      sample_order <- names(sort(condition))
      melted$Gene   <- factor(melted$Gene, levels = gene_order)
      melted$Sample <- factor(melted$Sample, levels = sample_order)
      
      p_heat <- ggplot(melted, aes(x = Sample, y = Gene, fill = Z_score)) +
        geom_tile(color = NA) +
        scale_fill_gradient2(low = COL_HEALTHY, mid = "white", high = COL_DCM,
                             midpoint = 0, name = "Z-score") +
        facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
        labs(title = paste0("Top DE Proteins (blood-excluded) — ", label),
             subtitle = sprintf("Top %d by FDR | square cells | Disease – Healthy",
                                nrow(de_for_heatmap)),
             x = NULL, y = NULL) +
        theme_minimal(base_size = 8) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
              axis.text.y = element_text(size = 6),
              strip.text  = element_text(face = "bold", size = 10),
              plot.title  = element_text(face = "bold", size = 12),
              aspect.ratio = nrow(de_for_heatmap) / ncol(norm_mat)) +
        THEME_TRANSPARENT +
        coord_fixed(ratio = 1)  # square cells
      
      ggsave(file.path(out_path, "heatmap_top_DE_no_blood.png"), p_heat,
             width = max(10, ncol(norm_mat) * 0.22 + 4),
             height = max(8, nrow(de_for_heatmap) * 0.22 + 3),
             dpi = 300, bg = "transparent")
      message(sprintf("  ✓ Heatmap saved (%d proteins, blood markers excluded)",
                      nrow(de_for_heatmap)))
    } else {
      message("  ℹ Too few non-blood DE proteins for heatmap (<5)")
    }
  }, error = function(e) message("  ⚠ Heatmap: ", e$message))
  
  # =============== ENRICHMENT ===============
  cp_results    <- NULL
  fgsea_results <- NULL
  
  if (exists("ENRICHMENT_AVAILABLE") && ENRICHMENT_AVAILABLE) {
    # 1. clusterProfiler (ORA: GO BP, KEGG, Reactome)
    cp_results <- run_clusterprofiler_enrichment(top_table, out_path, label)
    
    # 2. fgsea GSEA (Hallmarks, KEGG, Reactome — ranked)
    fgsea_results <- run_fgsea_enrichment(top_table, out_path, label)
    
    if (is.null(cp_results) && is.null(fgsea_results)) {
      cat(sprintf("  ⚠ ENRICHMENT PRODUCED NO RESULTS for %s\n", label))
      cat("    Possible causes: too few genes mapped, or gene symbols don't match MSigDB.\n")
      cat(sprintf("    Genes with symbols: %d / %d\n",
                  sum(!is.na(top_table$Gene) & top_table$Gene != "" &
                        top_table$Gene != top_table$Protein),
                  nrow(top_table)))
    }
  } else {
    cat("  ── Enrichment SKIPPED (packages not installed) ──\n")
  }
  
  # ---- Return ----
  return(list(
    top_table      = top_table,
    fit            = fit2,
    norm_mat       = norm_mat,
    avg_r          = avg_r,
    cp_enrichment  = cp_results,
    fgsea_enrichment = fgsea_results,
    nD = nD, nH = nH,
    n_total = n_total,
    n_fdr01 = n_fdr01, n_fdr05 = n_fdr05, n_fdr10 = n_fdr10, n_fdr20 = n_fdr20,
    n_pv01 = n_pv01,
    n_up_fdr05 = n_up_fdr05, n_down_fdr05 = n_down_fdr05,
    pct_sig = pct_sig
  ))
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\n╔══════════════════════════════════════════════════════════════════════════════╗")
message("║          HUMAN DCM PROTEOMICS BENCHMARKING — v3                            ║")
message("╚══════════════════════════════════════════════════════════════════════════════╝\n")

message("Loading: ", PROTEOMICS_DATA_PATH)
prot_raw <- read.csv(PROTEOMICS_DATA_PATH, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)

if (any(grepl("^X", colnames(prot_raw)))) {
  colnames(prot_raw) <- gsub("^X", "", colnames(prot_raw))
}

message(sprintf("  ✓ %d proteins x %d columns", nrow(prot_raw), ncol(prot_raw)))

# ---- Build intensity matrix ----
protein_ids <- prot_raw$PG.GroupLabel
meta_cols   <- c("rowname", "PG.GroupLabel")
intensity_col_names <- setdiff(colnames(prot_raw), meta_cols)

expected <- c(healthy_cols, dcm_filtered_cols)
missing  <- expected[!expected %in% intensity_col_names]
if (length(missing) > 0) {
  warning("MISSING expected columns: ", paste(missing, collapse = ", "))
} else {
  message("  ✓ All 27 core columns (14 Healthy + 13 DCM_filtered) verified")
}

d_present <- sum(healthy_cols %in% intensity_col_names)
v_present <- sum(all_vad_cols %in% intensity_col_names)
t_present <- sum(all_tx_cols  %in% intensity_col_names)
message(sprintf("  Column breakdown: %d -D (Healthy), %d -V (VAD), %d -T (Tx)",
                d_present, v_present, t_present))

data_mat <- as.matrix(prot_raw[, intensity_col_names])
rownames(data_mat) <- protein_ids
storage.mode(data_mat) <- "numeric"
data_mat[data_mat == 0] <- NA

message(sprintf("  ✓ Intensity matrix: %d proteins x %d samples",
                nrow(data_mat), ncol(data_mat)))
message(sprintf("  Missing values: %d (%.1f%%)",
                sum(is.na(data_mat)),
                100 * sum(is.na(data_mat)) / length(data_mat)))

# ---- Map UniProt → Gene symbols ----
message("\nMapping UniProt accessions to gene symbols...")
gene_symbols <- map_uniprot_to_gene(protein_ids)

dir.create(OUTPUT_BASE_DIR, showWarnings = FALSE)
write.csv(data.frame(UniProt = protein_ids, Gene = gene_symbols,
                     stringsAsFactors = FALSE),
          file.path(OUTPUT_BASE_DIR, "uniprot_gene_mapping.csv"), row.names = FALSE)

# ==============================================================================
# PACKAGE AVAILABILITY CHECK — enrichment + gene mapping
# ==============================================================================
# This runs BEFORE comparisons so you know immediately what's missing.

cat("\n")
cat("================================================================================\n")
cat("  PACKAGE AVAILABILITY CHECK\n")
cat("================================================================================\n")

ENRICHMENT_PACKAGES <- list(
  # Package name         = Install command
  "org.Hs.eg.db"    = 'BiocManager::install("org.Hs.eg.db")',
  "AnnotationDbi"   = 'BiocManager::install("AnnotationDbi")',
  "clusterProfiler" = 'BiocManager::install("clusterProfiler")',
  "enrichplot"      = 'BiocManager::install("enrichplot")',
  "ReactomePA"      = 'BiocManager::install("ReactomePA")',
  "fgsea"           = 'install.packages("fgsea")',
  "msigdbr"         = 'install.packages("msigdbr")'
)

pkg_status <- sapply(names(ENRICHMENT_PACKAGES), requireNamespace, quietly = TRUE)
ENRICHMENT_AVAILABLE <- all(pkg_status)

if (ENRICHMENT_AVAILABLE) {
  cat("  ✓ All enrichment packages installed — enrichment will run.\n")
} else {
  missing_pkgs <- names(pkg_status[!pkg_status])
  cat("\n")
  cat("  ╔════════════════════════════════════════════════════════════════════╗\n")
  cat("  ║  ⚠⚠⚠  ENRICHMENT PACKAGES MISSING — ENRICHMENT WILL BE SKIPPED  ║\n")
  cat("  ╚════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("  Missing packages:\n")
  for (pkg in missing_pkgs) {
    cat(sprintf("    ✗ %-20s  →  %s\n", pkg, ENRICHMENT_PACKAGES[[pkg]]))
  }
  cat("\n  Install all at once:\n")
  cat("    if (!requireNamespace('BiocManager', quietly = TRUE))\n")
  cat("      install.packages('BiocManager')\n")
  cat("    BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db',\n")
  cat("                           'enrichplot', 'ReactomePA'))\n")
  cat("    install.packages(c('fgsea', 'msigdbr'))\n")
  cat("\n  DE analysis, QC, PCA, and volcano plots will still run normally.\n")
  cat("  Re-run after installing to get enrichment results.\n")
}

# Also check if gene mapping worked (critical for enrichment)
n_mapped <- sum(!is.na(gene_symbols))
if (n_mapped == 0) {
  cat("\n  ⚠ Gene symbol mapping returned 0 matches!\n")
  cat("    Enrichment requires gene symbols, not UniProt IDs.\n")
  cat("    Install org.Hs.eg.db:  BiocManager::install('org.Hs.eg.db')\n")
  ENRICHMENT_AVAILABLE <- FALSE
} else {
  cat(sprintf("  ✓ Gene mapping: %d / %d UniProt → gene symbol\n",
              n_mapped, length(gene_symbols)))
}

cat("================================================================================\n\n")

# ==============================================================================
# RUN ALL 8 COMPARISONS
# ==============================================================================

all_results <- list()

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  result <- tryCatch(
    run_comparison(
      data_mat     = data_mat,
      gene_symbols = gene_symbols,
      disease_cols = comp$d,
      healthy_cols = comp$h,
      label        = comp$label,
      output_dir   = comp$dir
    ),
    error = function(e) {
      message(sprintf("  ✗ %s FAILED: %s", comp_name, e$message))
      return(NULL)
    }
  )
  if (!is.null(result)) {
    all_results[[comp_name]] <- result
  }
}

# ==============================================================================
# BENCHMARK SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  BENCHMARK SUMMARY — PROTEOMICS (Disease vs Healthy)\n")
cat("================================================================================\n")

bench_df <- do.call(rbind, lapply(names(all_results), function(comp_name) {
  r <- all_results[[comp_name]]
  data.frame(
    Comparison = comparisons[[comp_name]]$label,
    nD = r$nD, nH = r$nH, Total = r$n_total,
    FDR01 = r$n_fdr01, FDR05 = r$n_fdr05, FDR10 = r$n_fdr10, FDR20 = r$n_fdr20,
    pv01 = r$n_pv01,
    Up_FDR05 = r$n_up_fdr05, Down_FDR05 = r$n_down_fdr05,
    pct_sig = r$pct_sig,
    stringsAsFactors = FALSE
  )
}))

bench_df <- bench_df[order(-bench_df$FDR05), ]

cat(sprintf("\n%-65s %3s %3s %5s %5s %5s %5s %5s %5s %5s\n",
            "Comparison", "nD", "nH", "Total", "FDR01", "FDR05", "FDR10", "FDR20", "pv01", "%sig"))
cat(paste(rep("─", 120), collapse = ""), "\n")

for (i in seq_len(nrow(bench_df))) {
  cat(sprintf("%-65s %3d %3d %5d %5d %5d %5d %5d %5d %5.1f\n",
              bench_df$Comparison[i],
              bench_df$nD[i], bench_df$nH[i], bench_df$Total[i],
              bench_df$FDR01[i], bench_df$FDR05[i], bench_df$FDR10[i],
              bench_df$FDR20[i], bench_df$pv01[i], bench_df$pct_sig[i]))
}

write.csv(bench_df, file.path(OUTPUT_BASE_DIR, "benchmark_summary.csv"), row.names = FALSE)

# ==============================================================================
# BENCHMARK BARPLOTS
# ==============================================================================

tryCatch({
  plot_data <- bench_df %>%
    pivot_longer(cols = c(FDR01, FDR05, FDR10),
                 names_to = "Threshold", values_to = "N_sig") %>%
    mutate(
      Threshold  = factor(Threshold,
                          levels = c("FDR01", "FDR05", "FDR10"),
                          labels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.10")),
      Comparison = factor(Comparison, levels = rev(bench_df$Comparison))
    )
  
  p_bench <- ggplot(plot_data, aes(x = Comparison, y = N_sig, fill = Threshold)) +
    geom_col(position = "dodge", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = BENCH_COLORS) +
    labs(title = "Benchmark: DCM vs Healthy — significant proteins",
         subtitle = "Multiple grouping strategies | one-per-patient where applicable",
         x = NULL, y = "N significant proteins", fill = "Threshold") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom") +
    geom_text(aes(label = N_sig),
              position = position_dodge(width = 0.8),
              hjust = -0.1, size = 2.8)
  
  ggsave(file.path(OUTPUT_BASE_DIR, "benchmark_barplot.png"), p_bench,
         width = 14, height = 9, dpi = 300, bg = "transparent")
  message("  ✓ Benchmark barplot saved")
}, error = function(e) message("  ⚠ Barplot: ", e$message))

# Direction barplot
tryCatch({
  dir_data <- bench_df %>%
    pivot_longer(cols = c(Up_FDR05, Down_FDR05),
                 names_to = "Direction", values_to = "N") %>%
    mutate(
      Direction  = factor(Direction,
                          levels = c("Up_FDR05", "Down_FDR05"),
                          labels = c("Up in DCM", "Down in DCM")),
      Comparison = factor(Comparison, levels = rev(bench_df$Comparison))
    )
  
  p_dir <- ggplot(dir_data, aes(x = Comparison, y = N, fill = Direction)) +
    geom_col(position = "dodge", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = DIRECTION_COLORS) +
    labs(title = "Proteomics: Up vs Down regulation in DCM (FDR < 0.05)",
         subtitle = "Contrast: Disease – Healthy",
         x = NULL, y = "N significant proteins", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom") +
    geom_text(aes(label = N),
              position = position_dodge(width = 0.8),
              hjust = -0.1, size = 2.8)
  
  ggsave(file.path(OUTPUT_BASE_DIR, "benchmark_direction_barplot.png"), p_dir,
         width = 14, height = 9, dpi = 300, bg = "transparent")
  message("  ✓ Direction barplot saved")
}, error = function(e) message("  ⚠ Direction barplot: ", e$message))

# ==============================================================================
# CROSS-COMPARISON: logFC CORRELATION MATRIX + DENDROGRAM
# ==============================================================================
# How similar are the DE profiles across comparisons?
# Uses logFC vectors aligned by gene symbol (Pearson r).

tryCatch({
  message("\n  Cross-comparison correlation analysis...")
  
  # Build a gene × comparison logFC matrix (use Gene, keep best p per gene)
  logfc_list <- lapply(names(all_results), function(comp) {
    tt <- all_results[[comp]]$top_table %>%
      filter(!is.na(Gene) & Gene != "") %>%
      arrange(P.Value) %>%
      filter(!duplicated(Gene))
    setNames(tt$logFC, tt$Gene)
  })
  names(logfc_list) <- names(all_results)
  
  # Intersect genes present in ALL comparisons
  shared_genes <- Reduce(intersect, lapply(logfc_list, names))
  message(sprintf("  Shared genes across all comparisons: %d", length(shared_genes)))
  
  if (length(shared_genes) > 100) {
    logfc_mat <- sapply(logfc_list, function(v) v[shared_genes])
    rownames(logfc_mat) <- shared_genes
    
    # ---- Pearson correlation matrix ----
    cor_mat <- cor(logfc_mat, method = "pearson", use = "pairwise.complete.obs")
    
    write.csv(cor_mat, file.path(OUTPUT_BASE_DIR, "comparison_logFC_correlation_matrix.csv"),
              row.names = TRUE)
    
    # Correlation heatmap
    melted <- reshape2::melt(cor_mat)
    p_cor <- ggplot(melted, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "C", name = "Pearson r",
                           limits = c(min(melted$value, na.rm = TRUE), 1)) +
      geom_text(aes(label = round(value, 2)), size = 3, color = "white") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 13)) +
      labs(title = "Proteomics logFC Correlation Between Comparisons",
           subtitle = paste0(length(shared_genes), " shared genes"),
           x = NULL, y = NULL)
    
    ggsave(file.path(OUTPUT_BASE_DIR, "comparison_logFC_correlation_heatmap.png"),
           p_cor, width = 9, height = 8, dpi = 300, bg = "transparent")
    message("  ✓ Correlation heatmap saved")
    
    # ---- Dendrogram (Ward clustering by 1 - Pearson r) ----
    dist_mat <- as.dist(1 - cor_mat)
    hc <- hclust(dist_mat, method = "ward.D2")
    
    png(file.path(OUTPUT_BASE_DIR, "comparison_logFC_dendrogram.png"),
        width = 900, height = 600, res = 150)
    par(mar = c(5, 4, 4, 2))
    plot(hc, main = "Comparison Clustering (Ward, 1 - Pearson r on logFC)",
         xlab = "", sub = paste0(length(shared_genes), " shared genes"),
         hang = -1, cex = 0.9)
    dev.off()
    message("  ✓ Dendrogram saved")
  }
}, error = function(e) message("  ⚠ Cross-comparison correlation: ", e$message))

# ==============================================================================
# CROSS-COMPARISON: LINEAR REGRESSION SCATTER PLOTS
# ==============================================================================
# Individual scatter + regression per comparison vs C2 (primary),
# plus a combined faceted panel.

tryCatch({
  message("\n  Linear regression scatter plots...")
  
  # Use C2 as reference (primary comparison)
  ref_comp <- "C2"
  
  if (ref_comp %in% names(all_results)) {
    ref_tt <- all_results[[ref_comp]]$top_table %>%
      filter(!is.na(Gene) & Gene != "") %>%
      arrange(P.Value) %>%
      filter(!duplicated(Gene))
    ref_logfc <- setNames(ref_tt$logFC, ref_tt$Gene)
    ref_fdr   <- setNames(ref_tt$adj.P.Val, ref_tt$Gene)
    
    scatter_dir <- file.path(OUTPUT_BASE_DIR, "cross_comparison_scatters")
    dir.create(scatter_dir, showWarnings = FALSE)
    
    facet_data <- data.frame()
    
    other_comps <- setdiff(names(all_results), ref_comp)
    
    for (comp in other_comps) {
      tt <- all_results[[comp]]$top_table %>%
        filter(!is.na(Gene) & Gene != "") %>%
        arrange(P.Value) %>%
        filter(!duplicated(Gene))
      comp_logfc <- setNames(tt$logFC, tt$Gene)
      comp_fdr   <- setNames(tt$adj.P.Val, tt$Gene)
      
      shared <- intersect(names(ref_logfc), names(comp_logfc))
      if (length(shared) < 50) next
      
      df <- data.frame(
        Gene      = shared,
        ref_logFC = ref_logfc[shared],
        comp_logFC = comp_logfc[shared],
        ref_fdr   = ref_fdr[shared],
        comp_fdr  = comp_fdr[shared],
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          sig_both = ref_fdr < 0.05 & comp_fdr < 0.05,
          label = ifelse(sig_both & abs(ref_logFC) > quantile(abs(ref_logFC[sig_both]),
                                                              0.9, na.rm = TRUE),
                         Gene, "")
        )
      
      r_val <- cor(df$ref_logFC, df$comp_logFC, use = "complete.obs")
      comp_label <- comparisons[[comp]]$label
      
      # ---- Individual scatter ----
      p <- ggplot(df, aes(x = ref_logFC, y = comp_logFC)) +
        geom_point(aes(color = sig_both), alpha = 0.4, size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = COL_HEALTHY, linewidth = 0.8) +
        geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 20,
                        fontface = "bold", color = "black") +
        scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = COL_DCM),
                           labels = c("NS", "FDR<0.05 both"),
                           name = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        annotate("text", x = Inf, y = -Inf,
                 label = sprintf("r = %.3f\nn = %d genes", r_val, nrow(df)),
                 hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic",
                 color = "grey40") +
        labs(title = paste0(comp, " vs ", ref_comp, " — Protein logFC"),
             subtitle = comp_label,
             x = paste0("logFC (", ref_comp, ": ", comparisons[[ref_comp]]$label, ")"),
             y = paste0("logFC (", comp, ")")) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold", size = 12),
              plot.subtitle = element_text(size = 8, color = "grey40"),
              legend.position = "bottom")
      
      ggsave(file.path(scatter_dir, paste0("scatter_", comp, "_vs_", ref_comp, ".png")),
             p, width = 9, height = 8, dpi = 300, bg = "transparent")
      
      # Collect for faceted panel
      df$Comparison <- sprintf("%s (r=%.2f)", comp, r_val)
      facet_data <- rbind(facet_data, df)
    }
    
    message(sprintf("  ✓ %d individual scatter plots saved", length(other_comps)))
    
    # ---- Combined faceted panel ----
    if (nrow(facet_data) > 0) {
      p_facet <- ggplot(facet_data, aes(x = ref_logFC, y = comp_logFC)) +
        geom_point(aes(color = sig_both), alpha = 0.3, size = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = COL_HEALTHY, linewidth = 0.6) +
        scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = COL_DCM),
                           labels = c("NS", "FDR<0.05 both"), name = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
        facet_wrap(~ Comparison, ncol = 3, scales = "fixed") +
        labs(title = paste0("All Comparisons vs ", ref_comp, " — Protein logFC"),
             x = paste0("logFC (", ref_comp, ")"),
             y = "logFC (comparison)") +
        theme_minimal(base_size = 9) +
        theme(plot.title = element_text(face = "bold", size = 13),
              strip.text = element_text(face = "bold", size = 8),
              legend.position = "bottom")
      
      ggsave(file.path(OUTPUT_BASE_DIR, "comparison_ALL_vs_C2_faceted_scatter.png"),
             p_facet, width = 14, height = 12, dpi = 300, bg = "transparent")
      message("  ✓ Faceted scatter panel saved")
    }
  }
}, error = function(e) message("  ⚠ Scatter plots: ", e$message))

# ==============================================================================
# CORE DCM PROTEIN SIGNATURE
# ==============================================================================

tryCatch({
  sig_lists <- lapply(all_results[c("C1","C2","C3","C4")], function(r) {
    genes <- unique(r$top_table$Gene[r$top_table$adj.P.Val < 0.05 &
                                       !is.na(r$top_table$Gene) & r$top_table$Gene != ""])
  })
  
  core_genes <- Reduce(intersect, sig_lists)
  message(sprintf("\n  Core DCM proteins (FDR<0.05 in C1–C4): %d", length(core_genes)))
  
  if (length(core_genes) > 0) {
    c2_tt <- all_results[["C2"]]$top_table
    core_df <- c2_tt %>%
      filter(Gene %in% core_genes) %>%
      arrange(adj.P.Val) %>%
      filter(!duplicated(Gene)) %>%
      select(Gene, Protein, logFC, AveExpr, t, P.Value, adj.P.Val)
    
    write.csv(core_df, file.path(OUTPUT_BASE_DIR, "core_DCM_proteins.csv"),
              row.names = FALSE)
    
    cat("    ", paste(head(sort(core_genes), 40), collapse = ", "), "\n")
    if (length(core_genes) > 40) cat("    ... and", length(core_genes) - 40, "more\n")
    
    n_core_up   <- sum(core_df$logFC > 0)
    n_core_down <- sum(core_df$logFC < 0)
    cat(sprintf("    Direction in C2: %d up / %d down in DCM\n", n_core_up, n_core_down))
  }
  
  # Gene overlap matrix
  sig_all <- lapply(all_results, function(r) {
    unique(r$top_table$Gene[r$top_table$adj.P.Val < 0.05 &
                              !is.na(r$top_table$Gene) & r$top_table$Gene != ""])
  })
  all_genes <- unique(unlist(sig_all))
  if (length(all_genes) > 0 && length(sig_all) >= 2) {
    overlap <- sapply(sig_all, function(g) as.integer(all_genes %in% g))
    rownames(overlap) <- all_genes
    write.csv(overlap, file.path(OUTPUT_BASE_DIR, "benchmark_gene_overlap.csv"),
              row.names = TRUE)
    message("  ✓ Gene overlap matrix saved")
  }
  
}, error = function(e) message("  ⚠ Core gene analysis: ", e$message))

# ==============================================================================
# RECOMMENDATION
# ==============================================================================

best_idx   <- which.max(bench_df$FDR05)
best_label <- bench_df$Comparison[best_idx]

cat("\n================================================================================\n")
cat("  RECOMMENDATION\n")
cat("================================================================================\n")
cat(sprintf("  Best comparison at FDR < 0.05:\n"))
cat(sprintf("    → %s\n", best_label))
cat(sprintf("    → %d significant (%d up / %d down) = %.1f%%\n",
            bench_df$FDR05[best_idx], bench_df$Up_FDR05[best_idx],
            bench_df$Down_FDR05[best_idx], bench_df$pct_sig[best_idx]))
cat(sprintf("    → %d Disease vs %d Healthy\n",
            bench_df$nD[best_idx], bench_df$nH[best_idx]))
cat(sprintf("\n  Per-comparison outputs include:\n"))
cat(sprintf("    • QC: boxplots, density, correlation heatmap\n"))
cat(sprintf("    • PCA/MDS with outlier detection\n"))
cat(sprintf("    • Volcano plot (png + pdf)\n"))
cat(sprintf("    • DE tables at multiple thresholds\n"))
cat(sprintf("    • clusterProfiler: GO BP (all/up/down), KEGG, Reactome + dotplots\n"))
cat(sprintf("    • fgsea GSEA: Hallmarks, KEGG, Reactome\n"))
cat(sprintf("\n  Cross-comparison outputs:\n"))
cat(sprintf("    • comparison_logFC_correlation_matrix.csv  — Pearson r between all pairs\n"))
cat(sprintf("    • comparison_logFC_correlation_heatmap.png — correlation heatmap\n"))
cat(sprintf("    • comparison_logFC_dendrogram.png          — Ward clustering\n"))
cat(sprintf("    • cross_comparison_scatters/               — individual regression plots\n"))
cat(sprintf("    • comparison_ALL_vs_C2_faceted_scatter.png — combined faceted panel\n"))
cat(sprintf("\n  Next steps:\n"))
cat(sprintf("    1. Use proteomics for host protein normalization of phosphoproteomics\n"))
cat(sprintf("    2. Compare core DCM protein vs phosphosite signatures\n"))
cat(sprintf("    3. Integrate with Reitz et al. adult DCM proteomics\n"))
cat(sprintf("    4. Cross-reference with NRG1/ERBB in vitro timecourse\n"))
cat("================================================================================\n\n")

message("\n✓ Benchmarking complete. All outputs in: ", OUTPUT_BASE_DIR, "/")
