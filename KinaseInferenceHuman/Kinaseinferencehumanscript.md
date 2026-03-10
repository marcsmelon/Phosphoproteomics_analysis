################################################################################
# HUMAN DCM PHOSPHOPROTEOMICS BENCHMARK PIPELINE (v3 — DCM vs Healthy focus)
# Adapted from: Script_with_NRG1_ERBB_analysis_kinases_update_protein_normalisation14.R
#
# Purpose: Find the best significance for DCM vs Healthy comparison across
#          multiple grouping strategies, for integration with in vitro NRG1/ERBB data.
#
# Cohort: Pediatric DCM hearts — 66 intensity columns
#   - 14 Healthy donors       (-D suffix)
#   - 28 VAD disease samples  (-V suffix)
#   - 24 Tx disease samples   (-T suffix)
#   - 9 patients have BOTH -V and -T columns → need one-per-patient selection
#
# Pseudo-replication handling:
#   9 paired patients have two columns each (-V and -T). Since both represent
#   "disease" state of the SAME patient, including both would inflate n and
#   violate independence. For each comparison we select ONE column per patient.
#   The script tests both preferences (V-preferred, T-preferred) to see which
#   gives better results.
#
# Comparisons:
#   C1: DCM_filtered (13) vs Healthy — your curated metadata subset
#   C2: All DCM one-per-patient, prefer V (35 DCM) vs Healthy — max power
#   C3: All DCM one-per-patient, prefer T (35 DCM) vs Healthy — alternative
#   C4: All DCM naïve (44 columns, includes duplicates) vs Healthy — for reference
#   C5: Sarcomeric mutations vs Healthy
#   C6: Cytoskeletal/structural vs Healthy
#   C7: Idiopathic DCM vs Healthy
#   C8: All disease incl. myocarditis, one-per-patient (43) vs Healthy
#
# DATA NOTE: Phospho file >32 MB — run this locally pointing at the full file.
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PHOSPHO_DATA_PATH <- "/Users/marcos.sandemelon/Downloads/analysis human/fullCohortPhos_PivotTable_extra.tsv"       # <-- EDIT: your phospho export path
PROTEIN_DATA_PATH <- NULL                       # <-- proteinGroups.txt or NULL
setwd("/Users/marcos.sandemelon/Downloads/analysis human")       # <-- Uncomment & set

# ==============================================================================
# PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(scales)
  library(viridis)
  library(ggrepel)
})

# ==============================================================================
# VERIFIED COLUMN MAPPINGS
# ==============================================================================

# ---- 14 Healthy donor controls (-D) ----
healthy_cols <- c(
  "166-D", "246-D", "263-D", "270-D", "277-D", "285-D", "345-D",
  "382-D", "391-D", "450-D", "460-D", "480-D", "493-D", "508-D"
)

# ---- 13 DCM_filtered (curated from metadata) ----
dcm_filtered_cols <- c(
  "120-V",   # MCHTB120  | vad    | TPM1        | Sarcomeric
  "0091-V",  # 0091      | vad    | Idiopathic   | Idiopathic
  "63-T",    # 96/MCHTB63| paired | TNNT2       | Sarcomeric
  "272-V",   # MCHTB272  | paired | Idiopathic   | Idiopathic
  "309-V",   # MCHTB309  | vad    | TNNI3       | Sarcomeric
  "328-V",   # MCHTB328  | vad    | PPP1R13L    | Cytoskeletal
  "0062-V",  # 0062      | vad    | Idiopathic   | Idiopathic
  "0021-V",  # 0021      | vad    | NEXN (vus)  | Cytoskeletal
  "288-V",   # MCHTB288  | paired | Idiopathic   | Idiopathic
  "0019-V",  # 0019      | vad    | DSP         | Cytoskeletal
  "468-V",   # MCHTB468  | paired | Idiopathic   | Idiopathic
  "499-V",   # MCHTB499  | paired | DSP         | Cytoskeletal
  "545-T"    # MCHTB545  | tx     | TTNI3K      | Cytoskeletal
)

# ---- 9 patients with BOTH -V and -T ----
# For these, we pick ONE per patient depending on the comparison strategy.
paired_patients <- list(
  MCHTB499 = list(V = "499-V",  T = "499-T",  ds = "DCM", cause = "DSP",       group = "Cytoskeletal"),
  MCHTB468 = list(V = "468-V",  T = "468-T",  ds = "DCM", cause = "Idiopathic", group = "Idiopathic"),
  MCHTB288 = list(V = "288-V",  T = "288-T",  ds = "DCM", cause = "Idiopathic", group = "Idiopathic"),
  MCHTB426 = list(V = "426-V",  T = "426-T",  ds = "DCM", cause = "Idiopathic", group = "Idiopathic"),
  MCHTB53  = list(V = "53-V",   T = "53-T",   ds = "Myocarditis", cause = "",   group = "Myocarditis"),
  d0050    = list(V = "0050-V", T = "0050-T", ds = "DCM", cause = "DSP",        group = "Cytoskeletal"),
  MCHTB154 = list(V = "154-V",  T = "154-T",  ds = "DCM", cause = "ACTC1",      group = "Sarcomeric"),
  MCHTB172 = list(V = "172-V",  T = "172-T",  ds = "DCM", cause = "TAZ (vus)",  group = "Cytoskeletal"),
  MCHTB272 = list(V = "272-V",  T = "272-T",  ds = "DCM", cause = "Idiopathic", group = "Idiopathic")
)

# ---- VAD-only patients (no matching T column) ----
vad_only_cols <- c(
  "0019-V",  # DSP            | Cytoskeletal
  "0020-V",  # Idiopathic     | Idiopathic
  "0021-V",  # NEXN (vus)     | Cytoskeletal
  "0034-V",  # TTN            | Cytoskeletal
  "0039-V",  # Idiopathic     | Idiopathic
  "0054-V",  # Myocarditis
  "0062-V",  # Idiopathic     | Idiopathic
  "0077-V",  # Idiopathic     | Idiopathic
  "0091-V",  # Idiopathic     | Idiopathic
  "0096-V",  # Unknown
  "104-V",   # Idiopathic     | Idiopathic
  "120-V",   # TPM1           | Sarcomeric
  "148-V",   # Myocarditis
  "180-V",   # Myocarditis
  "309-V",   # TNNI3          | Sarcomeric
  "328-V",   # PPP1R13L       | Cytoskeletal
  "395-V",   # Idiopathic     | Idiopathic
  "72-V",    # Myocarditis
  "79-V"     # Idiopathic     | Idiopathic
)

# ---- Tx-only patients (no matching V column) ----
tx_only_cols <- c(
  "0025-T",  # Idiopathic     | Idiopathic
  "0029-T",  # Idiopathic     | Idiopathic
  "0031-T",  # Myocarditis
  "0053-T",  # Myocarditis    (=MCHTB53 but only if we didn't pick 53-V)
  "0058-T",  # Unknown
  "223-T",   # NEXN           | Cytoskeletal
  "262-T",   # Idiopathic     | Idiopathic
  "269-T",   # RBM20          | Cytoskeletal
  "344-T",   # TAZ            | Cytoskeletal
  "367-T",   # MYH7           | Sarcomeric
  "449-T",   # MYH7           | Sarcomeric
  "479-T",   # JPH2           | Cytoskeletal
  "490-T",   # TPM1 (vus 3A)  | Sarcomeric
  "545-T",   # TTNI3K         | Cytoskeletal
  "63-T"     # TNNT2          | Sarcomeric
)

# ==============================================================================
# BUILD ONE-PER-PATIENT SETS
# ==============================================================================

# Strategy A: prefer V (VAD = acute disease, before mechanical unloading)
one_per_patient_preferV <- c(
  vad_only_cols,
  tx_only_cols,
  sapply(paired_patients, function(p) p$V)  # pick V for all paired
)
# 0034-V is already in vad_only; 63-T already in tx_only
one_per_patient_preferV <- unique(one_per_patient_preferV)

# Strategy B: prefer T (transplant = end-stage)
one_per_patient_preferT <- c(
  vad_only_cols,
  tx_only_cols,
  sapply(paired_patients, function(p) p$T)  # pick T for all paired
)
one_per_patient_preferT <- unique(one_per_patient_preferT)

# ---- Filter to DCM-only (exclude myocarditis + unknown) ----
myocarditis_cols <- c("0054-V", "148-V", "180-V", "72-V",    # VAD-only myocarditis
                      "0031-T", "0053-T",                      # Tx-only myocarditis
                      "53-V", "53-T")                          # Paired myocarditis (MCHTB53)
unknown_cols     <- c("0058-T", "0096-V")

dcm_one_per_patient_V <- setdiff(one_per_patient_preferV,
                                  c(myocarditis_cols, unknown_cols))
dcm_one_per_patient_T <- setdiff(one_per_patient_preferT,
                                  c(myocarditis_cols, unknown_cols))

# ---- ALL DCM columns naïve (includes both V+T for paired patients) ----
all_vad <- c("0019-V","0020-V","0021-V","0034-V","0039-V","0050-V","0054-V",
             "0062-V","0077-V","0091-V","0096-V","104-V","120-V","148-V",
             "154-V","172-V","180-V","272-V","288-V","309-V","328-V","395-V",
             "426-V","468-V","499-V","53-V","72-V","79-V")
all_tx  <- c("0025-T","0029-T","0031-T","0050-T","0053-T","0058-T","154-T",
             "172-T","223-T","262-T","269-T","272-T","288-T","344-T","367-T",
             "426-T","449-T","468-T","479-T","490-T","499-T","53-T","545-T","63-T")
all_disease_naive <- c(all_vad, all_tx)
all_dcm_naive <- setdiff(all_disease_naive, c(myocarditis_cols, unknown_cols))

# ---- Genetic subgroups (one-per-patient, prefer V) ----
# Annotation: column -> gene_group
col_group <- c(
  # Sarcomeric: MYH7, TNNT2, ACTC1, TPM1, TNNI3, MYL2/MYH6
  "120-V"  = "Sarcomeric", "309-V"  = "Sarcomeric", "63-T"   = "Sarcomeric",
  "154-V"  = "Sarcomeric", "154-T"  = "Sarcomeric",
  "367-T"  = "Sarcomeric", "449-T"  = "Sarcomeric", "490-T"  = "Sarcomeric",
  # Cytoskeletal/structural: DSP, TTN, TAZ, RBM20, NEXN, JPH2, TTNI3K, PPP1R13L
  "0019-V" = "Cytoskeletal", "0021-V" = "Cytoskeletal", "0034-V" = "Cytoskeletal",
  "0050-V" = "Cytoskeletal", "0050-T" = "Cytoskeletal",
  "172-V"  = "Cytoskeletal", "172-T"  = "Cytoskeletal",
  "223-T"  = "Cytoskeletal", "269-T"  = "Cytoskeletal", "328-V"  = "Cytoskeletal",
  "344-T"  = "Cytoskeletal", "479-T"  = "Cytoskeletal", "499-V"  = "Cytoskeletal",
  "499-T"  = "Cytoskeletal", "545-T"  = "Cytoskeletal",
  # Idiopathic
  "0020-V" = "Idiopathic", "0039-V" = "Idiopathic", "0062-V" = "Idiopathic",
  "0077-V" = "Idiopathic", "0091-V" = "Idiopathic", "104-V"  = "Idiopathic",
  "272-V"  = "Idiopathic", "272-T"  = "Idiopathic", "288-V"  = "Idiopathic",
  "288-T"  = "Idiopathic", "395-V"  = "Idiopathic", "426-V"  = "Idiopathic",
  "426-T"  = "Idiopathic", "468-V"  = "Idiopathic", "468-T"  = "Idiopathic",
  "79-V"   = "Idiopathic", "0025-T" = "Idiopathic", "0029-T" = "Idiopathic",
  "262-T"  = "Idiopathic"
)

# For subgroups, use one-per-patient (prefer V)
sarcomeric_1pp   <- intersect(dcm_one_per_patient_V, names(col_group[col_group == "Sarcomeric"]))
cytoskeletal_1pp <- intersect(dcm_one_per_patient_V, names(col_group[col_group == "Cytoskeletal"]))
idiopathic_1pp   <- intersect(dcm_one_per_patient_V, names(col_group[col_group == "Idiopathic"]))

# ---- Print setup summary ----
cat("\n")
cat("================================================================================\n")
cat("  HUMAN DCM PHOSPHOPROTEOMICS BENCHMARK (v3 — DCM vs Healthy)\n")
cat("================================================================================\n\n")
cat(sprintf("  Healthy donors                     : %d\n", length(healthy_cols)))
cat(sprintf("  DCM_filtered (curated)             : %d\n", length(dcm_filtered_cols)))
cat(sprintf("  All DCM one-per-patient (prefer V) : %d\n", length(dcm_one_per_patient_V)))
cat(sprintf("  All DCM one-per-patient (prefer T) : %d\n", length(dcm_one_per_patient_T)))
cat(sprintf("  All DCM naive (with duplicates)    : %d\n", length(all_dcm_naive)))
cat(sprintf("  Sarcomeric (1pp)                   : %d\n", length(sarcomeric_1pp)))
cat(sprintf("  Cytoskeletal (1pp)                 : %d\n", length(cytoskeletal_1pp)))
cat(sprintf("  Idiopathic (1pp)                   : %d\n", length(idiopathic_1pp)))
cat(sprintf("  All disease incl. myocarditis (1pp): %d\n",
            length(setdiff(one_per_patient_preferV, unknown_cols))))

# ==============================================================================
# CORE COMPARISON FUNCTION
# ==============================================================================

run_comparison <- function(SN_data,
                           disease_col_names,
                           healthy_col_names,
                           comparison_label,
                           filtering_threshold = 3,
                           imputation_method = "min",
                           normalization_method = "median_center",
                           protein_data = NULL,
                           protein_norm_threshold = 4,
                           output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("benchmark_", gsub("[^A-Za-z0-9_]", "_", comparison_label))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  %-70s║\n", comparison_label))
  cat(sprintf("║  Disease: %-3d  |  Healthy: %-3d  |  Total: %-3d %23s║\n",
              length(disease_col_names), length(healthy_col_names),
              length(disease_col_names) + length(healthy_col_names), ""))
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")
  
  # Verify columns exist
  disease_col_names <- disease_col_names[disease_col_names %in% colnames(SN_data)]
  healthy_col_names <- healthy_col_names[healthy_col_names %in% colnames(SN_data)]
  
  if (length(disease_col_names) < 2 || length(healthy_col_names) < 2) {
    message("  ✗ Too few samples. Skipping.")
    return(NULL)
  }
  
  all_sample_cols <- c(disease_col_names, healthy_col_names)
  
  # ---- GenePhos identifier ----
  if ("GenePhos" %in% colnames(SN_data)) {
    genephos <- SN_data$GenePhos
  } else if (all(c("PG.Genes", "PTM.CollapseKey") %in% colnames(SN_data))) {
    phos <- gsub("^[^_]*", "", SN_data$PTM.CollapseKey)
    genephos <- paste0(SN_data$PG.Genes, phos)
  } else {
    genephos <- paste0("site_", seq_len(nrow(SN_data)))
  }
  
  # ---- Build matrix ----
  intensity_matrix <- as.matrix(SN_data[, all_sample_cols, drop = FALSE])
  intensity_matrix[intensity_matrix == "NaN" | intensity_matrix == ""] <- NA
  storage.mode(intensity_matrix) <- "numeric"
  rownames(intensity_matrix) <- genephos
  
  dups <- duplicated(rownames(intensity_matrix))
  if (any(dups)) {
    message(sprintf("  Removing %d duplicate GenePhos", sum(dups)))
    intensity_matrix <- intensity_matrix[!dups, , drop = FALSE]
  }
  
  condition <- factor(
    ifelse(colnames(intensity_matrix) %in% disease_col_names, "Disease", "Healthy"),
    levels = c("Healthy", "Disease")
  )
  
  # ---- Filtering ----
  if (is.numeric(filtering_threshold)) {
    is_valid <- rowSums(!is.na(intensity_matrix) & intensity_matrix > 0) >= filtering_threshold
  } else if (filtering_threshold == "by_condition") {
    h_idx <- which(condition == "Healthy")
    d_idx <- which(condition == "Disease")
    is_valid <- (rowSums(!is.na(intensity_matrix[, h_idx, drop=FALSE]) &
                           intensity_matrix[, h_idx, drop=FALSE] > 0) >= max(2, floor(length(h_idx)*0.5))) &
                (rowSums(!is.na(intensity_matrix[, d_idx, drop=FALSE]) &
                           intensity_matrix[, d_idx, drop=FALSE] > 0) >= max(2, floor(length(d_idx)*0.5)))
  } else {
    is_valid <- rep(TRUE, nrow(intensity_matrix))
  }
  
  intensity_matrix <- intensity_matrix[is_valid, , drop = FALSE]
  message(sprintf("  After filtering: %d phosphosites", nrow(intensity_matrix)))
  
  # ---- Imputation ----
  intensity_matrix[intensity_matrix == 0] <- NA
  if (imputation_method == "min") {
    for (i in seq_len(nrow(intensity_matrix))) {
      na_mask <- is.na(intensity_matrix[i, ])
      if (all(na_mask)) next
      intensity_matrix[i, na_mask] <- min(intensity_matrix[i, !na_mask]) / 2
    }
  }
  
  # ---- Log2 transform ----
  min_nz <- min(intensity_matrix[intensity_matrix > 0], na.rm = TRUE)
  offset <- ifelse(min_nz < 1, 1, 0)
  log2_mat <- log2(intensity_matrix + offset)
  
  # ---- Optional: host protein normalization ----
  if (!is.null(protein_data)) {
    tryCatch({
      prot_norm <- normalize_phospho_to_protein(
        log2_phospho_matrix    = log2_mat,
        protein_data           = protein_data,
        protein_norm_threshold = protein_norm_threshold
      )
      log2_mat <- prot_norm$normalized_matrix
      message(sprintf("  ✓ Protein norm: %d matched", prot_norm$n_matched))
    }, error = function(e) message("  ⚠ Protein norm failed: ", e$message))
  }
  
  # ---- Median centering ----
  if (normalization_method == "median_center") {
    norm_mat <- apply(log2_mat, 2, function(y) y - median(y, na.rm = TRUE))
  } else if (normalization_method == "quantile") {
    temp_elist <- new("EList", list(E = 2^log2_mat))
    temp_elist <- limma::normalizeBetweenArrays(temp_elist, method = "quantile")
    norm_mat <- log2(temp_elist$E + offset)
  } else {
    norm_mat <- log2_mat
  }
  
  # ---- limma: Disease vs Healthy ----
  design   <- model.matrix(~0 + condition)
  colnames(design) <- levels(condition)
  rownames(design) <- colnames(norm_mat)
  contrast <- limma::makeContrasts(Disease - Healthy, levels = design)
  
  fit  <- limma::lmFit(norm_mat, design = design)
  fit  <- limma::contrasts.fit(fit, contrast)
  fit  <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  tt <- limma::topTable(fit, n = Inf, coef = 1)
  tt$GenePhos <- rownames(tt)
  tt$Gene     <- gsub("_[STY].*$", "", tt$GenePhos)
  tt <- tt[, c("GenePhos", "Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  
  # ---- Significance tiers ----
  tiers <- data.frame(
    Comparison      = comparison_label,
    N_disease       = length(disease_col_names),
    N_healthy       = length(healthy_col_names),
    N_total_tested  = nrow(tt),
    FDR_0.01        = sum(tt$adj.P.Val < 0.01, na.rm = TRUE),
    FDR_0.05        = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
    FDR_0.1         = sum(tt$adj.P.Val < 0.1,  na.rm = TRUE),
    FDR_0.2         = sum(tt$adj.P.Val < 0.2,  na.rm = TRUE),
    rawP_0.01       = sum(tt$P.Value   < 0.01, na.rm = TRUE),
    rawP_0.05       = sum(tt$P.Value   < 0.05, na.rm = TRUE),
    Up_FDR0.05      = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0, na.rm = TRUE),
    Down_FDR0.05    = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0, na.rm = TRUE),
    Pct_FDR0.05     = round(100 * sum(tt$adj.P.Val < 0.05, na.rm=TRUE) / nrow(tt), 2),
    stringsAsFactors = FALSE
  )
  
  # ---- Save results ----
  write.csv(tt,    file.path(output_dir, "complete_statistics.csv"), row.names = FALSE)
  write.csv(tiers, file.path(output_dir, "significance_tiers.csv"), row.names = FALSE)
  write.csv(tt[tt$adj.P.Val < 0.05, ], file.path(output_dir, "sig_FDR0.05.csv"), row.names = FALSE)
  write.csv(tt[tt$adj.P.Val < 0.1, ],  file.path(output_dir, "sig_FDR0.1.csv"),  row.names = FALSE)
  write.csv(tt[tt$adj.P.Val < 0.2, ],  file.path(output_dir, "sig_FDR0.2.csv"),  row.names = FALSE)
  write.csv(tt[tt$P.Value < 0.01, ],   file.path(output_dir, "sig_rawP0.01.csv"), row.names = FALSE)
  
  # ---- Volcano plot ----
  tryCatch({
    tt_plot <- tt %>%
      dplyr::mutate(sig_group = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 1  ~ "Up (FDR<0.05, FC>2)",
        adj.P.Val < 0.05 & logFC < -1 ~ "Down (FDR<0.05, FC>2)",
        adj.P.Val < 0.05              ~ "Sig (FDR<0.05)",
        P.Value   < 0.01              ~ "Nominal (p<0.01)",
        TRUE                          ~ "NS"))
    
    top_lab <- tt_plot %>% dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = 25)
    
    p <- ggplot2::ggplot(tt_plot, ggplot2::aes(x = logFC, y = -log10(P.Value), color = sig_group)) +
      ggplot2::geom_point(alpha = 0.6, size = 1.2) +
      ggplot2::scale_color_manual(values = c("Up (FDR<0.05, FC>2)" = "#E63946",
                                              "Down (FDR<0.05, FC>2)" = "#2E86AB",
                                              "Sig (FDR<0.05)" = "#F4A261",
                                              "Nominal (p<0.01)" = "#A8DADC",
                                              "NS" = "grey80"), name = NULL) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = comparison_label,
                    subtitle = sprintf("n=%d | DCM:%d | Healthy:%d | FDR<0.05: %d (%d↑ %d↓)",
                                       nrow(tt), length(disease_col_names),
                                       length(healthy_col_names), tiers$FDR_0.05,
                                       tiers$Up_FDR0.05, tiers$Down_FDR0.05),
                    x = "log2 FC (DCM / Healthy)", y = "-log10(p-value)") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13),
                     plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
                     legend.position = "bottom")
    if (nrow(top_lab) > 0) {
      p <- p + ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = Gene),
                                         size = 2.8, max.overlaps = 15, segment.size = 0.3)
    }
    ggplot2::ggsave(file.path(output_dir, "volcano.png"), p, width = 10, height = 8, dpi = 300)
  }, error = function(e) message("  ⚠ Volcano failed: ", e$message))
  
  # ---- P-value histogram ----
  tryCatch({
    ph <- ggplot2::ggplot(tt, ggplot2::aes(x = P.Value)) +
      ggplot2::geom_histogram(bins = 50, fill = "#457B9D", color = "white", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste0("P-value distribution: ", comparison_label),
                    subtitle = "Enrichment near 0 = true differential signal",
                    x = "Raw p-value", y = "Count") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(file.path(output_dir, "pvalue_histogram.png"), ph, width = 8, height = 5, dpi = 300)
  }, error = function(e) NULL)
  
  # ---- Console summary ----
  cat(sprintf("\n  %s\n", comparison_label))
  cat(sprintf("  Tested: %d | FDR<0.01: %d | FDR<0.05: %d (%d↑ %d↓) | FDR<0.1: %d | rawP<0.01: %d\n\n",
              nrow(tt), tiers$FDR_0.01, tiers$FDR_0.05,
              tiers$Up_FDR0.05, tiers$Down_FDR0.05,
              tiers$FDR_0.1, tiers$rawP_0.01))
  
  return(list(top_table = tt, fit = fit, tiers = tiers,
              normalized_mat = norm_mat, condition = condition,
              comparison_label = comparison_label, output_dir = output_dir))
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\nLoading: ", PHOSPHO_DATA_PATH)
SN <- read.csv(PHOSPHO_DATA_PATH, header = TRUE, sep = "\t",
               check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(SN) > 0 && any(grepl("^X", colnames(SN)))) {
  colnames(SN) <- gsub("^X", "", colnames(SN))
}
message(sprintf("  ✓ %d rows x %d columns", nrow(SN), ncol(SN)))

# Quick column sanity check
expected <- c(healthy_cols, dcm_filtered_cols)
missing  <- expected[!expected %in% colnames(SN)]
if (length(missing) > 0) {
  warning("MISSING: ", paste(missing, collapse = ", "))
} else {
  message("  ✓ All 27 core columns verified")
}

# ==============================================================================
# DEFINE COMPARISONS
# ==============================================================================

comparisons <- list(
  
  # C1: Curated DCM_filtered from metadata (13 vs 14)
  C1 = list(d = dcm_filtered_cols, h = healthy_cols,
            label = "C1: DCM_filtered curated (13 vs 14)",
            dir   = "C1_DCM_filtered"),
  
  # C2: All DCM, one-per-patient, prefer V (35 vs 14) — MAX POWER, CLEAN
  C2 = list(d = dcm_one_per_patient_V, h = healthy_cols,
            label = "C2: All DCM 1-per-patient prefer-V (35 vs 14)",
            dir   = "C2_AllDCM_1pp_preferV"),
  
  # C3: All DCM, one-per-patient, prefer T (35 vs 14)
  C3 = list(d = dcm_one_per_patient_T, h = healthy_cols,
            label = "C3: All DCM 1-per-patient prefer-T (35 vs 14)",
            dir   = "C3_AllDCM_1pp_preferT"),
  
  # C4: All DCM naive — includes both V+T for 9 patients (44 vs 14)
  #     NOTE: Pseudo-replicated! Use as reference only.
  C4 = list(d = all_dcm_naive, h = healthy_cols,
            label = "C4: All DCM naive incl. duplicates (44 vs 14) [REFERENCE]",
            dir   = "C4_AllDCM_naive"),
  
  # C5: Sarcomeric mutations only
  C5 = list(d = sarcomeric_1pp, h = healthy_cols,
            label = "C5: Sarcomeric mutations (vs 14 Healthy)",
            dir   = "C5_Sarcomeric"),
  
  # C6: Cytoskeletal/structural mutations only
  C6 = list(d = cytoskeletal_1pp, h = healthy_cols,
            label = "C6: Cytoskeletal/structural (vs 14 Healthy)",
            dir   = "C6_Cytoskeletal"),
  
  # C7: Idiopathic DCM only
  C7 = list(d = idiopathic_1pp, h = healthy_cols,
            label = "C7: Idiopathic DCM (vs 14 Healthy)",
            dir   = "C7_Idiopathic"),
  
  # C8: All disease incl. myocarditis, one-per-patient (43 vs 14)
  C8 = list(d = setdiff(one_per_patient_preferV, unknown_cols), h = healthy_cols,
            label = "C8: All disease incl. myocarditis 1pp (43 vs 14)",
            dir   = "C8_AllDisease_1pp")
)

# Drop comparisons with < 3 disease samples
comparisons <- Filter(function(x) length(x$d) >= 3, comparisons)

cat(sprintf("\n  Running %d comparisons:\n\n", length(comparisons)))
for (nm in names(comparisons)) {
  pseudo_flag <- if (grepl("naive", comparisons[[nm]]$label)) " ⚠ pseudo-replicated" else ""
  cat(sprintf("  %-4s %-60s (%d vs %d)%s\n",
              nm, comparisons[[nm]]$label,
              length(comparisons[[nm]]$d), length(comparisons[[nm]]$h),
              pseudo_flag))
}

# ==============================================================================
# RUN ALL COMPARISONS
# ==============================================================================

all_results <- list()
all_tiers   <- list()

for (nm in names(comparisons)) {
  comp <- comparisons[[nm]]
  res <- tryCatch(
    run_comparison(SN, comp$d, comp$h, comp$label,
                   filtering_threshold = 3,
                   imputation_method = "min",
                   normalization_method = "median_center",
                   protein_data = PROTEIN_DATA_PATH,
                   output_dir = comp$dir),
    error = function(e) { message("  ✗ FAILED: ", nm, " -- ", e$message); NULL }
  )
  if (!is.null(res)) {
    all_results[[nm]] <- res
    all_tiers[[nm]]   <- res$tiers
  }
}

# ==============================================================================
# BENCHMARK SUMMARY
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                   BENCHMARK: DCM vs HEALTHY SIGNIFICANCE                       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════════╝\n\n")

if (length(all_tiers) > 0) {
  bench <- dplyr::bind_rows(all_tiers)
  write.csv(bench, "benchmark_summary.csv", row.names = FALSE)
  
  bench_sorted <- bench %>% dplyr::arrange(desc(FDR_0.05))
  
  # Print compact table
  cat(sprintf("%-62s %3s %3s %5s %5s %5s %5s %5s %5s  %5s\n",
              "Comparison", "nD", "nH", "Total", "FDR01", "FDR05", "FDR10", "FDR20", "pv01", "%sig"))
  cat(paste(rep("─", 115), collapse = ""), "\n")
  for (i in seq_len(nrow(bench_sorted))) {
    r <- bench_sorted[i, ]
    cat(sprintf("%-62s %3d %3d %5d %5d %5d %5d %5d %5d  %5.1f\n",
                r$Comparison, r$N_disease, r$N_healthy, r$N_total_tested,
                r$FDR_0.01, r$FDR_0.05, r$FDR_0.1, r$FDR_0.2, r$rawP_0.01, r$Pct_FDR0.05))
  }
  
  # ---- Benchmark barplot ----
  tryCatch({
    plot_df <- bench %>%
      tidyr::pivot_longer(cols = c(FDR_0.01, FDR_0.05, FDR_0.1, rawP_0.01),
                          names_to = "Threshold", values_to = "N_sig") %>%
      dplyr::mutate(
        Comp_short = gsub("^C\\d+: ", "", Comparison),
        Threshold = factor(Threshold,
                           levels = c("FDR_0.01", "FDR_0.05", "FDR_0.1", "rawP_0.01"),
                           labels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1", "raw p < 0.01")))
    
    p_bench <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(x = reorder(Comp_short, N_sig),
                                            y = N_sig, fill = Threshold)) +
      ggplot2::geom_col(position = "dodge", alpha = 0.85) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("FDR < 0.01" = "#264653",
                                             "FDR < 0.05" = "#2A9D8F",
                                             "FDR < 0.1"  = "#E9C46A",
                                             "raw p < 0.01" = "#E76F51")) +
      ggplot2::labs(title = "Benchmark: DCM vs Healthy — significant phosphosites",
                    subtitle = "Multiple grouping strategies | one-per-patient where applicable",
                    x = NULL, y = "N significant phosphosites", fill = "Threshold") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
                     legend.position = "bottom") +
      ggplot2::geom_text(ggplot2::aes(label = N_sig),
                         position = ggplot2::position_dodge(width = 0.9),
                         hjust = -0.1, size = 2.8)
    
    ggplot2::ggsave("benchmark_barplot.png", p_bench, width = 14, height = 9, dpi = 300)
    message("  ✓ Benchmark barplot saved")
  }, error = function(e) message("  ⚠ Barplot failed: ", e$message))
  
  # ---- Gene overlap (FDR < 0.05) ----
  tryCatch({
    sig_lists <- lapply(all_results, function(r)
      unique(r$top_table$Gene[r$top_table$adj.P.Val < 0.05]))
    all_genes <- unique(unlist(sig_lists))
    
    if (length(all_genes) > 0 && length(sig_lists) >= 2) {
      overlap <- sapply(sig_lists, function(g) as.integer(all_genes %in% g))
      rownames(overlap) <- all_genes
      write.csv(overlap, "benchmark_gene_overlap.csv", row.names = TRUE)
      
      # Core genes: significant in ALL DCM-only comparisons (C1–C4)
      dcm_comps <- intersect(c("C1","C2","C3","C4"), colnames(overlap))
      if (length(dcm_comps) >= 2) {
        core <- all_genes[rowSums(overlap[, dcm_comps, drop=FALSE]) == length(dcm_comps)]
        if (length(core) > 0) {
          cat(sprintf("\n  Core DCM genes (FDR<0.05 in C1–C4): %d\n", length(core)))
          cat(sprintf("    %s\n", paste(head(sort(core), 40), collapse = ", ")))
          write.csv(data.frame(Gene = sort(core)),
                    "core_DCM_genes.csv", row.names = FALSE)
        }
      }
    }
  }, error = function(e) NULL)
  
  # ---- Recommendation ----
  # Exclude the pseudo-replicated C4 from "best" selection
  bench_clean <- bench_sorted %>% dplyr::filter(!grepl("REFERENCE", Comparison))
  best <- bench_clean[1, ]
  
  cat("\n")
  cat("================================================================================\n")
  cat("  RECOMMENDATION\n")
  cat("================================================================================\n\n")
  cat(sprintf("  Best clean comparison at FDR < 0.05:\n"))
  cat(sprintf("    → %s\n", best$Comparison))
  cat(sprintf("    → %d significant (%d up / %d down) = %.1f%%\n",
              best$FDR_0.05, best$Up_FDR0.05, best$Down_FDR0.05, best$Pct_FDR0.05))
  cat(sprintf("    → %d DCM vs %d Healthy\n\n", best$N_disease, best$N_healthy))
  
  # Flag if naive C4 beats 1pp — indicates pseudo-replication boost
  if ("C4" %in% names(all_tiers)) {
    c4_sig <- all_tiers[["C4"]]$FDR_0.05
    if (c4_sig > best$FDR_0.05) {
      cat(sprintf("  ⚠ NOTE: Naïve C4 (%d sig) beats the best 1pp (%d sig).\n",
                  c4_sig, best$FDR_0.05))
      cat("    This is expected — pseudo-replication inflates significance.\n")
      cat("    Use the 1pp result for honest statistics.\n\n")
    }
  }
  
  cat("  Next steps for integration with NRG1/ERBB in vitro timecourse:\n")
  cat("    1. Source the main pipeline:  source('Script_with_NRG1_ERBB_..._14.R')\n")
  cat("    2. Run ERBB analysis:         analyze_erbb_phosphosites_final(fit2 = best$fit)\n")
  cat("    3. Run kinase activity:        run_kinase_activity(results = best_result)\n")
  cat("    4. Compare kinase profiles between human DCM and in vitro timecourse\n")
  cat("    5. Export for phosphonetworks: export_for_phosphonetworks.R\n\n")
  
} else {
  message("  No comparisons completed.")
}

cat("================================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("================================================================================\n\n")

# ==============================================================================
# DOWNSTREAM (uncomment after reviewing benchmark_summary.csv)
# ==============================================================================
#
# source("Script_with_NRG1_ERBB_analysis_kinases_update_protein_normalisation14.R")
#
# best_result <- all_results[["C2"]]  # or whichever won
#
# # Wrap in the format expected by downstream functions
# results_for_pipeline <- list(
#   top_table = best_result$top_table,
#   fit       = best_result$fit
# )
#
# # Kinase activity inference (5 methods)
# kinase_results <- run_kinase_activity(
#   results        = results_for_pipeline,
#   methods        = c("ksea", "mean", "median", "ulm", "wmean"),
#   min_substrates = 3, min_resources = 2,
#   top_n_kinases  = 25, create_plots = TRUE,
#   output_dir     = file.path(best_result$output_dir, "kinase_activity")
# )
#
# # Kinase pathway analysis (fgsea + PROGENy)
# pathway_results <- run_kinase_pathway_analysis(
#   kinase_results = kinase_results,
#   output_dir     = file.path(best_result$output_dir, "kinase_activity"),
#   fdr_threshold = 0.05, top_n_paths = 20, create_plots = TRUE
# )
#
# # ERBB-specific phosphosites
# erbb_results <- analyze_erbb_phosphosites_final(fit2 = best_result$fit, alpha = 0.05)
# nrg1_results <- check_nrg1_genes_in_de(fit2 = best_result$fit, alpha = 0.05)




#############

################################################################################
# KINASE ACTIVITY INFERENCE — ALL BENCHMARK COMPARISONS + CORE DCM GENES
#
# Run AFTER: Human_DCM_Phosphoproteomics_Benchmark_v3.R
# Requires:  all_results list from the benchmark (must still be in memory)
#            run_kinase_activity() from the v14 pipeline
#
# What it does:
#   1. Sources the v14 pipeline to load run_kinase_activity()
#   2. Runs 5-method kinase inference on each benchmark comparison (C1–C8)
#   3. Extracts the 821 core DCM genes from core_DCM_genes.csv
#   4. Runs kinase inference on core DCM genes using C2 statistics
#   5. Produces a cross-comparison kinase summary
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Path to the main pipeline (needed for run_kinase_activity function)
PIPELINE_PATH <- "Script_with_NRG1_ERBB_analysis_kinases_update_protein_normalisation14.R"

# Which methods to run (all 5 from the pipeline)
METHODS <- c("ksea", "mean", "median", "ulm", "wmean")

# Kinase filtering parameters
MIN_SUBSTRATES <- 3
MIN_RESOURCES  <- 2
TOP_N_KINASES  <- 25

# ==============================================================================
# SOURCE THE PIPELINE (to get run_kinase_activity)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  KINASE ACTIVITY INFERENCE — ALL COMPARISONS + CORE DCM\n")
cat("================================================================================\n\n")

message("Sourcing pipeline functions from: ", PIPELINE_PATH)
message("  (Loading function definitions only — skipping execution block)")
original_wd <- getwd()

# Read only the function definitions (lines 1–3513), skip the execution block
# which starts at line 3514 ("STEP 7: RUN THE ANALYSIS")
pipeline_lines <- readLines(PIPELINE_PATH)
exec_start <- grep("^# STEP 7: RUN THE ANALYSIS", pipeline_lines)
if (length(exec_start) > 0) {
  # Keep everything before the execution block
  func_lines <- pipeline_lines[1:(exec_start[1] - 2)]
  message(sprintf("  → Parsing lines 1–%d (skipping execution from line %d)",
                  exec_start[1] - 2, exec_start[1]))
} else {
  func_lines <- pipeline_lines
  message("  → Could not find execution block marker — loading entire file")
}
eval(parse(text = func_lines), envir = globalenv())

setwd(original_wd)  # Restore in case any setwd() was in the function defs
message("  ✓ Functions loaded: run_kinase_activity(), run_kinase_pathway_analysis(), etc.")
message("  ✓ Working directory preserved: ", getwd(), "\n")

# ==============================================================================
# VERIFY all_results EXIST
# ==============================================================================

if (!exists("all_results") || length(all_results) == 0) {
  stop("all_results not found in environment.\n",
       "  Run Human_DCM_Phosphoproteomics_Benchmark_v3.R first,\n",
       "  then run this script in the SAME R session.")
}

cat(sprintf("  Found %d benchmark comparisons in memory:\n", length(all_results)))
for (nm in names(all_results)) {
  cat(sprintf("    • %s\n", all_results[[nm]]$comparison_label))
}

# ==============================================================================
# RUN KINASE ACTIVITY ON EACH COMPARISON
# ==============================================================================

kinase_all_comparisons <- list()

for (nm in names(all_results)) {
  res <- all_results[[nm]]
  out_dir <- file.path(res$output_dir, "kinase_activity")
  
  cat("\n")
  cat(paste(rep("━", 78), collapse = ""), "\n")
  cat(sprintf("  KINASE INFERENCE: %s\n", res$comparison_label))
  cat(paste(rep("━", 78), collapse = ""), "\n")
  
  kinase_res <- tryCatch(
    run_kinase_activity(
      results        = list(top_table = res$top_table),
      methods        = METHODS,
      min_substrates = MIN_SUBSTRATES,
      min_resources  = MIN_RESOURCES,
      top_n_kinases  = TOP_N_KINASES,
      create_plots   = TRUE,
      output_dir     = out_dir
    ),
    error = function(e) {
      message("  ✗ FAILED: ", e$message)
      NULL
    }
  )
  
  if (!is.null(kinase_res)) {
    kinase_all_comparisons[[nm]] <- kinase_res
    
    # ---- Pathway analysis: fgsea (Hallmarks, KEGG, Reactome) + PROGENy ----
    tryCatch({
      pathway_res <- run_kinase_pathway_analysis(
        kinase_results = kinase_res,
        output_dir     = out_dir,
        fdr_threshold  = 0.05,
        top_n_paths    = 20,
        create_plots   = TRUE
      )
      message(sprintf("  ✓ Pathway analysis complete for %s", nm))
    }, error = function(e) {
      message(sprintf("  ⚠ Pathway analysis failed for %s: %s", nm, e$message))
    })
  }
}

# ==============================================================================
# CORE DCM GENES — KINASE ACTIVITY ON THE 821-GENE CONSENSUS SET
# ==============================================================================

cat("\n")
cat(paste(rep("━", 78), collapse = ""), "\n")
cat("  KINASE INFERENCE: CORE DCM GENES (821 genes, FDR<0.05 across C1–C4)\n")
cat(paste(rep("━", 78), collapse = ""), "\n\n")

# Load core gene list
core_genes_file <- "core_DCM_genes.csv"
if (file.exists(core_genes_file)) {
  core_genes <- read.csv(core_genes_file, stringsAsFactors = FALSE)$Gene
  message(sprintf("  ✓ Loaded %d core DCM genes from %s", length(core_genes), core_genes_file))
} else {
  # Rebuild from all_results if file not found
  message("  core_DCM_genes.csv not found — rebuilding from C1–C4 overlap...")
  dcm_comps <- intersect(c("C1", "C2", "C3", "C4"), names(all_results))
  if (length(dcm_comps) >= 2) {
    sig_gene_lists <- lapply(all_results[dcm_comps], function(r) {
      unique(r$top_table$Gene[r$top_table$adj.P.Val < 0.05])
    })
    all_genes <- unique(unlist(sig_gene_lists))
    overlap_mat <- sapply(sig_gene_lists, function(g) all_genes %in% g)
    core_genes <- all_genes[rowSums(overlap_mat) == length(dcm_comps)]
    message(sprintf("  ✓ Rebuilt: %d core DCM genes", length(core_genes)))
    write.csv(data.frame(Gene = sort(core_genes)), core_genes_file, row.names = FALSE)
  } else {
    core_genes <- character(0)
    message("  ✗ Could not rebuild core genes — fewer than 2 DCM comparisons available")
  }
}

if (length(core_genes) > 0) {
  # Use C2 (primary comparison) statistics, filtered to core genes
  c2_tt <- all_results[["C2"]]$top_table
  
  core_tt <- c2_tt[c2_tt$Gene %in% core_genes, ]
  message(sprintf("  ✓ Core DCM phosphosites from C2: %d (from %d genes)",
                  nrow(core_tt), length(unique(core_tt$Gene))))
  
  core_out_dir <- "Core_DCM/kinase_activity"
  dir.create("Core_DCM", showWarnings = FALSE)
  
  # Save the core gene top table for reference
  write.csv(core_tt, "Core_DCM/core_DCM_phosphosite_statistics.csv", row.names = FALSE)
  
  kinase_core <- tryCatch(
    run_kinase_activity(
      results        = list(top_table = core_tt),
      methods        = METHODS,
      min_substrates = MIN_SUBSTRATES,
      min_resources  = MIN_RESOURCES,
      top_n_kinases  = TOP_N_KINASES,
      create_plots   = TRUE,
      output_dir     = core_out_dir
    ),
    error = function(e) {
      message("  ✗ FAILED: ", e$message)
      NULL
    }
  )
  
  if (!is.null(kinase_core)) {
    kinase_all_comparisons[["Core_DCM"]] <- kinase_core
    
    # ---- Pathway analysis for Core DCM ----
    tryCatch({
      pathway_core <- run_kinase_pathway_analysis(
        kinase_results = kinase_core,
        output_dir     = core_out_dir,
        fdr_threshold  = 0.05,
        top_n_paths    = 20,
        create_plots   = TRUE
      )
      message("  ✓ Pathway analysis complete for Core DCM")
    }, error = function(e) {
      message("  ⚠ Pathway analysis failed for Core DCM: ", e$message)
    })
  }
}

# ==============================================================================
# CROSS-COMPARISON KINASE SUMMARY
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                 KINASE ACTIVITY: CROSS-COMPARISON SUMMARY                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════════╝\n\n")

# Build a master KSEA table: kinase × comparison
if (length(kinase_all_comparisons) > 0) {
  
  # Extract KSEA z-scores from each comparison
  ksea_cross <- lapply(names(kinase_all_comparisons), function(nm) {
    ks <- kinase_all_comparisons[[nm]]$ksea_results
    if (is.null(ks)) return(NULL)
    ks %>%
      dplyr::select(kinase, z_score, fdr, n_substrates) %>%
      dplyr::mutate(comparison = nm)
  })
  ksea_cross <- dplyr::bind_rows(ksea_cross)
  
  if (nrow(ksea_cross) > 0) {
    
    # ---- Wide format: kinase × comparison z-scores ----
    ksea_wide_z <- ksea_cross %>%
      dplyr::select(kinase, comparison, z_score) %>%
      tidyr::pivot_wider(names_from = comparison, values_from = z_score,
                         names_prefix = "zscore_")
    
    ksea_wide_fdr <- ksea_cross %>%
      dplyr::select(kinase, comparison, fdr) %>%
      tidyr::pivot_wider(names_from = comparison, values_from = fdr,
                         names_prefix = "FDR_")
    
    ksea_wide <- dplyr::left_join(ksea_wide_z, ksea_wide_fdr, by = "kinase")
    
    # Add substrate count and number of comparisons where FDR < 0.05
    ksea_summary <- ksea_cross %>%
      dplyr::group_by(kinase) %>%
      dplyr::summarise(
        n_comparisons_tested = dplyr::n(),
        n_sig_FDR05 = sum(fdr < 0.05, na.rm = TRUE),
        mean_z = mean(z_score, na.rm = TRUE),
        median_z = median(z_score, na.rm = TRUE),
        max_substrates = max(n_substrates, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(n_sig_FDR05), desc(abs(mean_z)))
    
    write.csv(ksea_wide,    "kinase_KSEA_cross_comparison.csv",  row.names = FALSE)
    write.csv(ksea_summary, "kinase_KSEA_summary.csv",           row.names = FALSE)
    
    # ---- Print top kinases ----
    cat("  Kinases significant (KSEA FDR<0.05) across most comparisons:\n\n")
    cat(sprintf("  %-15s %5s %8s %8s %4s\n",
                "Kinase", "nSig", "mean_z", "med_z", "nSub"))
    cat(paste(rep("─", 50), collapse = ""), "\n")
    
    top_summary <- head(ksea_summary, 30)
    for (i in seq_len(nrow(top_summary))) {
      r <- top_summary[i, ]
      sig_flag <- if (r$n_sig_FDR05 >= length(kinase_all_comparisons) - 1) "★" else ""
      cat(sprintf("  %-15s %2d/%d  %+7.2f %+7.2f  %3d %s\n",
                  r$kinase, r$n_sig_FDR05, r$n_comparisons_tested,
                  r$mean_z, r$median_z, r$max_substrates, sig_flag))
    }
    
    # ---- Consistently significant kinases (FDR<0.05 in ALL comparisons) ----
    n_total_comps <- length(kinase_all_comparisons)
    consistent_kinases <- ksea_summary %>%
      dplyr::filter(n_sig_FDR05 == n_total_comps)
    
    if (nrow(consistent_kinases) > 0) {
      cat(sprintf("\n\n  ★ Kinases FDR<0.05 in ALL %d comparisons: %d\n",
                  n_total_comps, nrow(consistent_kinases)))
      cat(sprintf("    Active (mean z > 0):   %s\n",
                  paste(consistent_kinases$kinase[consistent_kinases$mean_z > 0], collapse = ", ")))
      cat(sprintf("    Inactive (mean z < 0): %s\n",
                  paste(consistent_kinases$kinase[consistent_kinases$mean_z < 0], collapse = ", ")))
      
      write.csv(consistent_kinases, "kinase_consistent_across_all.csv", row.names = FALSE)
    }
    
    # ---- Cross-comparison heatmap ----
    tryCatch({
      # Select kinases significant in >= 2 comparisons
      plot_kinases <- ksea_summary %>%
        dplyr::filter(n_sig_FDR05 >= 2) %>%
        dplyr::arrange(mean_z) %>%
        dplyr::pull(kinase)
      
      if (length(plot_kinases) >= 3) {
        heat_df <- ksea_cross %>%
          dplyr::filter(kinase %in% plot_kinases) %>%
          dplyr::mutate(
            kinase = factor(kinase, levels = plot_kinases),
            comparison = gsub("^C", "", comparison),
            sig_star = ifelse(fdr < 0.05, "*", "")
          )
        
        n_kin  <- length(plot_kinases)
        n_comp <- length(unique(heat_df$comparison))
        
        p_heat <- ggplot2::ggplot(heat_df,
                                  ggplot2::aes(x = comparison, y = kinase, fill = z_score)) +
          ggplot2::geom_tile(color = "white", linewidth = 0.3) +
          ggplot2::geom_text(ggplot2::aes(label = sig_star), size = 4, fontface = "bold") +
          ggplot2::scale_fill_gradient2(
            low = "#2E86AB", mid = "white", high = "#E63946", midpoint = 0,
            name = "KSEA\nz-score"
          ) +
          ggplot2::labs(
            title = "Kinase Activity Across Comparisons (KSEA)",
            subtitle = paste0("Kinases FDR<0.05 in ≥2 comparisons | * = FDR<0.05 in that comparison"),
            x = "Comparison", y = NULL
          ) +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13),
            plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = ggplot2::element_text(size = 8),
            panel.grid = ggplot2::element_blank(),
            legend.key.height = ggplot2::unit(1.2, "cm")
          )
        
        ggplot2::ggsave("kinase_cross_comparison_heatmap.png", p_heat,
                        width = max(8, n_comp * 0.8 + 4),
                        height = max(6, n_kin * 0.35 + 3), dpi = 300)
        message("  ✓ Cross-comparison heatmap saved")
      }
    }, error = function(e) message("  ⚠ Heatmap failed: ", e$message))
    
    # ── Correlation matrix: how similar are kinase profiles across comparisons? ──
    tryCatch({
      message("\n  Building kinase profile correlation matrix...")
      
      # Pivot to wide: kinase × comparison, values = z_score
      z_wide <- ksea_cross %>%
        dplyr::select(kinase, comparison, z_score) %>%
        tidyr::pivot_wider(names_from = comparison, values_from = z_score) %>%
        tibble::column_to_rownames("kinase")
      
      # Only keep kinases present in at least 2 comparisons
      z_wide <- z_wide[rowSums(!is.na(z_wide)) >= 2, , drop = FALSE]
      
      if (ncol(z_wide) >= 2 && nrow(z_wide) >= 5) {
        # Pearson correlation between comparisons (column-wise)
        cor_mat <- cor(z_wide, use = "pairwise.complete.obs", method = "pearson")
        
        write.csv(cor_mat, "kinase_profile_correlation_matrix.csv", row.names = TRUE)
        
        # Print correlation matrix
        cat("\n\n  Pearson correlation of KSEA z-score profiles between comparisons:\n")
        cat(sprintf("  (based on %d kinases scored in ≥2 comparisons)\n\n", nrow(z_wide)))
        
        # Pretty print
        comp_names <- colnames(cor_mat)
        cat(sprintf("  %15s", ""))
        for (cn in comp_names) cat(sprintf(" %8s", cn))
        cat("\n")
        for (i in seq_along(comp_names)) {
          cat(sprintf("  %-15s", comp_names[i]))
          for (j in seq_along(comp_names)) {
            cat(sprintf(" %8.3f", cor_mat[i, j]))
          }
          cat("\n")
        }
        
        # ---- Correlation heatmap plot ----
        cor_melted <- reshape2::melt(cor_mat)
        colnames(cor_melted) <- c("Comp1", "Comp2", "Correlation")
        
        p_cor <- ggplot2::ggplot(cor_melted,
                                 ggplot2::aes(x = Comp1, y = Comp2, fill = Correlation)) +
          ggplot2::geom_tile(color = "white", linewidth = 0.5) +
          ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                             size = 3.5, color = "black", fontface = "bold") +
          ggplot2::scale_fill_gradient2(
            low = "#2E86AB", mid = "white", high = "#E63946", midpoint = 0.5,
            limits = c(min(cor_mat) - 0.05, 1),
            name = "Pearson r"
          ) +
          ggplot2::coord_fixed() +
          ggplot2::labs(
            title = "Kinase Profile Correlation Between Comparisons",
            subtitle = sprintf("Pearson r of KSEA z-scores | %d kinases | DCM vs Healthy benchmark",
                               nrow(z_wide)),
            x = NULL, y = NULL
          ) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13),
            plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = ggplot2::element_text(size = 9),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right"
          )
        
        ggplot2::ggsave("kinase_profile_correlation_heatmap.png", p_cor,
                        width = max(7, ncol(cor_mat) * 0.9 + 3),
                        height = max(6, nrow(cor_mat) * 0.9 + 3), dpi = 300)
        message("  ✓ Correlation heatmap saved")
        
        # ---- Hierarchical clustering of comparisons by kinase profile ----
        if (ncol(z_wide) >= 3) {
          tryCatch({
            dist_mat <- as.dist(1 - cor_mat)
            hc <- hclust(dist_mat, method = "ward.D2")
            
            png("kinase_profile_dendrogram.png", width = 900, height = 600, res = 150)
            par(mar = c(6, 4, 3, 1))
            plot(hc, main = "Clustering of Comparisons by Kinase Activity Profile",
                 sub = paste0("Ward.D2 | 1 - Pearson r | ", nrow(z_wide), " kinases"),
                 xlab = "", ylab = "Distance (1 - r)",
                 hang = -1, cex = 0.9)
            dev.off()
            message("  ✓ Dendrogram saved")
          }, error = function(e) NULL)
        }
        
        # ---- Scatter plots: each comparison vs Core_DCM kinase z-scores ----
        if ("Core_DCM" %in% colnames(z_wide)) {
          
          other_comps <- setdiff(colnames(z_wide), "Core_DCM")
          scatter_plots <- list()
          
          for (comp_nm in other_comps) {
            scatter_df <- data.frame(
              kinase = rownames(z_wide),
              comp   = z_wide[, comp_nm],
              Core   = z_wide[, "Core_DCM"],
              stringsAsFactors = FALSE
            ) %>% dplyr::filter(!is.na(comp) & !is.na(Core))
            
            if (nrow(scatter_df) < 5) next
            
            r_val <- cor(scatter_df$comp, scatter_df$Core, use = "complete.obs")
            
            # Label kinases significant in both
            ksea_comp <- kinase_all_comparisons[[comp_nm]]$ksea_results
            ksea_core <- kinase_all_comparisons[["Core_DCM"]]$ksea_results
            sig_both <- character(0)
            if (!is.null(ksea_comp) && !is.null(ksea_core)) {
              sig_both <- intersect(
                ksea_comp$kinase[ksea_comp$fdr < 0.05],
                ksea_core$kinase[ksea_core$fdr < 0.05]
              )
            }
            scatter_df$label <- ifelse(scatter_df$kinase %in% sig_both, scatter_df$kinase, "")
            
            # Get a clean comparison label
            comp_label <- if (!is.null(all_results[[comp_nm]])) {
              gsub("^C\\d+: ", "", all_results[[comp_nm]]$comparison_label)
            } else { comp_nm }
            
            p_sc <- ggplot2::ggplot(scatter_df,
                                    ggplot2::aes(x = comp, y = Core)) +
              ggplot2::geom_point(alpha = 0.6, size = 2, color = "#457B9D") +
              ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
              ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#E63946",
                                   linewidth = 0.8, alpha = 0.15) +
              ggrepel::geom_text_repel(ggplot2::aes(label = label),
                                       size = 2.8, max.overlaps = 20, segment.size = 0.3) +
              ggplot2::labs(
                title = paste0("Kinase Activity: ", comp_nm, " vs Core DCM"),
                subtitle = sprintf("Pearson r = %.3f | %d kinases | labeled = FDR<0.05 in both",
                                   r_val, nrow(scatter_df)),
                x = paste0("KSEA z-score — ", comp_label),
                y = "KSEA z-score — Core DCM (821 genes)"
              ) +
              ggplot2::theme_minimal(base_size = 11) +
              ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 13),
                plot.subtitle = ggplot2::element_text(size = 9, color = "grey40")
              )
            
            fname <- sprintf("kinase_%s_vs_CoreDCM_scatter.png", comp_nm)
            ggplot2::ggsave(fname, p_sc, width = 9, height = 8, dpi = 300)
            scatter_plots[[comp_nm]] <- p_sc
            message(sprintf("  ✓ %s vs Core_DCM scatter saved (r = %.3f)", comp_nm, r_val))
          }
          
          # ---- Combined faceted panel: all comparisons vs Core DCM ----
          if (length(other_comps) >= 2) {
            tryCatch({
              facet_df <- do.call(rbind, lapply(other_comps, function(comp_nm) {
                df <- data.frame(
                  kinase = rownames(z_wide),
                  comp   = z_wide[, comp_nm],
                  Core   = z_wide[, "Core_DCM"],
                  stringsAsFactors = FALSE
                ) %>% dplyr::filter(!is.na(comp) & !is.na(Core))
                
                if (nrow(df) == 0) return(NULL)
                
                r_val <- cor(df$comp, df$Core, use = "complete.obs")
                comp_label <- if (!is.null(all_results[[comp_nm]])) {
                  gsub("^C\\d+: ", "", all_results[[comp_nm]]$comparison_label)
                } else { comp_nm }
                
                df$panel <- sprintf("%s (r=%.2f)", comp_label, r_val)
                df
              }))
              
              if (!is.null(facet_df) && nrow(facet_df) > 0) {
                p_facet <- ggplot2::ggplot(facet_df,
                                           ggplot2::aes(x = comp, y = Core)) +
                  ggplot2::geom_point(alpha = 0.4, size = 1, color = "#457B9D") +
                  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                                       color = "grey50", linewidth = 0.4) +
                  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#E63946",
                                       linewidth = 0.7, alpha = 0.15) +
                  ggplot2::facet_wrap(~ panel, scales = "free_x", ncol = 3) +
                  ggplot2::labs(
                    title = "Kinase Activity: All Comparisons vs Core DCM Signature",
                    subtitle = "Each panel shows one comparison's KSEA z-scores against the 821-gene core DCM set",
                    x = "KSEA z-score (comparison)",
                    y = "KSEA z-score (Core DCM)"
                  ) +
                  ggplot2::theme_minimal(base_size = 10) +
                  ggplot2::theme(
                    plot.title = ggplot2::element_text(face = "bold", size = 14),
                    plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
                    strip.text = ggplot2::element_text(size = 8, face = "bold"),
                    panel.spacing = ggplot2::unit(0.8, "lines")
                  )
                
                n_panels <- length(unique(facet_df$panel))
                n_rows <- ceiling(n_panels / 3)
                ggplot2::ggsave("kinase_ALL_vs_CoreDCM_faceted.png", p_facet,
                                width = 14, height = 4.5 * n_rows + 1.5, dpi = 300)
                message("  ✓ Combined faceted scatter saved")
              }
            }, error = function(e) message("  ⚠ Faceted scatter failed: ", e$message))
          }
        }
      }
    }, error = function(e) message("  ⚠ Correlation analysis failed: ", e$message))
    
  }
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  KINASE INFERENCE COMPLETE\n")
cat("================================================================================\n\n")
cat("  Per-comparison outputs (in each C*/kinase_activity/ folder):\n")
cat("    • ksea_results.csv, mean_results.csv, median_results.csv\n")
cat("    • ulm_results.csv, wmean_results.csv\n")
cat("    • kinase_activity_all_methods.csv\n")
cat("    • ksea_barplot.png, ksea_volcano.png\n")
cat("    • kinase_method_comparison_heatmap.png\n")
cat("    • fgsea_hallmarks_results.csv + _barplot.png\n")
cat("    • fgsea_kegg_results.csv      + _barplot.png\n")
cat("    • fgsea_reactome_results.csv  + _barplot.png\n")
cat("    • gsea_plots_hallmarks/ gsea_plots_kegg/ gsea_plots_reactome/\n")
cat("    • progeny_pathway_activity.csv + .png\n\n")
cat("  Core DCM outputs (in Core_DCM/):\n")
cat("    • core_DCM_phosphosite_statistics.csv\n")
cat("    • kinase_activity/ (same structure as above)\n\n")
cat("  Cross-comparison outputs (in working directory):\n")
cat("    • kinase_KSEA_cross_comparison.csv  — z-scores + FDR per comparison\n")
cat("    • kinase_KSEA_summary.csv           — ranked by consistency\n")
cat("    • kinase_consistent_across_all.csv  — FDR<0.05 in every comparison\n")
cat("    • kinase_cross_comparison_heatmap.png  — z-score heatmap (kinase × comparison)\n")
cat("    • kinase_profile_correlation_matrix.csv — Pearson r between comparisons\n")
cat("    • kinase_profile_correlation_heatmap.png — correlation heatmap\n")
cat("    • kinase_profile_dendrogram.png     — hierarchical clustering of comparisons\n")
cat("    • kinase_C1_vs_CoreDCM_scatter.png  — individual scatter per comparison\n")
cat("    • kinase_C2_vs_CoreDCM_scatter.png    ... (one per comparison)\n")
cat("    • kinase_ALL_vs_CoreDCM_faceted.png — combined faceted panel\n")
cat("================================================================================\n\n")
