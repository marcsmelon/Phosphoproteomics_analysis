################################################################################
# iPSC-CM DCM PHOSPHOPROTEOMICS PIPELINE
#
# 4 genotypes (TTN, DSP, TNNT2, TPM1) × 3 timepoints (D15, D22, D29)
# Each genotype: isogenic control (Healthy) vs DCM mutation (Disease)
# 4 biological replicates per condition (with exceptions noted below)
#
# Analyses:
#   1.  Per-genotype per-timepoint DE (12 contrasts — DSP D22 skipped = 11)
#   2.  Host protein normalization (subtract log2 protein from log2 phosphosite)
#   3.  Kinase activity inference (KSEA, ULM, WMEAN) per contrast
#   4.  TWO KSEA volcano styles per contrast:
#       a) Original: scatter with size = n_substrates, labels on FDR < 0.05
#       b) Compact: theme_void, top 3 labels, n= counts, ~2.35×2.43 inches
#   5.  PCA & UMAP (global: by timepoint, by genotype)
#   6.  Global heatmaps per timepoint (all FDR < 0.05 sites)
#   7.  Top 50 heatmaps per genotype (25 up + 25 down)
#   8.  fgsea (Hallmarks, KEGG, Reactome) + PROGENy per contrast
#   9.  GO BP/MF/CC (all/up/down) with dotplots per contrast
#  10.  Cross-timepoint contrasts (D15 vs D22 vs D29 within genotype)
#  11.  Scatterplots + correlation matrices (within/across genotypes)
#  12.  Interaction model (Disease × Timepoint)
#  13.  Core DCM iPSC phospho signature per timepoint
#  14.  Maturation trajectory (PCA/UMAP by timepoint)
#
# Known data issues (handled):
#   - DSP D22 Healthy: only 1 replicate → D22 contrast skipped
#   - DSP D29 labels swapped: 522p3DSPKI = Disease, WT522p3 = Healthy (fixed)
#   - TNNT2 D22 Disease: only 2 replicates (run anyway)
#   - PB005_WT_Batch_1: batch control → excluded
#
# Colour scheme: gold (#F5CD6A) = up in DCM, purple (#3A2044) = down
################################################################################

# ==============================================================================
# SECTION 0: CONFIGURATION
# ==============================================================================

PHOSPHO_FILE   <- "Phospho_STY_.tsv"
PROTEIN_FILE   <- "General_proteome_WithIBAQ.tsv"
OUTPUT_BASE    <- "iPSC_Phosphoproteomics_Benchmark"

# v14 pipeline path — sourced for kinase inference + protein normalisation functions
V14_SCRIPT     <- "Script_with_NRG1_ERBB_analysis_kinases_update_protein_normalisation14.R"

# Filtering
FILTER_THRESHOLD  <- 3        # min non-NA per row (across all samples in a contrast)
PROTEIN_NORM_THRESHOLD <- 30  # min non-NA for protein normalisation

# Kinase inference
MIN_SUBSTRATES    <- 3
MIN_RESOURCES     <- 2

# Significance tiers for export
SIG_TIERS <- list(
  "FDR005"  = list(col = "adj.P.Val", thresh = 0.05),
  "FDR010"  = list(col = "adj.P.Val", thresh = 0.10),
  "FDR020"  = list(col = "adj.P.Val", thresh = 0.20),
  "rawP001" = list(col = "P.Value",   thresh = 0.01),
  "rawP005" = list(col = "P.Value",   thresh = 0.05)
)

# ==============================================================================
# SECTION 1: PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(decoupleR)
  library(OmnipathR)
  library(umap)
})

# ==============================================================================
# SECTION 2: COLOUR SCHEME (matches proteomics benchmark)
# ==============================================================================

COL_DCM     <- "#F5CD6A"    # gold — Disease / Up in DCM
COL_HEALTHY <- "#3A2044"    # dark purple — Healthy / Down in DCM
COL_NS      <- "#EAEAEA"    # grey — not significant

CONDITION_COLORS  <- c("Healthy" = COL_HEALTHY, "Disease" = COL_DCM)
DIRECTION_COLORS  <- c("Up in DCM" = COL_DCM, "Down in DCM" = COL_HEALTHY)

VOLCANO_COLORS <- c(
  "Up (FDR<0.05, FC>2)"   = COL_DCM,
  "Down (FDR<0.05, FC>2)" = COL_HEALTHY,
  "Sig (FDR<0.05)"        = "#C4A84E",
  "Nominal (p<0.01)"      = "#8B7BA5",
  "NS"                    = "grey80"
)

KINASE_COLORS <- c("Active" = COL_DCM, "Inactive" = COL_HEALTHY,
                   "UP" = COL_DCM, "DOWN" = COL_HEALTHY, "NS" = COL_NS)

# ---- Global plot dimensions (all plots except heatmaps) ----
# 169.0586 × 174.6366 pt ≈ 2.35 × 2.43 in ≈ 5.96 × 6.16 cm
PLOT_W <- 169.0586 / 72   # inches
PLOT_H <- 174.6366 / 72   # inches

# ---- Global plot theme (all plots except heatmaps) ----
# theme_void base + 1 pt axis lines + 5 pt text
THEME_PLOT <- theme_void() +
  theme(
    axis.title        = element_text(size = 5, color = "grey20"),
    axis.text         = element_text(size = 5, color = "grey30"),
    axis.line         = element_line(color = "black", linewidth = 1 / ggplot2::.pt),
    axis.ticks        = element_line(color = "black", linewidth = 1 / ggplot2::.pt),
    legend.text       = element_text(size = 5),
    legend.title      = element_text(size = 5),
    plot.title        = element_text(size = 5, face = "bold"),
    plot.subtitle     = element_text(size = 5, color = "grey40"),
    plot.margin       = margin(4, 4, 4, 4),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key.size   = unit(0.25, "cm"),
    strip.text        = element_text(size = 5, face = "bold")
  )

# ==============================================================================
# SECTION 3: SOURCE v14 FUNCTIONS (safe — skip execution block)
# ==============================================================================

cat("Sourcing kinase + normalisation functions from v14 pipeline...\n")
v14_lines <- readLines(V14_SCRIPT, warn = FALSE)
exec_line <- grep("^# STEP 7: RUN THE ANALYSIS", v14_lines)
if (length(exec_line) == 0) {
  exec_line <- grep("^################################################################################$", v14_lines)
  exec_line <- exec_line[exec_line > 3000]
  if (length(exec_line) > 0) exec_line <- exec_line[1] else exec_line <- length(v14_lines)
}
func_block <- paste(v14_lines[1:(exec_line[1] - 1)], collapse = "\n")
eval(parse(text = func_block))
cat("  ✓ Functions sourced (up to line", exec_line[1] - 1, ")\n")

# ==============================================================================
# SECTION 4: METADATA — SAMPLE MAPPING
# ==============================================================================

# Disease cell lines (DCM mutations)
DISEASE_LINES <- list(
  TTN   = "10p5TTNTrunc_TTN",
  DSP   = "522p3DSPKI_DSP",
  TNNT2 = "MCHTB63i8_TNNT2",
  TPM1  = "TPM1120i6_TPM1"
)

# Healthy cell lines (isogenic controls)
HEALTHY_LINES <- list(
  TTN   = "10p5_TTN",
  DSP   = "WT522p3_DSP",
  TNNT2 = "MCHTB63i8GC_TNNT2",
  TPM1  = "TPM19p1_TPM1"
)

TIMEPOINTS <- c("D15", "D22", "D29")
SKIP_CONTRASTS <- c("DSP_D22")  # only 1 healthy replicate
EXCLUDE_SAMPLES <- c("PB005_WT_Batch_1")

# DSP D29 label swap: in the raw file, 522p3DSPKI is labeled "Healthy"
# and WT522p3 is labeled "Disease" — this is REVERSED
DSP_D29_SWAP <- TRUE

# ==============================================================================
# SECTION 5: LOAD PHOSPHO DATA
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  LOADING PHOSPHOPROTEOMICS DATA\n")
cat("════════════════════════════════════════════════════════════════════════════════\n\n")

raw <- `Phospho(STY)`
cat(sprintf("  Loaded: %d rows × %d columns\n", nrow(raw), ncol(raw)))

# ---- Detect intensity columns ----
# Try multiple MaxQuant Phospho(STY) naming conventions
intensity_cols <- grep("^Intensity ", colnames(raw), value = TRUE)
if (length(intensity_cols) == 0) {
  intensity_cols <- grep("^LFQ\\.intensity\\.", colnames(raw), value = TRUE)
}
if (length(intensity_cols) == 0) {
  # Spectronaut-style: columns containing cell line names
  all_lines <- c(unlist(DISEASE_LINES), unlist(HEALTHY_LINES))
  intensity_cols <- colnames(raw)[sapply(colnames(raw), function(cn) {
    any(sapply(all_lines, function(cl) grepl(cl, cn, fixed = TRUE)))
  })]
}
cat(sprintf("  Detected %d intensity columns\n", length(intensity_cols)))
if (length(intensity_cols) == 0) stop("No intensity columns detected. Check file format.")

# ---- Exclude batch control samples ----
for (excl in EXCLUDE_SAMPLES) {
  drop_cols <- grep(excl, intensity_cols, value = TRUE, fixed = TRUE)
  if (length(drop_cols) > 0) {
    cat(sprintf("  Excluding %d columns matching '%s'\n", length(drop_cols), excl))
    intensity_cols <- setdiff(intensity_cols, drop_cols)
  }
}

# ---- Build GenePhos identifiers ----
# MaxQuant Phospho(STY) has: Gene names, Amino acid, Position, Multiplicity
if ("Gene.names" %in% colnames(raw) && "Amino.acid" %in% colnames(raw) &&
    "Position" %in% colnames(raw)) {
  # MaxQuant format
  raw$Gene <- sapply(raw$Gene.names, function(x) strsplit(as.character(x), ";")[[1]][1])
  if ("Multiplicity" %in% colnames(raw)) {
    raw$GenePhos <- paste0(raw$Gene, "_", raw$Amino.acid, raw$Position, "_M", raw$Multiplicity)
  } else {
    raw$GenePhos <- paste0(raw$Gene, "_", raw$Amino.acid, raw$Position)
  }
  cat("  Built GenePhos from MaxQuant columns (Gene.names + Amino.acid + Position)\n")
} else if ("GenePhos" %in% colnames(raw)) {
  cat("  GenePhos column already present\n")
  raw$Gene <- gsub("_[STY].*$", "", raw$GenePhos)
} else if ("EG.Gene" %in% colnames(raw) || "PG.Genes" %in% colnames(raw)) {
  # Spectronaut-like format
  gene_col <- if ("EG.Gene" %in% colnames(raw)) "EG.Gene" else "PG.Genes"
  if ("PTM.SiteAA" %in% colnames(raw) && "PTM.SiteLocation" %in% colnames(raw)) {
    raw$Gene <- raw[[gene_col]]
    raw$GenePhos <- paste0(raw$Gene, "_", raw$PTM.SiteAA, raw$PTM.SiteLocation)
    if ("PTM.Multiplicity" %in% colnames(raw)) {
      raw$GenePhos <- paste0(raw$GenePhos, "_M", raw$PTM.Multiplicity)
    }
    cat(sprintf("  Built GenePhos from Spectronaut columns (%s + PTM.SiteAA + PTM.SiteLocation)\n", gene_col))
  } else {
    stop("Cannot build GenePhos identifiers: missing site annotation columns.")
  }
} else {
  stop("Cannot identify gene/site columns. Expected MaxQuant or Spectronaut format.")
}

# Remove rows with NA gene names
raw <- raw[!is.na(raw$Gene) & raw$Gene != "" & raw$Gene != "NA", ]
cat(sprintf("  After gene filter: %d phosphosites\n", nrow(raw)))

# Remove reverse hits and contaminants
if ("Reverse" %in% colnames(raw)) {
  raw <- raw[raw$Reverse != "+" & (is.na(raw$Reverse) | raw$Reverse == ""), ]
}
if ("Potential.contaminant" %in% colnames(raw)) {
  raw <- raw[raw$Potential.contaminant != "+" &
               (is.na(raw$Potential.contaminant) | raw$Potential.contaminant == ""), ]
}
cat(sprintf("  After contaminant filter: %d phosphosites\n", nrow(raw)))

# ---- Extract intensity matrix ----
int_matrix <- as.matrix(raw[, intensity_cols])
rownames(int_matrix) <- raw$GenePhos

# Replace 0 with NA
int_matrix[int_matrix == 0] <- NA

# Log2 transform
log2_matrix <- log2(int_matrix)
cat(sprintf("  Intensity matrix: %d phosphosites × %d samples\n",
            nrow(log2_matrix), ncol(log2_matrix)))

# Strip MaxQuant prefix for cleaner column names
clean_names <- gsub("^Intensity |^LFQ\\.intensity\\.", "", colnames(log2_matrix))
colnames(log2_matrix) <- clean_names
intensity_cols_clean <- clean_names

# ==============================================================================
# SECTION 6: MAP SAMPLES TO CONDITIONS
# ==============================================================================

cat("\n  Mapping samples to genotype/timepoint/condition...\n")

sample_info <- data.frame(
  sample   = intensity_cols_clean,
  genotype = NA_character_,
  timepoint = NA_character_,
  condition = NA_character_,
  stringsAsFactors = FALSE
)

for (gene in names(DISEASE_LINES)) {
  dis_line  <- DISEASE_LINES[[gene]]
  heal_line <- HEALTHY_LINES[[gene]]
  
  for (tp in TIMEPOINTS) {
    # Disease samples
    dis_idx <- grep(paste0(dis_line, ".*", tp), sample_info$sample)
    sample_info$genotype[dis_idx]  <- gene
    sample_info$timepoint[dis_idx] <- tp
    sample_info$condition[dis_idx] <- "Disease"
    
    # Healthy samples
    heal_idx <- grep(paste0(heal_line, ".*", tp), sample_info$sample)
    sample_info$genotype[heal_idx]  <- gene
    sample_info$timepoint[heal_idx] <- tp
    sample_info$condition[heal_idx] <- "Healthy"
  }
}

# DSP D29 label swap fix
if (DSP_D29_SWAP) {
  dsp_d29 <- which(sample_info$genotype == "DSP" & sample_info$timepoint == "D29")
  for (i in dsp_d29) {
    if (grepl("522p3DSPKI", sample_info$sample[i])) {
      sample_info$condition[i] <- "Disease"
    } else if (grepl("WT522p3", sample_info$sample[i])) {
      sample_info$condition[i] <- "Healthy"
    }
  }
  cat("  ✓ DSP D29 label swap corrected\n")
}

# Report mapping
mapped <- sum(!is.na(sample_info$genotype))
cat(sprintf("  ✓ Mapped %d / %d samples\n", mapped, nrow(sample_info)))
unmapped <- sample_info$sample[is.na(sample_info$genotype)]
if (length(unmapped) > 0) {
  cat("  ⚠ Unmapped samples:\n")
  for (u in unmapped) cat("    -", u, "\n")
}

# Remove unmapped columns from matrix
keep_samples <- sample_info$sample[!is.na(sample_info$genotype)]
log2_matrix <- log2_matrix[, keep_samples, drop = FALSE]
sample_info <- sample_info[!is.na(sample_info$genotype), ]
cat(sprintf("  Working matrix: %d phosphosites × %d samples\n",
            nrow(log2_matrix), ncol(log2_matrix)))

# ==============================================================================
# SECTION 7: PROTEIN NORMALISATION
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  HOST PROTEIN NORMALISATION\n")
cat("════════════════════════════════════════════════════════════════════════════════\n\n")

# Load protein data and prepare for normalisation
tryCatch({
  prot_raw <- read.csv(PROTEIN_FILE, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded protein data: %d proteins\n", nrow(prot_raw)))
  
  # Detect protein intensity columns (same samples)
  prot_int_cols <- grep("^LFQ\\.intensity\\.|^Intensity ", colnames(prot_raw), value = TRUE)
  if (length(prot_int_cols) == 0) {
    # Try matching by sample names
    prot_int_cols <- colnames(prot_raw)[sapply(colnames(prot_raw), function(cn) {
      any(sapply(keep_samples, function(s) grepl(s, cn, fixed = TRUE)))
    })]
  }
  
  # Get gene column
  if ("Gene.names" %in% colnames(prot_raw)) {
    prot_raw$Gene <- sapply(prot_raw$Gene.names, function(x) strsplit(as.character(x), ";")[[1]][1])
  } else if ("PG.Genes" %in% colnames(prot_raw)) {
    prot_raw$Gene <- prot_raw$PG.Genes
  } else if ("PG.GroupLabel" %in% colnames(prot_raw)) {
    # UniProt accessions — map to gene symbols
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    up_ids <- prot_raw$PG.GroupLabel
    mapping <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys = up_ids,
                            columns = "SYMBOL", keytype = "UNIPROT"),
      error = function(e) NULL
    )
    if (!is.null(mapping)) {
      mapping <- mapping[!duplicated(mapping$UNIPROT), ]
      prot_raw$Gene <- mapping$SYMBOL[match(prot_raw$PG.GroupLabel, mapping$UNIPROT)]
    } else {
      prot_raw$Gene <- prot_raw$PG.GroupLabel
    }
  }
  
  # Extract protein intensity matrix
  prot_int_matrix <- as.matrix(prot_raw[, prot_int_cols])
  rownames(prot_int_matrix) <- prot_raw$Gene
  prot_int_matrix[prot_int_matrix == 0] <- NA
  
  # Clean column names to match phospho
  prot_clean_names <- gsub("^LFQ\\.intensity\\.|^Intensity ", "", colnames(prot_int_matrix))
  colnames(prot_int_matrix) <- prot_clean_names
  
  # Log2 transform
  log2_prot <- log2(prot_int_matrix)
  
  # Find shared samples
  shared_samples <- intersect(colnames(log2_matrix), colnames(log2_prot))
  cat(sprintf("  Shared samples with phospho: %d\n", length(shared_samples)))
  
  if (length(shared_samples) > 0) {
    # Remove proteins with < threshold non-NA values
    prot_valid <- rowSums(!is.na(log2_prot[, shared_samples, drop = FALSE]))
    log2_prot_filt <- log2_prot[prot_valid >= PROTEIN_NORM_THRESHOLD, shared_samples, drop = FALSE]
    cat(sprintf("  Proteins passing filter (>=%d values): %d / %d\n",
                PROTEIN_NORM_THRESHOLD, nrow(log2_prot_filt), nrow(log2_prot)))
    
    # Extract gene from phosphosite GenePhos
    phospho_genes <- gsub("_[STY].*$", "", rownames(log2_matrix))
    
    # Subtract protein log2 from phosphosite log2 per gene per sample
    log2_corrected <- log2_matrix
    n_corrected <- 0
    for (i in seq_len(nrow(log2_matrix))) {
      gene <- phospho_genes[i]
      if (gene %in% rownames(log2_prot_filt)) {
        for (samp in shared_samples) {
          if (!is.na(log2_matrix[i, samp]) && !is.na(log2_prot_filt[gene, samp])) {
            log2_corrected[i, samp] <- log2_matrix[i, samp] - log2_prot_filt[gene, samp]
          }
        }
        n_corrected <- n_corrected + 1
      }
    }
    cat(sprintf("  ✓ Normalised %d phosphosites to host protein level\n", n_corrected))
    cat("  Method: Müller-Dott et al., Nat Commun 2025\n")
    
    # Use corrected matrix going forward
    log2_matrix_uncorrected <- log2_matrix
    log2_matrix <- log2_corrected
    PROTEIN_NORM_DONE <- TRUE
  } else {
    cat("  ⚠ No shared samples — skipping protein normalisation\n")
    PROTEIN_NORM_DONE <- FALSE
  }
}, error = function(e) {
  cat(sprintf("  ⚠ Protein normalisation failed: %s\n", e$message))
  cat("  Continuing with uncorrected phosphosite data\n")
  PROTEIN_NORM_DONE <<- FALSE
})

# ==============================================================================
# SECTION 8: ENRICHMENT RESOURCES (load once)
# ==============================================================================

cat("\n  Loading enrichment databases...\n")

# MSigDB gene sets
msig_h  <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
msig_re <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

hallmarks_list <- split(msig_h$gene_symbol, msig_h$gs_name)
kegg_list      <- split(msig_c2$gene_symbol, msig_c2$gs_name)
reactome_list  <- split(msig_re$gene_symbol, msig_re$gs_name)

cat(sprintf("  ✓ Hallmarks: %d | KEGG: %d | Reactome: %d gene sets\n",
            length(hallmarks_list), length(kegg_list), length(reactome_list)))

# ==============================================================================
# SECTION 9: CORE ANALYSIS FUNCTIONS
# ==============================================================================

# ---- 9a: run_contrast — limma DE for one genotype × timepoint ----
run_contrast <- function(gene, tp, log2_mat, sinfo, output_dir) {
  
  contrast_name <- paste0(gene, "_", tp)
  cat(sprintf("\n  ── %s ──\n", contrast_name))
  
  # Select samples
  idx <- sinfo$genotype == gene & sinfo$timepoint == tp
  if (sum(idx) == 0) { cat("  ⚠ No samples\n"); return(NULL) }
  
  sub_info <- sinfo[idx, ]
  sub_mat  <- log2_mat[, sub_info$sample, drop = FALSE]
  
  n_dis  <- sum(sub_info$condition == "Disease")
  n_heal <- sum(sub_info$condition == "Healthy")
  cat(sprintf("  Disease: %d | Healthy: %d\n", n_dis, n_heal))
  
  if (n_dis < 2 || n_heal < 1) { cat("  ⚠ Insufficient samples\n"); return(NULL) }
  
  # Filter: keep rows with >= FILTER_THRESHOLD non-NA
  keep <- rowSums(!is.na(sub_mat)) >= FILTER_THRESHOLD
  sub_mat <- sub_mat[keep, , drop = FALSE]
  cat(sprintf("  After filter (>=%d values): %d phosphosites\n", FILTER_THRESHOLD, nrow(sub_mat)))
  
  # Median-centre normalisation
  col_medians <- apply(sub_mat, 2, median, na.rm = TRUE)
  sub_mat <- sweep(sub_mat, 2, col_medians - median(col_medians, na.rm = TRUE))
  
  # limma
  cond <- factor(sub_info$condition, levels = c("Healthy", "Disease"))
  design <- model.matrix(~ 0 + cond)
  colnames(design) <- levels(cond)
  
  fit <- lmFit(sub_mat, design)
  contrast_mat <- makeContrasts(Disease - Healthy, levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)
  
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt$GenePhos <- rownames(tt)
  tt$Gene     <- gsub("_[STY].*$", "", tt$GenePhos)
  
  # Create output dir
  cdir <- file.path(output_dir, contrast_name)
  if (!dir.exists(cdir)) dir.create(cdir, recursive = TRUE)
  
  # ---- Export significance tiers ----
  write.csv(tt, file.path(cdir, "top_table_all.csv"), row.names = FALSE)
  for (tier_name in names(SIG_TIERS)) {
    tier <- SIG_TIERS[[tier_name]]
    sig <- tt[tt[[tier$col]] < tier$thresh, ]
    write.csv(sig, file.path(cdir, paste0("DE_", tier_name, ".csv")), row.names = FALSE)
  }
  
  n_sig  <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_up   <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0, na.rm = TRUE)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0, na.rm = TRUE)
  cat(sprintf("  FDR < 0.05: %d total (%d up, %d down)\n", n_sig, n_up, n_down))
  
  # ---- PCA ----
  tryCatch({
    mat_complete <- sub_mat[complete.cases(sub_mat), ]
    if (nrow(mat_complete) > 10) {
      pca <- prcomp(t(mat_complete), scale. = TRUE)
      pca_df <- data.frame(
        PC1 = pca$x[, 1], PC2 = pca$x[, 2],
        Condition = sub_info$condition,
        Sample = sub_info$sample
      )
      var_exp <- summary(pca)$importance[2, 1:2] * 100
      
      p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
        geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
        geom_text_repel(aes(label = Sample), size = 1.5, max.overlaps = 15,
                        segment.linewidth = 0.2, min.segment.length = 0) +
        scale_color_manual(values = CONDITION_COLORS) +
        labs(title = paste0("PCA — ", contrast_name),
             x = sprintf("PC1 (%.1f%%)", var_exp[1]),
             y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
        THEME_PLOT +
        theme(legend.position = "bottom")
      
      ggsave(file.path(cdir, "pca.png"), p_pca, width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
      ggsave(file.path(cdir, "pca.pdf"), p_pca, width = PLOT_W, height = PLOT_H, units = "in")
    }
  }, error = function(e) cat(sprintf("  ⚠ PCA failed: %s\n", e$message)))
  
  # ---- Volcano plots (two styles) ----
  tryCatch({
    tt_plot <- tt %>%
      mutate(
        sig_group = case_when(
          adj.P.Val < 0.05 & logFC > 1   ~ "Up (FDR<0.05, FC>2)",
          adj.P.Val < 0.05 & logFC < -1  ~ "Down (FDR<0.05, FC>2)",
          adj.P.Val < 0.05               ~ "Sig (FDR<0.05)",
          P.Value < 0.01                  ~ "Nominal (p<0.01)",
          TRUE                           ~ "NS"
        ),
        neg_log10p = -log10(pmax(P.Value, 1e-300))
      )
    
    # Smart top-15 labels
    top_labels <- tt_plot %>%
      filter(adj.P.Val < 0.05) %>%
      arrange(P.Value) %>%
      slice_head(n = 15)
    
    x_lim <- max(abs(tt_plot$logFC), na.rm = TRUE) * 1.1
    
    # ---- Style A: Canonical (8×8, theme_minimal, legend, FC lines) ----
    p_volc_canon <- ggplot(tt_plot, aes(logFC, neg_log10p)) +
      geom_point(aes(color = sig_group), alpha = 0.6, size = 1.2, stroke = 0) +
      scale_color_manual(values = VOLCANO_COLORS, name = NULL) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_text_repel(data = top_labels, aes(label = Gene),
                      size = 2.5, fontface = "italic", max.overlaps = 20,
                      segment.linewidth = 0.2, min.segment.length = 0) +
      xlim(-x_lim, x_lim) +
      labs(title = paste0("Volcano — ", contrast_name),
           subtitle = sprintf("FDR < 0.05: %d up, %d down", n_up, n_down),
           x = "log2 Fold Change (Disease − Healthy)",
           y = "-log10(p-value)") +
      theme_minimal(base_size = 11) +
      theme(
        panel.background  = element_rect(fill = "transparent", color = NA),
        plot.background   = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom"
      )
    
    ggsave(file.path(cdir, "volcano.png"), p_volc_canon,
           width = 8, height = 8, dpi = 300, bg = "transparent")
    ggsave(file.path(cdir, "volcano.pdf"), p_volc_canon,
           width = 8, height = 8, bg = "transparent")
    
    # ---- Style B: Compact (169×174 pt, theme_void, top 3, n= counts) ----
    n_up_fc   <- sum(tt_plot$sig_group == "Up (FDR<0.05, FC>2)")
    n_down_fc <- sum(tt_plot$sig_group == "Down (FDR<0.05, FC>2)")
    
    top3_labels <- tt_plot %>%
      filter(adj.P.Val < 0.05) %>%
      mutate(dir = ifelse(logFC > 0, "UP", "DOWN")) %>%
      group_by(dir) %>%
      slice_min(order_by = P.Value, n = 3) %>%
      ungroup()
    
    p_volc_compact <- ggplot(tt_plot, aes(logFC, neg_log10p)) +
      geom_point(aes(color = sig_group), alpha = 0.9, size = 1.5, stroke = 0) +
      scale_color_manual(values = VOLCANO_COLORS) +
      geom_vline(xintercept = 0, color = "grey30", linewidth = 0.2) +
      geom_hline(yintercept = -log10(0.05), color = "grey30",
                 linetype = "dashed", linewidth = 0.2) +
      annotate("text",
               x = -x_lim * 0.8,
               y = max(tt_plot$neg_log10p, na.rm = TRUE) * 0.98,
               label = paste0("n=", n_down), color = COL_HEALTHY,
               size = 5 / .pt, fontface = "bold") +
      annotate("text",
               x = x_lim * 0.8,
               y = max(tt_plot$neg_log10p, na.rm = TRUE) * 0.98,
               label = paste0("n=", n_up), color = COL_DCM,
               size = 5 / .pt, fontface = "bold") +
      geom_text_repel(data = top3_labels, aes(label = Gene),
                      size = 5 / .pt, color = "black", fontface = "italic",
                      box.padding = 0.2, point.padding = 0.3,
                      segment.linewidth = 0.15, min.segment.length = 0) +
      xlim(-x_lim, x_lim) +
      labs(x = "log2 FC (Disease − Healthy)", y = "-log10(p-value)") +
      THEME_PLOT +
      theme(legend.position = "none")
    
    ggsave(file.path(cdir, "volcano_compact.png"), p_volc_compact,
           width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
    ggsave(file.path(cdir, "volcano_compact.pdf"), p_volc_compact,
           width = PLOT_W, height = PLOT_H, units = "in")
  }, error = function(e) cat(sprintf("  ⚠ Volcano failed: %s\n", e$message)))
  
  # ---- fgsea (Hallmarks, KEGG, Reactome) ----
  tryCatch({
    # Rank by signed -log10(p) using Gene (not GenePhos)
    rank_df <- tt %>%
      mutate(rank_score = sign(logFC) * (-log10(pmax(P.Value, 1e-300)))) %>%
      group_by(Gene) %>%
      slice_min(order_by = P.Value, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    ranks <- setNames(rank_df$rank_score, rank_df$Gene)
    ranks <- sort(ranks, decreasing = TRUE)
    
    for (gs_name in c("Hallmarks", "KEGG", "Reactome")) {
      gs_list <- switch(gs_name,
                        "Hallmarks" = hallmarks_list,
                        "KEGG"      = kegg_list,
                        "Reactome"  = reactome_list)
      
      fgsea_res <- fgsea(pathways = gs_list, stats = ranks, minSize = 10, maxSize = 500)
      fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
      fgsea_out <- as.data.frame(fgsea_res)
      fgsea_out$leadingEdge <- sapply(fgsea_out$leadingEdge, paste, collapse = ";")
      write.csv(fgsea_out, file.path(cdir, paste0("fgsea_", gs_name, ".csv")), row.names = FALSE)
      
      # Barplot for top pathways
      if (gs_name == "Hallmarks" && nrow(fgsea_res[fgsea_res$padj < 0.25, ]) > 0) {
        top_paths <- fgsea_res %>%
          filter(padj < 0.25) %>%
          arrange(NES) %>%
          mutate(
            pathway = gsub("^HALLMARK_", "", pathway),
            pathway = factor(pathway, levels = pathway),
            direction = ifelse(NES > 0, "Up in DCM", "Down in DCM")
          )
        
        p_bar <- ggplot(top_paths, aes(NES, pathway, fill = direction)) +
          geom_col(alpha = 0.85, width = 0.7) +
          scale_fill_manual(values = DIRECTION_COLORS) +
          geom_vline(xintercept = 0, linewidth = 0.2) +
          labs(title = paste0("Hallmarks — ", contrast_name),
               x = "NES", y = NULL, fill = NULL) +
          THEME_PLOT +
          theme(legend.position = "bottom")
        
        ggsave(file.path(cdir, "fgsea_hallmarks_barplot.png"), p_bar,
               width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
        ggsave(file.path(cdir, "fgsea_hallmarks_barplot.pdf"), p_bar,
               width = PLOT_W, height = PLOT_H, units = "in")
      }
    }
    cat("  ✓ fgsea (Hallmarks/KEGG/Reactome)\n")
  }, error = function(e) cat(sprintf("  ⚠ fgsea failed: %s\n", e$message)))
  
  # ---- GO enrichment (BP, MF, CC) ----
  tryCatch({
    sig_genes <- tt$Gene[tt$adj.P.Val < 0.05]
    up_genes  <- tt$Gene[tt$adj.P.Val < 0.05 & tt$logFC > 0]
    down_genes <- tt$Gene[tt$adj.P.Val < 0.05 & tt$logFC < 0]
    bg_genes  <- unique(tt$Gene)
    
    for (ont in c("BP", "MF", "CC")) {
      for (dir_label in c("all", "up", "down")) {
        gene_set <- switch(dir_label,
                           "all"  = sig_genes,
                           "up"   = up_genes,
                           "down" = down_genes)
        if (length(gene_set) < 5) next
        
        ego <- tryCatch(
          enrichGO(gene = unique(gene_set), universe = bg_genes,
                   OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = ont, pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2),
          error = function(e) NULL)
        
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          write.csv(as.data.frame(ego),
                    file.path(cdir, sprintf("GO_%s_%s.csv", ont, dir_label)),
                    row.names = FALSE)
        }
      }
    }
    cat("  ✓ GO enrichment (BP/MF/CC × all/up/down)\n")
  }, error = function(e) cat(sprintf("  ⚠ GO failed: %s\n", e$message)))
  
  # Return results for downstream
  list(
    contrast_name = contrast_name,
    gene = gene, tp = tp,
    top_table = tt,
    sub_mat = sub_mat,
    sub_info = sub_info,
    n_sig = n_sig, n_up = n_up, n_down = n_down,
    nD = n_dis, nH = n_heal
  )
}

# ---- 9b: KSEA Volcano — Original style (scatter, size = n_substrates) ----
plot_ksea_volcano_original <- function(ksea_df, contrast_name, output_dir) {
  
  df <- ksea_df %>%
    mutate(
      neg_log10_fdr = -log10(pmax(fdr, 1e-10)),
      label         = ifelse(fdr < 0.05, kinase, ""),
      direction     = ifelse(activity > 0, "Active", "Inactive")
    )
  
  p <- ggplot(df, aes(x = activity, y = neg_log10_fdr, color = direction, label = label)) +
    geom_point(aes(size = n_substrates), alpha = 0.7, stroke = 0) +
    geom_text_repel(size = 3, max.overlaps = 20, fontface = "bold",
                    segment.linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = c("Active" = COL_DCM, "Inactive" = COL_HEALTHY)) +
    scale_size_continuous(range = c(2, 8), name = "n substrates") +
    labs(
      title    = paste0("KSEA Volcano — ", contrast_name),
      subtitle = "Dashed line: FDR = 0.05",
      x        = "KSEA z-score (activity)",
      y        = "-log10(FDR)",
      color    = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      aspect.ratio = 1
    )
  
  ggsave(file.path(output_dir, "ksea_volcano_original.png"), p,
         width = 8, height = 8, dpi = 300, bg = "transparent")
  ggsave(file.path(output_dir, "ksea_volcano_original.pdf"), p,
         width = 8, height = 8, bg = "transparent")
  cat("  ✓ KSEA volcano (original style)\n")
}

# ---- 9c: KSEA Volcano — Compact style (theme_void, top 3, n= counts) ----
plot_ksea_volcano_compact <- function(ksea_df, contrast_name, output_dir) {
  
  df <- ksea_df %>%
    mutate(
      neg_log10_p = -log10(pmax(p_value, 1e-300)),
      group = case_when(
        p_value < 0.05 & z_score > 0 ~ "UP",
        p_value < 0.05 & z_score < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  
  n_up   <- sum(df$group == "UP")
  n_down <- sum(df$group == "DOWN")
  
  # Top 3 most significant per group
  top_labels <- df %>%
    filter(group != "NS") %>%
    group_by(group) %>%
    slice_min(order_by = p_value, n = 3) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = z_score, y = neg_log10_p)) +
    geom_point(aes(color = group), alpha = 0.9, size = 1.5, stroke = 0) +
    scale_color_manual(values = c("UP" = COL_DCM, "DOWN" = COL_HEALTHY, "NS" = COL_NS)) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "grey30",
               linetype = "dashed", linewidth = 0.2) +
    annotate("text",
             x = min(df$z_score, na.rm = TRUE) * 0.8,
             y = max(df$neg_log10_p, na.rm = TRUE) * 0.98,
             label = paste0("n=", n_down), color = COL_HEALTHY,
             size = 5 / .pt, fontface = "bold") +
    annotate("text",
             x = max(df$z_score, na.rm = TRUE) * 0.8,
             y = max(df$neg_log10_p, na.rm = TRUE) * 0.98,
             label = paste0("n=", n_up), color = COL_DCM,
             size = 5 / .pt, fontface = "bold") +
    geom_text_repel(data = top_labels, aes(label = kinase),
                    size = 5 / .pt, color = "black", fontface = "italic",
                    box.padding = 0.2, point.padding = 0.3,
                    segment.linewidth = 0.15, min.segment.length = 0) +
    labs(x = "Kinase Activity (z-score)", y = "-log10(p-value)") +
    THEME_PLOT +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "ksea_volcano_compact.pdf"), p,
         width = PLOT_W, height = PLOT_H, units = "in", device = "pdf")
  ggsave(file.path(output_dir, "ksea_volcano_compact.png"), p,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  cat("  ✓ KSEA volcano (compact style)\n")
}

# ---- 9d: run_kinase_for_contrast — wrapper calling v14 + both volcanos ----
run_kinase_for_contrast <- function(results, output_dir) {
  
  cdir <- file.path(output_dir, results$contrast_name, "kinase_activity")
  
  tryCatch({
    kinase_res <- run_kinase_activity(
      results      = results,
      min_substrates = MIN_SUBSTRATES,
      min_resources  = MIN_RESOURCES,
      methods      = c("ksea", "ulm", "wmean"),
      top_n_kinases = 25,
      create_plots  = TRUE,
      output_dir   = cdir
    )
    
    # Both KSEA volcano styles (if KSEA ran)
    if (!is.null(kinase_res) && "ksea" %in% names(kinase_res$per_method)) {
      ksea_df <- kinase_res$per_method[["ksea"]]
      
      # Original style
      plot_ksea_volcano_original(ksea_df, results$contrast_name, cdir)
      
      # Compact style — needs columns: z_score, p_value, kinase
      compact_df <- ksea_df %>%
        rename(z_score = activity, p_value = pval)
      if (!"p_value" %in% colnames(compact_df) && "fdr" %in% colnames(compact_df)) {
        # fall back to using raw p if pval not available
        compact_df$p_value <- compact_df$fdr
      }
      plot_ksea_volcano_compact(compact_df, results$contrast_name, cdir)
    }
    
    return(kinase_res)
  }, error = function(e) {
    cat(sprintf("  ⚠ Kinase inference failed for %s: %s\n",
                results$contrast_name, e$message))
    return(NULL)
  })
}

# ==============================================================================
# SECTION 10: RUN ALL 11 CONTRASTS
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  RUNNING PER-GENOTYPE PER-TIMEPOINT CONTRASTS\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

if (!dir.exists(OUTPUT_BASE)) dir.create(OUTPUT_BASE, recursive = TRUE)

disease_results <- list()
kinase_results  <- list()

for (gene in names(DISEASE_LINES)) {
  for (tp in TIMEPOINTS) {
    contrast_name <- paste0(gene, "_", tp)
    
    # Skip DSP D22
    if (contrast_name %in% SKIP_CONTRASTS) {
      cat(sprintf("\n  ── %s ── SKIPPED (insufficient healthy replicates)\n", contrast_name))
      next
    }
    
    # Run DE
    res <- run_contrast(gene, tp, log2_matrix, sample_info, OUTPUT_BASE)
    
    if (!is.null(res)) {
      disease_results[[contrast_name]] <- res
      
      # Run kinase inference + both volcano styles
      kinase_results[[contrast_name]] <- run_kinase_for_contrast(res, OUTPUT_BASE)
    }
  }
}

cat(sprintf("\n  ✓ Completed %d / 11 contrasts\n", length(disease_results)))

# ==============================================================================
# SECTION 11: GLOBAL HEATMAPS PER TIMEPOINT
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  GLOBAL HEATMAPS\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

heatmap_dir <- file.path(OUTPUT_BASE, "Heatmaps")
if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)

for (tp in TIMEPOINTS) {
  tryCatch({
    # Collect FDR < 0.05 phosphosites across all genotypes at this timepoint
    tp_contrasts <- grep(paste0("_", tp, "$"), names(disease_results), value = TRUE)
    if (length(tp_contrasts) == 0) next
    
    sig_sites <- unique(unlist(lapply(tp_contrasts, function(cn) {
      tt <- disease_results[[cn]]$top_table
      tt$GenePhos[tt$adj.P.Val < 0.05]
    })))
    
    if (length(sig_sites) < 5) {
      cat(sprintf("  %s: only %d sig sites — skipping heatmap\n", tp, length(sig_sites)))
      next
    }
    
    # Get samples for this timepoint
    tp_samples <- sample_info$sample[sample_info$timepoint == tp]
    tp_mat <- log2_matrix[sig_sites[sig_sites %in% rownames(log2_matrix)], tp_samples, drop = FALSE]
    
    # Impute NAs with row median for heatmap only
    tp_mat_imp <- tp_mat
    for (i in seq_len(nrow(tp_mat_imp))) {
      na_idx <- is.na(tp_mat_imp[i, ])
      if (any(na_idx) && !all(na_idx)) {
        tp_mat_imp[i, na_idx] <- median(tp_mat_imp[i, !na_idx])
      }
    }
    tp_mat_imp <- tp_mat_imp[complete.cases(tp_mat_imp), ]
    
    # Row scale
    tp_scaled <- t(scale(t(tp_mat_imp)))
    
    # Annotation
    tp_info <- sample_info[sample_info$timepoint == tp, ]
    tp_info <- tp_info[match(colnames(tp_scaled), tp_info$sample), ]
    
    # Order: genotype → condition (Healthy first) → sample
    col_order <- order(tp_info$genotype, tp_info$condition != "Healthy", tp_info$sample)
    tp_scaled <- tp_scaled[, col_order]
    tp_info   <- tp_info[col_order, ]
    
    anno_col <- data.frame(
      Condition = tp_info$condition,
      Genotype  = tp_info$genotype,
      row.names = tp_info$sample
    )
    
    anno_colors <- list(
      Condition = CONDITION_COLORS,
      Genotype  = c("TTN" = "#1b9e77", "DSP" = "#d95f02",
                    "TNNT2" = "#7570b3", "TPM1" = "#e7298a")
    )
    
    # Gaps between genotypes
    geno_rle <- rle(tp_info$genotype)
    gaps_col <- cumsum(geno_rle$lengths[-length(geno_rle$lengths)])
    
    # Cap height
    h <- min(49, max(6, nrow(tp_scaled) * 0.015 + 3))
    w <- max(6, ncol(tp_scaled) * 0.25 + 3)
    
    png(file.path(heatmap_dir, paste0("heatmap_global_", tp, ".png")),
        width = w, height = h, units = "in", res = 300, bg = "transparent")
    pheatmap(tp_scaled,
             cluster_rows = TRUE, cluster_cols = FALSE,
             show_rownames = FALSE, show_colnames = TRUE,
             annotation_col = anno_col, annotation_colors = anno_colors,
             gaps_col = gaps_col,
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             main = sprintf("Phosphoproteomics — %s (FDR < 0.05, %d sites)", tp, nrow(tp_scaled)),
             fontsize = 8, fontsize_col = 6)
    dev.off()
    
    pdf(file.path(heatmap_dir, paste0("heatmap_global_", tp, ".pdf")),
        width = w, height = h)
    pheatmap(tp_scaled,
             cluster_rows = TRUE, cluster_cols = FALSE,
             show_rownames = FALSE, show_colnames = TRUE,
             annotation_col = anno_col, annotation_colors = anno_colors,
             gaps_col = gaps_col,
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             main = sprintf("Phosphoproteomics — %s (FDR < 0.05, %d sites)", tp, nrow(tp_scaled)),
             fontsize = 8, fontsize_col = 6)
    dev.off()
    
    cat(sprintf("  ✓ %s: %d sig sites × %d samples\n", tp, nrow(tp_scaled), ncol(tp_scaled)))
  }, error = function(e) cat(sprintf("  ⚠ Heatmap %s failed: %s\n", tp, e$message)))
}

# ==============================================================================
# SECTION 12: TOP 50 HEATMAPS PER GENOTYPE
# ==============================================================================

cat("\n  Top 50 heatmaps (25 up + 25 down) per genotype...\n")

for (gene in names(DISEASE_LINES)) {
  tryCatch({
    # Use D29 as the primary timepoint (most mature)
    primary_cn <- paste0(gene, "_D29")
    if (!primary_cn %in% names(disease_results)) {
      # Fall back to D22 or D15
      primary_cn <- paste0(gene, "_D22")
      if (!primary_cn %in% names(disease_results)) primary_cn <- paste0(gene, "_D15")
    }
    if (!primary_cn %in% names(disease_results)) next
    
    tt <- disease_results[[primary_cn]]$top_table
    
    top_up   <- tt %>% filter(logFC > 0) %>% arrange(P.Value) %>% slice_head(n = 25)
    top_down <- tt %>% filter(logFC < 0) %>% arrange(P.Value) %>% slice_head(n = 25)
    top_sites <- c(top_up$GenePhos, top_down$GenePhos)
    
    if (length(top_sites) < 5) next
    
    # Get all timepoint samples for this genotype
    gene_samples <- sample_info$sample[sample_info$genotype == gene]
    gene_mat <- log2_matrix[top_sites[top_sites %in% rownames(log2_matrix)], gene_samples, drop = FALSE]
    
    # Impute + scale
    gene_mat_imp <- gene_mat
    for (i in seq_len(nrow(gene_mat_imp))) {
      na_idx <- is.na(gene_mat_imp[i, ])
      if (any(na_idx) && !all(na_idx)) {
        gene_mat_imp[i, na_idx] <- median(gene_mat_imp[i, !na_idx])
      }
    }
    gene_mat_imp <- gene_mat_imp[complete.cases(gene_mat_imp), ]
    gene_scaled <- t(scale(t(gene_mat_imp)))
    
    # Annotation
    gene_info <- sample_info[sample_info$genotype == gene, ]
    gene_info <- gene_info[match(colnames(gene_scaled), gene_info$sample), ]
    
    # Order: timepoint → condition
    col_order <- order(gene_info$timepoint, gene_info$condition != "Healthy")
    gene_scaled <- gene_scaled[, col_order]
    gene_info   <- gene_info[col_order, ]
    
    anno_col <- data.frame(
      Condition = gene_info$condition,
      Timepoint = gene_info$timepoint,
      row.names = gene_info$sample
    )
    anno_colors <- list(
      Condition = CONDITION_COLORS,
      Timepoint = c("D15" = "#66c2a5", "D22" = "#fc8d62", "D29" = "#8da0cb")
    )
    
    # Gaps between timepoints
    tp_rle <- rle(gene_info$timepoint)
    gaps_col <- cumsum(tp_rle$lengths[-length(tp_rle$lengths)])
    
    # Row labels: show Gene name
    row_labels <- gsub("_[STY].*$", "", rownames(gene_scaled))
    
    h <- min(49, max(6, nrow(gene_scaled) * 0.22 + 3))
    w <- max(6, ncol(gene_scaled) * 0.35 + 3)
    
    png(file.path(heatmap_dir, paste0("heatmap_top50_", gene, ".png")),
        width = w, height = h, units = "in", res = 300, bg = "transparent")
    pheatmap(gene_scaled,
             cluster_rows = TRUE, cluster_cols = FALSE,
             labels_row = row_labels, show_colnames = TRUE,
             annotation_col = anno_col, annotation_colors = anno_colors,
             gaps_col = gaps_col,
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             main = sprintf("%s — Top 50 phosphosites (25 up + 25 down)", gene),
             fontsize = 8, fontsize_row = 6, fontsize_col = 6)
    dev.off()
    
    pdf(file.path(heatmap_dir, paste0("heatmap_top50_", gene, ".pdf")),
        width = w, height = h)
    pheatmap(gene_scaled,
             cluster_rows = TRUE, cluster_cols = FALSE,
             labels_row = row_labels, show_colnames = TRUE,
             annotation_col = anno_col, annotation_colors = anno_colors,
             gaps_col = gaps_col,
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             main = sprintf("%s — Top 50 phosphosites (25 up + 25 down)", gene),
             fontsize = 8, fontsize_row = 6, fontsize_col = 6)
    dev.off()
    
    cat(sprintf("  ✓ %s: %d sites across %d samples\n", gene, nrow(gene_scaled), ncol(gene_scaled)))
  }, error = function(e) cat(sprintf("  ⚠ Top50 heatmap %s failed: %s\n", gene, e$message)))
}

# ==============================================================================
# SECTION 13: CROSS-TIMEPOINT CONTRASTS
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  CROSS-TIMEPOINT CONTRASTS\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

tp_contrast_dir <- file.path(OUTPUT_BASE, "Timepoint_Contrasts")
if (!dir.exists(tp_contrast_dir)) dir.create(tp_contrast_dir, recursive = TRUE)

for (gene in names(DISEASE_LINES)) {
  tryCatch({
    # Get Disease samples only for this genotype
    gene_dis_info <- sample_info[sample_info$genotype == gene &
                                   sample_info$condition == "Disease", ]
    
    if (nrow(gene_dis_info) < 4) next
    
    gene_dis_mat <- log2_matrix[, gene_dis_info$sample, drop = FALSE]
    
    # Filter
    keep <- rowSums(!is.na(gene_dis_mat)) >= FILTER_THRESHOLD
    gene_dis_mat <- gene_dis_mat[keep, , drop = FALSE]
    
    # Median centre
    col_meds <- apply(gene_dis_mat, 2, median, na.rm = TRUE)
    gene_dis_mat <- sweep(gene_dis_mat, 2, col_meds - median(col_meds, na.rm = TRUE))
    
    tp_factor <- factor(gene_dis_info$timepoint, levels = c("D15", "D22", "D29"))
    design <- model.matrix(~ 0 + tp_factor)
    colnames(design) <- levels(tp_factor)
    
    fit <- lmFit(gene_dis_mat, design)
    
    # Pairwise timepoint contrasts
    tp_pairs <- list(
      "D22_vs_D15" = "D22 - D15",
      "D29_vs_D15" = "D29 - D15",
      "D29_vs_D22" = "D29 - D22"
    )
    
    for (pair_name in names(tp_pairs)) {
      contrast_mat <- makeContrasts(contrasts = tp_pairs[[pair_name]], levels = design)
      fit2 <- contrasts.fit(fit, contrast_mat)
      fit2 <- eBayes(fit2)
      
      tt <- topTable(fit2, number = Inf, sort.by = "none")
      tt$GenePhos <- rownames(tt)
      tt$Gene     <- gsub("_[STY].*$", "", tt$GenePhos)
      
      out_dir <- file.path(tp_contrast_dir, paste0(gene, "_", pair_name))
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      write.csv(tt, file.path(out_dir, "top_table_all.csv"), row.names = FALSE)
      
      for (tier_name in names(SIG_TIERS)) {
        tier <- SIG_TIERS[[tier_name]]
        sig <- tt[tt[[tier$col]] < tier$thresh, ]
        write.csv(sig, file.path(out_dir, paste0("DE_", tier_name, ".csv")), row.names = FALSE)
      }
      
      cat(sprintf("  ✓ %s_%s: %d FDR < 0.05\n",
                  gene, pair_name, sum(tt$adj.P.Val < 0.05, na.rm = TRUE)))
    }
  }, error = function(e) cat(sprintf("  ⚠ Timepoint contrasts %s failed: %s\n", gene, e$message)))
}

# ==============================================================================
# SECTION 14: SCATTERPLOTS & CORRELATION MATRICES
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  SCATTERPLOTS & CORRELATIONS\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

scatter_dir <- file.path(OUTPUT_BASE, "Scatterplots")
if (!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)

# ---- 14a: Within-genotype: logFC correlations across timepoints ----
for (gene in names(DISEASE_LINES)) {
  tryCatch({
    gene_contrasts <- grep(paste0("^", gene, "_D"), names(disease_results), value = TRUE)
    if (length(gene_contrasts) < 2) next
    
    # Merge logFC tables
    merged <- NULL
    for (cn in gene_contrasts) {
      tt <- disease_results[[cn]]$top_table[, c("GenePhos", "logFC")]
      colnames(tt)[2] <- cn
      if (is.null(merged)) merged <- tt
      else merged <- merge(merged, tt, by = "GenePhos", all = TRUE)
    }
    
    # Pairwise scatterplots
    for (i in 1:(length(gene_contrasts) - 1)) {
      for (j in (i + 1):length(gene_contrasts)) {
        cn1 <- gene_contrasts[i]; cn2 <- gene_contrasts[j]
        df <- merged[complete.cases(merged[, c(cn1, cn2)]), ]
        if (nrow(df) < 10) next
        
        cor_val <- cor(df[[cn1]], df[[cn2]], method = "pearson")
        
        df$direction <- case_when(
          df[[cn1]] > 0 & df[[cn2]] > 0 ~ "Both Up",
          df[[cn1]] < 0 & df[[cn2]] < 0 ~ "Both Down",
          TRUE ~ "Discordant"
        )
        
        p <- ggplot(df, aes_string(x = cn1, y = cn2)) +
          geom_point(aes(color = direction), alpha = 0.4, size = 0.5, stroke = 0) +
          scale_color_manual(values = c("Both Up" = COL_DCM, "Both Down" = COL_HEALTHY,
                                        "Discordant" = "#CC3333")) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.2) +
          geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.3) +
          annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
                   label = sprintf("r = %.3f", cor_val), size = 5 / .pt, fontface = "bold") +
          labs(title = sprintf("%s: %s vs %s", gene, cn1, cn2),
               x = paste0("logFC (", cn1, ")"),
               y = paste0("logFC (", cn2, ")"),
               color = NULL) +
          THEME_PLOT +
          theme(legend.position = "none")
        
        fname <- paste0("scatter_", gene, "_", gsub(paste0(gene, "_"), "", cn1),
                        "_vs_", gsub(paste0(gene, "_"), "", cn2))
        ggsave(file.path(scatter_dir, paste0(fname, ".png")), p,
               width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
        ggsave(file.path(scatter_dir, paste0(fname, ".pdf")), p,
               width = PLOT_W, height = PLOT_H, units = "in")
      }
    }
    cat(sprintf("  ✓ %s within-genotype scatterplots\n", gene))
  }, error = function(e) cat(sprintf("  ⚠ Scatter %s failed: %s\n", gene, e$message)))
}

# ---- 14b: Across-genotype correlation heatmaps per timepoint ----
for (tp in TIMEPOINTS) {
  tryCatch({
    tp_contrasts <- grep(paste0("_", tp, "$"), names(disease_results), value = TRUE)
    if (length(tp_contrasts) < 2) next
    
    merged <- NULL
    for (cn in tp_contrasts) {
      tt <- disease_results[[cn]]$top_table[, c("GenePhos", "logFC")]
      colnames(tt)[2] <- cn
      if (is.null(merged)) merged <- tt
      else merged <- merge(merged, tt, by = "GenePhos", all = TRUE)
    }
    
    merged_complete <- merged[complete.cases(merged), ]
    if (nrow(merged_complete) < 20) next
    
    cor_mat <- cor(merged_complete[, -1], method = "pearson")
    
    png(file.path(scatter_dir, paste0("cor_heatmap_across_", tp, ".png")),
        width = 6, height = 5, units = "in", res = 300, bg = "transparent")
    pheatmap(cor_mat,
             display_numbers = TRUE, number_format = "%.2f",
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = paste0("Cross-genotype logFC correlation — ", tp),
             fontsize = 10)
    dev.off()
    
    pdf(file.path(scatter_dir, paste0("cor_heatmap_across_", tp, ".pdf")),
        width = 6, height = 5)
    pheatmap(cor_mat,
             display_numbers = TRUE, number_format = "%.2f",
             color = colorRampPalette(c(COL_HEALTHY, "white", COL_DCM))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = paste0("Cross-genotype logFC correlation — ", tp),
             fontsize = 10)
    dev.off()
    
    cat(sprintf("  ✓ %s across-genotype correlation heatmap\n", tp))
  }, error = function(e) cat(sprintf("  ⚠ Cross-genotype cor %s failed: %s\n", tp, e$message)))
}

# ==============================================================================
# SECTION 15: INTERACTION MODEL (Disease × Timepoint)
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  INTERACTION MODEL (Disease × Timepoint)\n")
cat("════════════════════════════════════════════════════════════════════════════════\n\n")

# Tests whether the DCM effect changes magnitude from D15 → D22 or D15 → D29.
# Model: ~ 0 + Condition:Timepoint  with contrasts extracting interaction terms.
# A significant interaction means the disease effect is timepoint-dependent.

interaction_dir <- file.path(OUTPUT_BASE, "Interaction_Model")
if (!dir.exists(interaction_dir)) dir.create(interaction_dir, recursive = TRUE)

for (gene in names(DISEASE_LINES)) {
  tryCatch({
    gene_info <- sample_info[sample_info$genotype == gene, ]
    
    # Need at least 2 timepoints with >= 2 disease + 1 healthy
    usable_tps <- c()
    for (tp in TIMEPOINTS) {
      tp_info <- gene_info[gene_info$timepoint == tp, ]
      if (sum(tp_info$condition == "Disease") >= 2 && sum(tp_info$condition == "Healthy") >= 1) {
        usable_tps <- c(usable_tps, tp)
      }
    }
    if (length(usable_tps) < 2) {
      cat(sprintf("  %s: only %d usable timepoints — skipping interaction\n", gene, length(usable_tps)))
      next
    }
    
    gene_info_sub <- gene_info[gene_info$timepoint %in% usable_tps, ]
    gene_mat <- log2_matrix[, gene_info_sub$sample, drop = FALSE]
    
    # Filter
    keep <- rowSums(!is.na(gene_mat)) >= FILTER_THRESHOLD
    gene_mat <- gene_mat[keep, , drop = FALSE]
    
    # Median centre
    col_meds <- apply(gene_mat, 2, median, na.rm = TRUE)
    gene_mat <- sweep(gene_mat, 2, col_meds - median(col_meds, na.rm = TRUE))
    
    cond <- factor(gene_info_sub$condition, levels = c("Healthy", "Disease"))
    tp   <- factor(gene_info_sub$timepoint, levels = usable_tps)
    
    group <- factor(paste0(cond, "_", tp))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    fit <- lmFit(gene_mat, design)
    
    # Interaction contrasts: (Disease_D22 - Healthy_D22) - (Disease_D15 - Healthy_D15)
    ref_tp <- usable_tps[1]  # D15
    for (test_tp in usable_tps[-1]) {
      contrast_str <- sprintf("(Disease_%s - Healthy_%s) - (Disease_%s - Healthy_%s)",
                              test_tp, test_tp, ref_tp, ref_tp)
      contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
      fit2 <- contrasts.fit(fit, contrast_mat)
      fit2 <- eBayes(fit2)
      
      tt <- topTable(fit2, number = Inf, sort.by = "none")
      tt$GenePhos <- rownames(tt)
      tt$Gene     <- gsub("_[STY].*$", "", tt$GenePhos)
      
      write.csv(tt, file.path(interaction_dir,
                              sprintf("interaction_%s_%s_vs_%s.csv", gene, test_tp, ref_tp)),
                row.names = FALSE)
      
      n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
      cat(sprintf("  ✓ %s: %s vs %s interaction — %d FDR < 0.05\n",
                  gene, test_tp, ref_tp, n_sig))
    }
  }, error = function(e) cat(sprintf("  ⚠ Interaction %s failed: %s\n", gene, e$message)))
}

# ==============================================================================
# SECTION 16: CORE DCM iPSC PHOSPHO SIGNATURE PER TIMEPOINT
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  CORE DCM iPSC PHOSPHO SIGNATURE\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

core_dir <- file.path(OUTPUT_BASE, "Core_Signatures")
if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)

for (tp in TIMEPOINTS) {
  tryCatch({
    tp_contrasts <- grep(paste0("_", tp, "$"), names(disease_results), value = TRUE)
    if (length(tp_contrasts) < 2) {
      cat(sprintf("  %s: only %d contrasts — need >= 2 for core signature\n", tp, length(tp_contrasts)))
      next
    }
    
    # Phosphosites significant (FDR < 0.05) in ALL genotypes at this timepoint
    sig_lists <- lapply(tp_contrasts, function(cn) {
      tt <- disease_results[[cn]]$top_table
      tt$GenePhos[tt$adj.P.Val < 0.05]
    })
    
    core_sites <- Reduce(intersect, sig_lists)
    cat(sprintf("  %s: %d sites significant in all %d genotypes\n",
                tp, length(core_sites), length(tp_contrasts)))
    
    if (length(core_sites) > 0) {
      # Extract statistics from each contrast
      core_df <- data.frame(GenePhos = core_sites, stringsAsFactors = FALSE)
      core_df$Gene <- gsub("_[STY].*$", "", core_df$GenePhos)
      
      for (cn in tp_contrasts) {
        tt <- disease_results[[cn]]$top_table
        tt_core <- tt[tt$GenePhos %in% core_sites, c("GenePhos", "logFC", "adj.P.Val")]
        colnames(tt_core)[2:3] <- paste0(c("logFC_", "FDR_"), cn)
        core_df <- merge(core_df, tt_core, by = "GenePhos", all.x = TRUE)
      }
      
      # Direction consistency
      logfc_cols <- grep("^logFC_", colnames(core_df), value = TRUE)
      core_df$mean_logFC <- rowMeans(core_df[, logfc_cols], na.rm = TRUE)
      core_df$direction  <- ifelse(core_df$mean_logFC > 0, "Up in DCM", "Down in DCM")
      
      core_df <- core_df[order(core_df$mean_logFC), ]
      write.csv(core_df, file.path(core_dir, sprintf("core_DCM_iPSC_phospho_%s.csv", tp)),
                row.names = FALSE)
    }
  }, error = function(e) cat(sprintf("  ⚠ Core signature %s failed: %s\n", tp, e$message)))
}

# ==============================================================================
# SECTION 17: MATURATION TRAJECTORY (PCA + UMAP)
# ==============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("  MATURATION TRAJECTORY\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

traj_dir <- file.path(OUTPUT_BASE, "Maturation_Trajectory")
if (!dir.exists(traj_dir)) dir.create(traj_dir, recursive = TRUE)

tryCatch({
  # Use all samples, complete cases only
  traj_mat <- log2_matrix[complete.cases(log2_matrix), ]
  if (nrow(traj_mat) < 50) {
    # Relax: allow up to 20% NA per row
    max_na <- ceiling(ncol(log2_matrix) * 0.2)
    keep <- rowSums(is.na(log2_matrix)) <= max_na
    traj_mat <- log2_matrix[keep, ]
    # Impute remaining NAs with row median
    for (i in seq_len(nrow(traj_mat))) {
      na_idx <- is.na(traj_mat[i, ])
      if (any(na_idx)) traj_mat[i, na_idx] <- median(traj_mat[i, !na_idx], na.rm = TRUE)
    }
  }
  cat(sprintf("  Trajectory matrix: %d phosphosites × %d samples\n",
              nrow(traj_mat), ncol(traj_mat)))
  
  # Colour palettes
  tp_colors   <- c("D15" = "#66c2a5", "D22" = "#fc8d62", "D29" = "#8da0cb")
  gene_colors <- c("TTN" = "#1b9e77", "DSP" = "#d95f02",
                   "TNNT2" = "#7570b3", "TPM1" = "#e7298a")
  
  traj_info <- sample_info[match(colnames(traj_mat), sample_info$sample), ]
  
  # ---- PCA ----
  pca <- prcomp(t(traj_mat), scale. = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    Timepoint = traj_info$timepoint,
    Genotype  = traj_info$genotype,
    Condition = traj_info$condition,
    Sample    = traj_info$sample
  )
  
  # PCA by timepoint
  p_pca_tp <- ggplot(pca_df, aes(PC1, PC2, color = Timepoint, shape = Condition)) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = tp_colors) +
    labs(title = "PCA — Phospho (Timepoint)",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "pca_by_timepoint.png"), p_pca_tp,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "pca_by_timepoint.pdf"), p_pca_tp,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  # PCA by genotype
  p_pca_gene <- ggplot(pca_df, aes(PC1, PC2, color = Genotype, shape = Condition)) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = gene_colors) +
    labs(title = "PCA — Phospho (Genotype)",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "pca_by_genotype.png"), p_pca_gene,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "pca_by_genotype.pdf"), p_pca_gene,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  # PCA faceted by genotype
  p_pca_facet <- ggplot(pca_df, aes(PC1, PC2, color = Timepoint, shape = Condition)) +
    geom_point(size = 1, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = tp_colors) +
    facet_wrap(~ Genotype) +
    labs(title = "PCA by Genotype",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "pca_faceted_by_genotype.png"), p_pca_facet,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "pca_faceted_by_genotype.pdf"), p_pca_facet,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  cat("  ✓ PCA (timepoint, genotype, faceted)\n")
  
  # ---- UMAP ----
  n_neighbors <- min(15, ncol(traj_mat) - 1)
  set.seed(42)
  umap_res <- umap(t(traj_mat), n_neighbors = n_neighbors)
  
  umap_df <- data.frame(
    UMAP1 = umap_res$layout[, 1], UMAP2 = umap_res$layout[, 2],
    Timepoint = traj_info$timepoint,
    Genotype  = traj_info$genotype,
    Condition = traj_info$condition,
    Sample    = traj_info$sample
  )
  
  # UMAP by timepoint
  p_umap_tp <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Timepoint, shape = Condition)) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = tp_colors) +
    labs(title = "UMAP — Phospho (Timepoint)",
         x = "UMAP1", y = "UMAP2") +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "umap_by_timepoint.png"), p_umap_tp,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "umap_by_timepoint.pdf"), p_umap_tp,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  # UMAP by genotype
  p_umap_gene <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Genotype, shape = Condition)) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = gene_colors) +
    labs(title = "UMAP — Phospho (Genotype)",
         x = "UMAP1", y = "UMAP2") +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "umap_by_genotype.png"), p_umap_gene,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "umap_by_genotype.pdf"), p_umap_gene,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  # UMAP faceted by genotype
  p_umap_facet <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Timepoint, shape = Condition)) +
    geom_point(size = 1, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = tp_colors) +
    facet_wrap(~ Genotype) +
    labs(title = "UMAP by Genotype",
         x = "UMAP1", y = "UMAP2") +
    THEME_PLOT + theme(legend.position = "bottom")
  
  ggsave(file.path(traj_dir, "umap_faceted_by_genotype.png"), p_umap_facet,
         width = PLOT_W, height = PLOT_H, units = "in", dpi = 300)
  ggsave(file.path(traj_dir, "umap_faceted_by_genotype.pdf"), p_umap_facet,
         width = PLOT_W, height = PLOT_H, units = "in")
  
  cat("  ✓ UMAP (timepoint, genotype, faceted)\n")
  
  # Export coordinates
  write.csv(pca_df, file.path(traj_dir, "pca_coordinates.csv"), row.names = FALSE)
  write.csv(umap_df, file.path(traj_dir, "umap_coordinates.csv"), row.names = FALSE)
  write.csv(sample_info, file.path(traj_dir, "sample_metadata.csv"), row.names = FALSE)
  
}, error = function(e) cat(sprintf("  ⚠ Maturation trajectory failed: %s\n", e$message)))

# ==============================================================================
# SECTION 18: SUMMARY
# ==============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")
cat("  iPSC-CM DCM PHOSPHOPROTEOMICS — COMPLETE\n")
cat("════════════════════════════════════════════════════════════════════════════════\n\n")

cat("  Disease vs Healthy summary:\n")
cat(sprintf("  %-20s %3s %3s %6s %6s %6s\n", "Contrast", "nD", "nH", "FDR05", "Up", "Down"))
cat(paste(rep("-", 62), collapse = ""), "\n")
for (cn in names(disease_results)) {
  r <- disease_results[[cn]]
  cat(sprintf("  %-20s %3d %3d %6d %6d %6d\n",
              cn, r$nD, r$nH, r$n_sig, r$n_up, r$n_down))
}

cat(sprintf("\n  Protein normalisation: %s\n",
            ifelse(exists("PROTEIN_NORM_DONE") && PROTEIN_NORM_DONE, "YES", "NO")))

cat(sprintf("\n  Output: %s/\n\n", OUTPUT_BASE))
cat("  Per-contrast (11 folders):\n")
cat("    • top_table_all.csv + DE_FDR005/010/020/rawP001/005.csv\n")
cat("    • pca.png/.pdf\n")
cat("    • volcano.png/.pdf  (canonical 8×8)\n")
cat("    • volcano_compact.png/.pdf  (169×174 pt, theme_void, top 3)\n")
cat("    • fgsea_Hallmarks/KEGG/Reactome.csv + barplot\n")
cat("    • GO_BP/MF/CC_all/up/down.csv\n")
cat("    • kinase_activity/\n")
cat("      - ksea_barplot.png, ksea_volcano_original.png/.pdf\n")
cat("      - ksea_volcano_compact.png/.pdf  (theme_void, top 3 labels)\n")
cat("      - combined_kinase_results.csv\n")
cat("  Timepoint_Contrasts/:\n")
cat("    • D22_vs_D15, D29_vs_D15, D29_vs_D22 per genotype\n")
cat("  Interaction_Model/:\n")
cat("    • interaction_{gene}_{D22/D29}_vs_D15.csv\n")
cat("  Heatmaps/:\n")
cat("    • heatmap_global_{D15/D22/D29}.png/.pdf\n")
cat("    • heatmap_top50_{gene}.png/.pdf\n")
cat("  Scatterplots/:\n")
cat("    • Within-genotype logFC correlations\n")
cat("    • cor_heatmap_across_{D15/D22/D29}.png/.pdf\n")
cat("  Core_Signatures/:\n")
cat("    • core_DCM_iPSC_phospho_{D15/D22/D29}.csv\n")
cat("  Maturation_Trajectory/:\n")
cat("    • PCA + UMAP by timepoint, genotype, faceted (.png/.pdf)\n")
cat("    • pca_coordinates.csv, umap_coordinates.csv, sample_metadata.csv\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")
