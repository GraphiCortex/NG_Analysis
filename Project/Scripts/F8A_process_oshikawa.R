#!/usr/bin/env Rscript

# ============================================================
# F8A — Full reprocessing of Oshikawa et al. 2013
# GSE37577 / GPL11202
#
# Goal:
#   Reprocess the 12 raw Agilent one-color arrays and estimate:
#
#   1. CAG-TKO  vs CAG-Control
#   2. MAP2-TKO vs MAP2-Control
#   3. Context interaction:
#
#      (MAP2_TKO - MAP2_Control)
#        -
#      (CAG_TKO - CAG_Control)
#
# The interaction contrast identifies genes whose response to
# Rb-family loss depends on developmental / differentiation
# context.
#
# Run from:
#   C:/Users/Seywe/Documents/GitHub/BIOL_295/Project
#
# Command:
#   Rscript .\Scripts\F8A_process_oshikawa.R
# ============================================================


suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
})


# ============================================================
# 0. Paths
# ============================================================

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

raw_dir <- file.path(
  project_root,
  "Data",
  "Raw",
  "microarray",
  "oshikawa2013",
  "Otxt"
)

out_dir <- file.path(
  project_root,
  "Data",
  "Processed",
  "Dynamics",
  "F8A_Oshikawa"
)

qc_dir <- file.path(out_dir, "QC")
de_dir <- file.path(out_dir, "DE")
annotation_dir <- file.path(out_dir, "Annotation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. Discover raw files
# ============================================================

files <- list.files(
  raw_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

files <- sort(files)

if (length(files) != 12) {
  stop(
    paste0(
      "Expected exactly 12 raw Oshikawa TXT files, found ",
      length(files),
      "."
    )
  )
}

file_names <- basename(files)

cat("\n")
cat("============================================================\n")
cat("F8A — OSHIKAWA RAW ARRAY REPROCESSING\n")
cat("============================================================\n")
cat("Raw arrays found:", length(files), "\n\n")

print(file_names)


# ============================================================
# 2. Build sample metadata from filenames
# ============================================================

get_group <- function(x) {

  if (grepl("_cCont", x)) {
    return("CAG_Control")
  }

  if (grepl("_cTKO", x)) {
    return("CAG_TKO")
  }

  if (grepl("_mCont", x)) {
    return("MAP2_Control")
  }

  if (grepl("_mTKO", x)) {
    return("MAP2_TKO")
  }

  return(NA_character_)
}


groups <- vapply(
  file_names,
  get_group,
  character(1)
)

if (any(is.na(groups))) {
  stop(
    "Could not determine group for: ",
    paste(file_names[is.na(groups)], collapse = ", ")
  )
}


sample_names <- sub("\\.txt$", "", file_names)

sample_info <- data.frame(
  Sample = sample_names,
  File = file_names,
  Group = groups,
  stringsAsFactors = FALSE
)

sample_info$Group <- factor(
  sample_info$Group,
  levels = c(
    "CAG_Control",
    "CAG_TKO",
    "MAP2_Control",
    "MAP2_TKO"
  )
)

write.csv(
  sample_info,
  file.path(out_dir, "sample_metadata.csv"),
  row.names = FALSE
)


cat("\nSample groups:\n")
print(table(sample_info$Group))


expected_counts <- c(
  CAG_Control = 3,
  CAG_TKO = 3,
  MAP2_Control = 3,
  MAP2_TKO = 3
)

observed_counts <- table(sample_info$Group)

if (!all(observed_counts[names(expected_counts)] == expected_counts)) {
  stop("Unexpected replicate counts.")
}


# ============================================================
# 3. Read raw one-color Agilent Feature Extraction files
# ============================================================
#
# source="agilent" uses:
#
#   foreground = gMedianSignal
#   background = gBGMedianSignal
#
# We additionally retain gIsWellAboveBG for expression filtering.
# ============================================================

cat("\nReading raw Agilent arrays...\n")

raw <- read.maimages(
  files = files,
  source = "agilent",
  green.only = TRUE,
  other.columns = "gIsWellAboveBG"
)

colnames(raw$E) <- sample_names
colnames(raw$Eb) <- sample_names

if (!is.null(raw$other$gIsWellAboveBG)) {
  colnames(raw$other$gIsWellAboveBG) <- sample_names
}

cat("Raw probe dimensions:", dim(raw$E), "\n")


# ============================================================
# 4. Basic raw QC
# ============================================================

raw_log2 <- log2(
  pmax(raw$E, 1)
)

png(
  file.path(qc_dir, "01_raw_signal_boxplot.png"),
  width = 1800,
  height = 1100,
  res = 160
)

par(mar = c(10, 5, 3, 1))

boxplot(
  raw_log2,
  las = 2,
  outline = FALSE,
  ylab = "log2 raw foreground intensity",
  main = "Oshikawa raw Agilent signals"
)

dev.off()


# ============================================================
# 5. Background correction
# ============================================================
#
# normexp is the standard limma workflow for Agilent arrays.
# ============================================================

bg <- backgroundCorrect(
  raw,
  method = "normexp",
  offset = 50
)

# ============================================================
# 6. Quantile normalization
# ============================================================

cat("Applying quantile normalization...\n")

norm <- normalizeBetweenArrays(
  bg,
  method = "quantile"
)

colnames(norm$E) <- sample_names


png(
  file.path(qc_dir, "02_normalized_signal_boxplot.png"),
  width = 1800,
  height = 1100,
  res = 160
)

par(mar = c(10, 5, 3, 1))

boxplot(
  norm$E,
  las = 2,
  outline = FALSE,
  ylab = "Normalized log2 expression",
  main = "Oshikawa quantile-normalized signals"
)

dev.off()


# ============================================================
# 7. Download / cache GPL11202 annotation
# ============================================================
#
# GPL11202 uses Agilent ProbeName as the platform ID.
# ============================================================

annotation_cache <- file.path(
  annotation_dir,
  "GPL11202_annotation.csv"
)

if (!file.exists(annotation_cache)) {

  cat("\nDownloading GPL11202 annotation from GEO...\n")

  gpl <- getGEO("GPL11202")

  annotation <- Table(gpl)

  write.csv(
    annotation,
    annotation_cache,
    row.names = FALSE
  )

} else {

  cat("\nUsing cached GPL11202 annotation...\n")

  annotation <- read.csv(
    annotation_cache,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


required_annotation_columns <- c(
  "ID",
  "GENE_SYMBOL",
  "GENE"
)

missing_ann <- setdiff(
  required_annotation_columns,
  colnames(annotation)
)

if (length(missing_ann) > 0) {

  stop(
    paste0(
      "GPL11202 annotation is missing expected columns: ",
      paste(missing_ann, collapse = ", ")
    )
  )
}


# ============================================================
# 8. Attach gene annotation
# ============================================================

annotation$ID <- trimws(
  as.character(annotation$ID)
)

annotation$GENE_SYMBOL <- trimws(
  as.character(annotation$GENE_SYMBOL)
)

annotation$GENE <- trimws(
  as.character(annotation$GENE)
)

probe_index <- match(
  norm$genes$ProbeName,
  annotation$ID
)

norm$genes$Symbol <- annotation$GENE_SYMBOL[
  probe_index
]

norm$genes$EntrezID <- annotation$GENE[
  probe_index
]


# Empty annotation fields -> NA

norm$genes$Symbol[
  norm$genes$Symbol %in% c("", "---", "NA", "null")
] <- NA

norm$genes$EntrezID[
  norm$genes$EntrezID %in% c("", "---", "NA", "null")
] <- NA


# Standardize gene symbols to uppercase because the existing
# BIOL_295 downstream pipeline uses uppercase symbols.

norm$genes$Symbol <- toupper(
  norm$genes$Symbol
)


# ============================================================
# 9. Probe filtering
# ============================================================
#
# Keep:
#
#   - non-control probes
#   - probes with gene symbol annotation
#   - probes above background in at least 3 arrays
#
# "3" corresponds to one complete biological replicate group.
# ============================================================

not_control <- (
  norm$genes$ControlType == 0
)

has_symbol <- (
  !is.na(norm$genes$Symbol)
  &
  nzchar(norm$genes$Symbol)
)

if (is.null(norm$other$gIsWellAboveBG)) {

  warning(
    "gIsWellAboveBG not available. ",
    "Skipping above-background expression filtering."
  )

  is_expressed <- rep(
    TRUE,
    nrow(norm$E)
  )

} else {

  is_expressed <- rowSums(
    norm$other$gIsWellAboveBG > 0,
    na.rm = TRUE
  ) >= 3
}


keep <- (
  not_control
  &
  has_symbol
  &
  is_expressed
)


cat("\nProbe filtering:\n")

cat(
  "  Total probes:       ",
  nrow(norm$E),
  "\n"
)

cat(
  "  Non-control:        ",
  sum(not_control, na.rm = TRUE),
  "\n"
)

cat(
  "  Annotated symbol:   ",
  sum(has_symbol, na.rm = TRUE),
  "\n"
)

cat(
  "  Above background:   ",
  sum(is_expressed, na.rm = TRUE),
  "\n"
)

cat(
  "  Final kept probes:  ",
  sum(keep, na.rm = TRUE),
  "\n"
)


filtered <- norm[
  keep,
]


# ============================================================
# 10. Remove ambiguous symbol mappings
# ============================================================
#
# If a probe is annotated to multiple symbols, it should not
# define a gene-level expression estimate.
# ============================================================

ambiguous <- grepl(
  "///|;|\\|",
  filtered$genes$Symbol
)

ambiguous[is.na(ambiguous)] <- FALSE

if (any(ambiguous)) {

  cat(
    "  Ambiguously annotated probes removed:",
    sum(ambiguous),
    "\n"
  )

  filtered <- filtered[
    !ambiguous,
  ]
}


# ============================================================
# 11. Collapse replicate probes to gene level
# ============================================================
#
# Multiple Agilent probes mapping to the same gene are averaged
# on the normalized log2 scale.
# ============================================================

cat("\nCollapsing probes to gene level...\n")

gene_expression <- avereps(
  filtered$E,
  ID = filtered$genes$Symbol
)

cat(
  "Gene-level expression matrix:",
  nrow(gene_expression),
  "genes x",
  ncol(gene_expression),
  "arrays\n"
)


gene_matrix_out <- data.frame(
  Gene = rownames(gene_expression),
  gene_expression,
  check.names = FALSE
)

write.csv(
  gene_matrix_out,
  file.path(
    out_dir,
    "oshikawa_gene_level_expression.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. Gene-level QC
# ============================================================

# ------------------------------------------------------------
# Correlation matrix
# ------------------------------------------------------------

cor_matrix <- cor(
  gene_expression,
  method = "pearson"
)

write.csv(
  cor_matrix,
  file.path(
    qc_dir,
    "sample_correlations.csv"
  )
)


png(
  file.path(
    qc_dir,
    "03_sample_correlation_heatmap.png"
  ),
  width = 1500,
  height = 1400,
  res = 160
)

heatmap(
  cor_matrix,
  symm = TRUE,
  margins = c(12, 12),
  main = "Sample correlation"
)

dev.off()


# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

pca <- prcomp(
  t(gene_expression),
  center = TRUE,
  scale. = FALSE
)

pca_var <- (
  pca$sdev^2 /
  sum(pca$sdev^2)
) * 100

pca_table <- data.frame(
  Sample = sample_names,
  Group = as.character(sample_info$Group),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

write.csv(
  pca_table,
  file.path(qc_dir, "PCA_coordinates.csv"),
  row.names = FALSE
)


group_num <- as.numeric(
  sample_info$Group
)

png(
  file.path(qc_dir, "04_PCA.png"),
  width = 1500,
  height = 1200,
  res = 160
)

plot(
  pca_table$PC1,
  pca_table$PC2,
  pch = group_num + 14,
  cex = 1.6,
  xlab = paste0(
    "PC1 (",
    round(pca_var[1], 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(pca_var[2], 1),
    "%)"
  ),
  main = "Oshikawa 2013 — PCA"
)

text(
  pca_table$PC1,
  pca_table$PC2,
  labels = seq_len(nrow(pca_table)),
  pos = 3,
  cex = 0.75
)

legend(
  "topright",
  legend = levels(sample_info$Group),
  pch = seq_along(levels(sample_info$Group)) + 14,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# MDS
# ------------------------------------------------------------

png(
  file.path(qc_dir, "05_MDS.png"),
  width = 1500,
  height = 1200,
  res = 160
)

plotMDS(
  gene_expression,
  labels = seq_len(ncol(gene_expression)),
  col = group_num,
  main = "Oshikawa 2013 — MDS"
)

legend(
  "topright",
  legend = levels(sample_info$Group),
  col = seq_along(levels(sample_info$Group)),
  pch = 16,
  bty = "n"
)

dev.off()


# ============================================================
# 13. Limma design matrix
# ============================================================

group <- sample_info$Group

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(group)

rownames(design) <- sample_names

write.csv(
  design,
  file.path(out_dir, "design_matrix.csv")
)


cat("\nDesign matrix:\n")
print(design)


# ============================================================
# 14. Fit gene-level linear models
# ============================================================

fit <- lmFit(
  gene_expression,
  design
)


# ============================================================
# 15. Define contrasts
# ============================================================
#
# Primary contrasts:
#
# CAG:
#   CAG_TKO - CAG_Control
#
# MAP2:
#   MAP2_TKO - MAP2_Control
#
# Interaction:
#
#   (MAP2_TKO - MAP2_Control)
#        -
#   (CAG_TKO - CAG_Control)
#
# Positive interaction logFC means the TKO effect is more
# positive in MAP2 than in CAG.
# ============================================================

contrast_matrix <- makeContrasts(

  CAG_TKO_vs_Control =
    CAG_TKO - CAG_Control,

  MAP2_TKO_vs_Control =
    MAP2_TKO - MAP2_Control,

  MAP2_vs_CAG_interaction =
    (MAP2_TKO - MAP2_Control)
    -
    (CAG_TKO - CAG_Control),

  levels = design
)


write.csv(
  contrast_matrix,
  file.path(out_dir, "contrast_matrix.csv")
)


cat("\nContrast matrix:\n")
print(contrast_matrix)


# ============================================================
# 16. Contrast fitting + empirical Bayes
# ============================================================

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(
  fit2,
  trend = TRUE,
  robust = TRUE
)


# ============================================================
# 17. Export helper
# ============================================================

export_contrast <- function(
  fit_object,
  coef_name,
  prefix,
  output_dir
) {

  tt <- topTable(
    fit_object,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )

  tt$Gene <- rownames(tt)

  tt <- tt[
    ,
    c(
      "Gene",
      setdiff(
        colnames(tt),
        "Gene"
      )
    )
  ]

  tt$Direction <- ifelse(
    tt$logFC > 0,
    "Up",
    ifelse(
      tt$logFC < 0,
      "Down",
      "Zero"
    )
  )

  tt$Significant_FDR05 <- (
    tt$adj.P.Val < 0.05
  )

  tt$Significant_LFC05 <- (
    tt$adj.P.Val < 0.05
    &
    abs(tt$logFC) >= 0.5
  )

  tt$Significant_LFC1 <- (
    tt$adj.P.Val < 0.05
    &
    abs(tt$logFC) >= 1
  )


  # --------------------------------------
  # Full results
  # --------------------------------------

  write.csv(
    tt,
    file.path(
      output_dir,
      paste0(prefix, "_full.csv")
    ),
    row.names = FALSE
  )


  # --------------------------------------
  # FDR-only results
  # --------------------------------------

  write.csv(
    tt[
      tt$Significant_FDR05,
    ],
    file.path(
      output_dir,
      paste0(prefix, "_FDR05.csv")
    ),
    row.names = FALSE
  )


  # --------------------------------------
  # |LFC| >= 0.5
  # --------------------------------------

  write.csv(
    tt[
      tt$Significant_LFC05,
    ],
    file.path(
      output_dir,
      paste0(prefix, "_FDR05_LFC05.csv")
    ),
    row.names = FALSE
  )


  # --------------------------------------
  # |LFC| >= 1
  # --------------------------------------

  write.csv(
    tt[
      tt$Significant_LFC1,
    ],
    file.path(
      output_dir,
      paste0(prefix, "_FDR05_LFC1.csv")
    ),
    row.names = FALSE
  )


  return(tt)
}


# ============================================================
# 18. Export all three contrasts
# ============================================================

cag <- export_contrast(
  fit2,
  "CAG_TKO_vs_Control",
  "CAG_TKO_vs_Control",
  de_dir
)

map2 <- export_contrast(
  fit2,
  "MAP2_TKO_vs_Control",
  "MAP2_TKO_vs_Control",
  de_dir
)

interaction <- export_contrast(
  fit2,
  "MAP2_vs_CAG_interaction",
  "MAP2_vs_CAG_interaction",
  de_dir
)


# ============================================================
# 19. Construct one combined Oshikawa gene-response table
# ============================================================

combined <- merge(
  cag[
    ,
    c(
      "Gene",
      "logFC",
      "P.Value",
      "adj.P.Val"
    )
  ],
  map2[
    ,
    c(
      "Gene",
      "logFC",
      "P.Value",
      "adj.P.Val"
    )
  ],
  by = "Gene",
  suffixes = c(
    "_CAG",
    "_MAP2"
  )
)


combined <- merge(
  combined,
  interaction[
    ,
    c(
      "Gene",
      "logFC",
      "P.Value",
      "adj.P.Val"
    )
  ],
  by = "Gene"
)


colnames(combined)[
  colnames(combined) == "logFC"
] <- "logFC_Interaction"

colnames(combined)[
  colnames(combined) == "P.Value"
] <- "P.Value_Interaction"

colnames(combined)[
  colnames(combined) == "adj.P.Val"
] <- "adj.P.Val_Interaction"


# ============================================================
# 20. Add qualitative context-response classes
# ============================================================

combined$CAG_direction <- ifelse(
  combined$logFC_CAG > 0,
  "Up",
  ifelse(
    combined$logFC_CAG < 0,
    "Down",
    "Zero"
  )
)

combined$MAP2_direction <- ifelse(
  combined$logFC_MAP2 > 0,
  "Up",
  ifelse(
    combined$logFC_MAP2 < 0,
    "Down",
    "Zero"
  )
)


combined$CAG_sig_LFC1 <- (
  combined$adj.P.Val_CAG < 0.05
  &
  abs(combined$logFC_CAG) >= 1
)

combined$MAP2_sig_LFC1 <- (
  combined$adj.P.Val_MAP2 < 0.05
  &
  abs(combined$logFC_MAP2) >= 1
)

combined$Interaction_sig <- (
  combined$adj.P.Val_Interaction < 0.05
)


combined$Context_class <- "No_strong_context_signal"


# Same-direction strict DE in both contexts

same_up <- (
  combined$CAG_sig_LFC1
  &
  combined$MAP2_sig_LFC1
  &
  combined$logFC_CAG > 0
  &
  combined$logFC_MAP2 > 0
)

same_down <- (
  combined$CAG_sig_LFC1
  &
  combined$MAP2_sig_LFC1
  &
  combined$logFC_CAG < 0
  &
  combined$logFC_MAP2 < 0
)

combined$Context_class[
  same_up
] <- "Conserved_up"

combined$Context_class[
  same_down
] <- "Conserved_down"


# Opposite signs and both strict-significant

reversed <- (
  combined$CAG_sig_LFC1
  &
  combined$MAP2_sig_LFC1
  &
  sign(combined$logFC_CAG)
  !=
  sign(combined$logFC_MAP2)
)

combined$Context_class[
  reversed
] <- "Direction_reversal"


# Strict in only one context

cag_only <- (
  combined$CAG_sig_LFC1
  &
  !combined$MAP2_sig_LFC1
)

map2_only <- (
  !combined$CAG_sig_LFC1
  &
  combined$MAP2_sig_LFC1
)

combined$Context_class[
  cag_only
] <- "CAG_specific"

combined$Context_class[
  map2_only
] <- "MAP2_specific"


# Significant formal interaction overrides less-specific labels,
# except true sign reversals which remain their own class.

interaction_only <- (
  combined$Interaction_sig
  &
  !reversed
)

combined$Context_class[
  interaction_only
] <- "Context_interaction"


write.csv(
  combined,
  file.path(
    out_dir,
    "oshikawa_gene_response_matrix.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. Export formal context-dependent genes
# ============================================================

context_dependent <- combined[
  combined$adj.P.Val_Interaction < 0.05,
]

context_dependent <- context_dependent[
  order(
    context_dependent$adj.P.Val_Interaction
  ),
]

write.csv(
  context_dependent,
  file.path(
    de_dir,
    "context_dependent_genes_FDR05.csv"
  ),
  row.names = FALSE
)


context_reversed <- combined[
  combined$Context_class == "Direction_reversal",
]

write.csv(
  context_reversed,
  file.path(
    de_dir,
    "context_direction_reversals.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. Differential-expression counts
# ============================================================

count_direction <- function(
  table,
  threshold = 1
) {

  sig <- (
    table$adj.P.Val < 0.05
    &
    abs(table$logFC) >= threshold
  )

  c(
    Total = sum(sig, na.rm = TRUE),
    Up = sum(
      sig & table$logFC > 0,
      na.rm = TRUE
    ),
    Down = sum(
      sig & table$logFC < 0,
      na.rm = TRUE
    )
  )
}


cag_counts <- count_direction(
  cag,
  threshold = 1
)

map2_counts <- count_direction(
  map2,
  threshold = 1
)

interaction_counts_fdr <- c(
  Total = sum(
    interaction$adj.P.Val < 0.05,
    na.rm = TRUE
  ),
  Positive = sum(
    interaction$adj.P.Val < 0.05
    &
    interaction$logFC > 0,
    na.rm = TRUE
  ),
  Negative = sum(
    interaction$adj.P.Val < 0.05
    &
    interaction$logFC < 0,
    na.rm = TRUE
  )
)


# ============================================================
# 23. Cross-context correlation
# ============================================================

cor_cag_map2 <- cor(
  combined$logFC_CAG,
  combined$logFC_MAP2,
  method = "pearson",
  use = "complete.obs"
)

spearman_cag_map2 <- cor(
  combined$logFC_CAG,
  combined$logFC_MAP2,
  method = "spearman",
  use = "complete.obs"
)


# ============================================================
# 24. CAG-vs-MAP2 fold-change scatter plot
# ============================================================

png(
  file.path(
    qc_dir,
    "06_CAG_vs_MAP2_logFC.png"
  ),
  width = 1400,
  height = 1400,
  res = 160
)

plot(
  combined$logFC_CAG,
  combined$logFC_MAP2,
  pch = 16,
  cex = 0.45,
  xlab = "CAG TKO vs Control log2FC",
  ylab = "MAP2 TKO vs Control log2FC",
  main = paste0(
    "CAG vs MAP2 Rb-TKO response\n",
    "Pearson r = ",
    round(cor_cag_map2, 3)
  )
)

abline(
  h = 0,
  v = 0,
  lty = 3
)

abline(
  0,
  1,
  lty = 2
)

dev.off()


# ============================================================
# 25. Summary file
# ============================================================

summary_file <- file.path(
  out_dir,
  "F8A_summary.txt"
)

sink(summary_file)

cat("F8A — OSHIKAWA 2013 FULL RAW REANALYSIS\n")
cat("========================================\n\n")

cat("Dataset\n")
cat("-------\n")
cat("Raw arrays:", length(files), "\n")
cat("Gene-level features:", nrow(gene_expression), "\n\n")


cat("Sample groups\n")
cat("-------------\n")
print(table(sample_info$Group))
cat("\n")


cat("Probe filtering\n")
cat("---------------\n")
cat("Total probes:", nrow(norm$E), "\n")
cat(
  "Non-control probes:",
  sum(not_control, na.rm = TRUE),
  "\n"
)
cat(
  "Annotated probes:",
  sum(has_symbol, na.rm = TRUE),
  "\n"
)
cat(
  "Above-background probes:",
  sum(is_expressed, na.rm = TRUE),
  "\n"
)
cat(
  "Final kept probes:",
  nrow(filtered$E),
  "\n\n"
)


cat("CAG TKO vs Control — FDR < 0.05 & |log2FC| >= 1\n")
cat("------------------------------------------------\n")
print(cag_counts)
cat("\n")


cat("MAP2 TKO vs Control — FDR < 0.05 & |log2FC| >= 1\n")
cat("-------------------------------------------------\n")
print(map2_counts)
cat("\n")


cat("MAP2-vs-CAG interaction — FDR < 0.05\n")
cat("-------------------------------------\n")
print(interaction_counts_fdr)
cat("\n")


cat("Cross-context similarity\n")
cat("------------------------\n")
cat(
  "Pearson correlation of CAG and MAP2 log2FC:",
  round(cor_cag_map2, 4),
  "\n"
)

cat(
  "Spearman correlation of CAG and MAP2 log2FC:",
  round(spearman_cag_map2, 4),
  "\n\n"
)


cat("Context classes\n")
cat("---------------\n")
print(
  sort(
    table(combined$Context_class),
    decreasing = TRUE
  )
)
cat("\n")


cat("Direction reversals:", nrow(context_reversed), "\n")
cat(
  "Formal interaction genes (FDR < 0.05):",
  nrow(context_dependent),
  "\n"
)


sink()


# ============================================================
# 26. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("F8A COMPLETE\n")
cat("============================================================\n\n")

cat("Gene-level features:", nrow(gene_expression), "\n\n")

cat(
  "CAG strict DEGs:",
  cag_counts["Total"],
  "(",
  cag_counts["Up"],
  "up /",
  cag_counts["Down"],
  "down )\n"
)

cat(
  "MAP2 strict DEGs:",
  map2_counts["Total"],
  "(",
  map2_counts["Up"],
  "up /",
  map2_counts["Down"],
  "down )\n"
)

cat(
  "Formal context-interaction genes:",
  nrow(context_dependent),
  "\n"
)

cat(
  "Strict direction reversals:",
  nrow(context_reversed),
  "\n"
)

cat(
  "CAG/MAP2 Pearson logFC correlation:",
  round(cor_cag_map2, 3),
  "\n"
)

cat(
  "CAG/MAP2 Spearman logFC correlation:",
  round(spearman_cag_map2, 3),
  "\n\n"
)


cat("Context classes:\n")
print(
  sort(
    table(combined$Context_class),
    decreasing = TRUE
  )
)


cat("\nOutputs written to:\n")
cat(out_dir, "\n\n")

cat(
  "[NEXT] Inspect QC, compare the new CAG/MAP2 contrasts ",
  "against the old Oshikawa L2FC1 table, then build the ",
  "Fong × E2F1/3 × CAG × MAP2 context matrix.\n"
)
