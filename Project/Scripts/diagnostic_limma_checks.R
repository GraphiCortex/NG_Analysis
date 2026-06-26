#!/usr/bin/env Rscript

# diagnostic_limma_checks.R
# Comprehensive diagnostics for your limma pipeline – no external Rdata required.

# 0. Load required package
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
library(limma)

# 1. Discover processed M-value files
files <- list.files(
  path = "Project/Data/Raw/ChipSeq/Julian2016/Processed/",
  pattern = "_Processed\\.txt$",
  full.names = TRUE
)
if (length(files) < 6) stop("Expected at least 6 processed files, found ", length(files))
cat("Found", length(files), "processed files.\n")

# 2. Inspect first file structure
cat("\n[1] Inspecting first file columns & head:\n")
d <- read.delim(files[1], stringsAsFactors = FALSE, nrows = 6)
print(colnames(d))
print(d)

# 3. Build M-value matrix from all files
cat("\n[2] Building M-value matrix...\n")
Mlist <- lapply(files, function(f) {
  df <- read.delim(f, stringsAsFactors = FALSE)
  if ("M" %in% names(df)) {
    vals <- df$M
  } else if (all(c("IP", "INPUT") %in% names(df))) {
    vals <- log2(df$IP / df$INPUT)
  } else {
    stop("File ", f, " missing M or IP/INPUT columns")
  }
  setNames(vals, df$ID_REF)
})
Mvals <- do.call(cbind, Mlist)
colnames(Mvals) <- sub("_Processed\\.txt$", "", basename(files))
cat("Dimensions (probes × samples):", dim(Mvals), "\n")
cat("Any missing values?", any(is.na(Mvals)), "\n")

# 4. Infer sample groups from filenames
cat("\n[3] Inferring sample groups from filenames...\n")
groups <- sapply(colnames(Mvals), function(nm) {
  if (grepl("E2f4", nm, ignore.case=TRUE)) return("E2F4")
  if (grepl("E2f3_wt", nm, ignore.case=TRUE)) return("E2F3_wt")
  if (grepl("E2f3_3amut", nm, ignore.case=TRUE)) return("E2F3_3amut")
  if (grepl("E2f3_3bmut", nm, ignore.case=TRUE)) return("E2F3_3bmut")
  stop("Cannot parse group from sample name: ", nm)
})
group <- factor(groups)
print(table(group))

# 5. Build design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cat("Design columns:", paste(colnames(design), collapse=", "), "\n")

# 6. Fit limma model and contrasts
cat("\n[4] Fitting model and computing contrasts...\n")
fit <- lmFit(Mvals, design)
cont <- makeContrasts(
  E2F4_vs_E2F3   = E2F4 - E2F3_wt,
  levels = design
)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

# 7. Summarize unfiltered statistics for E2F4_vs_E2F3
cat("\n[5] Unfiltered stats for E2F4_vs_E2F3:\n")
tab_all <- topTable(fit2, coef="E2F4_vs_E2F3", number=Inf, adjust.method="none")
cat("Total probes in topTable:", nrow(tab_all), "\n")
cat("logFC summary:\n"); print(summary(tab_all$logFC))
cat("P.Value summary:\n"); print(summary(tab_all$P.Value))
cat("Raw P<0.05 count:", sum(tab_all$P.Value < .05), "\n")

# 8. Summarize filtered results
cat("\n[6] Applying significance thresholds (FDR<0.05 & |logFC|>=0.5):\n")
sig_fdr <- topTable(
  fit2,
  coef          = "E2F4_vs_E2F3",
  number        = Inf,
  adjust.method = "fdr",
  p.value       = 0.05,
  lfc           = 0.5
)
cat("Significant probes count:", nrow(sig_fdr), "\n")

# End of diagnostics
cat("\nDiagnostics complete.\n")
