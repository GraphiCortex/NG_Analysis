#!/usr/bin/env Rscript

# ==============================================
# Pipeline_2 — Step 3: Differential binding with limma
# (non-promoter: DOWNSTREAM + INSIDE)
# ==============================================

# 0) Packages
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(stringr)
})

# 1) Load M-value matrix (Pipeline_2)
Mmat <- read_csv("Project/New_Data_2/Mvalue_matrix.csv", col_types = cols())

# 2) Subset to non-promoter probes (DOWNSTREAM + INSIDE)
map_np <- read_csv("Project/New_Data_2/GPL_nonpromoter_probe2gene.csv", col_types = cols())
Mmat_np <- Mmat %>% filter(SPOT_ID %in% map_np$SPOT_ID)

if (nrow(Mmat_np) == 0) stop("No non-promoter probes found in M matrix after filtering.")

# 3) Build expression matrix (rows = SPOT_ID, cols = samples)
rownames(Mmat_np) <- Mmat_np$SPOT_ID
expr <- as.matrix(Mmat_np[, -1, drop = FALSE])

# 4) Define sample-group factor from column names
samples <- colnames(expr)
groups <- sapply(samples, function(x) {
  parts <- strsplit(x, "_")[[1]]
  paste(parts[2], parts[3], sep = "_")
})
groups <- factor(
  groups,
  levels = c("E2f4_wt", "E2f3_wt", "E2f3_3amut", "E2f3_3bmut")
)

# 5) Design and contrast
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
cont.matrix <- makeContrasts(E2F4_vs_E2F3 = E2f4_wt - E2f3_wt, levels = design)

# 6) Fit model + moderated stats
fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, cont.matrix)
eb   <- eBayes(fit2)

# 7) Save full probe-level results
out_csv <- "Project/New_Data_2/E2F4_vs_E2F3_fullLimma.csv"
res_all <- topTable(eb, coef = "E2F4_vs_E2F3", number = Inf, sort.by = "none")
write_csv(res_all, out_csv)

message("✔ Step 3 (non-promoter) complete: wrote ",
        nrow(res_all), " probe-level results to ", out_csv)
