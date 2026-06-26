#!/usr/bin/env Rscript

# =========================================================
# Pipeline_2 — Step 4: TREAT-based significance (non-promoter)
# Change from original: TREAT lfc = 0.25 (was 0.5)
# Significant probes: raw P.Value < 0.05 (no FDR here)
# Outputs:
#   - Project/New_Data_2/E2F4_vs_E2F3_treat_full_nonpromoter.csv
#   - Project/New_Data_2/E2F4_vs_E2F3_significant_probes_nonpromoter.csv
# =========================================================

# 0) Packages
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(tibble)   # rownames_to_column
  library(stringr)
})

# 1) Load non-promoter M-matrix and mapping
Mmat   <- read_csv("Project/New_Data_2/Mvalue_matrix.csv", col_types = cols())
np_map <- read_csv("Project/New_Data_2/GPL_nonpromoter_probe2gene.csv", col_types = cols())

# keep only non-promoter probes
Mmat_np <- Mmat %>% filter(SPOT_ID %in% np_map$SPOT_ID)
stopifnot(nrow(Mmat_np) > 0)

# 2) Build matrix for limma
expr <- as.matrix(Mmat_np %>% dplyr::select(-SPOT_ID))  # avoid 'select' S4 clash
rownames(expr) <- Mmat_np$SPOT_ID

# 3) Design/contrast (same parsing as original)
samples <- colnames(expr)
groups <- sapply(samples, function(x) {
  parts <- strsplit(x, "_")[[1]]
  paste(parts[2], parts[3], sep = "_")
})
groups <- factor(groups, levels = c("E2f4_wt", "E2f3_wt", "E2f3_3amut", "E2f3_3bmut"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
cont.matrix <- makeContrasts(E2F4_vs_E2F3 = E2f4_wt - E2f3_wt, levels = design)

# 4) Fit + TREAT with lfc = 0.25
fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fitT <- treat(fit2, lfc = 0.25)

# 5) Full TREAT table & significant probes (raw P < 0.05)
res_treat <- topTreat(fitT, coef = "E2F4_vs_E2F3", number = Inf, sort.by = "none") %>%
  rownames_to_column(var = "SPOT_ID")

sig_probes <- res_treat %>% filter(P.Value < 0.05)

# 6) Write outputs (Pipeline_2)
out_dir <- "Project/New_Data_2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(res_treat, file.path(out_dir, "E2F4_vs_E2F3_treat_full_nonpromoter.csv"))
write_csv(sig_probes, file.path(out_dir, "E2F4_vs_E2F3_significant_probes_nonpromoter.csv"))

message("✔ Step 4 (non-promoter) complete at LFC=0.25: ",
        nrow(sig_probes),
        " significant probes (raw P < 0.05) written to New_Data_2/")
