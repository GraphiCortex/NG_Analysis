#!/usr/bin/env Rscript
# ==============================================
# Step 4 (PROMOTER): TREAT-based significance with adjustable LFC
# Default LFC = 0.25 (log2 units); override via CLI: Rscript step4_call_significance.R 0.20
# ==============================================

if (!requireNamespace("limma", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}
suppressMessages({ library(limma); library(readr); library(dplyr); library(tibble) })

# ---- Adjustable LFC ----
args <- commandArgs(trailingOnly = TRUE)
LFC  <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
if (is.na(LFC)) LFC <- 0.25

# ---- Inputs ----
in_dir  <- "Project/New_Data"
Mmat    <- read_csv(file.path(in_dir, "Mvalue_matrix.csv"), col_types = cols())
prommap <- read_csv(file.path(in_dir, "GPL_promoter_probe2gene.csv"), col_types = cols())

# ---- Promoter subset & expression matrix ----
Mmat_prom <- Mmat %>% filter(SPOT_ID %in% prommap$SPOT_ID)
expr <- as.matrix(Mmat_prom %>% select(-SPOT_ID))
rownames(expr) <- Mmat_prom$SPOT_ID

# ---- Design / contrasts (unchanged) ----
samples <- colnames(expr)
groups <- sapply(samples, function(x) {
  p <- strsplit(x, "_")[[1]]
  paste(p[2], p[3], sep = "_")
})
groups <- factor(groups, levels = c("E2f4_wt","E2f3_wt","E2f3_3amut","E2f3_3bmut"))
design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)
ct     <- makeContrasts(E2F4_vs_E2F3 = E2f4_wt - E2f3_wt, levels = design)

# ---- Fit + TREAT at relaxed LFC ----
fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, ct)
fitT <- treat(fit2, lfc = LFC)

res_treat <- topTreat(fitT, coef = "E2F4_vs_E2F3", number = Inf, sort.by = "none") %>%
  rownames_to_column("SPOT_ID")
sig_probes <- res_treat %>% filter(P.Value < 0.05)

# ---- Outputs (same filenames for downstream) ----
write_csv(res_treat, file.path(in_dir, "E2F4_vs_E2F3_treat_full.csv"))
write_csv(sig_probes, file.path(in_dir, "E2F4_vs_E2F3_significant_probes.csv"))

message(sprintf("✔ Step 4 (promoter) at LFC=%.2f → %d significant probes (P<0.05).", 
                LFC, nrow(sig_probes)))
