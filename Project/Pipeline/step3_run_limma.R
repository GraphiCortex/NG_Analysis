#!/usr/bin/env Rscript

# ==============================================
# Step 3: Differential binding with limma
# ==============================================

# 0. Ensure packages are present
if (!requireNamespace("limma", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
}
library(limma)
library(tidyverse)

# 1. Load M-value matrix
Mmat <- read_csv("Project/New_Data/Mvalue_matrix.csv", col_types = cols())

# 2. Subset to promoter probes
prom_map <- read_csv("Project/New_Data/GPL_promoter_probe2gene.csv", col_types = cols())
Mmat_prom <- Mmat %>% 
  filter(SPOT_ID %in% prom_map$SPOT_ID)

# 3. Build expr matrix (rows=SPOT_ID, cols=samples)
rownames(Mmat_prom) <- Mmat_prom$SPOT_ID
expr <- as.matrix(Mmat_prom[ , -1, drop=FALSE])

# 4. Define sample‐group factor
samples <- colnames(expr)
groups <- sapply(samples, function(x) {
  parts <- strsplit(x, "_")[[1]]
  paste(parts[2], parts[3], sep="_")
})
# ensure the levels order matches your experimental layout:
groups <- factor(groups, 
                 levels = c("E2f4_wt",
                            "E2f3_wt",
                            "E2f3_3amut",
                            "E2f3_3bmut"))

# 5. Create design and contrast matrices
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
cont.matrix <- makeContrasts(
  E2F4_vs_E2F3 = E2f4_wt - E2f3_wt,
  levels = design
)

# 6. Fit the linear model and compute moderated statistics
fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, cont.matrix)
eb   <- eBayes(fit2)

# 7. Save the full table of results (all probes)
res_all <- topTable(eb,
                    coef = "E2F4_vs_E2F3",
                    number = Inf,
                    sort.by = "none")
write_csv(res_all, "Project/New_Data/E2F4_vs_E2F3_fullLimma.csv")

message("✔ Step 3 complete: wrote ",
        nrow(res_all),
        " probe‐level results to Pipeline/E2F4_vs_E2F3_fullLimma.csv")
