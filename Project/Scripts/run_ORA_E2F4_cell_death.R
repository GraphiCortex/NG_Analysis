#!/usr/bin/env Rscript
# run_ORA_E2F4_cell_death.R
# Over‐representation test for cell‐death/​apoptosis in all E2F4-biased targets

# 1. Load libs (install if needed)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for(pkg in c("clusterProfiler","org.Mm.eg.db","dplyr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg)
}
suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
})

# 2. Read limma results
limma_file <- file.path("Project","Data","Processed","limma_results",
                        "limma_E2F4_vs_E2F3.csv")
res <- read.csv(limma_file, stringsAsFactors=FALSE, row.names=1)
res$ID_REF <- rownames(res)

# 3. Read probe→gene map
map_file <- file.path("Project","MetaData","GPL_probe2gene.csv")
annot   <- read.csv(map_file, stringsAsFactors=FALSE)

# 4. Merge to get SYMBOL for each probe
df <- res %>%
  inner_join(annot, by="ID_REF")       # now df has columns: ID_REF, logFC, SYMBOL, …

# 5. Define your gene set: those with stronger E2F4 binding (logFC > 0)
bound4 <- df %>% filter(logFC > 0) %>% pull(SYMBOL) %>% unique()

# sanity check
cat("Number of E2F4-biased genes:", length(bound4), "\n")

# 6. Run GO-BP over-representation
ora <- enrichGO(
  gene          = bound4,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ora_df <- as.data.frame(ora)

# 7. Filter for cell-death/apoptosis terms
cell_death_terms <- ora_df %>%
  filter(grepl("cell death|apoptosis", Description, ignore.case=TRUE))

# 8. Save results
outdir <- file.path("Project","Data","Processed","limma_results")
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
write.csv(
  cell_death_terms,
  file = file.path(outdir, "E2F4_cell_death_ORA.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("→ Wrote", nrow(cell_death_terms),
    "cell-death/apoptosis terms to", outdir, "\n")
