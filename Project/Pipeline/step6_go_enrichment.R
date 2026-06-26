#!/usr/bin/env Rscript
# Step 6: GO–BP enrichment with GO-ID filtering (5 families)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("clusterProfiler","org.Mm.eg.db","GO.db","readr","dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)   # switch to org.Hs.eg.db if human
  library(GO.db)
  library(readr); library(dplyr)
})

# 1) Significant genes (from Step 5)
sig_genes <- read_csv("Project/New_Data/E2F4_vs_E2F3_significant_genes.csv",
                      show_col_types = FALSE) %>%
  pull(GENE_SYMBOL) %>% unique()

# 2) Enrichment (BP)
ego_bp <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_df <- as.data.frame(ego_bp)

# 3) Helper to expand GO roots to all descendants
get_descendants <- function(goids) {
  goids <- as.character(goids)
  kids  <- unlist(lapply(goids, function(id) GOBPOFFSPRING[[id]]), use.names = FALSE)
  unique(c(goids, if (length(kids)) as.character(kids) else character(0)))
}

# 4) EXACTLY the five requested families
go_parents <- list(
  cell_cycle   = "GO:0007049",
  DNA_repair   = "GO:0006281",
  apoptosis    = "GO:0006915",
  autophagy    = "GO:0006914",
  necroptosis  = c("GO:0070266", "GO:0070265")  # necroptotic process, necrotic cell death
)

# 5) Build ID sets (roots + descendants)
go_sets <- lapply(go_parents, get_descendants)

# 6) Subsets by ID + keyword fallbacks (for robustness)
byID <- lapply(names(go_sets), function(nm) dplyr::filter(ego_df, ID %in% go_sets[[nm]]))
names(byID) <- names(go_sets)

kw_patterns <- list(
  cell_cycle   = "\\bcell cycle\\b",
  DNA_repair   = "DNA repair",
  apoptosis    = "apoptos",
  autophagy    = "autophag",
  necroptosis  = "necroptos|necrotic cell death"
)
byKW <- lapply(names(kw_patterns), function(nm)
  ego_df %>% filter(grepl(kw_patterns[[nm]], Description, ignore.case = TRUE)))
names(byKW) <- names(kw_patterns)

# 7) Write outputs
out_base <- "Project/New_Data"
write_csv(ego_df, file.path(out_base, "E2F4_vs_E2F3_GO_BP_full.csv"))
for (nm in names(go_sets)) {
  write_csv(byID[[nm]], file.path(out_base, sprintf("E2F4_vs_E2F3_GO_BP_%s_byID.csv", nm)))
  write_csv(byKW[[nm]], file.path(out_base, sprintf("E2F4_vs_E2F3_GO_BP_%s.csv", nm)))
}

# 8) Console summary (only the five families)
n <- function(x) if (is.null(x) || (is.data.frame(x) && nrow(x)==0)) 0L else nrow(x)
cat("✔ Enrichment (promoter) done:\n",
    sprintf(" • total BP terms: %d\n", nrow(ego_df)),
    sprintf(" • Cell cycle       (byID): %d\n", n(byID$cell_cycle)),
    sprintf(" • DNA repair       (byID): %d\n", n(byID$DNA_repair)),
    sprintf(" • Apoptosis        (byID): %d\n", n(byID$apoptosis)),
    sprintf(" • Autophagy        (byID): %d\n", n(byID$autophagy)),
    sprintf(" • Necroptosis      (byID): %d\n", n(byID$necroptosis)),
    "CSV files written under Project/New_Data/\n")
