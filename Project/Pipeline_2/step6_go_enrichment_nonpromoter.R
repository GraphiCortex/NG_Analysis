#!/usr/bin/env Rscript
# ==============================================
# Pipeline_2 — Step 6: GO–BP enrichment (non-promoter)
# Reads:  Project/New_Data_2/E2F4_vs_E2F3_significant_genes_nonpromoter.csv
# Writes: Project/New_Data_2/E2F4_vs_E2F3_GO_BP_*_nonpromoter.csv
# Keeps ONLY five families: cell_cycle, DNA_repair, apoptosis, autophagy, necrosis
# ==============================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("clusterProfiler","org.Mm.eg.db","GO.db","readr","dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)  # switch to org.Hs.eg.db if human
  library(GO.db)
  library(readr); library(dplyr)
})

in_genes <- "Project/New_Data_2/E2F4_vs_E2F3_significant_genes_nonpromoter.csv"
out_base <- "Project/New_Data_2"

# 1) Significant genes (SYMBOL)
sig_genes <- read_csv(in_genes, show_col_types = FALSE) %>%
  pull(GENE_SYMBOL) %>% unique() %>% na.omit()

if (length(sig_genes) < 5) stop("Too few non-promoter genes for GO enrichment (n < 5).")

# 2) Enrichment (BP; FDR < 0.05)
ego_bp <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Mm.eg.db,     # change to org.Hs.eg.db if needed
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_df <- as.data.frame(ego_bp)

# 3) Helper: parent + all children GO IDs
get_descendants <- function(goids) {
  goids <- as.character(goids)
  kids  <- unlist(lapply(goids, function(id) GOBPOFFSPRING[[id]]), use.names = FALSE)
  unique(c(goids, if (length(kids)) as.character(kids) else character(0)))
}

# 4) Category roots — FIVE families only
go_parents <- list(
  cell_cycle = "GO:0007049",
  DNA_repair = "GO:0006281",
  apoptosis  = "GO:0006915",
  autophagy  = "GO:0006914",
  necrosis   = c("GO:0070265","GO:0070266")  # necrotic/necroptotic processes
)

# 5) Build ID sets (roots + descendants)
go_sets <- lapply(go_parents, get_descendants)

# 6) Strict by-ID subsets
byID <- lapply(names(go_sets), function(nm) dplyr::filter(ego_df, ID %in% go_sets[[nm]]))
names(byID) <- names(go_sets)

# 7) (Optional) keyword fallbacks for quick browsing
kw_patterns <- list(
  cell_cycle = "\\bcell cycle\\b",
  DNA_repair = "DNA repair",
  apoptosis  = "apoptos",
  autophagy  = "autophag",
  necrosis   = "necros|necroptos"
)
byKW <- lapply(names(kw_patterns), function(nm)
  ego_df %>% filter(grepl(kw_patterns[[nm]], Description, ignore.case = TRUE)))
names(byKW) <- names(kw_patterns)

# 8) Write outputs (with _nonpromoter suffix)
write_csv(ego_df, file.path(out_base, "E2F4_vs_E2F3_GO_BP_full_nonpromoter.csv"))
for (nm in names(go_sets)) {
  write_csv(byID[[nm]], file.path(out_base, sprintf("E2F4_vs_E2F3_GO_BP_%s_byID_nonpromoter.csv", nm)))
  write_csv(byKW[[nm]], file.path(out_base, sprintf("E2F4_vs_E2F3_GO_BP_%s_nonpromoter.csv", nm)))
}

# 9) Console summary — five families only
n <- function(x) if (is.null(x) || (is.data.frame(x) && nrow(x)==0)) 0L else nrow(x)
cat("✔ Enrichment (non-promoter) done:\n",
    sprintf(" • total BP terms: %d\n", nrow(ego_df)),
    sprintf(" • Cell cycle   (byID): %d\n", n(byID$cell_cycle)),
    sprintf(" • DNA repair   (byID): %d\n", n(byID$DNA_repair)),
    sprintf(" • Apoptosis    (byID): %d\n", n(byID$apoptosis)),
    sprintf(" • Autophagy    (byID): %d\n", n(byID$autophagy)),
    sprintf(" • Necroptosis  (byID): %d\n", n(byID$necrosis)),
    "CSV files written under Project/New_Data_2/\n")
