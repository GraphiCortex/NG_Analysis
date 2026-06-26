#!/usr/bin/env Rscript
#
# run_all_ORA.R 
#   Perform GO–BP over-representation (ORA) for all three E2F4 gene-sets,
#   filter to the relevant terms (cell cycle / DNA repair / cell death),
#   and write out CSVs.

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
for (pkg in c("clusterProfiler","org.Mm.eg.db","readr","dplyr")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    BiocManager::install(pkg)
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(readr)
  library(dplyr)
})

# 1) Define your gene‐set files and the regex to pull out only the GO terms you care about
sets <- list(
  Cell_cycle   = list(file="E2F4_mitotic_genes.csv",
                      pattern="cell cycle|mitotic"),
  DNA_repair   = list(file="E2F4_DNArepair_genes.csv",
                      pattern="DNA repair|DNA recombination|DNA metabolic"),
  Cell_death   = list(file="E2F4_apoptosis_genes.csv",
                      pattern="apoptosis|cell death")
)

in_dir  <- file.path("Project","Data","Processed","gene_sets")
out_dir <- file.path("Project","Data","Processed","limma_results")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

for (nm in names(sets)) {
  # read the gene symbols
  gfile   <- file.path(in_dir, sets[[nm]]$file)
  genes   <- read_csv(gfile, col_types=cols(Gene="c"))$Gene

  # map SYMBOL → ENTREZID
  entrez  <- mapIds(org.Mm.eg.db,
                    keys      = genes,
                    column    = "ENTREZID",
                    keytype   = "SYMBOL",
                    multiVals = "first")  %>% na.omit()

  # run enrichGO (ORA)
  ego <- enrichGO(gene          = entrez,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

  # turn into a tibble, filter to only the terms matching our regex
  df <- as_tibble(ego) %>%
    filter(grepl(sets[[nm]]$pattern, Description, ignore.case=TRUE))

  # write CSV
  out_csv <- file.path(out_dir, paste0("E2F4_", nm, "_ORA.csv"))
  write_csv(df, out_csv)
  message("→ Wrote ", nrow(df), " ORA rows to ", out_csv)
}
