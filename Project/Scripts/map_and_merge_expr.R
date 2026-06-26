#!/usr/bin/env Rscript
# map_and_merge_expr.R
# 1) Map Ensembl → SYMBOL (strip version, use AnnotationDbi),
# 2) Read in your three gene-set CSVs,
# 3) Merge, flag, and write out expression-integrated tables.

# — install / load what we need
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("AnnotationDbi","org.Mm.eg.db","readr","dplyr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg)
}
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(readr)
  library(dplyr)
})

# 1. Load the raw Fong2022 microarray L2FC1 CSV
expr_raw <- read_csv(
  file.path("Project","Data","Raw","microarray","Fong2022","Fong 2022 L2FC1 microarray data.csv"),
  col_types = cols(.default="c")
)

# 2. Rename, strip Ensembl version, select only what we need
expr2 <- expr_raw %>%
  rename(
    ENSEMBL = Row.names,
    LFC      = log2FoldChange,
    adjP     = padj
  ) %>%
  mutate(
    # remove the “.###” version suffix
    ENSEMBL = sub("\\.\\d+$","",ENSEMBL),
    LFC      = as.numeric(LFC),
    adjP     = as.numeric(adjP)
  ) %>%
  select(ENSEMBL, LFC, adjP) %>%
  filter(!is.na(ENSEMBL))

# 3. Map each ENSEMBL to a SYMBOL
expr2$Gene <- mapIds(
  org.Mm.eg.db,
  keys      = expr2$ENSEMBL,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# drop anything that didn’t map
expr2 <- expr2 %>% filter(!is.na(Gene)) %>% select(Gene, LFC, adjP)

# 4. A helper to merge & flag significance
merge_flag <- function(gset_file, category) {
  gset <- read_csv(gset_file, col_types=cols(Gene="c"))
  df <- left_join(gset, expr2, by="Gene") %>%
    # keep every gene that actually has LFC & adjP
    filter(!is.na(LFC), !is.na(adjP)) %>%
    # now compute sig for every row
    mutate(
      sig      = (abs(LFC) >= 0.5) & (adjP < 0.10),
      category = category
    )
  return(df)
}




# 5. Paths to your three E2F4 gene sets
base           <- file.path("Project","Data","Processed","gene_sets")
mitotic_file   <- file.path(base,"E2F4_mitotic_genes.csv")
repair_file    <- file.path(base,"E2F4_DNArepair_genes.csv")
apoptosis_file <- file.path(base,"E2F4_apoptosis_genes.csv")

# 6. Merge & flag each
mitotic_df   <- merge_flag(mitotic_file,   "Cell-cycle")
repair_df    <- merge_flag(repair_file,    "DNA-repair")
apoptosis_df <- merge_flag(apoptosis_file, "Apoptosis")

# 7. Write out
outdir <- file.path("Project","Data","Processed","expression_integration")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

write_csv(mitotic_df,   file.path(outdir,"mitotic_expression.csv"))
write_csv(repair_df,    file.path(outdir,"DNArepair_expression.csv"))
write_csv(apoptosis_df, file.path(outdir,"apoptosis_expression.csv"))
write_csv(bind_rows(mitotic_df,repair_df,apoptosis_df),
          file.path(outdir,"all_three_categories.csv"))

cat("✅ Expression integration complete.\n",
    sprintf(" • Cell-cycle:   %d genes → %d significant\n",
            nrow(mitotic_df), sum(mitotic_df$sig, na.rm=TRUE)),
    sprintf(" • DNA-repair:   %d genes → %d significant\n",
            nrow(repair_df), sum(repair_df$sig, na.rm=TRUE)),
    sprintf(" • Apoptosis:     %d genes → %d significant\n",
            nrow(apoptosis_df), sum(apoptosis_df$sig, na.rm=TRUE)),
    sep = ""
)
