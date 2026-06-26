#!/usr/bin/env Rscript

# parse_AnnotGPL_promoter.R
# Read the manually downloaded AnnotGPL, extract promoter-only probes → gene symbols.

if (!requireNamespace("GEOquery", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)

# 1) Load the local AnnotGPL file
annot_file <- "Project/MetaData/GPL18280_family.soft"
message("Parsing AnnotGPL from ", annot_file)
gpl <- getGEO(filename = annot_file)
tbl <- Table(gpl)

# 2) Confirm column names
cols <- colnames(tbl)
message("Columns in AnnotGPL: ", paste(head(cols,20), collapse=", "))
# You should see ID_REF, Gene Symbol (or GENE_SYMBOL), DESCRIPTION, etc.

# 3) Subset to PROMOTER probes with non‐NA gene symbols
sym_col <- if("Gene Symbol" %in% cols) "Gene Symbol" else 
           if("GENE_SYMBOL" %in% cols) "GENE_SYMBOL" else
           stop("No gene‐symbol column in AnnotGPL: ", paste(cols, collapse=", "))
stopifnot("ID_REF" %in% cols, "DESCRIPTION" %in% cols)

prom <- subset(tbl,
               DESCRIPTION == "PROMOTER" &
               !is.na(tbl[[sym_col]]) &
               tbl[[sym_col]] != ""
)

message("→ Found ", nrow(prom), " promoter probes in AnnotGPL")

# 4) Write out the filtered probe→gene map
probe2gene <- data.frame(
  ID_REF = prom$ID_REF,
  SYMBOL = prom[[sym_col]],
  stringsAsFactors = FALSE
)
outf <- "Project/Data/MetaData/GPL_promoter_probe2gene.csv"
write.csv(probe2gene, outf, row.names=FALSE, quote=FALSE)
message("Wrote ", nrow(probe2gene), " mappings to ", outf)
