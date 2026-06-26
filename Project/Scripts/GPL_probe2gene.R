#!/usr/bin/env Rscript

# download_GPL_annotation.R
# Pull down Agilent GPL from GEO and extract probe → gene mapping

# 1) Install / load GEOquery and BiocManager if needed
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("GEOquery", quietly=TRUE)) {
    BiocManager::install("GEOquery")
}
library(GEOquery)

# 2) Fetch one of your GSMs to learn the platform ID
gsm <- getGEO("GSM1326858", GSEMatrix=FALSE)
plat_id <- Meta(gsm)$platform         # e.g. "GPL13534"

message("Detected platform: ", plat_id)

# 3) Download the platform (annotation) itself
gpl <- getGEO(plat_id, AnnotGPL=TRUE)

# 4) Pull out the table and select only the columns you need
tbl <- Table(gpl)
# Inspect names(tbl) if “Gene Symbol” is different
if (!"Gene Symbol" %in% colnames(tbl)) {
    stop("Couldn't find a column named 'Gene Symbol' in the GPL table. Columns are:\n",
         paste(colnames(tbl), collapse=", "))
}

probe2gene <- data.frame(
    ID_REF = tbl$ID_REF,
    SYMBOL = tbl$`Gene Symbol`,
    stringsAsFactors = FALSE
)

# 5) Write to your MetaData folder
outdir <- file.path("Project","Data","MetaData")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
outfile <- file.path(outdir, "GPL_probe2gene.csv")

write.csv(probe2gene, file=outfile, row.names=FALSE, quote=FALSE)
message("Wrote ", nrow(probe2gene), " mappings to ", outfile)
