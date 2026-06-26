#!/usr/bin/env Rscript
# filter_promoter_probe2gene_from_GPLquery.R
library(GEOquery)

# 1) Load the GEOquery GPL you just downloaded
gsm   <- getGEO("GSM1326858", GSEMatrix=FALSE)
plat  <- Meta(gsm)$platform
gpl   <- getGEO(plat, AnnotGPL=TRUE)
tbl   <- Table(gpl)

# 2) Subset to promoter probes and non‐NA gene symbols
promoters <- subset(tbl, DESCRIPTION=="PROMOTER" & !is.na(`Gene Symbol`))
probe2gene <- data.frame(
  ID_REF = promoters$ID_REF,
  SYMBOL = promoters$`Gene Symbol`,
  stringsAsFactors = FALSE
)

# 3) Write out only the ~24 K promoter mappings
outf <- "Project/Data/MetaData/GPL_promoter_probe2gene.csv"
write.csv(probe2gene, outf, row.names=FALSE, quote=FALSE)
message("Wrote ", nrow(probe2gene), " promoter mappings to ", outf)
