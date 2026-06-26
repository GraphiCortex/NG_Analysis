#!/usr/bin/env Rscript

# ==============================================
# Step 1: Probe → Gene mapping for promoters
# ==============================================

# 1. Read the full tiling array annotation, skipping comment lines:
gpl <- read.delim(
  "Project/MetaData/GPL18280-19916.txt",
  sep        = "\t",
  header     = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

# 2. Filter to promoter probes with non-blank gene symbols:
promoter_probes <- subset(
  gpl,
  DESCRIPTION == "PROMOTER" &
    GENE_SYMBOL != ""
)

# 3. Select just the probe ID → gene symbol columns:
probe2gene <- promoter_probes[, c("SPOT_ID", "GENE_SYMBOL")]

# 4. Write out the ∼24K promoter probe2gene mappings:
write.csv(
  probe2gene,
  file      = "Project/New_Data/GPL_promoter_probe2gene.csv",
  row.names = FALSE,
  quote     = FALSE
)
