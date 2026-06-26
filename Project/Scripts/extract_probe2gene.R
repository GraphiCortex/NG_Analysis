#!/usr/bin/env Rscript

# extract_probe2gene.R
# Read GPL18280-19916.txt, filter to promoter probes, and pull SPOT_ID → GENE_SYMBOL

# 1. Load the raw GPL annotation
gpl_raw <- read.delim(
  file.path("Project","MetaData","GPL18280-19916.txt"),
  comment.char    = "#",
  sep             = "\t",
  stringsAsFactors = FALSE,
  check.names     = FALSE
)

# 2. Filter to PROMOTER probes only
stopifnot(all(c("SPOT_ID","DESCRIPTION","GENE_SYMBOL") %in% colnames(gpl_raw)))
promoter_map <- subset(
  gpl_raw,
  DESCRIPTION == "PROMOTER" &
    !is.na(GENE_SYMBOL) &
    GENE_SYMBOL != ""
)
message(
  "Filtering to PROMOTER probes: found ",
  nrow(promoter_map),
  " entries"
)

# 3. Build the two-column mapping
probe2gene <- data.frame(
  ID_REF           = promoter_map$SPOT_ID,
  SYMBOL           = promoter_map$GENE_SYMBOL,
  stringsAsFactors = FALSE
)

# 4. Write out filtered mapping
outdir <- file.path("Project","MetaData")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write.csv(
  probe2gene,
  file      = file.path(outdir, "GPL_probe2gene.csv"),
  row.names = FALSE,
  quote     = FALSE
)

cat(
  "→ Wrote", nrow(probe2gene),
  "promoter probe→gene mappings to",
  file.path(outdir, "GPL_probe2gene.csv"),
  "\n"
)
