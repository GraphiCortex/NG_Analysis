#!/usr/bin/env Rscript

# full_pipeline_localGPL_E2F4vsE2F3.R
# Reproduce the original E2F4_vs_E2F3 promoter-array analysis:
# 1) Parse local promoter annotation (~24K promoters)
# 2) Read processed M-values and subset to promoter probes
# 3) limma + TREAT (|logFC|>=0.5) + BH-FDR<0.05 -> ~1,256 significant probes
# 4) Collapse to ~1,045 unique genes
# 5) GO-BP ORA (BH-FDR<0.05) -> 312 enriched terms

# 0. Install & load dependencies
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("limma","clusterProfiler","org.Mm.eg.db","readr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg)
}
suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(readr)
})

# 1. Load promoter annotation file
#    This should be the GPL18280-19916.txt containing only promoter flags
ann_file <- "Project/MetaData/GPL18280-19916.txt"
message("Reading promoter‐probe annotation from ", ann_file)
ann_tbl <- read.delim(
  ann_file,
  stringsAsFactors = FALSE,
  comment.char = "#"
)
stopifnot(all(c("SPOT_ID","DESCRIPTION","GENE_SYMBOL") %in% colnames(ann_tbl)))
# Filter to PROMOTER
promoters <- subset(ann_tbl, DESCRIPTION == "PROMOTER")
message("Total promoter probes in annotation: ", nrow(promoters))
# Build probe -> gene map
ann_map <- subset(promoters, GENE_SYMBOL != "")
ann_map <- ann_map[, c("SPOT_ID","GENE_SYMBOL")]
colnames(ann_map) <- c("ID_REF","GeneSymbol")
message("Probe->gene mappings retained: ", nrow(ann_map))

# 2. Read processed M-values and subset to promoter probes
processed_dir <- "Project/Data/Raw/ChipSeq/Julian2016/Processed/"
files <- list.files(
  processed_dir,
  pattern = "_Processed\\.txt$",
  full.names = TRUE
)
Mlist <- lapply(files, function(f) {
  df <- read.delim(f, stringsAsFactors = FALSE)
  mvals <- if ("M" %in% colnames(df)) df$M else log2(df$IP / df$INPUT)
  setNames(mvals, df$ID_REF)
})
Mvals <- do.call(cbind, Mlist)
colnames(Mvals) <- sub("_Processed\\.txt$", "", basename(files))
# Subset to promoter probes
keep_probes <- intersect(rownames(Mvals), ann_map$ID_REF)
Mvals <- Mvals[keep_probes, , drop = FALSE]
message("Promoter‐probe universe: ", nrow(Mvals), " probes")

# 3. limma fit & contrasts
groups <- sapply(colnames(Mvals), function(nm) {
  if      (grepl("E2f4",       nm, ignore.case=TRUE)) "E2F4"
  else if (grepl("E2f3_wt",    nm, ignore.case=TRUE)) "E2F3_wt"
  else if (grepl("E2f3_3amut", nm, ignore.case=TRUE)) "E2F3_3amut"
  else if (grepl("E2f3_3bmut", nm, ignore.case=TRUE)) "E2F3_3bmut"
  else stop("Cannot parse sample name: ", nm)
})
group  <- factor(groups)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit  <- lmFit(Mvals, design)
fit2 <- contrasts.fit(fit,
          makeContrasts(E2F4_vs_E2F3 = E2F4 - E2F3_wt,
                       levels = design))
fit2 <- eBayes(fit2)

# 4. TREAT + BH-FDR<0.05 -> significant probes (~1,256)
fit_tr     <- treat(fit2, lfc = 0.5)
sig_probes <- topTreat(
  fit_tr,
  coef          = "E2F4_vs_E2F3",
  number        = Inf,
  adjust.method = "fdr",
  p.value       = 0.05
)
message("Significant promoter probes: ", nrow(sig_probes))
write.csv(
  sig_probes,
  file = file.path("Project/Data/Processed/limma_results",
                    "limma_E2F4_vs_E2F3_promoter_probes.csv"),
  quote = FALSE
)

# 5. Collapse to genes (~1,045)
sig_genes <- unique(ann_map$GeneSymbol[match(rownames(sig_probes), ann_map$ID_REF)])
message("Unique significant genes: ", length(sig_genes))
write.csv(
  data.frame(GeneSymbol = sig_genes),
  file = file.path("Project/Data/Processed/limma_results",
                   "E2F4_vs_E2F3_sig_genes.csv"),
  quote = FALSE,
  row.names = FALSE
)

# 6. GO-BP ORA (BH-FDR<0.05) -> 312 enriched terms
ego <- enrichGO(
  gene          = sig_genes,
  universe      = unique(ann_map$GeneSymbol),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
message("Enriched GO-BP terms: ", nrow(ego))
write.csv(
  as.data.frame(ego),
  file = file.path("Project/Data/Processed/limma_results",
                   "E2F4_vs_E2F3_GO_BP_enrichment.csv"),
  quote = FALSE,
  row.names = FALSE
)

message("✅ Pipeline complete. Check Project/Data/Processed/limma_results for outputs.")
