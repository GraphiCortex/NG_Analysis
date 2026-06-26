#!/usr/bin/env Rscript

# limma_promoter_rawP_lfc_analysis.R
# Differential binding analysis using raw P-value (P<0.05) & |logFC|>=0.5
# then collapsing significant probes to genes for promoter-mapped array.

# 0. Load packages
if (!requireNamespace("limma", quietly=TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
library(limma)

# 1. File paths
processed_dir   <- "Project/Data/Raw/ChipSeq/Julian2016/Processed/"
ann_file        <- "Project/MetaData/GPL_probe2gene.csv"
gpl_file        <- "Project/MetaData/GPL18280-19916.txt"
outdir          <- "Project/Data/Processed/limma_results"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# 2. Load annotation
ann <- read.csv(ann_file, stringsAsFactors=FALSE)
probe_col <- if ("ProbeID" %in% colnames(ann)) {
  "ProbeID"
} else {
  "ID_REF"
}

# Properly chained if/else for gene column
if ("SYMBOL" %in% colnames(ann)) {
  gene_col <- "SYMBOL"
} else if ("GeneSymbol" %in% colnames(ann)) {
  gene_col <- "GeneSymbol"
} else if ("Gene.Symbol" %in% colnames(ann)) {
  gene_col <- "Gene.Symbol"
} else {
  stop("No gene symbol column in annotation.")
}

# Standardize
ann$ID_REF     <- ann[[probe_col]]
ann$GeneSymbol <- ann[[gene_col]]
ann <- ann[!is.na(ann$GeneSymbol) & ann$GeneSymbol!="", ]

# 3. Read M-values
files <- list.files(processed_dir, pattern="_Processed\\.txt$",
                    full.names=TRUE)
Mlist <- lapply(files, function(f) {
  df <- read.delim(f, stringsAsFactors=FALSE)
  if ("M" %in% colnames(df)) {
    vals <- df$M
  } else {
    vals <- log2(df$IP / df$INPUT)
  }
  setNames(vals, df$ID_REF)
})
Mvals <- do.call(cbind, Mlist)
colnames(Mvals) <- sub("_Processed\\.txt$", "", basename(files))

# 4. Restrict to true promoter probes via GPL18280 DESCRIPTION
cat("\nFiltering to GPL18280 PROMOTER probes…\n")
gpl <- read.delim(gpl_file, stringsAsFactors=FALSE, comment.char="#")
stopifnot(all(c("SPOT_ID","DESCRIPTION") %in% colnames(gpl)))
promoter_ids    <- gpl$SPOT_ID[gpl$DESCRIPTION == "PROMOTER"]
cat("  → GPL flags", length(promoter_ids), "probes as PROMOTER\n")
common_all      <- intersect(rownames(Mvals), promoter_ids)
cat("  → Promoter probes in your data:", length(common_all), "\n")
common_annotated <- intersect(common_all, ann$ID_REF)
Mvals            <- Mvals[common_annotated, , drop=FALSE]
cat("  → Final probe universe:", nrow(Mvals), "probes\n\n")

# 5. Design matrix
groups <- sapply(colnames(Mvals), function(nm) {
  if      (grepl("E2f4",      nm, ignore.case=TRUE)) "E2F4"
  else if (grepl("E2f3_wt",   nm, ignore.case=TRUE)) "E2F3_wt"
  else if (grepl("E2f3_3amut",nm, ignore.case=TRUE)) "E2F3_3amut"
  else if (grepl("E2f3_3bmut",nm, ignore.case=TRUE)) "E2F3_3bmut"
  else stop("Cannot parse sample name: ", nm)
})
group  <- factor(groups)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 6. Fit limma & contrasts
fit   <- lmFit(Mvals, design)
cont  <- makeContrasts(E2F4_vs_E2F3 = E2F4 - E2F3_wt, levels=design)
fit2  <- contrasts.fit(fit, cont)
fit2  <- eBayes(fit2)

# 7. Apply TREAT & extract significant probes
fit_treat <- treat(fit2, lfc = 0.5)
tab_sig   <- topTreat(
  fit_treat,
  coef          = "E2F4_vs_E2F3",
  number        = Inf,
  adjust.method = "none",
  p.value       = 0.05
)
cat("Significant probes (TREAT raw P<0.05 & |logFC|>=0.5):", nrow(tab_sig), "\n")

# 8. Collapse to genes
gene_hits <- unique(ann$GeneSymbol[ match(rownames(tab_sig), ann$ID_REF) ])
cat("Unique significant genes:", length(gene_hits), "\n")

# 9. Write files
write.csv(tab_sig,
          file = file.path(outdir, "limma_E2F4_vs_E2F3_promoter_probes.csv"),
          quote=FALSE)
write.csv(data.frame(GeneSymbol=gene_hits),
          file = file.path(outdir, "E2F4_vs_E2F3_sig_genes.csv"),
          quote=FALSE, row.names=FALSE)

cat("Done. Outputs in", outdir, "\n")

#
# 10. GO–BP over-representation analysis (BH FDR<0.05)
library(clusterProfiler)
library(org.Mm.eg.db)

# Foreground: your 1 523 genes
fg <- gene_hits

# Background: all promoter-filtered genes
bg <- unique(ann$GeneSymbol[ann$ID_REF %in% common_annotated])

ego <- enrichGO(
  gene          = fg,
  universe      = bg,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("Enriched GO–BP terms (FDR<0.05):", nrow(ego), "\n")
write.csv(
  as.data.frame(ego),
  file = file.path(outdir, "E2F4_vs_E2F3_GO_BP_enrichment.csv"),
  quote = FALSE, row.names = FALSE
)

cat("Done. All outputs in", outdir, "\n")