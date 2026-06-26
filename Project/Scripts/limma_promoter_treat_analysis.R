#!/usr/bin/env Rscript

# limma_promoter_treat_analysis.R
# Differential binding analysis restricted to promoter-mapped probes
# using limma's TREAT method (|logFC| ≥ 0.5) and FDR < 0.05.

# 0. Load required packages
if (!requireNamespace("limma", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("limma")
}
library(limma)

# 1. Paths & annotation
processed_dir <- "Project/Data/Raw/ChipSeq/Julian2016/Processed/"
annotation_file <- "Project/MetaData/GPL_probe2gene.csv"
outdir <- "Project/Data/Processed/limma_results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 2. Read and inspect probe-to-gene annotation
ann <- read.csv(annotation_file, stringsAsFactors = FALSE)
cat("Annotation columns detected:\n")
print(colnames(ann))
cat("First rows of annotation:\n")
print(head(ann, 3))

# Attempt to identify key columns
probe_col <- if ("ProbeID" %in% colnames(ann)) {
    "ProbeID"
} else if ("ID_REF" %in% colnames(ann)) {
    "ID_REF"
} else {
    stop("Cannot find 'ProbeID' or 'ID_REF' column in annotation file.")
}

gene_col <- if ("GeneSymbol" %in% colnames(ann)) {
    "GeneSymbol"
} else if ("Gene.Symbol" %in% colnames(ann)) {
    "Gene.Symbol"
} else if ("Symbol" %in% colnames(ann)) {
    "Symbol"
} else if ("SYMBOL" %in% colnames(ann)) {
    "SYMBOL"
} else {
    stop("Cannot find a gene-symbol column ('GeneSymbol', 'Gene.Symbol', 'Symbol', or 'SYMBOL') in annotation file.")
}

# Standardize names
ann$ID_REF <- ann[[probe_col]]
ann$GeneSymbol <- ann[[gene_col]]
# Filter out missing symbols
ann <- ann[!is.na(ann$GeneSymbol) & ann$GeneSymbol != "", ]
cat("Retained", nrow(ann), "rows with valid probe→gene mappings.\n")

# 3. Load processed M-values and subset to promoter probes
files <- list.files(processed_dir, pattern = "_Processed\\.txt$", full.names = TRUE)
Mlist <- lapply(files, function(f) {
    d <- read.delim(f, stringsAsFactors = FALSE)
    if ("M" %in% colnames(d)) {
        mvals <- d$M
    } else if (all(c("IP","INPUT") %in% colnames(d))) {
        mvals <- log2(d$IP / d$INPUT)
    } else {
        stop("File ", f, " missing M or IP/INPUT columns.")
    }
    setNames(mvals, d$ID_REF)
})
Mvals_full <- do.call(cbind, Mlist)
colnames(Mvals_full) <- sub("_Processed\\.txt$", "", basename(files))

common_probes <- intersect(rownames(Mvals_full), ann$ID_REF)
Mvals <- Mvals_full[common_probes, , drop = FALSE]
cat("Total processed probes:", nrow(Mvals_full), "\n")
cat("Promoter-mapped probes:", nrow(Mvals), "\n")

# 4. Define experimental design
groups <- sapply(colnames(Mvals), function(nm) {
    if (grepl("E2f4", nm, ignore.case = TRUE)) return("E2F4")
    if (grepl("E2f3_wt", nm, ignore.case = TRUE)) return("E2F3_wt")
    if (grepl("E2f3_3amut", nm, ignore.case = TRUE)) return("E2F3_3amut")
    if (grepl("E2f3_3bmut", nm, ignore.case = TRUE)) return("E2F3_3bmut")
    stop("Cannot parse group from sample name: ", nm)
})
group <- factor(groups)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 5. Fit linear model & contrasts
fit <- lmFit(Mvals, design)
cont <- makeContrasts(
    E2F4_vs_E2F3 = E2F4 - E2F3_wt,
    levels = design
)
fit2 <- contrasts.fit(fit, cont)

# 6. Apply TREAT for fold-change threshold
fit_treat <- treat(fit2, lfc = 0.5)

# 7. Extract significant probes with FDR < 0.05
tab_sig <- topTreat(
    fit_treat,
    coef          = "E2F4_vs_E2F3",
    number        = Inf,
    adjust.method = "fdr",
    p.value       = 0.05
)
cat("Significant promoter-mapped probes (FDR<0.05 & |logFC|>=0.5):", nrow(tab_sig), "\n")

# 8. Map to genes and collapse
gene_hits <- unique(ann$GeneSymbol[match(rownames(tab_sig), ann$ID_REF)])
cat("Unique significant genes:", length(gene_hits), "\n")

# 9. Write outputs
write.csv(
    tab_sig,
    file = file.path(outdir, "limma_E2F4_vs_E2F3_sig_promoter_probes.csv"),
    quote = FALSE,
    row.names = TRUE
)
write.csv(
    data.frame(GeneSymbol = gene_hits),
    file = file.path(outdir, "E2F4_vs_E2F3_sig_genes.csv"),
    quote = FALSE,
    row.names = FALSE
)

cat("Finished promoter-TREAT analysis. Files written to:", outdir, "\n")
