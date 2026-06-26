#!/usr/bin/env Rscript

# limma_analysis.R
# Differential analysis of Chip–Chip M-values using limma

#– Load libraries
if (!requireNamespace("limma", quietly=TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("limma")
}
library(limma)

#– 1. Read in your processed M-value text files
files <- c(
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326858_E2f4_wt_replicate1_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326859_E2f4_wt_replicate2_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326860_E2f4_wt_replicate3_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326861_E2f3_wt_replicate1_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326862_E2f3_wt_replicate2_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326863_E2f3_wt_replicate3_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326864_E2f3_3amut_replicate1_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326865_E2f3_3amut_replicate2_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326866_E2f3_3amut_replicate3_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326867_E2f3_3bmut_replicate1_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326868_E2f3_3bmut_replicate2_Processed.txt",
    "Project/Data/Raw/ChipSeq/Julian2016/Processed/GSM1326869_E2f3_3bmut_replicate3_Processed.txt"
)

# Read each file and pull out the M column, indexed by probe ID
Mlist <- lapply(files, function(f){
  d <- read.delim(f, stringsAsFactors=FALSE)
  # compute the log2 enrichment
  d$M <- log2(d$IP / d$INPUT)
  setNames(d$M, d$ID_REF)
})

# Combine into one matrix
Mvals <- do.call(cbind, Mlist)
colnames(Mvals) <- c(
  paste0("E2F4_rep", 1:3),
  paste0("E2F3_wt_rep", 1:3),
  paste0("E2F3_3amut_rep", 1:3),
  paste0("E2F3_3bmut_rep", 1:3)
)
rownames(Mvals) <- names(Mlist[[1]])

#– 2. Build your design matrix
group <- factor(rep(c("E2F4","E2F3_wt","E2F3_3amut","E2F3_3bmut"), each=3))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#– 3. Fit the linear model & specify your contrasts
fit <- lmFit(Mvals, design)
cont <- makeContrasts(
  E2F4_vs_E2F3   = E2F4      - E2F3_wt,
  E2F4_vs_3amut  = E2F4      - E2F3_3amut,
  E2F4_vs_3bmut  = E2F4      - E2F3_3bmut,
  E2F3_wt_vs_3amut = E2F3_wt - E2F3_3amut,
  E2F3_wt_vs_3bmut = E2F3_wt - E2F3_3bmut,
  levels=design
)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

#– 4. Save topTables & MA-plots
outdir <- "Project/Data/Processed/limma_results"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Write results
for(co in colnames(cont)) {
  tab <- topTable(fit2, coef=co, number=Inf, adjust.method="fdr")
  write.csv(tab, file=file.path(outdir, paste0("limma_", co, ".csv")), quote=FALSE)
}

# MA-plots into a single PDF
pdf(file.path(outdir, "limma_MA_plots.pdf"))
for(co in colnames(cont)) {
  plotMA(fit2, coef=co, main=co)
}
dev.off()

cat("Limma analysis complete. Results in", outdir, "\n")
