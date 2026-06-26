#!/usr/bin/env Rscript

# new_limma_sig_probes.R
# Performs limma differential binding analysis and exports significant probes
# (FDR < 0.05 & |logFC| ≥ 0.5) for each contrast.

# -- 0. Setup ------------------------------------------------------------------
if (!requireNamespace("limma", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("limma")
}
library(limma)

# -- 1. Specify input files --------------------------------------------------
# Adjust paths as needed to your processed M-value files
files <- list.files(
    path = "Project/Data/Raw/ChipSeq/Julian2016/Processed/",
    pattern = "_Processed\\.txt$", full.names = TRUE
)

# -- 2. Read & assemble M-value matrix ----------------------------------------
Mlist <- lapply(files, function(f) {
    d <- read.delim(f, stringsAsFactors = FALSE)
    if ("M" %in% names(d)) {
        vals <- d$M
    } else if (all(c("IP", "INPUT") %in% names(d))) {
        vals <- log2(d$IP / d$INPUT)
    } else {
        stop("File ", f, " must contain either M or IP & INPUT columns.")
    }
    setNames(vals, d$ID_REF)
})
Mvals <- do.call(cbind, Mlist)
colnames(Mvals) <- sub("_Processed\\.txt$", "", basename(files))

# -- 3. Define experimental design -------------------------------------------
groups <- ifelse(grepl("E2f4", colnames(Mvals), ignore.case = TRUE), "E2F4",
            ifelse(grepl("E2f3_wt", colnames(Mvals), ignore.case = TRUE), "E2F3_wt",
            ifelse(grepl("E2f3_3amut", colnames(Mvals), ignore.case = TRUE), "E2F3_3amut",
            ifelse(grepl("E2f3_3bmut", colnames(Mvals), ignore.case = TRUE), "E2F3_3bmut",
                   NA))))
if (any(is.na(groups))) stop("Could not parse group names from file names.")
group <- factor(groups)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# -- 4. Fit linear model & contrasts -----------------------------------------
fit <- lmFit(Mvals, design)
cont <- makeContrasts(
    E2F4_vs_E2F3   = E2F4 - E2F3_wt,
    E2F4_vs_3amut  = E2F4 - E2F3_3amut,
    E2F4_vs_3bmut  = E2F4 - E2F3_3bmut,
    E2F3_wt_vs_3amut = E2F3_wt - E2F3_3amut,
    E2F3_wt_vs_3bmut = E2F3_wt - E2F3_3bmut,
    levels = design
)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

# -- 5. Extract & save significant probes ------------------------------------
outdir <- "Project/Data/Processed/limma_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (co in colnames(cont)) {
    # retrieve all results then filter
    tab <- topTable(
        fit2,
        coef          = co,
        number        = Inf,
        adjust.method = "fdr"
    )
    sig <- tab[ tab$adj.P.Val < 0.05 & abs(tab$logFC) >= 0.5, ]

    # report counts
    message(sprintf("%s → %d significant probes", co, nrow(sig)))

    # write CSV
    write.csv(
        sig,
        file = file.path(outdir, paste0("limma_", co, "_sig_probes.csv")),
        quote = FALSE,
        row.names = TRUE
    )
}

message("Finished exporting significant probe lists.")
