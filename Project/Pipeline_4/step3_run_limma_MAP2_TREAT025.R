#!/usr/bin/env Rscript
# Pipeline_4/step3_run_limma_MAP2_TREAT025.R
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma); library(edgeR)
  library(dplyr); library(tibble); library(stringr); library(readr)
  library(AnnotationDbi); library(org.Mm.eg.db)
  library(janitor)
})

# ---- SETTINGS (edit if grouping guess needs tweaking) ------------------------
ACC          <- "GSE37577"
OUT_DIR      <- "Project/New_Data_4"
LFC_TREAT    <- 0.25               # TREAT threshold
PVAL_THRESH  <- 0.05               # nominal p-value cutoff for “significant” set
# Regex used to detect groups from pData. Left = WT/Control, Right = Mutant.
REGEX_WT     <- "(control|wt)"
REGEX_MUT    <- "(mutant|map2|ko|knockout|rb)"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Download series matrix + pheno --------------------------------------
message("Fetching GEO series: ", ACC)
gset <- getGEO(ACC, GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset) > 1) {
  # pick the largest by samples
  lens <- sapply(gset, ncol)
  eset <- gset[[which.max(lens)]]
} else {
  eset <- gset[[1]]
}

# Expression matrix (log-ish). Keep as numeric matrix.
expr <- Biobase::exprs(eset)
# Clean probe IDs
probe_ids <- rownames(expr)

# ---- 2) Build group factor (WT vs MUT) --------------------------------------
pd <- Biobase::pData(eset)
pd_str <- apply(pd, 1, function(r) paste(r, collapse = " | "))
grp <- ifelse(grepl(REGEX_MUT, pd_str, ignore.case = TRUE), "MUT",
         ifelse(grepl(REGEX_WT,  pd_str, ignore.case = TRUE), "WT", NA))

if (any(is.na(grp))) {
  stop("Some samples could not be assigned to WT/MUT by regex. ",
       "Open the script and adjust REGEX_WT/REGEX_MUT. Offending samples:\n",
       paste(colnames(expr)[is.na(grp)], collapse = ", "))
}
grp <- factor(grp, levels = c("WT","MUT"))
message("Group counts: WT=", sum(grp=="WT"), " | MUT=", sum(grp=="MUT"))

# ---- 3) limma (assume already log-scale; use limma directly) -----------------
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
fit <- lmFit(expr, design)
cont <- makeContrasts(MUT_vs_WT = MUT - WT, levels = design)
fit2 <- contrasts.fit(fit, cont)
fitT <- treat(fit2, lfc = LFC_TREAT)

# Full table (all rows), unsorted to preserve order
tt_full <- limma::topTreat(fitT, coef = "MUT_vs_WT", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column(var = "ID")
# ---- 4) Map probe -> Gene Symbol using platform annotation -------------------
fd <- Biobase::fData(eset)

# try to find a symbol-like column on the platform / series
sym_col <- NA_character_
candidates <- names(fd)
sym_hits <- grep("(^|_)gene( |_)?symbol|^symbol$|gene.symbol|geneSymbol|SYMBOL",
                 candidates, ignore.case = TRUE, value = TRUE)
if (length(sym_hits) > 0) sym_col <- sym_hits[1]

if (!is.na(sym_col)) {
  # Use platform-provided symbols
  map <- data.frame(
    ID     = rownames(fd),
    SYMBOL = as.character(fd[[sym_col]]),
    stringsAsFactors = FALSE
  )

  # Clean up multi-symbol cells like "Cdkn1a /// P21"
  map$SYMBOL <- gsub("\\s*///\\s*", ";", map$SYMBOL)
  map$SYMBOL <- sapply(strsplit(map$SYMBOL, "[;,/ ]+"), function(x) x[which(nzchar(x))[1]])
  map$SYMBOL[map$SYMBOL == "" | is.na(map$SYMBOL)] <- NA_character_

} else if (any(grepl("^ENSMUSG", rownames(expr)))) {
  # Fallback only if rows are Ensembl gene IDs
  suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Mm.eg.db) })
  map <- AnnotationDbi::select(org.Mm.eg.db,
                               keys    = rownames(expr),
                               keytype = "ENSEMBL",
                               columns = c("SYMBOL"))
  map <- dplyr::distinct(map, ENSEMBL, .keep_all = TRUE)
  names(map)[1] <- "ID"
} else {
  # No mapping available; keep NA symbols and use probe IDs later
  map <- data.frame(ID = rownames(expr), SYMBOL = NA_character_)
}

# join onto the full table
map <- dplyr::distinct(map, ID, .keep_all = TRUE)
tt_anno <- dplyr::left_join(tt_full, map, by = "ID")

# Prefer SYMBOL when available; else keep probe ID
tt_anno <- tt_anno %>%
  dplyr::mutate(Gene = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", ID, SYMBOL))

# ---- 5) Collapse to gene-level & export -------------------------------------
# Prefer SYMBOL when available; else keep ID
tt_anno <- tt_anno |>
  mutate(Gene = if_else(is.na(SYMBOL) | SYMBOL=="", ID, SYMBOL))

# Standard column names for downstream Python
out_full <- tt_anno |>
  transmute(
    Gene           = Gene,
    log2FoldChange = logFC,
    pvalue         = P.Value,
    padj           = adj.P.Val,
    baseMean       = AveExpr
  )

# Significant set: TREAT p<0.05 (already respects lfc via treat)
out_sig <- out_full |> filter(pvalue < PVAL_THRESH)

# Sanity summary
cat(sprintf("Rows: full = %d | TREAT p<%.2f = %d\n",
            nrow(out_full), PVAL_THRESH, nrow(out_sig)))
qs <- quantile(abs(out_sig$log2FoldChange), probs = c(.5,.9,.95,.99), na.rm = TRUE)
cat(sprintf("abs(logFC) percentiles (sig): 50%%=%.3f, 90%%=%.3f, 95%%=%.3f, 99%%=%.3f\n",
            qs[1], qs[2], qs[3], qs[4]))

# Write
write_csv(out_full, file.path(OUT_DIR, "Oshikawa_MAP2_DE_TREAT025_full.csv"))
write_csv(out_sig,  file.path(OUT_DIR, "Oshikawa_MAP2_DE_TREAT025_significant.csv"))
cat("→ Wrote:\n  ",
    file.path(OUT_DIR, "Oshikawa_MAP2_DE_TREAT025_full.csv"), "\n  ",
    file.path(OUT_DIR, "Oshikawa_MAP2_DE_TREAT025_significant.csv"), "\n")
