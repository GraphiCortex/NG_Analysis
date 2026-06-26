#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(tidyverse)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

gse_id  <- "GSE37577"
outdir  <- "Project/New_Data_5"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Fetching GEO series: ", gse_id)
gsets <- getGEO(gse_id, GSEMatrix = TRUE)
if (length(gsets) != 1) stop("Unexpected number of ExpressionSets in series.")
eset_all <- gsets[[1]]

pd <- pData(eset_all)

# ---- pick CAG samples (3 WT + 3 TKO) from titles ----
is_cag  <- grepl("^CAG-", pd$title) | grepl("\\bCAG\\b", pd$source_name_ch1 %||% "")
eset    <- eset_all[, is_cag]
pd_cag  <- pData(eset)

grp <- ifelse(grepl("CAG-Cont", pd_cag$title), "WT",
              ifelse(grepl("CAG-TKO", pd_cag$title), "TKO", NA))
if (any(is.na(grp))) stop("Could not parse WT/TKO for some CAG samples.")
grp <- factor(grp, levels = c("WT","TKO"))

message(sprintf("Group counts (CAG): WT=%d | TKO=%d",
                sum(grp=="WT"), sum(grp=="TKO")))

# ---- expression matrix & design ----
E <- exprs(eset)
design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)

# ---- limma + TREAT (|log2FC| >= 0.25) ----
fit  <- lmFit(E, design)
fit2 <- contrasts.fit(fit, makeContrasts(TKO - WT, levels = design))
fitT <- treat(fit2, lfc = 0.25)

tt_full <- topTreat(fitT, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("PROBE")

# ---- map probe -> SYMBOL from feature data, with a few common header aliases ----
fd   <- fData(eset)
sym_col <- c("Gene Symbol","GENE_SYMBOL","GeneSymbol","Symbol","SYMBOL")
sym_col <- sym_col[sym_col %in% colnames(fd)][1]
if (is.na(sym_col)) sym_col <- colnames(fd)[1]  # last resort

probe2sym <- tibble(PROBE = rownames(fd),
                    SYMBOL = as.character(fd[[sym_col]])) %>%
  distinct(PROBE, .keep_all = TRUE)

out_full <- tt_full %>%
  left_join(probe2sym, by = "PROBE") %>%
  transmute(
    Gene            = ifelse(is.na(SYMBOL) | SYMBOL=="", PROBE, SYMBOL),
    log2FoldChange  = logFC,
    pvalue          = P.Value,
    padj            = adj.P.Val,
    baseMean        = AveExpr
  )

sig <- out_full %>% filter(pvalue < 0.05, abs(log2FoldChange) >= 0.25)

message(sprintf("Rows: full = %d | TREAT p<0.05 = %d",
                nrow(out_full), sum(out_full$pvalue < 0.05)))
q <- quantile(abs(sig$log2FoldChange), c(.50,.90,.95,.99), na.rm = TRUE)
message(sprintf("abs(logFC) percentiles (sig): 50%%=%.3f, 90%%=%.3f, 95%%=%.3f, 99%%=%.3f",
                q[1], q[2], q[3], q[4]))

write.csv(out_full, file.path(outdir, "Oshikawa_CAG_DE_TREAT025_full.csv"), row.names = FALSE)
write.csv(sig,      file.path(outdir, "Oshikawa_CAG_DE_TREAT025_significant.csv"), row.names = FALSE)

# helpful bookkeeping
write.table(pd_cag[,c("title","geo_accession")],
            file.path(outdir, "CAG_sample_manifest.tsv"),
            sep="\t", row.names = FALSE, quote = FALSE)

message("→ Wrote:\n  ", file.path(outdir, "Oshikawa_CAG_DE_TREAT025_full.csv"),
        "\n  ",      file.path(outdir, "Oshikawa_CAG_DE_TREAT025_significant.csv"),
        "\n  ",      file.path(outdir, "CAG_sample_manifest.tsv"))
