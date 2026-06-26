#!/usr/bin/env Rscript
# Fong_pipeline/fong_step_limma_TREAT025.R  (fixed: log2 + between-array norm)

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(janitor)
  library(limma)
  library(AnnotationDbi); library(org.Mm.eg.db)
  library(tibble)
})
select  <- dplyr::select
rename  <- dplyr::rename
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

in_path <- "Project/Data/Raw/microarray/Fong2022/Fong2022_raw/FDSC Batch IHW.xlsx"
sheet   <- 1

# main “clean” outputs your pipeline consumes
out_full_csv <- "Project/New_Data/Fong2022_DE_TREAT025_full.csv"
out_sig_csv  <- "Project/New_Data/Fong2022_DE_TREAT025_significant.csv"

# -------------------------------
# 1) Load + basic cleaning
# -------------------------------
raw <- read_excel(in_path, sheet = sheet, .name_repair = "minimal") |>
  remove_empty(c("rows","cols")) |>
  clean_names()

# choose an ID column heuristically
id_col_candidates <- c("row_names","ensembl_id","gene","symbol","gene_name","id")
id_col <- intersect(id_col_candidates, names(raw))[1]
if (is.na(id_col)) id_col <- names(raw)[1]

# sample columns: RB_THC_* (WT) and RB_TKO_* (TKO)
sample_cols <- names(raw)[grepl("^rb_(thc|wt)_\\d+$|^rb_tko_\\d+$", names(raw), ignore.case = TRUE)]
if (length(sample_cols) < 4) {
  # fallback: numeric-looking columns excluding id
  num_cols <- names(raw)[sapply(raw, function(z) suppressWarnings(!all(is.na(as.numeric(z)))))]
  sample_cols <- setdiff(num_cols, id_col)
}
stopifnot(length(sample_cols) >= 4)

# -------------------------------
# 2) Build matrix, log2, normalize
# -------------------------------
expr_df <- raw |>
  select(all_of(c(id_col, sample_cols))) |>
  mutate(across(all_of(sample_cols), \(x) suppressWarnings(as.numeric(x)))) |>
  distinct(.data[[id_col]], .keep_all = TRUE)

gene_id <- as.character(expr_df[[id_col]])
X <- as.matrix(expr_df[, sample_cols])
mode(X) <- "numeric"
rownames(X) <- gene_id
colnames(X) <- sample_cols

# log2 transform (avoid log2(0))
X <- log2(pmax(X, 1))

# between-array normalization (quantile) — standard for microarrays
X <- normalizeBetweenArrays(X, method = "quantile")

# -------------------------------
# 3) Design/contrast (WT vs TKO)
# -------------------------------
group <- ifelse(grepl("tko", colnames(X), ignore.case = TRUE), "TKO", "WT")
group <- factor(group, levels = c("WT","TKO"))
design <- model.matrix(~ 0 + group); colnames(design) <- c("WT","TKO")

# fit limma on log2-normalized intensities
fit  <- lmFit(X, design)
fit2 <- contrasts.fit(fit, makeContrasts(TKO_vs_WT = TKO - WT, levels = design))

# TREAT with relaxed threshold (log2 scale)
fitT <- treat(fit2, lfc = 0.25)

# -------------------------------
# 4) Tables + quick sanity prints
# -------------------------------
tt_full <- limma::topTreat(fitT, coef = "TKO_vs_WT", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column(var = "ID")

USE_STRICT <- TRUE  # TRUE => keep TREAT p<0.05 AND |logFC|>=0.25
tt_sig <- dplyr::filter(tt_full, P.Value < 0.05)
if (USE_STRICT) tt_sig <- dplyr::filter(tt_sig, abs(logFC) >= 0.25)

cat(sprintf("Rows: full = %d | TREAT p<0.05 = %d%s\n",
            nrow(tt_full),
            sum(tt_full$P.Value < 0.05),
            if (USE_STRICT) sprintf(" | p<0.05 & |logFC|≥0.25 = %d", nrow(tt_sig)) else ""))

# sanity: logFC distribution should be modest (not hundreds/thousands)
q <- quantile(tt_full$logFC, c(.5,.9,.95,.99), na.rm=TRUE)
cat(sprintf("logFC summary (median/90/95/99%%): %.3f / %.3f / %.3f / %.3f\n",
            median(tt_full$logFC, na.rm=TRUE), q[2], q[3], q[4]))

# -------------------------------
## --- REPLACE your Step 5 mapping block with this robust version ---
## 5) Map IDs → SYMBOL (mouse), supporting mixed Ensembl/Alias keys ----------
ids <- as.character(tt_sig$ID)
ens_keys <- ids[grepl("^ENSMUSG", ids)]
ali_keys <- setdiff(ids, ens_keys)

map_ens <- if (length(ens_keys)) {
  suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db, keys = ens_keys, keytype = "ENSEMBL", columns = "SYMBOL"
  )) |> dplyr::rename(ID = ENSEMBL)
} else tibble::tibble(ID = character(), SYMBOL = character())

map_ali <- if (length(ali_keys)) {
  suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db, keys = ali_keys, keytype = "ALIAS", columns = "SYMBOL"
  )) |> dplyr::rename(ID = ALIAS)
} else tibble::tibble(ID = character(), SYMBOL = character())

map_sym <- dplyr::bind_rows(map_ens, map_ali) |>
  dplyr::mutate(ID = as.character(ID)) |>
  dplyr::distinct(ID, .keep_all = TRUE)

tt <- dplyr::left_join(tt_sig |> dplyr::mutate(ID = as.character(ID)),
                       map_sym, by = "ID")

# ensure SYMBOL exists even if nothing mapped (rare but safe)
if (!"SYMBOL" %in% names(tt)) tt$SYMBOL <- NA_character_

## 6) Pipeline-compatible columns (no 'B' in topTreat output) -----------------
out <- tt |>
  dplyr::transmute(
    Gene           = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", ID, SYMBOL),
    log2FoldChange = logFC,
    pvalue         = P.Value,
    padj           = adj.P.Val,
    baseMean       = AveExpr
  )

# also produce “significant” view if you still want a p<0.05 CSV
out_sig <- tt |>
  dplyr::filter(P.Value < 0.05) |>
  dplyr::transmute(
    Gene           = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", ID, SYMBOL),
    log2FoldChange = logFC,
    pvalue         = P.Value,
    padj           = adj.P.Val,
    baseMean       = AveExpr
  )

readr::write_csv(out,     "Project/New_Data/Fong2022_DE_TREAT025_full.csv")
readr::write_csv(out_sig, "Project/New_Data/Fong2022_DE_TREAT025_significant.csv")
cat("→ Wrote full + significant:\n",
    "   Project/New_Data/Fong2022_DE_TREAT025_full.csv\n",
    "   Project/New_Data/Fong2022_DE_TREAT025_significant.csv\n")
