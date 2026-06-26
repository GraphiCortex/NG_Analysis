#!/usr/bin/env Rscript
# Step 7: Write GO results to Excel (5 families; robust to column names)

library(readr); library(dplyr); library(purrr); library(openxlsx)

in_base <- "Project/New_Data"
probe_map_csv <- file.path(in_base,"E2F4_vs_E2F3_sigProbe2Gene.csv")
sig_probes_csv <- file.path(in_base,"E2F4_vs_E2F3_significant_probes.csv")
prom_map_csv   <- file.path(in_base,"GPL_promoter_probe2gene.csv")

full_csv <- file.path(in_base,"E2F4_vs_E2F3_GO_BP_full.csv")
out_xlsx <- file.path(in_base,"GO_results.xlsx")

# ---------- helpers ----------
pick_col <- function(df, candidates) {
  hits <- intersect(candidates, names(df))
  if (length(hits)) hits[[1]] else NA_character_
}

load_category <- function(cat){
  byid <- file.path(in_base, sprintf("E2F4_vs_E2F3_GO_BP_%s_byID.csv", cat))
  kw   <- file.path(in_base, sprintf("E2F4_vs_E2F3_GO_BP_%s.csv",    cat))
  if (file.exists(byid)) read_csv(byid, col_types = cols())
  else if (file.exists(kw)) read_csv(kw, col_types = cols())
  else tibble()
}

make_genes_df_from_sigmap <- function(sig_map) {
  gcol <- pick_col(sig_map, c("GENE_SYMBOL","Gene","gene_symbol","SYMBOL","gene"))
  lcol <- pick_col(sig_map, c("logFC","LogFC","LOGFC"))
  if (is.na(gcol) || is.na(lcol)) return(NULL)
  out <- sig_map %>%
    filter(!is.na(.data[[gcol]])) %>%
    group_by(across(all_of(gcol))) %>%
    summarize(mean_logFC = mean(.data[[lcol]], na.rm = TRUE), .groups = "drop")
  names(out)[names(out) == gcol] <- "Gene"
  arrange(out, desc(mean_logFC))
}

make_genes_df_fallback <- function() {
  if (!file.exists(sig_probes_csv) || !file.exists(prom_map_csv)) return(NULL)
  sp <- read_csv(sig_probes_csv, col_types = cols())
  pm <- read_csv(prom_map_csv,   col_types = cols())
  if (!("SPOT_ID" %in% names(sp)) || !("SPOT_ID" %in% names(pm))) return(NULL)
  df <- left_join(sp, pm, by = "SPOT_ID")
  make_genes_df_from_sigmap(df)
}

# ---------- 1) build genes_df robustly ----------
sig_map <- if (file.exists(probe_map_csv)) read_csv(probe_map_csv, col_types = cols()) else tibble()
genes_df <- NULL
if (nrow(sig_map)) genes_df <- make_genes_df_from_sigmap(sig_map)
if (is.null(genes_df)) {
  message("Rebuilding genes_df from significant_probes + promoter map (fallback)…")
  genes_df <- make_genes_df_fallback()
}
if (is.null(genes_df)) {
  stop("Step 7: Could not determine gene and/or logFC columns. ",
       "Checked: ", basename(probe_map_csv), ", ",
       basename(sig_probes_csv), ", ", basename(prom_map_csv))
}

# ---------- 2) full GO table ----------
full_df <- read_csv(full_csv, col_types = cols())

# ---------- 3) load the five categories ----------
cc_df <- load_category("cell_cycle")
dr_df <- load_category("DNA_repair")
ap_df <- load_category("apoptosis")
au_df <- load_category("autophagy")
ne_df <- load_category("necroptosis")  # NOTE: necroptosis (not necrosis)

# ---------- 4) add Genes_logFC to category sheets ----------
logfc_map <- set_names(genes_df$mean_logFC, genes_df$Gene)
fmt_genes <- function(ids){
  if (is.null(ids) || is.na(ids)) return("")
  g <- strsplit(ids, "/")[[1]]
  vals_chr <- ifelse(is.na(logfc_map[g]), "NA", sprintf("%.3f", logfc_map[g]))
  paste0(g, " (", vals_chr, ")", collapse = "; ")
}
add_glfc <- function(df){
  if (nrow(df) == 0) df else transform(df, Genes_logFC = purrr::map_chr(geneID, fmt_genes))
}
cc_df <- add_glfc(cc_df); dr_df <- add_glfc(dr_df); ap_df <- add_glfc(ap_df)
au_df <- add_glfc(au_df); ne_df <- add_glfc(ne_df)

# ---------- 5) annotate full_df with Category labels ----------
full_df <- full_df %>%
  mutate(Category = case_when(
    ID %in% cc_df$ID ~ "Cell cycle",
    ID %in% dr_df$ID ~ "DNA repair",
    ID %in% ap_df$ID ~ "Apoptosis",
    ID %in% au_df$ID ~ "Autophagy",
    ID %in% ne_df$ID ~ "Necroptosis",
    TRUE ~ "Other"
  ))

# ---------- 6) workbook ----------
wb <- createWorkbook()
for (s in paste0("Sh_", 1:7)) addWorksheet(wb, s)

writeData(wb,"Sh_1", genes_df)
writeData(wb,"Sh_2", full_df)
writeData(wb,"Sh_3", cc_df)
writeData(wb,"Sh_4", dr_df)
writeData(wb,"Sh_5", ap_df)
writeData(wb,"Sh_6", au_df)
writeData(wb,"Sh_7", ne_df)

# shading on Sh_2
styles <- list(
  "Cell cycle"  = createStyle(fgFill="#FFF2CC"),
  "DNA repair"  = createStyle(fgFill="#D9EAD3"),
  "Apoptosis"   = createStyle(fgFill="#F4CCCC"),
  "Autophagy"   = createStyle(fgFill="#C6EFCE"),
  "Necroptosis" = createStyle(fgFill="#E6E6FA")
)
for (i in seq_len(nrow(full_df))) {
  st <- styles[[ full_df$Category[i] ]]
  if (!is.null(st)) addStyle(wb,"Sh_2", st, rows = i+1, cols = seq_len(ncol(full_df)), gridExpand = TRUE)
}
freezePane(wb,"Sh_2", firstActiveRow = 2)
for (s in paste0("Sh_", 1:7)) setColWidths(wb, s, 1:50, "auto")

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("→ Excel workbook written to", out_xlsx, "\n")
