#!/usr/bin/env Rscript
# ==============================================
# Pipeline_2 — Step 7: Write GO results (non-promoter)
# Sheets: genes, full, Cell cycle, DNA repair, Apoptosis, Autophagy, Necroptosis
# ==============================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(openxlsx)
})

in_base   <- "Project/New_Data_2"
probe_csv <- file.path(in_base, "E2F4_vs_E2F3_sigProbe2Gene_nonpromoter.csv")
full_csv  <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_full_nonpromoter.csv")

cc_csv <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_cell_cycle_byID_nonpromoter.csv")
dr_csv <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_DNA_repair_byID_nonpromoter.csv")
ap_csv <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_apoptosis_byID_nonpromoter.csv")
au_csv <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_autophagy_byID_nonpromoter.csv")
ne_csv <- file.path(in_base, "E2F4_vs_E2F3_GO_BP_necrosis_byID_nonpromoter.csv")

out_xlsx  <- file.path(in_base, "GO_results_nonpromoter.xlsx")

safe_read <- function(p) if (file.exists(p)) read_csv(p, col_types = cols()) else tibble()

# --- 1) Build genes_df (Gene + mean_logFC), robust to column names ---
sig_map <- safe_read(probe_csv)
if (!nrow(sig_map)) stop("File not found or empty: ", probe_csv)

gene_candidates <- c("GENE_SYMBOL","Gene","SYMBOL","gene","GeneSymbol","GENE","GENE.SYMBOL","Gene_Symbol")
gene_col <- NA_character_
for (nm in names(sig_map)) {
  if (toupper(nm) %in% toupper(gene_candidates)) { gene_col <- nm; break }
}
if (is.na(gene_col)) stop("Could not find a gene symbol column in: ", probe_csv)

if ("logFC" %in% names(sig_map)) {
  genes_df <- sig_map %>%
    group_by(.data[[gene_col]]) %>%
    summarize(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_logFC))
} else {
  genes_df <- sig_map %>%
    distinct(.data[[gene_col]]) %>%
    mutate(mean_logFC = NA_real_) %>%
    arrange(.data[[gene_col]])
}
colnames(genes_df)[1] <- "Gene"

# --- 2) Read enrichment tables (five categories) ---
full_df <- safe_read(full_csv)
cc_df   <- safe_read(cc_csv)
dr_df   <- safe_read(dr_csv)
ap_df   <- safe_read(ap_csv)
au_df   <- safe_read(au_csv)
ne_df   <- safe_read(ne_csv)

# --- 3) Add “Genes_logFC” columns to category sheets ---
logfc_map <- setNames(genes_df$mean_logFC, genes_df$Gene)
fmt_genes <- function(ids){
  if (is.null(ids) || is.na(ids) || ids == "") return("")
  g <- strsplit(ids, "/")[[1]]
  v <- logfc_map[g]
  parts <- ifelse(is.na(v), paste0(g, " (NA)"),
                           paste0(g, " (", sprintf("%.3f", v), ")"))
  paste(parts, collapse = "; ")
}
add_glfc <- function(df){
  if (!nrow(df) || !"geneID" %in% names(df)) return(df)
  df$Genes_logFC <- map_chr(df$geneID, fmt_genes)
  df
}
cc_df <- add_glfc(cc_df)
dr_df <- add_glfc(dr_df)
ap_df <- add_glfc(ap_df)
au_df <- add_glfc(au_df)
ne_df <- add_glfc(ne_df)

# --- 4) Annotate full_df with the five categories ---
if (nrow(full_df)) {
  full_df$Category <- "Other"
  if (nrow(cc_df)) full_df$Category[full_df$ID %in% cc_df$ID] <- "Cell cycle"
  if (nrow(dr_df)) full_df$Category[full_df$ID %in% dr_df$ID] <- "DNA repair"
  if (nrow(ap_df)) full_df$Category[full_df$ID %in% ap_df$ID] <- "Apoptosis"
  if (nrow(au_df)) full_df$Category[full_df$ID %in% au_df$ID] <- "Autophagy"
  if (nrow(ne_df)) full_df$Category[full_df$ID %in% ne_df$ID] <- "Necroptosis"
}

# --- 5) Workbook with 7 sheets ---
wb <- createWorkbook()
addWorksheet(wb, "Sh_1")  # genes
addWorksheet(wb, "Sh_2")  # full
addWorksheet(wb, "Sh_3")  # Cell cycle
addWorksheet(wb, "Sh_4")  # DNA repair
addWorksheet(wb, "Sh_5")  # Apoptosis
addWorksheet(wb, "Sh_6")  # Autophagy
addWorksheet(wb, "Sh_7")  # Necroptosis

writeData(wb, "Sh_1", genes_df)
writeData(wb, "Sh_2", full_df)
writeData(wb, "Sh_3", cc_df)
writeData(wb, "Sh_4", dr_df)
writeData(wb, "Sh_5", ap_df)
writeData(wb, "Sh_6", au_df)
writeData(wb, "Sh_7", ne_df)

# --- 6) Shade Sh_2 rows by category ---
styles <- list(
  "Cell cycle"   = createStyle(fgFill = "#FFF2CC"),  # yellow
  "DNA repair"   = createStyle(fgFill = "#D9EAD3"),  # green
  "Apoptosis"    = createStyle(fgFill = "#F4CCCC"),  # red
  "Autophagy"    = createStyle(fgFill = "#C6EFCE"),  # light green
  "Necroptosis"  = createStyle(fgFill = "#E6E6FA")   # lavender
)
if (nrow(full_df)) {
  for (i in seq_len(nrow(full_df))) {
    st <- styles[[ full_df$Category[i] ]]
    if (!is.null(st)) addStyle(wb, "Sh_2", st, rows = i + 1,
                               cols = seq_len(ncol(full_df)), gridExpand = TRUE)
  }
  freezePane(wb, "Sh_2", firstActiveRow = 2)
}

for (s in paste0("Sh_", 1:7)) setColWidths(wb, s, 1:50, "auto")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("→ Excel workbook written to", out_xlsx, "\n")
