#!/usr/bin/env Rscript

# ==============================================
# Pipeline_2 — Step 5: Collapse significant probes → genes (non-promoter)
# (Exact mirror of old Step 5, just using New_Data_2 + non-promoter map)
# ==============================================

library(readr)
library(dplyr)

# 1) Significant non-promoter probes (with SPOT_ID)
sig_probes <- read_csv(
  "Project/New_Data_2/E2F4_vs_E2F3_significant_probes_nonpromoter.csv",
  col_types = cols()
)

# 2) Non-promoter mapping (probe → gene)
np_map <- read_csv(
  "Project/New_Data_2/GPL_nonpromoter_probe2gene.csv",
  col_types = cols()
)

# 3) Join probes to gene symbols
sig_map <- sig_probes %>%
  left_join(np_map, by = "SPOT_ID") %>%
  filter(!is.na(GENE_SYMBOL))

# 4) Write probe→gene table
write_csv(
  sig_map,
  "Project/New_Data_2/E2F4_vs_E2F3_sigProbe2Gene_nonpromoter.csv"
)

# 5) Unique gene list
sig_genes <- sig_map %>%
  distinct(GENE_SYMBOL) %>%
  arrange(GENE_SYMBOL)

write_csv(
  sig_genes,
  "Project/New_Data_2/E2F4_vs_E2F3_significant_genes_nonpromoter.csv"
)

message("✔ Step 5 (non-promoter) complete: collapsed ",
        nrow(sig_map), " probes into ",
        nrow(sig_genes), " unique genes.")
