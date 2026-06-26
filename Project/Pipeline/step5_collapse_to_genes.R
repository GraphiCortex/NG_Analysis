#!/usr/bin/env Rscript

# ==============================================
# Step 5: Collapse significant probes → genes
# ==============================================

library(readr)
library(dplyr)

# 1. Read in your significant‐probe list (with SPOT_ID column)
sig_probes <- read_csv(
  "Project/New_Data/E2F4_vs_E2F3_significant_probes.csv",
  col_types = cols()
)

# 2. Read the promoter mapping (probe → gene)
prom_map <- read_csv(
  "Project/New_Data/GPL_promoter_probe2gene.csv",
  col_types = cols()
)

# 3. Join probes to their gene symbols
sig_map <- sig_probes %>%
  left_join(prom_map, by = "SPOT_ID") %>%
  filter(!is.na(GENE_SYMBOL))

# 4. Write out the probe→gene table
write_csv(
  sig_map,
  "Project/New_Data/E2F4_vs_E2F3_sigProbe2Gene.csv"
)

# 5. Extract & write the unique gene list
sig_genes <- sig_map %>%
  distinct(GENE_SYMBOL) %>%
  arrange(GENE_SYMBOL)

write_csv(
  sig_genes,
  "Project/New_Data/E2F4_vs_E2F3_significant_genes.csv"
)

message("✔ Step 5 complete: collapsed ",
        nrow(sig_map), " probes into ",
        nrow(sig_genes), " unique genes.")
