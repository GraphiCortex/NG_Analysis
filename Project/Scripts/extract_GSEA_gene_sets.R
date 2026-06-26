#!/usr/bin/env Rscript
# extract_GSEA_gene_sets.R
# Pull out mitotic, DNA-repair & apoptosis gene sets from GSEA + ORA results

# 1) Load libs
if (!requireNamespace("readr", quietly=TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
suppressMessages({
  library(readr)
  library(dplyr)
  library(stringr)    # ← add this line
})


# 2) Read in the preranked GSEA table
gsea <- read_csv("Project/Data/Processed/limma_results/E2F4_vs_E2F3_GSEA_GO_BP.csv",
                 col_types = cols(.default="c"))  # read all as character to preserve lists

# 3) Identify the column with the slash-separated gene lists
#    (clusterProfiler calls this "core_enrichment")
stopifnot("core_enrichment" %in% colnames(gsea))

# 4) Extract mitotic / cell-cycle gene set
mitotic_terms <- gsea %>%
  filter(str_detect(Description, regex("mitotic|cell cycle|chromatid|chromosome segregation", ignore_case=TRUE)))
mitotic_genes <- mitotic_terms$core_enrichment %>%
  str_split("/") %>% unlist() %>% unique() %>% sort()
mitotic_df <- tibble(Gene=mitotic_genes)

# 5) Extract DNA-repair / recombination gene set
repair_terms <- gsea %>%
  filter(str_detect(Description, regex("DNA repair|recombination|DNA metabolic", ignore_case=TRUE)))
repair_genes <- repair_terms$core_enrichment %>%
  str_split("/") %>% unlist() %>% unique() %>% sort()
repair_df <- tibble(Gene=repair_genes)

# 6) Read in the apoptosis ORA results we already generated
ora <- read_csv("Project/Data/Processed/limma_results/E2F4_cell_death_ORA.csv",
                col_types = cols(.default="c"))
apoptosis_genes <- ora$geneID %>%
  str_split("/") %>% unlist() %>% unique() %>% sort()
apoptosis_df <- tibble(Gene=apoptosis_genes)

# 7) Write each out for downstream merging
outdir <- "Project/Data/Processed/gene_sets"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write_csv(mitotic_df,   file.path(outdir, "E2F4_mitotic_genes.csv"))
write_csv(repair_df,    file.path(outdir, "E2F4_DNArepair_genes.csv"))
write_csv(apoptosis_df, file.path(outdir, "E2F4_apoptosis_genes.csv"))

cat("Gene sets extracted:\n",
    " • Mitotic / Cell-cycle: ", length(mitotic_genes), "genes\n",
    " • DNA repair / recombination: ", length(repair_genes), "genes\n",
    " • Apoptosis / cell-death: ",  length(apoptosis_genes), "genes\n",
    "CSV files in ", outdir, "\n", sep="")
