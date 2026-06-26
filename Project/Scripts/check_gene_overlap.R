#!/usr/bin/env Rscript
# check_gene_overlap.R

library(readr)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Load expression SYMBOL list
expr2 <- read_csv("Project/Data/Processed/expression_expr2.csv")  # this is your expr2 from before
genes_expr <- expr2$Gene

# Function to report overlap
report_overlap <- function(path) {
  gs <- read_csv(path, col_types=cols(Gene="c"))$Gene
  overlap     <- sum(gs %in% genes_expr)
  total       <- length(gs)
  missing     <- setdiff(gs, genes_expr)
  cat(
    basename(path),"\n",
    " • Total genes:     ", total, "\n",
    " • Overlap:         ", overlap, "\n",
    " • Missing names:   ", length(missing), "\n",
    " • Examples missing:", paste(head(missing, 10), collapse=", "),"\n\n"
  )
}

# Run for each set
paths <- list.files("Project/Data/Processed/gene_sets", pattern="\\.csv$", full.names=TRUE)
for(p in paths) report_overlap(p)
