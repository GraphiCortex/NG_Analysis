#!/usr/bin/env Rscript
#
# run_preranked_GSEA.R  (updated to use dotplot instead of barplot)

# 0. Install & load required packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("clusterProfiler","org.Mm.eg.db","dplyr","ggplot2")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    BiocManager::install(pkg)
  }
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
})

# 1. Read in limma results
contrast   <- "E2F4_vs_E2F3"
limma_file <- file.path("Project","Data","Processed","limma_results",
                        paste0("limma_", contrast, ".csv"))
res <- read.csv(limma_file, row.names=1, stringsAsFactors=FALSE)

# 2. Read in probe→gene mapping
map_file <- file.path("Project","MetaData","GPL_probe2gene.csv")
annot    <- read.csv(map_file, stringsAsFactors=FALSE)

# 3. Merge & collapse probes to genes (max |logFC|)
df <- res %>%
  mutate(ID_REF = rownames(res)) %>%
  select(ID_REF, logFC) %>%
  inner_join(annot, by="ID_REF") %>%
  group_by(SYMBOL) %>%
  summarize(logFC = logFC[which.max(abs(logFC))]) %>%
  ungroup()

# 4. Build a named, sorted vector
geneList       <- df$logFC
names(geneList) <- df$SYMBOL
geneList       <- sort(geneList, decreasing=TRUE)

# 5. Run preranked GSEA on GO Biological Process
set.seed(2025)
gseaRes <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 20,
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  eps          = 0,
  verbose      = FALSE
)

# 6. Write out full GSEA table
outdir <- file.path("Project","Data","Processed","limma_results")
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
write.csv(
  as.data.frame(gseaRes),
  file = file.path(outdir, paste0(contrast,"_GSEA_GO_BP.csv")),
  row.names = FALSE,
  quote = FALSE
)
# 7. Create a dotplot of the GSEA results
library(ggplot2)
library(stringr)

# 7b. Clean dotplot with wrapped labels
p <- dotplot(
  gseaRes,
  showCategory = 20,
  title        = paste0(contrast, " GSEA: GO Biological Process")
) +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust=0.5),
    axis.text.y      = element_text(size=9),
    axis.title.y     = element_blank(),
    plot.margin      = margin(t=5, r=5, b=5, l=20, unit="pt")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25))

# Save it
ggsave(
  filename = file.path(outdir, paste0(contrast, "_GSEA_dotplot_wrapped.png")),
  plot     = p,
  width    = 10,
  height   = 8
)

message("🏁 GSEA complete. Results in: ", outdir)
