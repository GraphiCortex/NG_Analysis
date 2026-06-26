#!/usr/bin/env Rscript
# extract_limma_significant.R
# Pull out all significant probes plus top-50 summaries

# 1. Load
contrast <- "E2F4_vs_E2F3"       # change to the contrast you want
infile   <- file.path("Project","Data","Processed","limma_results",
                      paste0("limma_", contrast, ".csv"))
res      <- read.csv(infile, row.names=1, stringsAsFactors=FALSE)

# 2. Filter on FDR & fold-change
sig_all <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 0.5)

# 3. Write all significant
outdir <- file.path("Project","Data","Processed","limma_results")
write.csv(sig_all,
          file=file.path(outdir, paste0(contrast, "_significant.csv")),
          quote=FALSE)

# 4. Also write Top 50 up & Top 50 down
sig_up   <- head(sig_all[order(-sig_all$logFC), ], 50)
sig_down <- head(sig_all[order( sig_all$logFC), ], 50)

write.csv(sig_up,
          file=file.path(outdir, paste0(contrast, "_top50_up.csv")),
          quote=FALSE)
write.csv(sig_down,
          file=file.path(outdir, paste0(contrast, "_top50_down.csv")),
          quote=FALSE)

cat("Wrote", nrow(sig_all), "significant probes.\n",
    "Top 50 up/down in:", outdir, "\n")
