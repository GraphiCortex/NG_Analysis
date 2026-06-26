#!/usr/bin/env Rscript
# debug_limma_results.R

res <- read.csv("Project/Data/Processed/limma_results/limma_E2F4_vs_E2F3.csv",
                row.names=1, stringsAsFactors=FALSE)

# 1. How many total probes?
cat("Total probes:\t", nrow(res), "\n")

# 2. Distribution of adjusted p-values
summary(res$adj.P.Val)
hist(res$adj.P.Val, breaks=50, main="adj.P.Val distribution", xlab="adj.P.Val")

# 3. Distribution of logFC
summary(res$logFC)
hist(res$logFC, breaks=50, main="logFC distribution", xlab="logFC")

# 4. Count how many meet each threshold separately
cat("FDR<0.05 only:\t", sum(res$adj.P.Val < 0.05), "\n")
cat("|logFC|>0.5 only:\t", sum(abs(res$logFC) > 0.5), "\n")
cat("Both criteria:\t", sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5), "\n")
