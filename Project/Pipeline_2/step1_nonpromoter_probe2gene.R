#!/usr/bin/env Rscript
# Step 1 (Pipeline_2): Probe → Gene mapping (non-promoter: DOWNSTREAM + INSIDE)

gpl <- read.delim("Project/MetaData/GPL18280-19916.txt",
                  sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)
stopifnot(all(c("SPOT_ID","GENE_SYMBOL","DESCRIPTION") %in% names(gpl)))

# Use only the non-promoter labels available on this GPL
use_desc <- intersect(unique(gpl$DESCRIPTION), c("DOWNSTREAM","INSIDE"))
if (length(use_desc) == 0L) stop("No non-promoter labels (DOWNSTREAM/INSIDE) found.")

message("Using non-promoter DESCRIPTION levels: ", paste(sort(use_desc), collapse=", "))

nonpromoter_probes <- subset(gpl, DESCRIPTION %in% use_desc & GENE_SYMBOL != "")
probe2gene <- nonpromoter_probes[, c("SPOT_ID","GENE_SYMBOL","DESCRIPTION")]

out_dir <- "Project/New_Data_2"; dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(probe2gene, file.path(out_dir, "GPL_nonpromoter_probe2gene.csv"),
          row.names=FALSE, quote=FALSE)

summ <- as.data.frame(table(nonpromoter_probes$DESCRIPTION)); names(summ) <- c("DESCRIPTION","n_probes")
write.csv(summ, file.path(out_dir, "nonpromoter_label_counts.csv"), row.names=FALSE, quote=FALSE)

cat("→ Wrote ", nrow(probe2gene), " non-promoter mappings (", paste(sort(use_desc), collapse=", "),
    ") to ", file.path(out_dir, "GPL_nonpromoter_probe2gene.csv"), "\n", sep="")

