# 1) Load the GPL file and peek at its structure
gpl <- read.delim("Project/MetaData/GPL18280-19916.txt",
                  stringsAsFactors=FALSE,
                  comment.char="#")
cat("Columns in GPL:\n"); print(colnames(gpl))

# 2) Show the first few rows
cat("\nFirst 6 rows of GPL:\n")
print(head(gpl[, 1:6]))

# 3) What values appear in DESCRIPTION?
cat("\nUnique DESCRIPTION values:\n")
print(unique(gpl$DESCRIPTION))

# 4) How many probes have each DESCRIPTION?
cat("\nCounts per DESCRIPTION:\n")
print(table(gpl$DESCRIPTION))
