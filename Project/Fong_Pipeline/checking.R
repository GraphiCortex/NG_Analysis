full  <- read.csv("Project/New_Data/Fong2022_DE_TREAT025_full.csv")
sig   <- read.csv("Project/New_Data/Fong2022_DE_TREAT025_significant.csv")
nrow(full); nrow(sig)           # expect 30195 vs 6036
setdiff(names(full), names(sig))# should be empty (same columns)
