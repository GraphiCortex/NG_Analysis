#!/usr/bin/env Rscript

# ==============================================
# Step 2: Build M-value matrix from processed files
# ==============================================

install.packages("tidyverse")
library(tidyverse)


# 1. Locate all processed ChIP–chip output files
proc_dir <- "Project/Data/Raw/ChipSeq/Julian2016/Processed/"
files <- list.files(proc_dir, pattern = "_Processed\\.txt$", full.names = TRUE)

# 2. Read + compute M for each sample
m_list <- map(files, function(fp) {
  # derive a clean sample name from the filename
  samp <- basename(fp) %>% str_remove("_Processed\\.txt$")
  
  df <- read_delim(fp,
                   delim = "\t",
                   col_types = cols(
                     ID_REF = col_character(),
                     IP    = col_double(),
                     INPUT = col_double()
                   ))
  
  df %>%
    # compute M-value
    mutate(M = log2(IP / INPUT)) %>%
    # rename probe ID for downstream joins
    rename(SPOT_ID = ID_REF) %>%
    # keep only SPOT_ID + this sample’s M
    select(SPOT_ID, M) %>%
    # name the M column by sample
    rename(!!samp := M)
})

# 3. Merge all samples into one matrix
Mmat <- reduce(m_list, left_join, by = "SPOT_ID")

# 4. Write out the combined M-value matrix
write_csv(Mmat, "Project/New_Data/Mvalue_matrix.csv")

message("✔ Wrote M-value matrix with ",
        nrow(Mmat), " probes × ", ncol(Mmat) - 1,
        " samples to Pipeline/Mvalue_matrix.csv")
