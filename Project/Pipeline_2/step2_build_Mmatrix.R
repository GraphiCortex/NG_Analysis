#!/usr/bin/env Rscript
# Pipeline_2 (optional): rebuild M-value matrix from processed files

suppressPackageStartupMessages({library(readr); library(dplyr); library(stringr); library(purrr)})

proc_dir <- "Project/Data/Raw/ChipSeq/Julian2016/Processed/"
files <- list.files(proc_dir, pattern = "_Processed\\.txt$", full.names = TRUE)

m_list <- map(files, function(fp){
  samp <- basename(fp) %>% str_remove("_Processed\\.txt$")
  read_delim(fp, delim="\t",
    col_types = cols(ID_REF=col_character(), IP=col_double(), INPUT=col_double())
  ) %>%
    mutate(M = log2(IP/INPUT)) %>%
    rename(SPOT_ID = ID_REF) %>%
    transmute(SPOT_ID, !!samp := M)
})

Mmat <- reduce(m_list, left_join, by="SPOT_ID")
dir.create("Project/New_Data_2", recursive=TRUE, showWarnings=FALSE)
write_csv(Mmat, "Project/New_Data_2/Mvalue_matrix.csv")
message("✔ Wrote M-value matrix with ", nrow(Mmat), " probes × ", ncol(Mmat)-1, " samples to Project/New_Data_2/Mvalue_matrix.csv")
