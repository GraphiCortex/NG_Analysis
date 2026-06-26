library(readxl)
library(dplyr)
library(tidyr)
library(janitor)

load_fong_sheet <- function(path, sheet = 1) {
  x <- read_excel(path, sheet = sheet, .name_repair = "minimal")
  
  # If first row looks like headers, promote it
  looks_like_header <- function(v) {
    # many non-NA, mostly character, not all numeric
    ch <- sum(!is.na(v) & grepl("[A-Za-z]", as.character(v)))
    ch >= max(3, floor(length(v) * 0.4))
  }
  if (looks_like_header(colnames(x)) == FALSE && looks_like_header(x[1, ])) {
    nm <- as.character(unlist(x[1, ]))
    colnames(x) <- nm
    x <- x[-1, , drop = FALSE]
  }
  
  # Clean names, drop completely empty rows/cols
  x <- x %>%
    remove_empty(which = c("rows", "cols")) %>%
    clean_names()

  # Trim whitespace in the first column (often gene/probe id)
  if (ncol(x) > 0) x[[1]] <- trimws(as.character(x[[1]]))

  # Heuristic: expression matrix if there are ≥6 numeric columns
  num_cols <- sapply(x, function(col) suppressWarnings(is.numeric(as.numeric(col))))
  numeric_count <- sum(num_cols, na.rm = TRUE)

  kind <- if (numeric_count >= 6) "expression_matrix" else "deg_table"
  message(sprintf("→ Detected: %s (%d numeric columns)", kind, numeric_count))

  # Quick peek
  print(head(x, 5))
  return(invisible(list(data = x, kind = kind)))
}

p_raw  <- "Project/Data/Raw/microarray/Fong2022/Fong2022_raw/FDSC Batch IHW.xlsx"
p_sort <- "Project/Data/Raw/microarray/Fong2022/Fong2022_raw/Sorted.xlsx"

a1 <- load_fong_sheet(p_raw,  sheet = 1)
a2 <- load_fong_sheet(p_sort, sheet = 1)
