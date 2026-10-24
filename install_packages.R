# Install the CRAN packages used by the R scripts (tested with R 4.5.2).
# orchaRd is not on CRAN and is not needed; orchard plots are drawn with ggplot2.
pkgs <- c("metafor", "clubSandwich", "readxl", "dplyr", "tidyr", "stringr", "readr",
          "purrr", "ggplot2", "maps")
need <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))
cat("All packages available:", paste(pkgs, collapse = ", "), "\n")
