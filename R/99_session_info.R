# 99_session_info.R: record R and package versions used for this run
source("R/00_setup.R")
library(maps)
writeLines(c(capture.output(sessionInfo()), "",
             paste("metafor", packageVersion("metafor")),
             paste("clubSandwich", packageVersion("clubSandwich"))),
           "session_info.txt")
cat("Wrote session_info.txt\n")
