# 00_setup.R
# Shared paths, packages and helper functions. Every numbered script sources
# this file first. Run all scripts from the project root, for example:
#   Rscript R/01_load_clean.R

suppressPackageStartupMessages({
  library(maps)      # loaded first so purrr::map is not masked
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(metafor)
  library(clubSandwich)
})

if (!file.exists("R/00_setup.R")) {
  stop("Run the scripts from the project root (the folder that contains R/ and data_raw/).")
}

PATHS <- list(
  raw_xlsx  = "data_raw/SoilHealthDB_V1.xlsx",
  koppen    = "data_raw/Koeppen-Geiger-ASCII.txt",
  processed = "data_processed",
  tables    = "outputs/tables",
  figures   = "outputs/figures",
  logs      = "outputs/logs"
)
for (p in PATHS[c("processed", "tables", "figures", "logs")]) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# Minimum data for a practice x indicator cell to be pooled
MIN_K       <- 10   # effect sizes
MIN_STUDIES <- 5    # independent studies (papers)

# Order and labels used in tables and figures
PRACTICE_LEVELS  <- c("Cover crop", "No-till / reduced tillage", "Organic management",
                      "Residue retention", "Agroforestry", "Crop rotation")
INDICATOR_LEVELS <- c("SOC concentration", "SOC stock", "Bulk density",
                      "Aggregate stability", "Microbial biomass C", "Infiltration rate")

# Percent change from a log response ratio
pct <- function(x) (exp(x) - 1) * 100

# Parse the first number in a free-text field ("3", "3-4", "about 5 yr")
first_num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

# Three-level random-effects model: effect sizes nested in experiments
# (sites/trials) nested in studies (papers). REML, t-based inference.
fit_ml <- function(dat, mods = ~ 1) {
  # mods is a one-sided formula, e.g. ~ texture (default ~ 1 = intercept only).
  # If the default optimizer (nlminb) fails to converge, two other
  # optimizers are tried before giving up.
  go <- function(ctrl) rma.mv(yi, vi, mods = mods, data = dat,
                              random = ~ 1 | study_id / exp_id / es_id,
                              method = "REML", test = "t", dfs = "contain", sparse = TRUE,
                              control = ctrl)
  f <- try(go(list()), silent = TRUE)
  if (inherits(f, "try-error")) f <- try(go(list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = 5000)), silent = TRUE)
  if (inherits(f, "try-error")) f <- go(list(optimizer = "bobyqa"))
  f
}

# Multilevel I2 (Nakagawa & Santos 2012) using the "typical" sampling variance
i2_ml <- function(fit) {
  W <- diag(1 / fit$vi)
  X <- model.matrix(fit)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (fit$k - fit$p) / sum(diag(P))
  tot <- sum(fit$sigma2) + typical_v
  c(I2_total = 100 * sum(fit$sigma2) / tot,
    I2_study = 100 * fit$sigma2[1] / tot,
    I2_experiment = 100 * fit$sigma2[2] / tot,
    I2_within = 100 * fit$sigma2[3] / tot)
}

# Cluster-robust (CR2, small-sample corrected) test of the intercept, by study
robust_ci <- function(fit, cluster) {
  r <- robust(fit, cluster = cluster, clubSandwich = TRUE)
  c(est = unname(r$beta[1]), lo = r$ci.lb[1], hi = r$ci.ub[1], p = r$pval[1])
}

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92", colour = NA)))

save_png <- function(plot, file, w = 7, h = 5) {
  ggsave(file.path(PATHS$figures, file), plot, width = w, height = h, dpi = 300)
}

log_line <- function(..., file = "run_log.txt") {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(PATHS$logs, file), append = TRUE)
}
