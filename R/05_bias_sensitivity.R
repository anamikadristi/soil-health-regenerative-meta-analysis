# 05_bias_sensitivity.R
# Publication bias and sensitivity checks for every pooled cell.
#
# 1. Multilevel Egger test: add sqrt(1/n_eff) as a moderator to the
#    three-level model, n_eff = n_t n_c / (n_t + n_c). Effective sample size
#    is used instead of the SE because, with lnRR, the SE is computed from
#    the means and is therefore correlated with the effect itself, and
#    because most SEs here are imputed (Nakagawa et al. 2022, Methods Ecol
#    Evol 13:4).
# 2. Time-lag test: publication year as a moderator (do early studies report
#    larger effects?).
# 3. Classic Egger regression test and trim-and-fill on study-level
#    averages (one effect per study, aggregate() with rho = 0.5), because
#    these methods assume independent effect sizes.
# 4. Leave-one-study-out: refit the three-level model without each study.
# 5. Variance strategies A to D from script 02.
# 6. Refit without extreme ratios (|lnRR| > ln 3, a 3-fold change).

source("R/00_setup.R")
log_line("== 05_bias_sensitivity ==")

es <- read_csv(file.path(PATHS$processed, "effect_sizes.csv"), show_col_types = FALSE) %>%
  mutate(practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS),
         vi = vi_A, n_eff = n_t * n_c / (n_t + n_c), inv_sqrt_neff = sqrt(1 / n_eff))
pooled <- read_csv(file.path(PATHS$tables, "pooled_effects.csv"), show_col_types = FALSE)

agg_study <- function(dat) {
  e <- escalc(measure = "ROM", yi = yi, vi = vi, data = as.data.frame(dat))
  aggregate(e, cluster = study_id, rho = 0.5, struct = "CS")
}

# REML with step halving; falls back to the Paule-Mandel estimator if REML
# does not converge (recorded in the output)
rma_safe <- function(ag) {
  f <- try(rma(yi, vi, data = ag, method = "REML", control = list(stepadj = 0.5, maxiter = 1000)), silent = TRUE)
  if (inherits(f, "try-error")) f <- rma(yi, vi, data = ag, method = "PM")
  f
}

# trim-and-fill refits the model internally; use PM if REML fails there
trimfill_safe <- function(fr, ag) {
  tf <- try(trimfill(fr), silent = TRUE)
  if (inherits(tf, "try-error")) tf <- trimfill(rma(yi, vi, data = ag, method = "PM"))
  tf
}

bias <- list(); loo_all <- list(); sens <- list(); aggs <- list()
for (i in seq_len(nrow(pooled))) {
  ind <- pooled$indicator[i]; pr <- pooled$practice[i]
  dat <- es %>% filter(indicator == ind, practice == pr)
  key <- paste(ind, pr, sep = " | ")

  # 1-2. multilevel Egger and time-lag
  fe <- fit_ml(dat, ~ inv_sqrt_neff)
  ft <- try(fit_ml(dat %>% filter(!is.na(pub_year)), ~ pub_year), silent = TRUE)
  # 3. study-level Egger and trim-and-fill
  ag <- agg_study(dat); aggs[[key]] <- ag
  fr <- rma_safe(ag)
  eg <- if (nrow(ag) >= 10) regtest(fr) else NULL
  tf <- trimfill_safe(fr, ag)
  bias[[key]] <- tibble(
    indicator = ind, practice = pr, k = nrow(dat), studies = nrow(ag),
    ml_egger_slope = fe$b[2], ml_egger_p = fe$pval[2],
    timelag_slope_per_yr = if (inherits(ft, "try-error")) NA else ft$b[2],
    timelag_p = if (inherits(ft, "try-error")) NA else ft$pval[2],
    study_level_tau2_method = fr$method, study_level_pct = pct(fr$b[1]),
    egger_z = if (is.null(eg)) NA else eg$zval, egger_p = if (is.null(eg)) NA else eg$pval,
    trimfill_method = tf$method, trimfill_missing = tf$k0, trimfill_side = tf$side,
    trimfill_pct = pct(tf$b[1]), trimfill_pct_lb = pct(tf$ci.lb), trimfill_pct_ub = pct(tf$ci.ub))

  # 4. leave-one-study-out
  base <- fit_ml(dat)
  loo <- map_dfr(unique(dat$study_id), function(s) {
    f <- try(fit_ml(dat %>% filter(study_id != s)), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    tibble(left_out = s, pct = pct(f$b[1]), pct_lb = pct(f$ci.lb), pct_ub = pct(f$ci.ub), p = f$pval)
  })
  loo_all[[key]] <- loo %>% mutate(indicator = ind, practice = pr)

  # 5-6. variance strategies and extreme ratios
  strat <- function(label, d2) {
    if (nrow(d2) < MIN_K || n_distinct(d2$study_id) < MIN_STUDIES)
      return(tibble(analysis = label, k = nrow(d2), studies = n_distinct(d2$study_id),
                    pct = NA, pct_lb = NA, pct_ub = NA, p = NA, I2_total = NA))
    f <- fit_ml(d2)
    tibble(analysis = label, k = f$k, studies = n_distinct(d2$study_id), pct = pct(f$b[1]),
           pct_lb = pct(f$ci.lb), pct_ub = pct(f$ci.ub), p = f$pval, I2_total = i2_ml(f)[1])
  }
  sens[[key]] <- bind_rows(
    strat("A: reported SD, else mean-CV imputation (primary)", dat),
    strat("B: mean-CV imputation for all rows", dat %>% mutate(vi = vi_B)),
    strat("C: complete cases (reported SD only)", dat %>% filter(!is.na(vi_C)) %>% mutate(vi = vi_C)),
    strat("D: imputed CV doubled (variance x4)", dat %>% mutate(vi = vi_D)),
    strat("E: without |lnRR| > ln(3)", dat %>% filter(!extreme_flag)),
    strat("F: only rows with reported n", dat %>% filter(n_source == "reported"))
  ) %>% mutate(indicator = ind, practice = pr, .before = 1)

  loo_summary <- tibble(indicator = ind, practice = pr,
                        loo_min_pct = min(loo$pct), loo_max_pct = max(loo$pct),
                        loo_n_nonsig = sum(loo$p >= 0.05), loo_n_sig = sum(loo$p < 0.05),
                        base_p = base$pval)
  bias[[key]] <- bias[[key]] %>% left_join(loo_summary, by = c("indicator", "practice"))
}

bias_tab <- bind_rows(bias)
write_csv(bias_tab %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(PATHS$tables, "publication_bias.csv"))
write_csv(bind_rows(loo_all) %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(PATHS$tables, "leave_one_study_out.csv"))
sens_tab <- bind_rows(sens)
write_csv(sens_tab %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(PATHS$tables, "sensitivity_variance_and_outliers.csv"))
saveRDS(aggs, file.path(PATHS$processed, "study_level_aggregates.rds"))

print(as.data.frame(bias_tab %>% select(indicator, practice, studies, ml_egger_p, timelag_p, egger_p,
                                        trimfill_missing, trimfill_pct, loo_min_pct, loo_max_pct, loo_n_nonsig) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))))
print(as.data.frame(sens_tab %>% mutate(across(where(is.double), ~ signif(.x, 3)))))

# Funnel plots (study-level averages, with trim-and-fill imputed studies as open points)
main <- c("SOC concentration | Cover crop", "SOC concentration | Organic management",
          "SOC stock | Cover crop", "Bulk density | Cover crop",
          "Microbial biomass C | Cover crop", "Aggregate stability | Cover crop")
png(file.path(PATHS$figures, "fig05_funnel_plots.png"), width = 10, height = 6.5, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
for (m in main) {
  tf <- trimfill_safe(rma_safe(aggs[[m]]), aggs[[m]])
  funnel(tf, xlab = "lnRR (study average)", legend = FALSE, cex = 0.8)
  title(main = sub(" \\| ", "\n", m), cex.main = 0.95, font.main = 1)
}
dev.off()
log_line("Funnel plots written")
