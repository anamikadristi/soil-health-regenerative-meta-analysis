# 03_random_effects.R
# Pooled effect of each practice on each indicator.
#
# One study (paper) often contributes many effect sizes: several sites,
# years, depths and cover crop species, often against the same control.
# These are not independent. The model therefore has three random effects:
#   study_id / exp_id / es_id
#   - between studies (papers)
#   - between experiments (sites or trials) within a study
#   - between effect sizes within an experiment (depth, year, species ...)
# Inference uses t-distributions with "contain" degrees of freedom, and a
# cluster-robust (CR2) confidence interval by study is reported as a check
# against the remaining dependence among sampling errors (shared controls).

source("R/00_setup.R")
log_line("== 03_random_effects ==")

es <- read_csv(file.path(PATHS$processed, "effect_sizes.csv"), show_col_types = FALSE) %>%
  mutate(practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS),
         vi = vi_A)

cells <- es %>% group_by(indicator, practice) %>%
  summarise(k = n(), studies = n_distinct(study_id), .groups = "drop") %>%
  filter(k >= MIN_K, studies >= MIN_STUDIES) %>% arrange(indicator, practice)

fits <- list()
res <- pmap_dfr(cells %>% select(indicator, practice), function(indicator, practice) {
  dat <- es %>% filter(indicator == !!indicator, practice == !!practice)
  fit <- fit_ml(dat)
  key <- paste(indicator, practice, sep = " | ")
  fits[[key]] <<- fit
  pr <- predict(fit)
  rb <- robust_ci(fit, dat$study_id)
  i2 <- i2_ml(fit)
  tibble(indicator = as.character(indicator), practice = as.character(practice),
         k = fit$k, studies = n_distinct(dat$study_id), experiments = n_distinct(dat$exp_id),
         pct_sd_reported = round(100 * mean(dat$sd_source != "imputed"), 1),
         lnRR = fit$b[1], se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p = fit$pval,
         pct = pct(fit$b[1]), pct_lb = pct(fit$ci.lb), pct_ub = pct(fit$ci.ub),
         pi_pct_lb = pct(pr$pi.lb), pi_pct_ub = pct(pr$pi.ub),
         robust_pct_lb = pct(rb["lo"]), robust_pct_ub = pct(rb["hi"]), robust_p = rb["p"],
         sigma2_study = fit$sigma2[1], sigma2_experiment = fit$sigma2[2], sigma2_within = fit$sigma2[3],
         tau2_total = sum(fit$sigma2),
         I2_total = i2["I2_total"], I2_study = i2["I2_study"], I2_experiment = i2["I2_experiment"],
         I2_within = i2["I2_within"], Q = fit$QE, Q_p = fit$QEp)
})

res <- res %>% mutate(exploratory = indicator == "Infiltration rate" | studies < 8)
write_csv(res %>% mutate(across(where(is.double), ~ signif(.x, 4))),
          file.path(PATHS$tables, "pooled_effects.csv"))
saveRDS(fits, file.path(PATHS$processed, "fits_pooled.rds"))

show <- res %>% transmute(indicator, practice, k, studies,
                          `% change [95% CI]` = sprintf("%.1f [%.1f, %.1f]", pct, pct_lb, pct_ub),
                          `robust CI` = sprintf("[%.1f, %.1f]", robust_pct_lb, robust_pct_ub),
                          p = signif(p, 2), I2 = round(I2_total, 1), tau2 = signif(tau2_total, 3))
print(as.data.frame(show))
capture.output(print(as.data.frame(show)), file = file.path(PATHS$logs, "pooled_effects_print.txt"))

# Full model output for the main cell, for teaching
sink(file.path(PATHS$logs, "model_SOC_covercrop_full.txt"))
print(fits[["SOC concentration | Cover crop"]])
sink()
log_line("Cells pooled: ", nrow(res))
