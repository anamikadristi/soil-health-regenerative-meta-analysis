# 04_moderators.R
# Multilevel meta-regression: does the effect depend on climate, soil
# texture, years under the practice, sampling depth or land use?
#
# Each moderator is tested on its own (one moderator per model) because
# moderators are unevenly reported and partly confounded (e.g., climate and
# country). For a categorical moderator, levels with fewer than 5 effect sizes
# or 3 studies are dropped before fitting. The omnibus test (QM, F-test) asks
# whether the moderator explains any heterogeneity. Pseudo-R2 is the
# proportional drop in total tau2 relative to an intercept-only model fitted
# to the same rows. Holm-adjusted p-values are given because many tests are
# run; treat moderator results as exploratory.

source("R/00_setup.R")
log_line("== 04_moderators ==")

es <- read_csv(file.path(PATHS$processed, "effect_sizes.csv"), show_col_types = FALSE) %>%
  mutate(vi = vi_A, map100 = map_mm / 100,
         years_class = factor(years_class, c("1-3 yr", "4-10 yr", "> 10 yr")),
         depth_class = factor(depth_class, c("Surface (bottom <= 15 cm)", "0-30 cm layers",
                                             "Deeper (> 30 cm or starts >= 15 cm)")),
         texture = factor(texture, c("Coarse", "Medium", "Fine")))

cells <- es %>% group_by(indicator, practice) %>%
  summarise(k = n(), studies = n_distinct(study_id), .groups = "drop") %>%
  filter(k >= 50, studies >= 15)

mods_cat  <- c(climate_zone = "Koppen climate zone", texture = "Soil texture",
               years_class = "Years under practice (class)", depth_class = "Sampling depth (class)",
               land_use = "Land use", tillage_type = "No-till vs reduced tillage")
mods_cont <- c(map100 = "Mean annual precipitation (per 100 mm)", mat_c = "Mean annual temperature (per 1 C)",
               years = "Years under practice (per year)", depth_mid = "Depth midpoint (per cm)")

run_mod <- function(dat, mod, label, indicator, practice) {
  dat <- dat %>% filter(!is.na(.data[[mod]]))
  is_cat <- mod %in% names(mods_cat)
  if (is_cat) {
    keep <- dat %>% group_by(.data[[mod]]) %>%
      summarise(k = n(), s = n_distinct(study_id), .groups = "drop") %>% filter(k >= 5, s >= 3)
    dat <- dat %>% filter(.data[[mod]] %in% keep[[mod]]) %>% mutate(m = droplevels(factor(.data[[mod]])))
    if (nlevels(dat$m) < 2) return(NULL)
  } else {
    if (n_distinct(dat$study_id) < 10) return(NULL)
    dat$m <- dat[[mod]]
  }
  f0 <- try(fit_ml(dat), silent = TRUE)
  f1 <- try(fit_ml(dat, ~ m), silent = TRUE)
  if (inherits(f0, "try-error") || inherits(f1, "try-error")) return(NULL)
  r2 <- max(0, (sum(f0$sigma2) - sum(f1$sigma2)) / sum(f0$sigma2)) * 100
  base <- tibble(indicator, practice, moderator = label, k = f1$k, studies = n_distinct(dat$study_id),
                 QM = f1$QM, QM_df1 = f1$QMdf[1], QM_df2 = f1$QMdf[2], QM_p = f1$QMp, pseudo_R2 = r2)
  if (is_cat) {
    f2 <- fit_ml(dat, ~ m - 1)
    lv <- dat %>% group_by(m) %>% summarise(k_level = n(), studies_level = n_distinct(study_id), .groups = "drop")
    out <- tibble(level = levels(dat$m), lnRR = f2$b[, 1], ci_lb = f2$ci.lb, ci_ub = f2$ci.ub) %>%
      left_join(lv %>% mutate(level = as.character(m)) %>% select(-m), by = "level")
  } else {
    out <- tibble(level = "slope", lnRR = f1$b[2], ci_lb = f1$ci.lb[2], ci_ub = f1$ci.ub[2],
                  k_level = f1$k, studies_level = n_distinct(dat$study_id))
  }
  bind_cols(base[rep(1, nrow(out)), ], out)
}

res <- pmap_dfr(cells %>% select(indicator, practice), function(indicator, practice) {
  dat <- es %>% filter(indicator == !!indicator, practice == !!practice)
  map_dfr(names(c(mods_cat, mods_cont)), function(m) {
    if (m == "tillage_type" && practice != "No-till / reduced tillage") return(NULL)
    run_mod(dat, m, c(mods_cat, mods_cont)[[m]], indicator, practice)
  })
})

# Holm adjustment over the omnibus tests (one per indicator x practice x moderator)
omni <- res %>% distinct(indicator, practice, moderator, QM_p) %>%
  mutate(QM_p_holm = p.adjust(QM_p, "holm"))
res <- res %>% left_join(omni, by = c("indicator", "practice", "moderator", "QM_p")) %>%
  mutate(pct = ifelse(level == "slope", NA, pct(lnRR)), pct_lb = ifelse(level == "slope", NA, pct(ci_lb)),
         pct_ub = ifelse(level == "slope", NA, pct(ci_ub)),
         slope_pct_per_unit = ifelse(level == "slope", pct(lnRR), NA))
write_csv(res %>% mutate(across(where(is.double), ~ signif(.x, 4))),
          file.path(PATHS$tables, "moderator_results.csv"))

omni_tab <- res %>% distinct(indicator, practice, moderator, k, studies, QM, QM_df1, QM_df2, QM_p, QM_p_holm, pseudo_R2)
write_csv(omni_tab %>% mutate(across(where(is.double), ~ signif(.x, 4))),
          file.path(PATHS$tables, "moderator_omnibus_tests.csv"))
print(as.data.frame(omni_tab %>% mutate(across(where(is.double), ~ signif(.x, 3)))))

# Do practices differ for SOC concentration? (all pooled practices, one model)
soc <- es %>% filter(indicator == "SOC concentration",
                     practice %in% c("Cover crop", "No-till / reduced tillage", "Organic management",
                                     "Residue retention", "Agroforestry")) %>% mutate(practice = droplevels(factor(practice, PRACTICE_LEVELS)))
fp <- fit_ml(soc, ~ practice)
log_line("SOC concentration, practice as moderator: F(", fp$QMdf[1], ",", round(fp$QMdf[2], 1), ") = ",
         round(fp$QM, 2), ", p = ", signif(fp$QMp, 3), ", k = ", fp$k)
capture.output(print(fp), file = file.path(PATHS$logs, "model_SOC_practice_moderator.txt"))
