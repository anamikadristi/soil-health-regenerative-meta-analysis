# 02_effect_sizes.R
# Log response ratio, lnRR = ln(mean_treatment / mean_control), and its
# sampling variance.
#
# Most rows in SoilHealthDB have no SD (for cover crop x SOC, about 9% do).
# Dropping them would throw away most of the evidence and could bias the
# result toward studies that report more statistics. Instead, missing
# variances are imputed from the average coefficient of variation (CV) of
# the rows that do report SDs, following Nakagawa et al. (2023, Methods in
# Ecology and Evolution 14:2072) and Lajeunesse (2013):
#
#   v_i = CVbar_t^2 / n_t + CVbar_c^2 / n_c
#
# which is the delta-method variance of lnRR with the CVs replaced by their
# average. Four variance versions are stored so script 05 can test whether
# the conclusions depend on this choice:
#   vi_A  primary: reported SD where available, imputed CV otherwise ("missing-cases")
#   vi_B  imputed CV for every row ("all-cases"; ignores reported SDs)
#   vi_C  reported SD only (complete-case; rows without SD get NA)
#   vi_D  as vi_A, but the imputed CV is doubled (variance x 4) for imputed rows

source("R/00_setup.R")
log_line("== 02_effect_sizes ==")

d <- read_csv(file.path(PATHS$processed, "soil_health_long.csv"), show_col_types = FALSE)

# Clean SDs: non-positive values are treated as missing. If only one arm has
# an SD, that arm's CV is used for both arms.
d <- d %>% mutate(
  sd_c = ifelse(!is.na(sd_c) & sd_c > 0, sd_c, NA),
  sd_t = ifelse(!is.na(sd_t) & sd_t > 0, sd_t, NA),
  cv_c = sd_c / m_c, cv_t = sd_t / m_t,
  sd_source = case_when(!is.na(cv_c) & !is.na(cv_t) ~ "reported",
                        !is.na(cv_c) | !is.na(cv_t) ~ "one arm reported",
                        TRUE ~ "imputed"),
  cv_c = coalesce(cv_c, cv_t), cv_t = coalesce(cv_t, cv_c))

# lnRR and variance from reported SDs (metafor::escalc, measure = "ROM")
d <- escalc(measure = "ROM", m1i = m_t, sd1i = cv_t * m_t, n1i = n_t,
            m2i = m_c, sd2i = cv_c * m_c, n2i = n_c, data = d, var.names = c("yi", "vi_reported"))
d$yi <- log(d$m_t / d$m_c)   # identical to escalc's yi; also defined where SDs are missing

# Average CV^2 per indicator: average within each study first (so studies
# with many rows do not dominate), then across studies.
# Two documented exceptions:
#  * Bulk density CVs above 0.5 are physically implausible (BD varies by a
#    few percent within a field) and most likely SE/SD or entry errors. They
#    are left out of the imputation constant only; the rows stay in the data.
#  * If fewer than 5 studies report SDs for an indicator (infiltration has 1),
#    its own average is not trustworthy. The largest average CV^2 among the
#    other indicators is borrowed instead, which gives conservative (wide)
#    variances. Results for that indicator are labelled exploratory.
cv_src <- d %>% filter(sd_source != "imputed",
                       !(indicator == "Bulk density" & pmax(cv_c, cv_t) > 0.5))
cv_tab <- cv_src %>%
  group_by(indicator, study_id) %>%
  summarise(cv2_c = mean(cv_c^2), cv2_t = mean(cv_t^2), .groups = "drop") %>%
  group_by(indicator) %>%
  summarise(studies_with_sd = n(), cv2_c = mean(cv2_c), cv2_t = mean(cv2_t),
            .groups = "drop") %>%
  mutate(cv_rule = ifelse(studies_with_sd >= 5, "own indicator", "borrowed (fewer than 5 studies with SD)"))
ok <- cv_tab$studies_with_sd >= 5
max_c <- max(cv_tab$cv2_c[ok]); max_t <- max(cv_tab$cv2_t[ok])
cv_tab <- cv_tab %>%
  mutate(cv2_c = ifelse(studies_with_sd >= 5, cv2_c, max_c),
         cv2_t = ifelse(studies_with_sd >= 5, cv2_t, max_t),
         mean_cv_c = sqrt(cv2_c), mean_cv_t = sqrt(cv2_t))
# Indicators with no reported SDs at all would be missing here; add them
miss_ind <- setdiff(unique(d$indicator), cv_tab$indicator)
if (length(miss_ind)) cv_tab <- bind_rows(cv_tab, tibble(indicator = miss_ind, studies_with_sd = 0,
  cv2_c = max_c, cv2_t = max_t, cv_rule = "borrowed (no study with SD)",
  mean_cv_c = sqrt(max_c), mean_cv_t = sqrt(max_t)))
write_csv(cv_tab, file.path(PATHS$tables, "cv_imputation_values.csv"))
print(as.data.frame(cv_tab))

d <- d %>% left_join(cv_tab %>% select(indicator, cv2_c, cv2_t), by = "indicator") %>%
  mutate(vi_imp = cv2_t / n_t + cv2_c / n_c,
         vi_A = ifelse(sd_source == "imputed", vi_imp, vi_reported),
         vi_B = vi_imp,
         vi_C = ifelse(sd_source == "imputed", NA, vi_reported),
         vi_D = ifelse(sd_source == "imputed", 4 * vi_imp, vi_reported),
         pct_change = pct(yi))

# Extreme ratios are kept (the database authors checked units) but flagged
d$extreme_flag <- abs(d$yi) > log(3)
log_line("Effect sizes: ", nrow(d), "; with reported SD (both or one arm): ",
         sum(d$sd_source != "imputed"), " (", round(100 * mean(d$sd_source != "imputed"), 1), "%)",
         "; |lnRR| > ln(3): ", sum(d$extreme_flag))

write_csv(as.data.frame(d), file.path(PATHS$processed, "effect_sizes.csv"))

summ <- d %>% group_by(indicator, practice) %>%
  summarise(k = n(), studies = n_distinct(study_id),
            pct_sd_reported = round(100 * mean(sd_source != "imputed"), 1),
            pct_n_imputed = round(100 * mean(n_source == "imputed median"), 1),
            median_lnRR = round(median(yi), 3), raw_mean_pct = round(pct(mean(yi)), 1),
            .groups = "drop")
write_csv(summ, file.path(PATHS$tables, "effect_size_summary_raw.csv"))
print(as.data.frame(summ))
