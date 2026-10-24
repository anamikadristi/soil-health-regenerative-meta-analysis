# 07_great_plains_subset.R
# Focused look at evidence that resembles the southern Great Plains (SGP).
#
# Two definitions (set in 01_load_clean.R):
#  * Geographic subset (sgp_geo): U.S. sites at 30 to 40.5 N, 106 to 96 W.
#  * Climate-analog subset (sgp_climate): any site with Koppen BSk/BSh
#    (semi-arid steppe), or Cfa/Dfa with reported MAP <= 900 mm.
#
# For each subset: number of rows, papers and sites; pooled effects where
# the data allow (geographic: at least 5 effect sizes from 3 papers, which is
# very little; climate analog: the usual 10 and 5); and a test of whether
# the analog sites differ from all other sites (subset as a moderator).

source("R/00_setup.R")
log_line("== 07_great_plains_subset ==")

es <- read_csv(file.path(PATHS$processed, "effect_sizes.csv"), show_col_types = FALSE) %>%
  mutate(practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS), vi = vi_A)

# ---- site list for the geographic subset
geo_sites <- es %>% filter(sgp_geo) %>%
  group_by(study_id, author, pub_year, site_info, lat, lon, koppen, map_mm, practice) %>%
  summarise(indicators = paste(sort(unique(as.character(indicator))), collapse = "; "), k = n(), .groups = "drop") %>%
  arrange(lat)
write_csv(geo_sites, file.path(PATHS$tables, "sgp_geographic_site_list.csv"))
print(as.data.frame(geo_sites %>% select(study_id, author, site_info, koppen, map_mm, practice, k)))

# ---- pooled estimates within each subset
pool_subset <- function(flag, label, min_k, min_s) {
  es %>% filter(.data[[flag]]) %>% group_by(indicator, practice) %>% group_split() %>%
    map_dfr(function(dat) {
      out <- tibble(subset = label, indicator = as.character(dat$indicator[1]), practice = as.character(dat$practice[1]),
                    k = nrow(dat), studies = n_distinct(dat$study_id), sites = n_distinct(paste(dat$lat, dat$lon)),
                    median_raw_pct = pct(median(dat$yi)))
      if (nrow(dat) < min_k || n_distinct(dat$study_id) < min_s) return(out %>% mutate(pooled = FALSE))
      f <- fit_ml(dat)
      out %>% mutate(pooled = TRUE, pct = pct(f$b[1]), pct_lb = pct(f$ci.lb), pct_ub = pct(f$ci.ub),
                     p = f$pval, df = f$ddf[1], I2_total = i2_ml(f)[1])
    })
}
sgp <- bind_rows(pool_subset("sgp_geo", "Geographic (U.S. southern/central Great Plains)", 5, 3),
                 pool_subset("sgp_climate", "Climate analog (semi-arid / subhumid, MAP <= 900 mm)", MIN_K, MIN_STUDIES))
write_csv(sgp %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(PATHS$tables, "sgp_subset_pooled.csv"))
print(as.data.frame(sgp %>% mutate(across(where(is.double), ~ signif(.x, 3)))))

# ---- do analog sites differ from the rest? (within each pooled global cell)
pooled <- read_csv(file.path(PATHS$tables, "pooled_effects.csv"), show_col_types = FALSE)
cmp <- pmap_dfr(pooled %>% select(indicator, practice), function(indicator, practice) {
  dat <- es %>% filter(indicator == !!indicator, practice == !!practice) %>%
    mutate(grp = factor(ifelse(sgp_climate, "analog", "other"), c("other", "analog")))
  s_an <- n_distinct(dat$study_id[dat$grp == "analog"])
  if (sum(dat$grp == "analog") < 5 || s_an < 3) return(NULL)
  f <- fit_ml(dat, ~ grp)
  tibble(indicator, practice, k_analog = sum(dat$grp == "analog"), studies_analog = s_an,
         k_other = sum(dat$grp == "other"),
         diff_lnRR = f$b[2], diff_lb = f$ci.lb[2], diff_ub = f$ci.ub[2], p = f$pval[2])
})
write_csv(cmp %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(PATHS$tables, "sgp_analog_vs_other_test.csv"))
print(as.data.frame(cmp %>% mutate(across(where(is.double), ~ signif(.x, 3)))))

# ---- representativeness: compare site conditions
rep_tab <- es %>% mutate(group = case_when(sgp_geo ~ "Geographic SGP", sgp_climate ~ "Climate analog (non-SGP)", TRUE ~ "All other")) %>%
  group_by(group) %>%
  summarise(rows = n(), papers = n_distinct(study_id), sites = n_distinct(paste(lat, lon)),
            median_MAP_mm = median(map_mm, na.rm = TRUE), median_MAT_C = median(mat_c, na.rm = TRUE),
            pct_coarse_texture = round(100 * mean(texture == "Coarse", na.rm = TRUE), 1),
            median_years = median(years, na.rm = TRUE),
            pct_grassland = round(100 * mean(land_use != "Cropland"), 1),
            pct_cover_crop = round(100 * mean(practice == "Cover crop"), 1), .groups = "drop")
write_csv(rep_tab, file.path(PATHS$tables, "sgp_representativeness.csv"))
print(as.data.frame(rep_tab))

# ---- figures: site map and subset forest plot
world <- map_data("world")  # needs the maps package
sites <- es %>% distinct(lat, lon, practice, sgp_geo, sgp_climate) %>% filter(!is.na(lat)) %>%
  mutate(group = case_when(sgp_geo ~ "SGP geographic", sgp_climate ~ "Climate analog", TRUE ~ "Other"))
p_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey90", colour = "grey70", linewidth = 0.15) +
  geom_point(data = sites %>% arrange(group != "Other"), aes(lon, lat, colour = group), size = 1.3, alpha = 0.8) +
  annotate("rect", xmin = -106, xmax = -96, ymin = 30, ymax = 40.5, fill = NA, colour = "black", linewidth = 0.4) +
  scale_colour_manual(values = c("SGP geographic" = "#C0392B", "Climate analog" = "#E69F00", "Other" = "#4C72B0"), name = NULL) +
  coord_quickmap(ylim = c(-50, 70)) + labs(x = NULL, y = NULL,
  title = "Sites in the analysed subset of SoilHealthDB",
  subtitle = "Box: U.S. southern/central Great Plains window (30-40.5 N, 106-96 W)") +
  theme(legend.position = "bottom")
save_png(p_map, "fig07a_site_map_sgp.png", w = 9, h = 5.2)

fp <- sgp %>% filter(pooled) %>% mutate(lab = paste0(indicator, " | ", practice, "  (k=", k, ", n=", studies, ")"))
p_sgp <- ggplot(fp, aes(pct, reorder(lab, pct), colour = subset)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = pct_lb, xmax = pct_ub), size = 0.35) +
  facet_wrap(~ subset, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#C0392B", "#E69F00"), guide = "none") +
  labs(x = "Change relative to control (%), 95% CI", y = NULL,
       title = "Southern Great Plains subsets",
       subtitle = "Geographic subset is very small; treat as descriptive")
save_png(p_sgp, "fig07b_sgp_subset_forest.png", w = 8, h = 6)
log_line("SGP geographic rows: ", sum(es$sgp_geo), " from ", n_distinct(es$study_id[es$sgp_geo]),
         " papers; climate-analog rows: ", sum(es$sgp_climate), " from ", n_distinct(es$study_id[es$sgp_climate]), " papers")
