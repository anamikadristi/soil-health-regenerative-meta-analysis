# 06_figures.R
# Figures (PNG, 300 dpi) from the tables written by scripts 02 to 05.
# The orchard plots are drawn with ggplot2 in the style of Nakagawa et al.
# (2021); the orchaRd package is not on CRAN and is not required.

source("R/00_setup.R")
log_line("== 06_figures ==")

es <- read_csv(file.path(PATHS$processed, "effect_sizes.csv"), show_col_types = FALSE) %>%
  mutate(practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS), vi = vi_A)
pooled <- read_csv(file.path(PATHS$tables, "pooled_effects.csv"), show_col_types = FALSE) %>%
  mutate(practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS))
fits <- readRDS(file.path(PATHS$processed, "fits_pooled.rds"))
pal <- c("Cover crop" = "#4C72B0", "No-till / reduced tillage" = "#DD8452", "Organic management" = "#55A868",
         "Residue retention" = "#8172B3", "Agroforestry" = "#937860", "Crop rotation" = "grey50")
pct_axis <- sec_axis(~ pct(.), name = "Change relative to control (%)", breaks = c(-50, -25, 0, 25, 50, 100, 200))

# ---- 1. summary forest plot
fp <- pooled %>% mutate(lab = sprintf("%s  (k=%d, %d studies)", practice, k, studies),
                        pi_lb = log(1 + pi_pct_lb / 100), pi_ub = log(1 + pi_pct_ub / 100))
p1 <- ggplot(fp, aes(lnRR, reorder(lab, as.numeric(practice) * -1), colour = practice)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbar(orientation = "y", aes(xmin = pi_lb, xmax = pi_ub), width = 0, linewidth = 0.4, alpha = 0.45) +
  geom_errorbar(orientation = "y", aes(xmin = ci_lb, xmax = ci_ub), width = 0, linewidth = 1.1) +
  geom_point(size = 2.4, shape = 21, fill = "white", stroke = 1.1) +
  facet_grid(indicator ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(sec.axis = pct_axis) +
  labs(x = "Pooled lnRR (thick: 95% CI; thin: 95% prediction interval)", y = NULL,
       title = "Pooled effects of conservation practices on soil health indicators",
       subtitle = "Three-level random-effects models; re-analysis of SoilHealthDB (Jian et al. 2020)") +
  theme(strip.text.y = element_text(angle = 0))
save_png(p1, "fig06a_forest_pooled_effects.png", w = 9, h = 8)

# ---- 2. orchard plots
orchard <- function(dat, group, pool, title, file, w = 8, h = 5) {
  p <- ggplot(dat, aes(yi, .data[[group]])) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_jitter(aes(size = 1 / sqrt(vi), colour = practice), height = 0.25, alpha = 0.25, shape = 16) +
    geom_errorbar(orientation = "y", data = pool, aes(x = lnRR, xmin = pi_lb, xmax = pi_ub, y = .data[[group]]), width = 0, linewidth = 0.6, inherit.aes = FALSE) +
    geom_errorbar(orientation = "y", data = pool, aes(x = lnRR, xmin = ci_lb, xmax = ci_ub, y = .data[[group]]), width = 0, linewidth = 2, inherit.aes = FALSE) +
    geom_point(data = pool, aes(x = lnRR, y = .data[[group]]), size = 3.5, shape = 21, fill = "white", stroke = 1.2, inherit.aes = FALSE) +
    geom_text(data = pool, aes(x = Inf, y = .data[[group]], label = paste0("k=", k, " (", studies, ")")),
              hjust = 1.05, vjust = -1.2, size = 3, inherit.aes = FALSE) +
    scale_colour_manual(values = pal, guide = "none") + scale_size_continuous(range = c(0.4, 4), guide = "none") +
    scale_x_continuous(sec.axis = pct_axis) + coord_cartesian(xlim = c(-1.2, 1.6)) +
    labs(x = "lnRR (points: effect sizes sized by precision; thick bar: 95% CI; thin bar: 95% PI)", y = NULL, title = title)
  save_png(p, file, w, h)
}
pp <- fp %>% filter(indicator == "SOC concentration")
orchard(es %>% filter(indicator == "SOC concentration", practice %in% pp$practice), "practice", pp,
        "Orchard plot: SOC concentration by practice", "fig06b_orchard_SOC_by_practice.png")
pc <- fp %>% filter(practice == "Cover crop")
orchard(es %>% filter(practice == "Cover crop", indicator %in% pc$indicator), "indicator", pc,
        "Orchard plot: cover crop effects across indicators", "fig06c_orchard_covercrop_by_indicator.png", h = 5.5)

# ---- 3. caterpillar plot, SOC concentration x cover crop
cat_d <- es %>% filter(indicator == "SOC concentration", practice == "Cover crop") %>%
  mutate(lo = yi - 1.96 * sqrt(vi), hi = yi + 1.96 * sqrt(vi)) %>% arrange(yi) %>% mutate(rank = row_number())
f0 <- fits[["SOC concentration | Cover crop"]]
p3 <- ggplot(cat_d, aes(yi, rank)) +
  geom_errorbar(orientation = "y", aes(xmin = lo, xmax = hi, colour = sd_source != "imputed"), width = 0, linewidth = 0.25, alpha = 0.6) +
  geom_point(size = 0.5) +
  annotate("rect", xmin = f0$ci.lb, xmax = f0$ci.ub, ymin = -Inf, ymax = Inf, alpha = 0.25, fill = "#C0392B") +
  geom_vline(xintercept = f0$b[1], colour = "#C0392B") + geom_vline(xintercept = 0, linetype = 2) +
  scale_colour_manual(values = c(`TRUE` = "#4C72B0", `FALSE` = "grey65"), labels = c(`TRUE` = "SD reported", `FALSE` = "SD imputed"), name = NULL) +
  scale_x_continuous(sec.axis = pct_axis) + coord_cartesian(xlim = c(-1, 1.3)) +
  labs(x = "lnRR with 95% CI", y = "Effect sizes, ranked",
       title = sprintf("Caterpillar plot: cover crops and SOC concentration (k = %d, %d studies)", f0$k, nrow(distinct(cat_d, study_id))),
       subtitle = sprintf("Red line and band: pooled estimate %.1f%% [%.1f, %.1f]", pct(f0$b[1]), pct(f0$ci.lb), pct(f0$ci.ub))) +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_png(p3, "fig06d_caterpillar_SOC_covercrop.png", w = 7, h = 8)

# ---- 4. bubble plots for continuous moderators (SOC concentration x cover crop)
bubble <- function(var, xlab, file) {
  d <- es %>% filter(indicator == "SOC concentration", practice == "Cover crop", !is.na(.data[[var]])) %>% mutate(m = .data[[var]])
  f <- fit_ml(d, ~ m)
  nd <- tibble(m = seq(min(d$m), max(d$m), length.out = 100))
  pr <- predict(f, newmods = nd$m)
  nd <- nd %>% mutate(fit = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub)
  p <- ggplot(d, aes(m, yi)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(size = 1 / sqrt(vi)), alpha = 0.25, colour = "#4C72B0") +
    geom_ribbon(data = nd, aes(x = m, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.25) +
    geom_line(data = nd, aes(m, fit), inherit.aes = FALSE, linewidth = 1) +
    scale_size_continuous(range = c(0.4, 4), guide = "none") + coord_cartesian(ylim = c(-1, 1.3)) +
    scale_y_continuous(sec.axis = sec_axis(~ pct(.), name = "Change (%)")) +
    labs(x = xlab, y = "lnRR", title = paste("Cover crops, SOC concentration:", tolower(xlab)),
         subtitle = sprintf("Slope = %.4f lnRR per unit, p = %.3f, k = %d, %d studies", f$b[2], f$pval[2], f$k, n_distinct(d$study_id)))
  save_png(p, file, w = 7, h = 4.5)
}
bubble("years", "Years under practice", "fig06e_bubble_SOC_covercrop_years.png")
bubble("depth_mid", "Sampling depth midpoint (cm)", "fig06f_bubble_SOC_covercrop_depth.png")

# ---- 5. sensitivity of pooled estimates to variance strategy and outliers
sens <- read_csv(file.path(PATHS$tables, "sensitivity_variance_and_outliers.csv"), show_col_types = FALSE) %>%
  filter(!is.na(pct)) %>% mutate(cell = paste(indicator, practice, sep = "\n"),
                                 analysis = factor(analysis, rev(unique(analysis))))
main <- c("SOC concentration\nCover crop", "SOC concentration\nOrganic management", "SOC stock\nCover crop",
          "Bulk density\nCover crop", "Aggregate stability\nCover crop", "Microbial biomass C\nCover crop")
p5 <- ggplot(sens %>% filter(cell %in% main) %>% mutate(cell = factor(cell, main)), aes(pct, analysis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = pct_lb, xmax = pct_ub), size = 0.3, colour = "#4C72B0") +
  geom_text(aes(x = Inf, label = paste0("k=", k)), hjust = 1.1, size = 2.6) +
  facet_wrap(~ cell, scales = "free_x", ncol = 3) +
  labs(x = "Pooled change (%) with 95% CI", y = NULL, title = "Sensitivity of pooled estimates to variance handling and extreme ratios")
save_png(p5, "fig06g_sensitivity_variance_strategies.png", w = 11, h = 6)

# ---- 6. leave-one-study-out, SOC concentration x cover crop
loo <- read_csv(file.path(PATHS$tables, "leave_one_study_out.csv"), show_col_types = FALSE) %>%
  filter(indicator == "SOC concentration", practice == "Cover crop") %>% arrange(pct) %>% mutate(r = row_number())
p6 <- ggplot(loo, aes(pct, r)) + geom_vline(xintercept = pct(f0$b[1]), colour = "#C0392B") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = pct_lb, xmax = pct_ub), size = 0.15) +
  labs(x = "Pooled change in SOC (%) with one study removed", y = "Study removed (ranked)",
       title = "Leave-one-study-out: cover crops and SOC concentration") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_png(p6, "fig06h_leave_one_out_SOC_covercrop.png", w = 6.5, h = 6)
log_line("Figures written to ", PATHS$figures)
