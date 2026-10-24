# Results summary

This is a re-analysis of the public SoilHealthDB (Jian, Du & Stewart 2020, *Scientific Data* 7:16). No new field data were collected. All numbers below come from the tables in `outputs/tables/` and can be reproduced with `run_all.sh`.

**Notation.** lnRR = ln(treatment mean / control mean). Percent change = (exp(lnRR) - 1) x 100. k = number of effect sizes; studies = papers; experiments = sites or trials within a paper. 95% CI = confidence interval of the mean effect. 95% PI = prediction interval, the range where the effect at a new site would be expected to fall. I2 = share of total variance that is heterogeneity rather than sampling error (multilevel version, Nakagawa & Santos 2012). tau2 = sum of the three variance components (study + experiment + within experiment), on the lnRR scale.

## 1. Data used
- Database: 5,907 rows from 321 papers.
- After keeping single-practice contrasts and the six indicators: **3,344 effect sizes from 186 papers**.
- Measurement screening removed 83 rows (`audit_measurement_screening.csv`):
  - SOC rows that were total carbon;
  - aggregate rows that were single size classes or the wind-erodible fraction;
  - microbial biomass rows based on PLFA.
- Combined contrasts were left out. For example, cover crop plus a tillage change against a conventional control (`audit_practice_assignment.csv`).
- Only **22% of effect sizes have reported SDs**. For cover crop x SOC concentration the figure is 9%. Missing variances were imputed from the mean coefficient of variation (CV) of rows with SDs (`cv_imputation_values.csv`). For infiltration, only 1 paper reported SDs, so a conservative CV was borrowed from the indicator with the largest CV (aggregate stability).
- Crop rotation had only 3 usable SOC comparisons from 1 paper. It could not be analysed.
- Cells with fewer than 10 effect sizes or 5 papers were not pooled (`data_availability.csv`).

## 2. Pooled effects (three-level random-effects models)
From `pooled_effects.csv`. CR2 cluster-corrected CIs by study were almost identical to the model CIs in every cell.

| Indicator | Practice | k | Studies | % with SD | % change [95% CI] | p | 95% PI (%) | I2 total (%) | tau2 |
|---|---|---|---|---|---|---|---|---|---|
| SOC concentration | Cover crop | 905 | 62 | 9 | **12.1 [7.8, 16.5]** | <0.001 | -16.3 to 49.9 | 69.9 | 0.021 |
| SOC concentration | No-till / reduced tillage | 35 | 6 | 0 | 13.0 [-11.2, 43.7] | 0.25 | -54.3 to 179 | 81.0 | 0.115 |
| SOC concentration | Organic management | 316 | 74 | 56 | **17.8 [12.1, 23.8]** | <0.001 | -28.2 to 93.2 | 99.5 | 0.061 |
| SOC concentration | Residue retention | 34 | 5 | 0 | 18.2 [-13.6, 61.7] | 0.21 | -55.4 to 213 | 81.5 | 0.110 |
| SOC concentration | Agroforestry | 99 | 12 | 16 | **28.4 [9.3, 50.9]** | 0.006 | -28.4 to 130 | 83.4 | 0.065 |
| SOC stock | Cover crop | 448 | 21 | 32 | **8.4 [5.5, 11.5]** | <0.001 | -2.8 to 21.0 | 23.4 | 0.003 |
| SOC stock | Organic management | 148 | 42 | 60 | **16.5 [8.9, 24.6]** | <0.001 | -27.6 to 87.4 | 99.6 | 0.054 |
| Bulk density | Cover crop | 565 | 39 | 6 | **-1.8 [-3.0, -0.6]** | 0.003 | -7.2 to 3.8 | 16.0 | 0.001 |
| Bulk density | No-till / reduced tillage | 68 | 9 | 79 | 1.8 [-2.0, 5.8] | 0.32 | -13.9 to 20.3 | 90.2 | 0.005 |
| Bulk density | Organic management | 109 | 28 | 31 | **-3.9 [-5.4, -2.3]** | <0.001 | -13.6 to 6.9 | 52.8 | 0.003 |
| Aggregate stability | Cover crop | 123 | 16 | 7 | 8.6 [0.4, 17.4] | 0.041 | -10.5 to 31.7 | 8.7 | 0.007 |
| Aggregate stability | Organic management | 62 | 9 | 32 | 22.2 [-4.9, 57.1] | 0.10 | -48.2 to 189 | 97.8 | 0.127 |
| Microbial biomass C | Cover crop | 176 | 20 | 5 | **22.1 [8.4, 37.4]** | 0.002 | -37.1 to 137 | 72.8 | 0.097 |
| Microbial biomass C | Organic management | 121 | 20 | 55 | **18.2 [3.7, 34.8]** | 0.015 | -39.7 to 132 | 87.7 | 0.099 |
| Infiltration rate (exploratory) | Cover crop | 30 | 7 | 0 | 67.4 [29.5, 116] | 0.003 | -35.7 to 336 | 58.9 | 0.142 |

**Where the heterogeneity sits.** For cover crop x SOC concentration, almost all heterogeneity is between studies (I2 study = 68.7%; I2 experiment and I2 within-experiment are each under 1%). The effect of cover crops on SOC depends mainly on the study context, not on the depth or year within a study. For organic management, heterogeneity is split between studies (50%), within experiments (34%) and experiments (15%).

**Do practices differ for SOC concentration?** A single model with practice as a moderator (k = 1,389, 156 studies) found no significant difference: F(4, 151) = 1.81, p = 0.13. Agroforestry was higher than cover crops (difference in lnRR 0.13 [0.01, 0.25], p = 0.035), but that is one contrast out of four.

## 3. Moderators (one moderator per model; `moderator_omnibus_tests.csv`, `moderator_results.csv`)
Tested for cells with at least 50 effect sizes and 15 studies. Sixty-six omnibus tests were run, so Holm-adjusted p-values are given.

**Cover crops, SOC concentration (k up to 905)**
- **Depth:**
  - Surface (bottom <= 15 cm): +12.2% [7.9, 16.7].
  - 0-30 cm layers: +11.4% [6.9, 16.1].
  - Deeper layers: +6.9% [2.1, 12.0].
  - Tests: class p = 0.013; continuous slope -0.09% per cm, p = 0.001 (Holm p = 0.076).
- **Years under practice:** +0.50% per additional year (p = 0.0085, Holm 0.51). By class, 1-3 yr gave +9.7%, 4-10 yr +11.6% and > 10 yr +19.1% (p = 0.079).
- **Mean annual precipitation:** -2.5% per 100 mm (p = 0.004, Holm 0.25). This used only 284 effect sizes from 16 studies, because MAP is missing for most rows.
- **No detectable effect:**
  - Koppen climate zone (p = 0.80): temperate +13.2%, cold/continental +10.4%, tropical +9.5%.
  - Soil texture (p = 0.24): coarse +17.3%, medium +11.1%, fine +10.6%.
  - Land use (p = 0.21).

**Organic management, SOC concentration**
- **Depth** is the strongest moderator in the project.
  - Surface: +20.2% [12.8, 28.0].
  - 0-30 cm: +20.3% [11.3, 30.0].
  - Deeper layers: -5.9% [-16.0, 5.3].
  - Tests: Holm p = 0.005; continuous slope -1.1% per cm, Holm p = 0.0002.
- Gains in organic systems are concentrated in the topsoil.

**SOC stock**
- **Organic management:** stocks rose more after more than 10 years (+28.9% [18.1, 40.6]) than after 1-3 years (+10.1%). The continuous slope was +0.82% per year (p = 0.0015, Holm 0.09).
- **Cover crops:** a similar but weaker pattern (p = 0.019, Holm 1.0).

**Other cells**
- No moderator was significant for bulk density, aggregate stability or microbial biomass C. The one exception was MAT for organic bulk density (p = 0.009, Holm 0.51).

**Land use**
- Only 172 effect sizes (5%) come from grassland or pasture. For cover crop SOC these are 43 effect sizes from 3 studies (-0.8% [-18.4, 20.6]).
- There are **no rangeland grazing studies** in this database.

## 4. Publication bias and sensitivity (`publication_bias.csv`, `sensitivity_variance_and_outliers.csv`, `leave_one_study_out.csv`)

**Small-study effects**
- The multilevel Egger test (sqrt(1/n_eff) as moderator) was not significant in any cell (all p > 0.08).
- The classic Egger test on study averages was significant for:
  - agroforestry SOC (p = 0.020);
  - cover crop SOC stock (p = 0.048).
- Trim-and-fill on study averages did not change the cover crop SOC concentration estimate (0 studies added).
  - For cover crop SOC stock it added 4 studies and lowered the study-level estimate from +7.0% to +5.6% [1.9, 9.5]. The effect stays positive.
  - For organic SOC concentration it added 20 studies on the positive side. The asymmetry there is not the pattern expected from suppressed null results.
- **Time-lag:** no decline with publication year, except for cover crop aggregate stability (p = 0.03, k = 123, 16 studies).

**Leave-one-study-out (cover crops, SOC concentration)**
- Dropping any single study left the estimate between +10.5% and +12.4%, and it stayed significant every time.
- The same held for cover crop SOC stock, bulk density and microbial biomass C, and for organic SOC.

**Fragile results**
- **Cover crop aggregate stability** became non-significant when any of 10 of the 16 studies was dropped. It was also non-significant under the all-cases variance strategy (B), when rows without reported n were removed (F), and after trim-and-fill (+0.5% [-9.2, 11.3], 4 studies added).
- Treat the +8.6% as **weak evidence**.

**Variance strategies (cover crops, SOC concentration)**

| Strategy | Estimate [95% CI] |
|---|---|
| A, primary | +12.1% [7.8, 16.5] |
| B, imputed CV for all rows | +11.7% [7.6, 15.9] |
| C, complete cases only (k = 80, 12 studies) | +13.4% [1.2, 27.1] |
| D, imputed variance x 4 | +10.7% [7.1, 14.4] |
| E, without 3-fold ratios | +10.5% [7.3, 13.8] |

The direction and size are stable.

**Other cells under the variance strategies**
- Strategies A, B, D, E and F gave the same conclusions for:
  - cover crop SOC stock, bulk density and microbial biomass C;
  - organic SOC concentration, SOC stock, bulk density and microbial biomass C.
- The complete-case analysis (C) is much smaller, and two cells lost significance under it:
  - cover crop bulk density: -3.0% [-11.2, 6.0], k = 36, 5 studies;
  - organic microbial biomass C: +15.1% [-5.4, 40.0], k = 67, 13 studies.
- For cover crop SOC stock, aggregate stability and microbial biomass C, too few studies reported SDs for a complete-case analysis.

## 5. Southern Great Plains subsets (`sgp_*.csv`, figures `fig07a`, `fig07b`)

**Geographic subset** (U.S. sites at 30 to 40.5 N, 106 to 96 W)
- Size: **94 effect sizes from 5 papers at 6 sites**:
  - Bushland, TX;
  - Garden City, KS;
  - Hesston, KS;
  - Akron, CO;
  - two sites in south-central Nebraska.
- 89% of the effect sizes are cover crop comparisons; all are cropland.
- Only two cells could be pooled even with a relaxed threshold (5 effect sizes, 3 papers):
  - Cover crops, SOC concentration: **+17.7% [6.4, 30.2]**, k = 28, 3 papers, df = 2.
  - Cover crops, bulk density: -3.3% [-9.6, 3.6], k = 28, 3 papers.
- With 3 papers, these numbers are descriptive. They do not describe the region as a whole.

**Climate-analog subset** (Koppen BSk/BSh, or Cfa/Dfa with MAP <= 900 mm, anywhere)
- Size: 492 effect sizes from 26 papers, mainly from:
  - the U.S. Great Plains and Midwest;
  - Argentina's Pampas;
  - northern China;
  - Australia.
- Pooled results:
  - Cover crops, SOC concentration: +12.0% [4.2, 20.4] (k = 101, 5 papers).
  - No-till / reduced tillage, SOC concentration: +12.8% [-14.7, 49.1] (k = 34, 5 papers).
  - Organic management, SOC concentration: +17.5% [5.2, 31.3] (k = 47, 11 papers).
  - Organic management, SOC stock: +22.8% [-8.4, 64.7].
  - Organic management, bulk density: -3.5% [-8.0, 1.1].
- **Analog versus other sites:** analog sites did not differ from the rest in any cell. For cover crop SOC the lnRR difference was +0.04 [-0.02, 0.10], p = 0.16.

**Is it representative?**
- No, not fully.
- The geographic subset has a median MAP of 489 mm and no coarse-textured soils.
- Most of the data come from a few long-term experiment stations.
- Grazed rangeland, the other half of the target landscape, is not represented at all.

## 6. Interpretation
- **Cover crops and SOC.** The evidence that cover crops raise SOC in cropland is consistent and holds up under the checks run here. The typical gain is about 8 to 12%, depending on whether concentration or stock is used. The prediction interval is wide (-16% to +50% for concentration). A single field can show no gain, or even a loss, even though the average is positive.
- **Stock versus concentration.** The stock estimate is smaller than the concentration estimate. Stocks depend on bulk density and depth, and cover crops slightly lowered bulk density, so this gap makes sense.
- **Microbial biomass and bulk density.** Microbial biomass C responds faster and more strongly (+22%) than SOC, which fits its role as an early indicator. Bulk density changes are small (-2 to -4%) and probably below what most field sampling can detect at a single site.
- **Organic management.** It increases topsoil SOC, but the very high heterogeneity (I2 near 100%) and the depth pattern suggest the gain mostly reflects added organic inputs at the surface. A gain in total profile carbon is less clear.
- **Tillage.** This database is not a good source for tillage effects (6 to 9 papers). The no-till estimates here should not replace meta-analyses designed for tillage.
- **Southern Great Plains.** Semi-arid analog sites show effects similar to the global average. The direct evidence from the region itself is very thin. That is a data gap, and it is the kind of gap a regional synthesis would need to fill.

## 7. Limitations
1. **Re-analysis of one database.** Study selection, data extraction and unit conversion were done by the database authors. I did not re-check the original papers. The database covers papers published up to about 2019.
2. **Missing variances.** Most variances are imputed. The checks in Section 4 show that the main conclusions do not depend on this. Still, I2 and tau2 values depend partly on the imputed CVs. Cells with large imputed CVs (aggregate stability, infiltration) show low I2 partly for that reason.
3. **Non-independence.** The three-level model handles clustering of effects within experiments and studies. It does not model the correlation of sampling errors when many treatments share one control. The cluster-corrected (CR2) CIs, which make no assumption about that correlation, were nearly identical to the model CIs.
4. **Mixed measures.**
   - SOC concentration mixes SOC and SOM. This is fine for a ratio.
   - Aggregate stability mixes MWD and water-stable aggregate percentages, after screening.
   - SOC stocks were not corrected to an equivalent soil mass.
5. **Moderators.**
   - They are tested one at a time, are often missing (MAP is missing for most cover crop rows) and are confounded with country and study design.
   - Results are observational associations across studies, not causal effects.
   - Many tests were run; after Holm correction only the organic-management depth effect remains clearly significant.
   - The "contain" degrees-of-freedom method gives large df to some practice contrasts in the single SOC model. The practice comparison should be read with that in mind.
6. **Design differences.** Some cover crop contrasts also differ in fertilizer (Fert_Diff = "Yes"). These were kept, as in the original database, and could inflate or dampen effects.
7. **Rangeland.** The database is almost entirely cropland. It says nothing direct about grazing management, rangeland seeding or prescribed fire.
8. **Climate data.** Koppen classes come from a 0.5 degree grid and may misclassify sites near climate boundaries.
