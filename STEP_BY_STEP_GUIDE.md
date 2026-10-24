# Step-by-step guide

This guide explains how to re-run the project and what each step does, so you can explain every number. Run commands from the project folder (`05_MetaAnalysis_Regenerative_Practices_Soil_Health/`).

## Step 0. Set up and download the data
```bash
Rscript install_packages.R            # metafor, clubSandwich, readxl, dplyr, tidyr, stringr, readr, purrr, ggplot2, maps
python3 -I python/00_fetch_soilhealthdb.py data_raw
```
**What happens**
1. The script asks the figshare API for record 8292176 (SoilHealthDB, version 8, the version cited in the paper).
2. The record holds one 128 MB zip. The script reads only three files out of it with HTTP range requests, so the download is about 5 MB:
   - `SoilHealthDB_V1.xlsx`
   - `SoilHealthDB_V1.csv`
   - `Koeppen-Geiger-ASCII.txt`
3. `zipfile` checks each file's CRC. The SHA-256 checksums are written to `data_raw/download_log.txt`.

**Check:** `data_raw/` should contain the three files and the log. If the figshare record changes, compare the checksums before going on.

## Step 1. Load and clean (`R/01_load_clean.R`)
```bash
Rscript R/01_load_clean.R
```
**What it does**
1. Reads Sheet1 of the xlsx: 5,907 rows, 270 columns. Each row is one treatment-versus-control comparison.
2. **Practice.** Assigns a practice from `Conservation_Type`:

   | Group | Codes |
   |---|---|
   | Cover crop | CC, CC-NT, CC-RT, only when tillage is the same in both arms |
   | No-till / reduced tillage | NT, RT |
   | Organic management | OF, OF-CC, OF-NT |
   | Residue retention | SR |
   | Agroforestry | AF |
   | Crop rotation | RO |

   A row where the cover crop arm also changed tillage is a combined contrast, so it is dropped. See `outputs/tables/audit_practice_assignment.csv`.
3. **Indicators.** Builds one long table with one row per effect size. The six indicators are:
   - SOC concentration (`OC_C`, `OC_T`);
   - SOC stock (`OC_Con_*`);
   - bulk density (`BD_*`);
   - aggregate stability (`Aggre_*`);
   - microbial biomass C (`MBC_*`);
   - infiltration rate (`Infiltration_*`).
4. **Unit harmonization.**
   - lnRR is a ratio, so units cancel as long as the treatment and control in a row share a unit, which they do in the database. No conversion between g/kg and % is needed.
   - What must not be mixed is different *quantities*. The comment fields are therefore screened:
     - total carbon is dropped from SOC;
     - single aggregate size classes and the wind-erodible fraction (< 0.84 mm, where more is worse) are dropped;
     - PLFA-based "microbial biomass" is dropped.
   - Every decision is listed in `audit_measurement_screening.csv`.
5. **Moderators.** Adds:
   - Koppen climate class, from the nearest 0.5 degree grid cell to each site's coordinates;
   - texture class (coarse, medium or fine);
   - depth class and depth midpoint, parsed from text such as "0-to-15";
   - years under the practice;
   - land use (cropland or grassland/pasture, from crop names);
   - two southern Great Plains flags (see Step 7).
6. **Sample size.** Uses the number of plots per arm (`No_C`, `No_T`), else `NoSubsample`, else the median for that indicator. The source is recorded in `n_source`.

**Check:** the log should report 3,344 comparisons from 186 papers. Open `data_availability.csv` to see k and the number of studies for each practice x indicator.

## Step 2. Effect sizes (`R/02_effect_sizes.R`)
**Effect size.** lnRR = ln(mean_T / mean_C), computed with `escalc(measure = "ROM")`.
- A value of 0.10 means the treatment is exp(0.10) = 1.105 times the control, a 10.5% increase.

**Sampling variance.** v = SD_T^2 / (n_T mean_T^2) + SD_C^2 / (n_C mean_C^2). This is the delta-method variance, with CV = SD / mean.

**Missing SDs.** Only 22% of effect sizes have SDs. Rather than dropping the other 78%, the average CV^2 of rows that report SDs is used, by indicator. It is averaged within each study first, so one large study does not dominate. The imputed variance is v = CVbar_T^2 / n_T + CVbar_C^2 / n_C (Nakagawa et al. 2023). Two documented exceptions:
- Bulk density CVs above 0.5 are physically implausible, so they are left out of the average.
- Infiltration had SDs in only 1 paper, so it borrows the largest average CV (aggregate stability), which gives wide, conservative variances.

**Four variance versions** are stored and compared in Step 5:

| Version | Rule |
|---|---|
| A (primary) | Reported SD where available, otherwise imputed |
| B | Imputed for every row |
| C | Reported SDs only (complete cases) |
| D | Imputed variance multiplied by 4 |

**Check:** `cv_imputation_values.csv` shows the CVs used. For example, SOC concentration has CV of about 0.17 from 65 studies.

## Step 3. Pooled random-effects models (`R/03_random_effects.R`)
**Model.** For each practice x indicator with at least 10 effect sizes and 5 papers:
```r
rma.mv(yi, vi, random = ~ 1 | study_id / exp_id / es_id, method = "REML", test = "t", dfs = "contain")
```
**Why three levels.** One paper often reports many effect sizes (several sites, years, depths and cover crop species), and they share a control and a location. Treating them as independent would make the CI too narrow. The three random effects are:
- between papers;
- between experiments (sites) within a paper;
- between effect sizes within an experiment.

**Outputs in `pooled_effects.csv`**
- Pooled lnRR and % change, with 95% CI.
- 95% prediction interval: where a new site's effect is expected to fall.
- A cluster-corrected (CR2, `clubSandwich`) CI by paper, as a check on remaining dependence.
- The three variance components (sigma2), their sum (tau2), and multilevel I2, also split by level.

**How to read one row.** Cover crop x SOC concentration:
- lnRR = 0.114, which is +12.1% [7.8, 16.5], from 905 effect sizes in 62 papers.
- I2 = 70%: most of the spread in observed effects is real variation among studies, not sampling noise.
- PI = -16% to +50%: individual sites vary widely.

The full printout is in `outputs/logs/model_SOC_covercrop_full.txt`.

## Step 4. Moderators (`R/04_moderators.R`)
**What it tests.** For cells with at least 50 effect sizes and 15 papers, each moderator is added one at a time:
- Koppen zone;
- MAP;
- MAT;
- texture;
- years (as a class and continuous);
- depth (as a class and continuous);
- land use;
- NT versus RT.

Levels with fewer than 5 effect sizes or 3 papers are dropped.

**What each test reports**
- **QM (F-test):** does the moderator explain heterogeneity?
- **Pseudo-R2:** the drop in tau2.
- **Estimates:** one per level, or a slope for a continuous moderator.
- **Holm-adjusted p:** about 66 tests are run, so some p < 0.05 results are expected by chance.

The script also fits one SOC model with practice as the moderator, to test whether practices differ.

**Main findings:** depth matters for organic SOC (Holm p = 0.005). The cover crop SOC effect weakens with depth and grows with years, but these do not survive Holm correction.

## Step 5. Publication bias and sensitivity (`R/05_bias_sensitivity.R`)
**Bias tests**
1. **Multilevel Egger test.** sqrt(1/n_eff) is used as a moderator, with n_eff = n_T n_C / (n_T + n_C). The SE is not used for lnRR, because the SE is computed from the means and is mostly imputed here.
2. **Time-lag test.** Publication year as a moderator.
3. **Classic Egger test and trim-and-fill.** Both are run on **study averages** (`aggregate()` with rho = 0.5), because they assume independent effect sizes. Funnel plots are in `fig05_funnel_plots.png`.

**Sensitivity checks**
4. **Leave-one-study-out.** The three-level model is refit without each paper in turn.
5. **Variance strategies** A to D, plus E (without 3-fold ratios) and F (only rows with reported n).

**How to read the outputs**
- A result is "sturdy" if it keeps its sign and significance across these checks. Cover crop SOC, SOC stock, microbial biomass and bulk density are.
- Cover crop aggregate stability is not. Dropping single studies often removes its significance.

This step takes about 9 minutes.

## Step 6. Figures (`R/06_figures.R`)
| File | What it shows |
|---|---|
| `fig06a_forest_pooled_effects.png` | Every pooled cell, with CI (thick) and PI (thin) |
| `fig06b_orchard_SOC_by_practice.png`, `fig06c_orchard_covercrop_by_indicator.png` | Orchard plots: every effect size (sized by precision) plus the pooled CI and PI (Nakagawa et al. 2021 style; drawn with ggplot2 because orchaRd is not on CRAN) |
| `fig06d_caterpillar_SOC_covercrop.png` | All 905 cover crop SOC effects ranked, colored by whether the SD was reported |
| `fig06e`, `fig06f` | Bubble plots for years and depth |
| `fig06g_sensitivity_variance_strategies.png` | Estimates under each variance strategy |
| `fig06h_leave_one_out_SOC_covercrop.png` | Leave-one-study-out estimates |

## Step 7. Southern Great Plains subset (`R/07_great_plains_subset.R`)
**Two subsets**
- **Geographic.** U.S. sites at 30 to 40.5 N and 106 to 96 W. This covers:
  - the Texas Panhandle and Rolling Plains;
  - Oklahoma;
  - Kansas;
  - eastern Colorado;
  - the southern edge of Nebraska.
- **Climate analog.** Any site in a semi-arid steppe climate (Koppen BSk/BSh), or a humid subtropical or hot-summer continental climate (Cfa/Dfa) with reported MAP <= 900 mm. This matches the precipitation gradient from the Texas High Plains (about 450 mm) to central Oklahoma (about 900 mm).

**What the script does**
- Lists the sites (`sgp_geographic_site_list.csv`).
- Pools effects where possible (`sgp_subset_pooled.csv`).
- Tests analog versus other sites (`sgp_analog_vs_other_test.csv`).
- Compares site conditions (`sgp_representativeness.csv`).
- Draws a site map (`fig07a`).

**Main finding:** 94 effect sizes from 5 papers in the geographic window, and none from rangeland. Say this plainly when presenting.

## Step 8. Session info
`Rscript R/99_session_info.R` writes `session_info.txt` with R and package versions.

## Run everything
```bash
PY=python3 RS=Rscript ./run_all.sh        # about 12 minutes
```
Step 7 runs before step 6 in `run_all.sh`. The figure script does not depend on step 7, so the order makes no difference.

## Things to try next, to learn more
- Add a variance-covariance matrix for shared controls with `metafor::vcalc()` and compare the results.
- Fit the small-sample bias-corrected lnRR (Lajeunesse 2015) and compare.
- Add moderators jointly (for example depth + years) for cover crop SOC, and check collinearity.
- Add newer studies (2019 onward) for the Great Plains to see whether the regional gap can be filled.
