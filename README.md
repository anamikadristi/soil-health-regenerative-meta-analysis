# Meta-analysis: conservation and regenerative practices and soil health indicators

**Author:** Anamika Dristi, PhD (Louisiana State University, 2026)
**Status:** Self-directed training project done after the PhD. It is a **re-analysis of a published, public database**, not new field data, and it is not published or peer reviewed.

## Question
Across published field experiments, how much do cover crops, no-till or reduced tillage, organic management, residue retention and agroforestry change soil organic carbon (SOC) and other soil health indicators, compared with conventional controls? How much do the effects vary, what explains the variation, and what does the evidence say about semi-arid and subhumid farmland like the southern Great Plains?

| Step | Script | What it does |
|---|---|---|
| 0 | `python/00_fetch_soilhealthdb.py` | Downloads the database files from figshare (about 5 MB; reads only the needed files out of the 128 MB zip) |
| 1 | `R/01_load_clean.R` | Selects practices and indicators, screens measurement types, adds moderators (Koppen climate, texture, depth, years, land use, Great Plains flags) |
| 2 | `R/02_effect_sizes.R` | Log response ratio (lnRR) with `metafor::escalc(measure = "ROM")`; imputes missing variances from the mean coefficient of variation |
| 3 | `R/03_random_effects.R` | Three-level random-effects models (effect sizes within experiments within studies) per practice x indicator |
| 4 | `R/04_moderators.R` | Meta-regression on climate, precipitation, temperature, texture, years, depth, land use |
| 5 | `R/05_bias_sensitivity.R` | Funnel plots, Egger tests (multilevel and study-level), trim-and-fill, time-lag, leave-one-study-out, variance-imputation sensitivity |
| 6 | `R/06_figures.R` | Forest, orchard, caterpillar, bubble and sensitivity plots |
| 7 | `R/07_great_plains_subset.R` | Southern Great Plains subsets (geographic and climate analog), site map |

## Data source
- **Paper:** Jian, J., Du, X., & Stewart, R. D. (2020). A database for global soil health assessment. *Scientific Data*, 7, 16. https://doi.org/10.1038/s41597-020-0356-3
- **Data:** Jian, J., Stewart, R. D., & Du, X. SoilHealthDB (version 8). figshare. https://doi.org/10.6084/m9.figshare.8292176.v8 (file `SoilHealthDB.zip`, md5 68719dd87d60ddd7b97e2901a83821bc). The record is released under **CC0 1.0**. The USDA Ag Data Commons copy of the same dataset lists CC BY 4.0. Either way, cite the paper above.
- **Files used:** `SoilHealthDB_V1.xlsx` (5,907 rows x 270 columns, 321 papers; the version described in the paper), the CSV copy, and `Koeppen-Geiger-ASCII.txt` (0.5 degree Koppen-Geiger grid, Kottek et al. 2006, shipped in the same zip). SHA-256 checksums are in `data_raw/download_log.txt`.
- Each database row is one treatment-versus-control comparison with means, SDs (when the paper reported them) and site information.

## How to run
```bash
# R 4.3 or newer (tested with R 4.5.2)
Rscript install_packages.R
# full pipeline from the project folder; step 0 needs internet
PY=python3 RS=Rscript ./run_all.sh
# if data_raw/ is already filled:
SKIP_DOWNLOAD=1 ./run_all.sh
```
The Python step uses only the standard library. Run all R scripts from the project root (for example `Rscript R/03_random_effects.R`). The full run takes about 12 minutes; most of it is the leave-one-study-out refits in step 5. Package versions are in `session_info.txt`. See `STEP_BY_STEP_GUIDE.md` for a walk-through.

## Outputs
- `outputs/tables/`: pooled effects, moderator tests, bias and sensitivity tables, Great Plains subset, and audit tables that show every data-screening decision (CSV)
- `outputs/figures/`: PNG, 300 dpi
- `outputs/logs/`: run log and full model printouts
- `data_processed/`: cleaned comparison table and effect sizes (CSV)
- Git tracks the scripts, documents and outputs. It does not track `data_raw/` (except the download log with checksums) or `data_processed/`, because the scripts re-create both.

## Key results (details in `RESULTS_SUMMARY.md`)
Percent change = (exp(lnRR) - 1) x 100, with 95% CI; k = effect sizes; "studies" = papers.
1. **Cover crops** raised SOC concentration by **12.1% [7.8, 16.5]** (k = 905, 62 studies) and SOC stock by **8.4% [5.5, 11.5]** (k = 448, 21 studies). Microbial biomass C rose 22.1% [8.4, 37.4] (k = 176, 20 studies). Bulk density fell slightly, -1.8% [-3.0, -0.6]. Aggregate stability rose 8.6% [0.4, 17.4], but this result did not hold in several sensitivity checks.
2. **Organic management** raised SOC concentration by 17.8% [12.1, 23.8] (k = 316, 74 studies), with very high heterogeneity (I2 = 99.5%).
3. **No-till / reduced tillage** and **residue retention** had only 5 to 9 studies each in this database. Their SOC estimates (+13.0% and +18.2%) have confidence intervals that include zero.
4. Practices did not differ significantly from each other for SOC concentration (F(4, 151) = 1.81, p = 0.13).
5. The cover crop SOC effect was smaller in deeper soil layers and grew slowly with years under the practice (about +0.5% per year). Climate zone, texture and land use did not explain the variation. None of the moderator tests for cover crops stayed significant after Holm correction.
6. **Southern Great Plains:** only 94 comparisons from 5 papers (6 sites in TX, KS, CO and NE) fall inside the U.S. Great Plains window. They suggest a cover crop SOC gain of 17.7% [6.4, 30.2], but that estimate rests on 3 papers. A wider climate-analog subset (semi-arid or subhumid, MAP <= 900 mm; 26 papers worldwide) gave cover crop SOC +12.0% [4.2, 20.4], close to the global estimate. **No rangeland studies** are in the database.
7. Imputing missing variances was necessary (only about 22% of effect sizes have reported SDs). The SOC, bulk density and microbial biomass results changed little across four variance strategies.

## Documentation
- `RESULTS_SUMMARY.md`: all numbers, interpretation and limitations
- `STEP_BY_STEP_GUIDE.md`: how to run each step and what it means

## License

Code is released under the MIT License (see LICENSE). The Soil Health Database is CC0 and keeps its original terms.
