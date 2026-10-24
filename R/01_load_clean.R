# 01_load_clean.R
# Read SoilHealthDB_V1.xlsx, keep the practices and indicators of interest,
# reshape to one row per treatment-control comparison, screen measurement
# types, and attach moderators (climate, texture, depth, years, land use,
# southern Great Plains flags).
#
# Output: data_processed/soil_health_long.csv and several audit tables.

source("R/00_setup.R")
unlink(file.path(PATHS$logs, "run_log.txt"))
log_line("== 01_load_clean ==  ", format(Sys.time()))

raw <- read_excel(PATHS$raw_xlsx, sheet = "Sheet1", guess_max = 10000)
log_line("Rows in SoilHealthDB_V1.xlsx Sheet1: ", nrow(raw), "; columns: ", ncol(raw),
         "; papers (StudyID): ", n_distinct(raw$StudyID))

# ---------------------------------------------------------------- practices
# Conservation_Type codes are defined in the database's ConservationType sheet.
# A contrast is kept only when it isolates one practice. Rows where tillage
# also differs between treatment and control (Tillage_Diff = "Yes") are
# combined contrasts and are dropped from the cover crop group.
raw <- raw %>%
  mutate(row_id = row_number(),
         ct = str_trim(Conservation_Type),
         till_diff = toupper(str_trim(as.character(Tillage_Diff))),
         practice = case_when(
           ct %in% c("CC", "CC-NT", "CC-RT") & (is.na(till_diff) | till_diff == "NO") ~ "Cover crop",
           ct %in% c("NT", "RT") ~ "No-till / reduced tillage",
           ct %in% c("OF", "OF-CC", "OF-NT") ~ "Organic management",
           ct == "SR" ~ "Residue retention",
           ct == "AF" ~ "Agroforestry",
           ct == "RO" ~ "Crop rotation",
           TRUE ~ NA_character_),
         tillage_type = case_when(ct == "NT" ~ "No-till", ct == "RT" ~ "Reduced tillage"))

practice_audit <- raw %>% count(Conservation_Type, Tillage_Diff, practice, name = "rows") %>%
  arrange(desc(rows))
write_csv(practice_audit, file.path(PATHS$tables, "audit_practice_assignment.csv"))
log_line("Rows assigned to a practice: ", sum(!is.na(raw$practice)), " of ", nrow(raw))

# ---------------------------------------------------------------- indicators
# Each indicator: treatment and control means, their SDs. SDs in the
# database are already converted from SE or CI where papers reported those
# (see the database's SD_SE sheet).
ind_map <- tribble(
  ~indicator,            ~mc,              ~mt,              ~sc,           ~st,           ~cm,
  "SOC concentration",   "OC_C",           "OC_T",           "OC_C_SD%",    "OC_T_SD%",    "4.OC_Comments",
  "SOC stock",           "OC_Con_C",       "OC_Con_T",       "OC_C_SD",     "OC_T_SD",     "4.2OC_Comments",
  "Bulk density",        "BD_C",           "BD_T",           "BD_C_SD",     "BD_T_SD",     "3.BD_Comments",
  "Aggregate stability", "Aggre_C",        "Aggre_T",        "Aggre_C_SD",  "Aggre_T_SD",  "12.Aggre_Comments",
  "Microbial biomass C", "MBC_C",          "MBC_T",          "MBC_C_SD",    "MBC_T_SD",    "37.MBC_Comments",
  "Infiltration rate",   "Infiltration_C", "Infiltration_T", "Infiltration_C_SD", "Infiltration_T_SD", "15.Infiltration_Comments"
)

num <- function(x) suppressWarnings(as.numeric(x))
long <- pmap_dfr(ind_map, function(indicator, mc, mt, sc, st, cm) {
  raw %>% filter(!is.na(practice)) %>%
    transmute(row_id, indicator = indicator,
              m_c = num(.data[[mc]]), m_t = num(.data[[mt]]),
              sd_c = num(.data[[sc]]), sd_t = num(.data[[st]]),
              comment = as.character(.data[[cm]]))
}) %>% filter(!is.na(m_c), !is.na(m_t))
log_line("Candidate comparisons (practice assigned, both means present): ", nrow(long))

# ---------------------------------------------------------------- measurement screening
# lnRR is unit-free, so mixing g/kg and % across studies is fine as long as
# treatment and control share a unit within a row (they do in the database).
# What is not fine is mixing different quantities, so the free-text comments
# are screened:
#  * SOC concentration: drop rows flagged as total carbon (TC), which includes
#    carbonate C. SOM rows are kept: SOM = factor x SOC, the factor cancels
#    in the ratio.
#  * Aggregate stability: keep stability indices (MWD, GMD, water-stable
#    aggregates, macroaggregates, stability %); drop single size fractions
#    and the wind-erodible fraction (<0.84 mm, where higher is worse).
#  * Microbial biomass C: drop PLFA-based "total microbial biomass" and rows
#    that are actually rates.
classify_measure <- function(indicator, comment) {
  cmt <- tolower(ifelse(is.na(comment), "", comment))
  case_when(
    indicator == "SOC concentration" & str_detect(cmt, "\\btc\\b|total c") ~ "drop: total carbon (TC)",
    indicator == "Aggregate stability" & str_detect(cmt, "0\\.84|pass sieve|dmp|> ?10 ?mm|19\\.2|0\\.5-0\\.25|4\\.75-2|< ?2 but|>1mm|> ?2 ?mm|wsa\\) > 2") ~
      "drop: single size class or erodible fraction",
    indicator == "Aggregate stability" & str_detect(cmt, "mwd|weight diameter|gmd|stab|stable|macro|0\\.25|250|sieving index|^g/kg$|^g/g$|^all$") ~
      "keep: stability index",
    indicator == "Aggregate stability" & cmt == "" ~ "keep: unspecified (database default g/kg WSA)",
    indicator == "Aggregate stability" ~ "drop: unclear aggregate measure",
    indicator == "Microbial biomass C" & str_detect(cmt, "plfa|fatty|phospholipid|per 10 day") ~ "drop: not chloroform-fumigation MBC",
    TRUE ~ "keep")
}
long <- long %>% mutate(measure_decision = classify_measure(indicator, comment))
write_csv(long %>% count(indicator, comment, measure_decision, name = "rows") %>% arrange(indicator, desc(rows)),
          file.path(PATHS$tables, "audit_measurement_screening.csv"))

excl <- long %>% filter(str_starts(measure_decision, "drop")) %>% count(indicator, measure_decision, name = "rows")
log_line("Rows dropped by measurement screening: ", sum(excl$rows))
long <- long %>% filter(!str_starts(measure_decision, "drop"))

# Means must be positive for a log ratio
n_nonpos <- sum(long$m_c <= 0 | long$m_t <= 0)
log_line("Rows dropped because a mean was <= 0: ", n_nonpos)
long <- long %>% filter(m_c > 0, m_t > 0)

# ---------------------------------------------------------------- join site and design data
site <- raw %>% transmute(
  row_id, practice, tillage_type, conservation_type = ct,
  study_id = StudyID, exp_id = paste(StudyID, ExperimentID, sep = "_"),
  author = Author_F, pub_year = first_num(YearPublication),
  country = Country, site_info = SiteInfor,
  lat = num(Latitude), lon = num(Longitude),
  map_mm = first_num(MAP), mat_c = first_num(MAT),
  texture_raw = Texture, clay = num(ClayPerc), sand = num(SandPerc),
  depth_txt = as.character(SamplingDepth),
  years_raw = first_num(TimeAfterCoverCrop), duration = first_num(Duration),
  n_c_raw = first_num(No_C), n_t_raw = first_num(No_T), n_sub = first_num(NoSubsample),
  grain_crop = GrainCrop, fert_diff = Fert_Diff)

long <- long %>% left_join(site, by = "row_id")

# ---------------------------------------------------------------- sample size
# n = number of plots (replicates) per arm. If missing, use NoSubsample,
# else the median n of that indicator (flagged).
long <- long %>% mutate(
  n_c = coalesce(n_c_raw, n_t_raw, n_sub), n_t = coalesce(n_t_raw, n_c_raw, n_sub),
  n_source = ifelse(is.na(n_c_raw) & is.na(n_t_raw), ifelse(is.na(n_sub), "imputed median", "NoSubsample"), "reported"))
med_n <- long %>% group_by(indicator) %>% summarise(med = median(n_c, na.rm = TRUE))
long <- long %>% left_join(med_n, by = "indicator") %>%
  mutate(n_c = ifelse(is.na(n_c) | n_c < 2, pmax(med, 2), n_c),
         n_t = ifelse(is.na(n_t) | n_t < 2, pmax(med, 2), n_t)) %>% select(-med)

# ---------------------------------------------------------------- climate (Koppen-Geiger)
# Kottek et al. (2006) 0.5 degree grid shipped with SoilHealthDB. Each site
# takes the class of the nearest land cell (within 1 degree).
kg <- read_table(PATHS$koppen, col_types = "ddc")
names(kg) <- c("klat", "klon", "koppen")
sites_u <- long %>% distinct(lat, lon) %>% filter(!is.na(lat), !is.na(lon))
nearest_kg <- function(la, lo) {
  cand <- kg[abs(kg$klat - la) <= 1 & abs(kg$klon - lo) <= 1, ]
  if (!nrow(cand)) return(NA_character_)
  cand$koppen[which.min((cand$klat - la)^2 + ((cand$klon - lo) * cos(la * pi / 180))^2)]
}
sites_u$koppen <- map2_chr(sites_u$lat, sites_u$lon, nearest_kg)
long <- long %>% left_join(sites_u, by = c("lat", "lon")) %>%
  mutate(climate_zone = recode(substr(koppen, 1, 1), A = "Tropical (A)", B = "Arid/semi-arid (B)",
                               C = "Temperate (C)", D = "Cold/continental (D)", .default = NA_character_))

# ---------------------------------------------------------------- texture
tex <- tolower(str_trim(long$texture_raw))
long$texture <- case_when(
  tex %in% c("sand", "sandy", "loamy sand", "sandy loam") ~ "Coarse",
  tex %in% c("loam", "silt loam", "silt", "sandy clay loam") ~ "Medium",
  tex %in% c("clay", "clay loam", "silty clay", "silty clay loam", "sandy clay") ~ "Fine",
  !is.na(long$clay) & long$clay >= 35 ~ "Fine",
  !is.na(long$sand) & !is.na(long$clay) & long$sand >= 52 & long$clay < 20 ~ "Coarse",
  !is.na(long$clay) ~ "Medium",
  TRUE ~ NA_character_)

# ---------------------------------------------------------------- depth and years
dnum <- str_extract_all(long$depth_txt, "[0-9.]+")
long$depth_top <- map_dbl(dnum, ~ if (length(.x) >= 2) as.numeric(.x[1]) else NA_real_)
long$depth_bot <- map_dbl(dnum, ~ if (length(.x) >= 2) as.numeric(.x[2]) else NA_real_)
long <- long %>% mutate(
  depth_mid = (depth_top + depth_bot) / 2,
  depth_class = case_when(depth_bot <= 15 ~ "Surface (bottom <= 15 cm)",
                          depth_top < 15 & depth_bot <= 30 ~ "0-30 cm layers",
                          depth_top >= 15 | depth_bot > 30 ~ "Deeper (> 30 cm or starts >= 15 cm)"),
  years = coalesce(years_raw, duration),
  years_source = ifelse(!is.na(years_raw), "TimeAfterCoverCrop", ifelse(!is.na(duration), "Duration", NA)),
  years_class = cut(years, c(0, 3, 10, Inf), labels = c("1-3 yr", "4-10 yr", "> 10 yr"), right = TRUE))

# ---------------------------------------------------------------- land use
long$land_use <- ifelse(str_detect(tolower(paste(long$grain_crop, long$site_info)),
                                   "grass|pasture|range|prairie"), "Grassland/pasture", "Cropland")

# ---------------------------------------------------------------- southern Great Plains flags
# sgp_geo: U.S. sites inside a box that covers the southern and central Great
#   Plains (lat 30 to 40.5 N, lon 106 to 96 W: Texas Panhandle and Rolling
#   Plains, Oklahoma, Kansas, eastern Colorado, far southern Nebraska).
# sgp_climate: a climate analog anywhere in the world: semi-arid steppe
#   (Koppen BSk/BSh), or humid subtropical / hot-summer continental
#   (Cfa/Dfa) with reported MAP <= 900 mm. Sites with Koppen B and a
#   reported MAP > 900 mm are excluded.
long <- long %>% mutate(
  sgp_geo = country == "USA" & !is.na(lat) & lat >= 30 & lat <= 40.5 & lon >= -106 & lon <= -96,
  sgp_climate = (koppen %in% c("BSk", "BSh") & (is.na(map_mm) | map_mm <= 900)) |
                (koppen %in% c("Cfa", "Dfa") & !is.na(map_mm) & map_mm <= 900),
  sgp_climate = coalesce(sgp_climate, FALSE))

long <- long %>% arrange(indicator, practice, study_id, row_id) %>%
  mutate(es_id = row_number(),
         practice = factor(practice, PRACTICE_LEVELS), indicator = factor(indicator, INDICATOR_LEVELS))

write_csv(long, file.path(PATHS$processed, "soil_health_long.csv"))

avail <- long %>% group_by(indicator, practice) %>%
  summarise(k = n(), studies = n_distinct(study_id), experiments = n_distinct(exp_id),
            pct_with_sd = round(100 * mean(!is.na(sd_c) & !is.na(sd_t) & sd_c > 0 & sd_t > 0), 1),
            sgp_geo_k = sum(sgp_geo), sgp_climate_k = sum(sgp_climate), .groups = "drop") %>%
  mutate(analysed = k >= MIN_K & studies >= MIN_STUDIES)
write_csv(avail, file.path(PATHS$tables, "data_availability.csv"))
print(as.data.frame(avail))
log_line("Final comparisons: ", nrow(long), " from ", n_distinct(long$study_id), " papers; ",
         "grassland/pasture rows: ", sum(long$land_use != "Cropland"),
         "; sgp_geo rows: ", sum(long$sgp_geo), "; sgp_climate rows: ", sum(long$sgp_climate))
log_line("Koppen classes assigned: ", sum(!is.na(long$koppen)), " of ", nrow(long))
