### Cleaning data for SDM ----------------------------------------

# Packages
library(tidyr)
library(dplyr)
library(PerformanceAnalytics)
library(spatialsample)
library(sf)
library(purrr)

# -------------------------------------------------------------------------
# Read datasets
# -------------------------------------------------------------------------
mapbiomas <- read.csv("data/a_chamek_ter_mammals_lulc_Oct2022.csv")

climate <- read.csv("data/a_chamek_ter_mammals_climate_Oct2022.csv") %>%
  dplyr::select(row_code, bio13_precip_wettest_month, cmi_min)

amazon_basin_pnts <- read.csv("data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

# -------------------------------------------------------------------------
# Relabel MapBiomas classes + drop undefined classes
# -------------------------------------------------------------------------
# Inspect available classes if needed:
# unique(mapbiomas$class)

mapbiomas_clean <- mapbiomas %>%
  mutate(
    class = as.integer(class),
    class = recode(
      class,
      `3`  = "forest_formation",
      `14` = "farming",
      .default = "other"
    )
  ) %>%
  filter(class != "other")

# -------------------------------------------------------------------------
# Wide format: one row per row_code-year, columns for each class area
# -------------------------------------------------------------------------
mapbiomas_wide <- mapbiomas_clean %>%
  dplyr::select(row_code, year, class, area) %>%
  group_by(row_code, year, class) %>%
  summarise(area = sum(area, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = class,
    values_from = area,
    values_fill = 0
  )

# -------------------------------------------------------------------------
# Summaries per point across years + difference over study period
# -------------------------------------------------------------------------
mapbiomas_mean_diff <- mapbiomas_wide %>%
  group_by(row_code) %>%
  summarise(
    # mean area per class
    mean_forest  = mean(forest_formation, na.rm = TRUE),
    mean_farming = mean(farming, na.rm = TRUE),
    
    # difference over study period (requires both years)
    diff_forest_formation =
      forest_formation[year == 2020] - forest_formation[year == 2001],
    
    .groups = "drop"
  ) %>%
  # drop rows where either year is missing
  filter(lengths(diff_forest_formation) == 1) %>%
  mutate(diff_forest_formation = as.numeric(diff_forest_formation))

# Save cleaned LULC covariates
write.csv(
  mapbiomas_mean_diff,
  "data/a_chamek_ter_mammals_lulc_cleaned_Oct2022.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Merge response + climate + LULC covariates
# -------------------------------------------------------------------------
covariates <- inner_join(climate, mapbiomas_mean_diff, by = "row_code")

data0 <- inner_join(amazon_basin_pnts, covariates, by = "row_code")

# -------------------------------------------------------------------------
# Multicollinearity (Pearson correlation)
# -------------------------------------------------------------------------
pred_df <- data0 %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-row_code, -lon, -lat, -presence)

corr <- abs(cor(pred_df, use = "pairwise.complete.obs"))
corr

chart.Correlation(
  pred_df,
  histogram = TRUE,
  method = "pearson"
)

# -------------------------------------------------------------------------
# Spatial cross-validation folds
# -------------------------------------------------------------------------
n_folds <- 3

data0_sf <- st_as_sf(
  data0,
  coords = c("lon", "lat"),
  crs = 4326
)

set.seed(99)
clusters <- spatial_block_cv(
  data0_sf,
  method = "random",
  n = 30,
  relevant_only = TRUE,
  v = n_folds
)

splits_df <- map_dfr(seq_len(n_folds), function(i) {
  assessment(clusters$splits[[i]]) %>%
    st_drop_geometry() %>%
    dplyr::select(row_code) %>%
    mutate(fold = i)
})

analysis_data <- inner_join(data0, splits_df, by = "row_code")

# Sanity check: fold balance
table(analysis_data$fold, analysis_data$presence)

# Save final analysis dataset
write.csv(
  analysis_data,
  "data/a_chamek_ter_mammals_finalData_Oct22.csv",
  row.names = FALSE
)
