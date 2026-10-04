# ============================================================
# 04. MATCHED LAND-USE EFFECTS
# Mozambique protected areas
# ============================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(MatchIt)

matching_data_dir <- "data/processed/matching"
matching_output_dir <- "outputs/matching"
output_dir <- "outputs/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Raw outcome data
hansen_treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
hansen_treecover_file <- "Mozambique_treecover2000_percent.tif"

hansen_loss_file <- "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"
hansen_loss_file <- "Mozambique_forest_loss_year_2001_2025.tif"

crop_2015_file <- "data/raw/glad_cropland/Mozambique_GLAD_cropland_2015.tif"
crop_2015_file <- "glad_cropland/Mozambique_GLAD_cropland_2015.tif"

crop_2024_file <- "data/raw/glad_cropland/Mozambique_GLAD_cropland_2024.tif"
crop_2024_file <- "glad_cropland/Mozambique_GLAD_cropland_2024.tif"

built_dir <- "data/raw/ghsl"
built_dir <- "ghsl"

# Final balance rules used in the manuscript.
main_balance_threshold <- 0.25
strict_balance_threshold <- 0.10

# ------------------------------------------------------------
# Load grid, matches and balance tables
# ------------------------------------------------------------

grid_1000 <- readRDS(file.path(matching_data_dir, "matching_grid_1000m.rds"))

m_tree <- readRDS(file.path(matching_output_dir, "matching_tree_1000m_final.rds"))
m_crop <- readRDS(file.path(matching_output_dir, "matching_crop_1000m_final.rds"))
m_built <- readRDS(file.path(matching_output_dir, "matching_built_1000m_final.rds"))

balance_tree <- read.csv(file.path(matching_output_dir, "matching_tree_balance_final.csv"))
balance_crop <- read.csv(file.path(matching_output_dir, "matching_crop_balance_final.csv"))
balance_built <- read.csv(file.path(matching_output_dir, "matching_built_balance_final.csv"))

add_balance_flags <- function(x) {
  x %>% mutate(
    balance_025 = round(max_abs_SMD, 2) <= main_balance_threshold,
    balance_010 = max_abs_SMD <= strict_balance_threshold
  )
}

# ------------------------------------------------------------
# Common approximately 1-km raster/grid geometry
# ------------------------------------------------------------

res_1000 <- 1000 / 111320
r_1000 <- rast(
  ext(vect(grid_1000)),
  resolution = res_1000,
  crs = crs(vect(grid_1000))
)

grid_centres <- st_centroid(grid_1000)

# ------------------------------------------------------------
# Helper: weighted matched means and inside - control difference
# ------------------------------------------------------------

calculate_matched_effects <- function(matches, outcome_table, outcome_col, balance_table,
                                      inside_name, control_name) {

  effects <- bind_rows(lapply(names(matches), function(park) {

    md <- MatchIt::match.data(matches[[park]]) %>%
      left_join(outcome_table %>% select(cell_id, all_of(outcome_col)), by = "cell_id")

    y <- md[[outcome_col]]
    if (all(is.na(y))) return(NULL)

    inside <- weighted.mean(
      y[md$treated == 1],
      md$weights[md$treated == 1],
      na.rm = TRUE
    )

    control <- weighted.mean(
      y[md$treated == 0],
      md$weights[md$treated == 0],
      na.rm = TRUE
    )

    out <- data.frame(
      NAME_ENG = park,
      inside = inside,
      control = control,
      effect_pp = inside - control
    )

    names(out)[names(out) == "inside"] <- inside_name
    names(out)[names(out) == "control"] <- control_name
    out
  }))

  effects %>%
    left_join(
      balance_table %>% select(NAME_ENG, mean_abs_SMD, max_abs_SMD),
      by = "NAME_ENG"
    ) %>%
    add_balance_flags()
}

# ============================================================
# TREE-COVER LOSS, 2000-2025
# ============================================================

treecover <- rast(hansen_treecover_file)
lossyear <- rast(hansen_loss_file)

# Proportion of baseline tree cover lost by 2025.
tree_loss_prop <- ifel(lossyear > 0, treecover / 100, 0)
tree_loss_1000 <- resample(tree_loss_prop, r_1000, method = "average")
tree_vals <- terra::extract(tree_loss_1000, vect(grid_centres))

tree_cell_values <- st_drop_geometry(grid_1000) %>%
  select(cell_id, NAME_ENG, treated) %>%
  mutate(tree_loss_pp = tree_vals[, 2] * 100)

tree_effects <- calculate_matched_effects(
  m_tree,
  tree_cell_values,
  "tree_loss_pp",
  balance_tree,
  "loss_inside_pp",
  "loss_control_pp"
)

# ============================================================
# CROPLAND CHANGE, 2015-2024
# ============================================================

crop_2015 <- rast(crop_2015_file)
crop_2024 <- rast(crop_2024_file)
crop_change <- crop_2024 - crop_2015
crop_change_1000 <- resample(crop_change, r_1000, method = "average")
crop_vals <- terra::extract(crop_change_1000, vect(grid_centres))

crop_cell_values <- st_drop_geometry(grid_1000) %>%
  select(cell_id, NAME_ENG, treated) %>%
  mutate(crop_change_pp = crop_vals[, 2])

crop_effects <- calculate_matched_effects(
  m_crop,
  crop_cell_values,
  "crop_change_pp",
  balance_crop,
  "crop_inside_pp",
  "crop_control_pp"
)

# ============================================================
# BUILT-UP CHANGE, 2000-2025
# ============================================================

built_files <- sort(list.files(
  built_dir,
  pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
  full.names = TRUE
))


built_2000 <- rast(built_files[1])[[1]]
built_2025 <- rast(built_files[length(built_files)])[[1]]

# Remove invalid values. Native GHSL cells are 100 x 100 m,
# so valid built surface is 0-10,000 m2 per cell.
built_2000 <- ifel(built_2000 > 10000, NA, built_2000)
built_2025 <- ifel(built_2025 > 10000, NA, built_2025)

built_2000_pct <- built_2000 / 10000 * 100
built_2025_pct <- built_2025 / 10000 * 100
built_change_100m <- built_2025_pct - built_2000_pct

# Aggregate the native 100-m product to approximately 1 km first,
# then project to the WGS84 matching raster used for extraction.
built_change_1km_moll <- aggregate(
  built_change_100m,
  fact = 10,
  fun = "mean",
  na.rm = TRUE
)

built_change_1km_match <- project(
  built_change_1km_moll,
  r_1000,
  method = "bilinear"
)

built_vals <- terra::extract(built_change_1km_match, vect(grid_centres))

built_cell_values <- st_drop_geometry(grid_1000) %>%
  select(cell_id, NAME_ENG, treated) %>%
  mutate(built_change_pp = built_vals[, 2])

built_effects <- calculate_matched_effects(
  m_built,
  built_cell_values,
  "built_change_pp",
  balance_built,
  "built_inside_pp",
  "built_control_pp"
)

# ------------------------------------------------------------
# Save pressure-specific effects and extracted cell outcomes
# ------------------------------------------------------------

write.csv(tree_effects, file.path(output_dir, "PA_matched_effects_treecover_2000_2025.csv"), row.names = FALSE)
write.csv(crop_effects, file.path(output_dir, "PA_matched_effects_cropland_2015_2024.csv"), row.names = FALSE)
write.csv(built_effects, file.path(output_dir, "PA_matched_effects_built_2000_2025.csv"), row.names = FALSE)

saveRDS(tree_cell_values, file.path(output_dir, "Hansen_tree_loss_2000_2025_matching_cell_values.rds"))
saveRDS(crop_cell_values, file.path(output_dir, "GLAD_crop_change_2015_2024_matching_cell_values.rds"))
saveRDS(built_cell_values, file.path(output_dir, "GHSL_built_change_2000_2025_matching_cell_values.rds"))

# ------------------------------------------------------------
# One combined table for Figure 3 and supplementary material
# All 13 PAs are retained; unsuccessful matches are represented by NA.
# ------------------------------------------------------------

all_parks <- c(
  "Banhine", "Chimanimani", "Gilé", "Gorongosa", "Limpopo",
  "Magoe", "Maputo", "Marromeu", "Niassa", "Pomene", "Quirimbas",
  "Serra da Gorongosa", "Zinave"
)

format_effects <- function(x, pressure) {
  x %>% transmute(
    NAME_ENG,
    pressure = pressure,
    effect_pp,
    mean_abs_SMD,
    max_abs_SMD,
    adequate_balance = balance_025,
    strict_balance = balance_010
  )
}

effects_all <- bind_rows(
  format_effects(tree_effects, "a. Tree-cover loss"),
  format_effects(crop_effects, "b. Cropland expansion"),
  format_effects(built_effects, "c. Built-up expansion")
) %>%
  complete(
    NAME_ENG = all_parks,
    pressure = c(
      "a. Tree-cover loss",
      "b. Cropland expansion",
      "c. Built-up expansion"
    )
  ) %>%
  mutate(
    status = case_when(
      is.na(effect_pp) ~ "Matching unsuccessful",
      adequate_balance ~ "Adequate balance",
      TRUE ~ "Inadequate balance"
    )
  )

write.csv(
  effects_all,
  file.path(output_dir, "PA_matching_effects_all_pressures.csv"),
  row.names = FALSE
)

saveRDS(
  effects_all,
  file.path(output_dir, "PA_matching_effects_all_pressures.rds")
)


