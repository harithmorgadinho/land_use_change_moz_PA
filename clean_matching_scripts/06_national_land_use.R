# ============================================================
# 06. NATIONAL LAND-USE CHANGE
# Mozambique: tree cover, cropland and built-up surface
# ============================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)

terraOptions(progress = 1, memfrac = 0.7)

boundary_file <- "data/processed/mozambique_boundary.gpkg"
boundary_file <- "/Users/gdt366/Library/CloudStorage/Dropbox/gee_deforestation_mozambique/moz_admin_boundaries.shp/moz_admin0.shp"

treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
lossyear_file  <- "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"

crop_dir  <- "data/raw/glad_cropland"
built_dir <- "data/raw/ghsl"

output_dir <- "outputs/national"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Mozambique boundary and land area
# ------------------------------------------------------------

moz <- st_read(boundary_file, quiet = TRUE) %>%
  st_make_valid()

# Equal-area calculation used only for the national land denominator.
moz_equal_area <- st_transform(moz, "EPSG:6933")

mozambique_land_ha <- sum(as.numeric(st_area(moz_equal_area))) / 10000


raster_sum <- function(x) {
  as.numeric(global(x, "sum", na.rm = TRUE)[1, 1])
}

mask_to_mozambique <- function(r, moz_sf) {
  moz_v <- vect(st_transform(moz_sf, crs(r)))
  mask(crop(r, moz_v), moz_v)
}

find_year_file <- function(files, year) {
  f <- files[grepl(as.character(year), basename(files))]
  if (length(f) != 1) {
    stop("Expected exactly one file for ", year, "; found ", length(f))
  }
  f
}

# ============================================================
# TREE COVER, 2000-2025
# ============================================================

treecover <- rast(treecover_file)
lossyear  <- rast(lossyear_file)

treecover <- mask_to_mozambique(treecover, moz)
lossyear  <- mask_to_mozambique(lossyear, moz)

cell_area_ha <- cellSize(treecover, unit = "ha")

# Continuous tree-cover-equivalent area in 2000.
tree_area_2000 <- cell_area_ha * (treecover / 100)
names(tree_area_2000) <- "tree_area_ha"

tree_baseline_ha <- raster_sum(tree_area_2000)

tree_annual <- bind_rows(lapply(2001:2025, function(y) {

  loss_max <- global(lossyear, "max", na.rm = TRUE)[1, 1]
  code_y <- if (loss_max <= 25) y - 2000 else y

  loss_y <- ifel(lossyear == code_y, tree_area_2000, 0)

  tibble(
    year = y,
    annual_loss_ha = raster_sum(loss_y)
  )
})) %>%
  arrange(year) %>%
  mutate(cumulative_loss_ha = cumsum(annual_loss_ha))

tree_loss_ha <- sum(tree_annual$annual_loss_ha)
tree_final_ha <- tree_baseline_ha - tree_loss_ha
tree_loss_percent <- 100 * tree_loss_ha / tree_baseline_ha

write.csv(
  tree_annual,
  file.path(output_dir, "Mozambique_treecover_annual_loss_2001_2025.csv"),
  row.names = FALSE
)

# ============================================================
# CROPLAND, 2015-2024
# ============================================================

crop_years <- 2015:2024

crop_files <- list.files(
  crop_dir,
  pattern = "^Mozambique_GLAD_cropland_.*\\.tif$",
  full.names = TRUE
)

invisible(lapply(crop_years, function(y) find_year_file(crop_files, y)))

crop_annual <- bind_rows(lapply(crop_years, function(y) {

  message("Cropland: ", y)

  r <- rast(find_year_file(crop_files, y))
  r <- mask_to_mozambique(r, moz)

  rmax <- global(r, "max", na.rm = TRUE)[1, 1]

  crop_binary <- if (rmax <= 1.5) {
    ifel(r > 0, 1, 0)
  } else {
    ifel(r == 100, 1, 0)
  }

  crop_area_ha <- crop_binary * cellSize(r, unit = "ha")

  tibble(
    year = y,
    cropland_ha = raster_sum(crop_area_ha)
  )
})) %>%
  arrange(year) %>%
  mutate(
    change_from_2015_ha = cropland_ha - cropland_ha[year == 2015],
    change_from_2015_percent =
      100 * change_from_2015_ha / cropland_ha[year == 2015],
    percent_land = 100 * cropland_ha / mozambique_land_ha,
    change_from_2015_pp = percent_land - percent_land[year == 2015]
  )

crop_2015_ha <- crop_annual$cropland_ha[crop_annual$year == 2015]
crop_2024_ha <- crop_annual$cropland_ha[crop_annual$year == 2024]
crop_change_ha <- crop_2024_ha - crop_2015_ha
crop_relative_change <- 100 * crop_change_ha / crop_2015_ha

write.csv(
  crop_annual,
  file.path(output_dir, "Mozambique_cropland_annual_2015_2024.csv"),
  row.names = FALSE
)

# ============================================================
# BUILT-UP SURFACE, 2000-2025
# ============================================================

built_years <- c(2000, 2005, 2010, 2015, 2020, 2025)

built_files <- list.files(
  built_dir,
  pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
  full.names = TRUE
)

invisible(lapply(built_years, function(y) find_year_file(built_files, y)))

built_temporal <- bind_rows(lapply(built_years, function(y) {

  message("Built-up: ", y)

  r <- rast(find_year_file(built_files, y))[["built_surface"]]
  r <- mask_to_mozambique(r, moz)

  r <- ifel(r >= 0 & r <= 10000, r, NA)

  built_m2 <- raster_sum(r)
  built_ha <- built_m2 / 10000

  tibble(
    year = y,
    built_surface_m2 = built_m2,
    built_surface_ha = built_ha,
    percent_land = 100 * built_ha / mozambique_land_ha
  )
})) %>%
  arrange(year) %>%
  mutate(
    change_from_2000_ha =
      built_surface_ha - built_surface_ha[year == 2000],
    change_from_2000_percent =
      100 * change_from_2000_ha / built_surface_ha[year == 2000],
    change_from_2000_pp =
      percent_land - percent_land[year == 2000]
  )

built_2000_ha <- built_temporal$built_surface_ha[built_temporal$year == 2000]
built_2025_ha <- built_temporal$built_surface_ha[built_temporal$year == 2025]
built_change_ha <- built_2025_ha - built_2000_ha
built_relative_change <- 100 * built_change_ha / built_2000_ha

write.csv(
  built_temporal,
  file.path(output_dir, "Mozambique_builtup_temporal_2000_2025.csv"),
  row.names = FALSE
)

# ============================================================
# TABLE S1: NATIONAL LAND-USE CHANGE SUMMARY
# ============================================================

summary_table <- tibble(
  Threat = c("Tree cover", "Cropland", "Built-up surface"),
  `Baseline year` = c(2000, 2015, 2000),
  `Final year` = c(2025, 2024, 2025),

  `Baseline extent (ha)` = c(
    tree_baseline_ha,
    crop_2015_ha,
    built_2000_ha
  ),

  `Final extent (ha)` = c(
    tree_final_ha,
    crop_2024_ha,
    built_2025_ha
  ),

  `Loss / expansion (ha)` = c(
    tree_loss_ha,
    crop_change_ha,
    built_change_ha
  ),

  `Loss / expansion (%)` = c(
    tree_loss_percent,
    crop_relative_change,
    built_relative_change
  )
) %>%
  mutate(
    `Baseline extent (Mha)` = `Baseline extent (ha)` / 1e6,
    `Final extent (Mha)` = `Final extent (ha)` / 1e6,
    `Loss / expansion (Mha)` = `Loss / expansion (ha)` / 1e6,
    `Baseline land share (%)` =
      100 * `Baseline extent (ha)` / mozambique_land_ha,
    `Final land share (%)` =
      100 * `Final extent (ha)` / mozambique_land_ha,
    `Change in land share (pp)` =
      `Final land share (%)` - `Baseline land share (%)`
  ) %>%
  transmute(
    Threat,
    `Baseline year`,
    `Final year`,
    `Baseline extent (ha)` = round(`Baseline extent (ha)`, 0),
    `Baseline extent (Mha)` = round(`Baseline extent (Mha)`, 3),
    `Final extent (ha)` = round(`Final extent (ha)`, 0),
    `Final extent (Mha)` = round(`Final extent (Mha)`, 3),
    `Loss / expansion (ha)` = round(`Loss / expansion (ha)`, 0),
    `Loss / expansion (Mha)` = round(`Loss / expansion (Mha)`, 3),
    `Loss / expansion (%)` = round(`Loss / expansion (%)`, 2),
    `Baseline land share (%)` = round(`Baseline land share (%)`, 3),
    `Final land share (%)` = round(`Final land share (%)`, 3),
    `Change in land share (pp)` = round(`Change in land share (pp)`, 3)
  )

write.csv(
  summary_table,
  file.path(output_dir, "Table_S1_national_land_use_change_summary.csv"),
  row.names = FALSE
)

write.csv(
  check_table,
  file.path(output_dir, "national_reproducibility_checks.csv"),
  row.names = FALSE
)

