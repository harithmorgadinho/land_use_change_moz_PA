# ============================================================
# 05. TEMPORAL LAND-USE CHANGE WITHIN AND AROUND PROTECTED AREAS
# Mozambique protected areas
# ============================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)

load('parks_core.Rdata')
parks_file <- "data/processed/parks_core.Rdata"
rings_file <- "data/processed/Mozambique_PA_equal_area_surroundings_nationally_clipped.shp"

rings_file = 'Mozambique_PA_equal_area_surroundings_nationally_clipped.shp'

treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
lossyear_file  <- "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"

crop_dir  <- "data/raw/glad_cropland"
built_dir <- "data/raw/ghsl"

output_dir <- "outputs/temporal"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Protected areas and adjacent equal-area comparison landscapes
# ------------------------------------------------------------

load(parks_file) # loads object: parks_core

parks_core <- parks_core %>%
  filter(
    DESIG %in% c(
      "Parque Nacional",
      "National Park",
      "Reserva Nacional",
      "Reserva Especial"
    )
  ) %>%
  st_make_valid()

rings_all <- st_read(rings_file, quiet = TRUE) %>%
  rename(NAME_ENG = park_name) %>%
  st_make_valid()

# Retain only rings corresponding to the protected-area units analysed.
rings_all <- rings_all %>%
  filter(NAME_ENG %in% parks_core$NAME_ENG)

stopifnot(
  nrow(parks_core) == 13,
  nrow(rings_all) == 13,
  setequal(parks_core$NAME_ENG, rings_all$NAME_ENG)
)

# ------------------------------------------------------------
# Region areas
# ------------------------------------------------------------

moll_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

park_areas <- parks_core %>%
  st_transform(moll_crs) %>%
  mutate(
    region_area_ha = as.numeric(st_area(geometry)) / 10000,
    location = "Inside"
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(NAME_ENG, location, region_area_ha)

ring_areas <- rings_all %>%
  st_transform(moll_crs) %>%
  mutate(
    region_area_ha = as.numeric(st_area(geometry)) / 10000,
    location = "Outside"
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(NAME_ENG, location, region_area_ha)

region_areas <- bind_rows(park_areas, ring_areas)

get_region_area <- function(park_name, location_name) {
  x <- region_areas %>%
    filter(NAME_ENG == park_name, location == location_name) %>%
    pull(region_area_ha)

  stopifnot(length(x) == 1)
  x
}

# ============================================================
# TREE-COVER LOSS, 2000-2025
# ============================================================

treecover <- rast(treecover_file)
lossyear  <- rast(lossyear_file)

# Convert percentage tree cover in 2000 to tree-cover-equivalent hectares.
cell_area_ha <- cellSize(treecover, unit = "ha")
tree_area_2000 <- cell_area_ha * (treecover / 100)
names(tree_area_2000) <- "tree_area_ha"

hansen <- c(tree_area_2000, lossyear)
names(hansen) <- c("tree_area_ha", "loss_year")

tree_loss_series <- function(polygon, park_name, location_name, region_area_ha) {

  polygon_v <- vect(polygon)
  hansen_crop <- crop(hansen, polygon_v)

  values <- terra::extract(
    hansen_crop,
    polygon_v,
    weights = TRUE
  ) %>%
    mutate(tree_area_weighted_ha = tree_area_ha * weight)

  annual_loss <- values %>%
    filter(loss_year >= 2001, loss_year <= 2025) %>%
    group_by(loss_year) %>%
    summarise(
      annual_loss_ha = sum(tree_area_weighted_ha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(year = loss_year)

  annual_loss <- tibble(year = 2001:2025) %>%
    left_join(annual_loss, by = "year") %>%
    mutate(annual_loss_ha = replace_na(annual_loss_ha, 0))

  bind_rows(
    tibble(year = 2000, annual_loss_ha = 0, cumulative_loss_ha = 0),
    annual_loss %>%
      arrange(year) %>%
      mutate(cumulative_loss_ha = cumsum(annual_loss_ha))
  ) %>%
    mutate(
      NAME_ENG = park_name,
      location = location_name,
      region_area_ha = region_area_ha,
      change_pp = cumulative_loss_ha / region_area_ha * 100,
      threat = "Tree-cover loss"
    ) %>%
    dplyr::select(
      NAME_ENG, location, year, region_area_ha,
      annual_loss_ha, cumulative_loss_ha, change_pp, threat
    )
}

run_for_regions <- function(polygons, location_name, fun) {
  bind_rows(lapply(seq_len(nrow(polygons)), function(i) {
    park_name <- polygons$NAME_ENG[i]
    message(location_name, ": ", park_name)

    fun(
      polygon = polygons[i, ],
      park_name = park_name,
      location_name = location_name,
      region_area_ha = get_region_area(park_name, location_name)
    )
  }))
}

tree_inside  <- run_for_regions(parks_core, "Inside", tree_loss_series)
tree_outside <- run_for_regions(rings_all, "Outside", tree_loss_series)

tree_ts <- bind_rows(tree_inside, tree_outside) %>%
  arrange(NAME_ENG, location, year)

saveRDS(
  tree_ts,
  file.path(output_dir, "PA_treecover_temporal_inside_outside.rds")
)
write.csv(
  tree_ts,
  file.path(output_dir, "PA_treecover_temporal_inside_outside.csv"),
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

find_year_file <- function(files, year) {
  f <- files[grepl(as.character(year), basename(files))]
  if (length(f) != 1) {
    stop("Expected exactly one file for ", year, "; found ", length(f))
  }
  f
}

invisible(lapply(crop_years, function(y) find_year_file(crop_files, y)))

cropland_series <- function(polygon, park_name, location_name, region_area_ha) {

  polygon_v <- vect(polygon)

  result <- bind_rows(lapply(crop_years, function(y) {

    message(park_name, " | ", location_name, " | cropland ", y)

    # Load raster
    crop_y <- rast(find_year_file(crop_files, y))

    # FIRST crop to the protected area / surrounding polygon
    crop_y <- crop(crop_y, polygon_v)

    # Calculate cell area only for this small cropped raster
    cell_area_y <- cellSize(crop_y, unit = "ha")

    # Cropland area per cell
    crop_area_y <- cell_area_y * (crop_y / 100)
    names(crop_area_y) <- "crop_area_ha"

    # Weighted extraction for boundary cells
    values_y <- terra::extract(
      crop_area_y,
      polygon_v,
      weights = TRUE
    )

    total_crop_ha <- sum(
      values_y$crop_area_ha * values_y$weight,
      na.rm = TRUE
    )

    # Remove intermediate rasters before moving to next year
    rm(crop_y, cell_area_y, crop_area_y, values_y)
    gc()

    tibble(
      year = y,
      crop_area_ha = total_crop_ha
    )
  }))

  baseline_2015 <- result %>%
    filter(year == 2015) %>%
    pull(crop_area_ha)

  result %>%
    mutate(
      crop_change_ha = crop_area_ha - baseline_2015,
      change_pp = crop_change_ha / region_area_ha * 100,
      NAME_ENG = park_name,
      location = location_name,
      region_area_ha = region_area_ha,
      threat = "Cropland"
    ) %>%
    dplyr::select(
      NAME_ENG, location, year, region_area_ha,
      crop_area_ha, crop_change_ha, change_pp, threat
    )
}

crop_inside  <- run_for_regions(parks_core, "Inside", cropland_series)
crop_outside <- run_for_regions(rings_all, "Outside", cropland_series)

crop_ts <- bind_rows(crop_inside, crop_outside) %>%
  arrange(NAME_ENG, location, year)

saveRDS(
  crop_ts,
  file.path(output_dir, "PA_cropland_temporal_inside_outside.rds")
)
write.csv(
  crop_ts,
  file.path(output_dir, "PA_cropland_temporal_inside_outside.csv"),
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

builtup_series <- function(polygon, park_name, location_name, region_area_ha) {

  result <- bind_rows(lapply(built_years, function(y) {

    message(park_name, " | ", location_name, " | built-up ", y)

    built_y <- rast(find_year_file(built_files, y))[["built_surface"]]

    polygon_ghsl <- st_transform(polygon, crs(built_y))
    polygon_v <- vect(polygon_ghsl)

    values_y <- terra::extract(
      crop(built_y, polygon_v),
      polygon_v,
      weights = TRUE
    ) %>%
      # Retain the same fill-value treatment used in the final analysis.
      filter(is.na(built_surface) | built_surface <= 10000)

    tibble(
      year = y,
      built_area_ha = sum(
        values_y$built_surface * values_y$weight,
        na.rm = TRUE
      ) / 10000
    )
  }))

  baseline_2000 <- result %>%
    filter(year == 2000) %>%
    pull(built_area_ha)

  result %>%
    mutate(
      built_change_ha = built_area_ha - baseline_2000,
      change_pp = built_change_ha / region_area_ha * 100,
      NAME_ENG = park_name,
      location = location_name,
      region_area_ha = region_area_ha,
      threat = "Built-up"
    ) %>%
    dplyr::select(
      NAME_ENG, location, year, region_area_ha,
      built_area_ha, built_change_ha, change_pp, threat
    )
}

built_inside  <- run_for_regions(parks_core, "Inside", builtup_series)
built_outside <- run_for_regions(rings_all, "Outside", builtup_series)

built_ts <- bind_rows(built_inside, built_outside) %>%
  arrange(NAME_ENG, location, year)

saveRDS(
  built_ts,
  file.path(output_dir, "PA_builtup_temporal_inside_outside.rds")
)
write.csv(
  built_ts,
  file.path(output_dir, "PA_builtup_temporal_inside_outside.csv"),
  row.names = FALSE
)

