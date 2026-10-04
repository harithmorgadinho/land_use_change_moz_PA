# ============================================================
# 02. PREPARE MATCHING DATA
# Mozambique protected areas
# Final resolution:  1 x 1 km
# ============================================================

library(sf)
library(terra)
library(dplyr)

parks_file <- "data/processed/parks_core.Rdata"
#load('parks_core.Rdata')

elevation_file <- "data/raw/topography/elevation_1KMmd_GMTEDmd.tif"
elevation_file = '/Users/gdt366/Library/CloudStorage/Dropbox/Vertebrates_mozambique_project/Topographic layers/elevation_1KMmd_GMTEDmd.tif'
slope_file <- "data/raw/topography/slope_1KMmd_GMTEDmd.tif"
slope_file = '/Users/gdt366/Library/CloudStorage/Dropbox/Vertebrates_mozambique_project/Topographic layers/slope_1KMmd_GMTEDmd.tif'

treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
treecover_file <- "Mozambique_treecover2000_percent.tif"

settlement_distance_file <- "data/raw/accessibility/Mozambique_distance_to_settlements_2000.tif"
settlement_distance_file <- "Mozambique_distance_to_settlements_2000.tif"

road_distance_file <- "data/raw/accessibility/Mozambique_distance_to_roads.tif"
road_distance_file <- "Mozambique_distance_to_roads.tif"

crop_dir <- "data/raw/glad_cropland"
built_dir <- "data/raw/ghsl"

crop_dir <- "glad_cropland"
built_dir <- "ghsl"

output_dir <- "data/processed/matching"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Protected areas used in the analysis
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
  st_make_valid() %>%
  st_transform(4326)

# ------------------------------------------------------------
# 50-km comparison landscape around each protected area
# and divide it into approximately 1-km cells.
# Treatment is assigned from the cell centroid.
# ------------------------------------------------------------

moll_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
cellsize_m <- 1000
buffer_km <- 50
cellsize_deg <- cellsize_m / 111320

make_matching_grid <- function(parks, cellsize_deg, buffer_km) {

  grids <- vector("list", nrow(parks))

  for (i in seq_len(nrow(parks))) {

    park <- parks[i, ]
    message("Preparing grid: ", park$NAME_ENG)

    # Buffer in an equal-area projected CRS so 50 km is measured in metres.
    study_area <- park %>%
      st_transform(moll_crs) %>%
      st_buffer(dist = buffer_km * 1000) %>%
      st_transform(4326)

    grid <- st_make_grid(
      study_area,
      cellsize = cellsize_deg,
      square = TRUE
    )

    centres <- st_centroid(grid)
    keep <- lengths(st_intersects(centres, study_area)) > 0

    grid <- grid[keep]
    centres <- centres[keep]

    treated <- lengths(st_intersects(centres, park)) > 0

    grids[[i]] <- st_sf(
      NAME_ENG = park$NAME_ENG,
      treated = treated,
      resolution_m = cellsize_m,
      geometry = grid
    )
  }

  do.call(rbind, grids)
}

grid_1000 <- make_matching_grid(
  parks = parks_core,
  cellsize_deg = cellsize_deg,
  buffer_km = buffer_km
)

# Permanent spatial-unit id.
grid_1000$cell_id <- seq_len(nrow(grid_1000))

# ------------------------------------------------------------
# Matching covariates
# ------------------------------------------------------------

elevation <- rast(elevation_file)
slope <- rast(slope_file)
baseline_treecover <- rast(treecover_file)
distance_settlements <- rast(settlement_distance_file)
distance_roads <- rast(road_distance_file)

crop_files <- sort(list.files(
  crop_dir,
  pattern = "^Mozambique_GLAD_cropland_.*\\.tif$",
  full.names = TRUE
))

built_files <- sort(list.files(
  built_dir,
  pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
  full.names = TRUE
))

# The first temporal layer is the baseline used in matching:
# cropland = 2015; built-up = 2000.
baseline_cropland <- rast(crop_files[1])
baseline_built <- rast(built_files[1])

baseline_built = baseline_built[[1]]
# Match all covariates to the Hansen tree-cover template before stacking.
# Continuous covariates are interpolated bilinearly.
template <- baseline_treecover

project_to_template <- function(x) {
  project(x, template, method = "bilinear")
}

elevation <- project_to_template(elevation)
slope <- project_to_template(slope)
distance_roads <- project_to_template(distance_roads)
distance_settlements <- project_to_template(distance_settlements)
baseline_cropland <- project_to_template(baseline_cropland)
baseline_built <- project_to_template(baseline_built)

covariates <- c(
  elevation,
  slope,
  distance_roads,
  distance_settlements,
  baseline_treecover,
  baseline_cropland,
  baseline_built
)

names(covariates) <- c(
  "elevation",
  "slope",
  "dist_roads",
  "dist_settlements",
  "treecover_2000",
  "cropland_baseline",
  "baseline_built"
)

# Resample all covariates to the final approximately 1-km analytical grid.
res_1000 <- 1000 / 111320
r1000 <- rast(
  ext(vect(grid_1000)),
  resolution = res_1000,
  crs = crs(vect(grid_1000))
)

cov_1000 <- resample(covariates, r1000, method = "bilinear")

centres_1000 <- st_centroid(grid_1000)
vals_1000 <- terra::extract(cov_1000, vect(centres_1000))

dat_1000 <- cbind(
  data.frame(
    cell_id = grid_1000$cell_id,
    NAME_ENG = grid_1000$NAME_ENG,
    treated = grid_1000$treated,
    resolution_m = grid_1000$resolution_m
  ),
  vals_1000[, -1, drop = FALSE]
)


saveRDS(dat_1000, file.path(output_dir, "matching_input_1000m.rds"))
saveRDS(grid_1000, file.path(output_dir, "matching_grid_1000m.rds"))

