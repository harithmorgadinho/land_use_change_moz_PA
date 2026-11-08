# 02. Matching grid and covariates - 500 m sensitivity analysis

library(sf)
library(terra)
library(dplyr)

parks_file <- "data/processed/parks_core.Rdata"
wdpa_dir <- "data/raw/wdpa"
veg_file <- "/Users/gdt366/Library/CloudStorage/Dropbox/gee_deforestation_mozambique/Mozambique_historical_vegetation_2022.gpkg"
elevation_file <- "/Users/gdt366/Library/CloudStorage/Dropbox/Vertebrates_mozambique_project/Topographic layers/elevation_1KMmd_GMTEDmd.tif"
slope_file <- "/Users/gdt366/Library/CloudStorage/Dropbox/Vertebrates_mozambique_project/Topographic layers/slope_1KMmd_GMTEDmd.tif"
treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
settlement_distance_file <- "Mozambique_distance_to_settlements_2000.tif"
road_distance_file <- "Mozambique_distance_to_roads.tif"
built_dir <- "ghsl"
gaced_file <- "data/raw/GACED30/Mozambique_GACED30_cropland_2000_2024.tif"

out_dir <- "data/processed/matching_GACED30/revised_500m_forest_proportion_exclusion_exception"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "terra_tmp"), recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = file.path(out_dir, "terra_tmp"), progress = 1)

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
cellsize_m <- 500
inner_m <- 20000
outer100_m <- 100000
tree_threshold <- 30

# Protected areas

e <- new.env(parent = emptyenv())
load(parks_file, envir = e)

parks <- e$parks_core |>
  filter(DESIG %in% c("Parque Nacional", "National Park", "Reserva Nacional", "Reserva Especial")) |>
  st_make_valid() |>
  st_transform(moll_crs)

wdpa_files <- list.files(wdpa_dir, pattern = "\\.shp$", full.names = TRUE)

wdpa <- bind_rows(lapply(wdpa_files, function(f) st_read(f, quiet = TRUE))) |>
  filter(SITE_TYPE == "PA", STATUS == "Designated", REALM == "Terrestrial") |>
  st_make_valid() |>
  st_transform(moll_crs)

# Exclude Zambezi Delta and Tchuma Tchato from the WDPA control mask

exception_names <- c("Zambezi Delta", "Tchuma Tchato")

normalise_name <- function(x) {
  tolower(trimws(gsub("[[:space:]]+", " ", as.character(x))))
}

exception_rows <- normalise_name(wdpa$NAME_ENG) %in% normalise_name(exception_names)
exception_columns <- intersect(
  c("NAME_ENG", "DESIG", "STATUS", "IUCN_CAT", "WDPAID"), names(wdpa)
)

write.csv(
  st_drop_geometry(wdpa[exception_rows, exception_columns, drop = FALSE]),
  file.path(out_dir, "wdpa_exclusion_exceptions.csv"),
  row.names = FALSE
)

wdpa <- wdpa[!exception_rows, ]
no_control_zone <- st_buffer(st_union(st_geometry(wdpa)), inner_m)

# Historical vegetation

veg <- st_read(veg_file, quiet = TRUE) |>
  st_make_valid() |>
  st_transform(moll_crs)

veg$Regional_Ecosystem <- trimws(as.character(veg$Regional_Ecosystem))
veg <- veg[!is.na(veg$Regional_Ecosystem) & nzchar(veg$Regional_Ecosystem), ]
veg_classes <- veg$Regional_Ecosystem

# 500 m grid

study_100km <- st_buffer(st_union(st_geometry(parks)), outer100_m)
all_grid <- st_make_grid(study_100km, cellsize = cellsize_m, square = TRUE)
all_centres <- st_centroid(all_grid)

in_study <- lengths(st_intersects(all_centres, study_100km)) > 0
all_grid <- all_grid[in_study]
all_centres <- all_centres[in_study]
global_id <- seq_along(all_centres)

excluded <- lengths(st_intersects(all_centres, no_control_zone)) > 0
candidate_idx <- which(!excluded)

candidate_points <- st_sf(
  global_id = candidate_idx,
  geometry = all_centres[candidate_idx]
)

control_veg_hits <- st_intersects(candidate_points, veg)

# Protected and control cells

assignments <- vector("list", nrow(parks))

for (i in seq_len(nrow(parks))) {

  park <- parks[i, ]
  nm <- as.character(park$NAME_ENG[[1]])

  park_classes <- unique(veg_classes[st_intersects(park, veg)[[1]]])
  treated_idx <- which(lengths(st_intersects(all_centres, park)) > 0)

  distances <- as.numeric(st_distance(candidate_points, st_geometry(park)))
  within100 <- which(distances <= outer100_m)

  ecological <- vapply(control_veg_hits[within100], function(h) {
    length(h) > 0 && any(veg_classes[h] %in% park_classes)
  }, logical(1))

  eligible_rows <- within100[ecological]
  control_idx <- candidate_idx[eligible_rows]
  control_distance <- distances[eligible_rows] / 1000

  assignments[[i]] <- bind_rows(
    data.frame(
      NAME_ENG = nm,
      global_id = treated_idx,
      treated = TRUE,
      distance_km = 0,
      pool_20_50 = FALSE,
      pool_20_100 = FALSE
    ),
    data.frame(
      NAME_ENG = nm,
      global_id = control_idx,
      treated = FALSE,
      distance_km = control_distance,
      pool_20_50 = control_distance <= 50,
      pool_20_100 = TRUE
    )
  )
}

assignment <- bind_rows(assignments) |>
  mutate(resolution_m = cellsize_m) |>
  select(NAME_ENG, global_id, treated, resolution_m,
         distance_km, pool_20_50, pool_20_100)

grid_sf <- st_sf(global_id = global_id, geometry = all_grid)

saveRDS(grid_sf, file.path(out_dir, "grid_500m_mollweide.rds"))
saveRDS(assignment, file.path(out_dir, "park_cell_assignments.rds"))
write.csv(assignment, file.path(out_dir, "park_cell_assignments.csv"), row.names = FALSE)

# Covariates

needed_idx <- sort(unique(assignment$global_id))
needed_points <- st_sf(global_id = needed_idx, geometry = all_centres[needed_idx])

extract_one <- function(r, method = "bilinear") {
  pts <- st_transform(needed_points, st_crs(crs(r)))
  values <- terra::extract(r, terra::vect(pts), method = method, ID = FALSE)
  as.numeric(values[[1]])
}

cov_data <- data.frame(global_id = needed_idx)

cov_data$elevation <- extract_one(rast(elevation_file))
cov_data$slope <- extract_one(rast(slope_file))
cov_data$dist_roads <- extract_one(rast(road_distance_file))
cov_data$dist_settlements <- extract_one(rast(settlement_distance_file))

# Forest cover in 2000 (Hansen, >=30% canopy cover)

tree_raw <- rast(treecover_file)

forest_2000 <- terra::ifel(
  is.na(tree_raw), NA,
  terra::ifel(tree_raw >= tree_threshold, 1, 0),
  filename = file.path(out_dir, "hansen_forest_binary_30pct.tif"),
  overwrite = FALSE
)

bb <- st_bbox(grid_sf)

target_grid <- terra::rast(
  xmin = as.numeric(bb["xmin"]),
  xmax = as.numeric(bb["xmax"]),
  ymin = as.numeric(bb["ymin"]),
  ymax = as.numeric(bb["ymax"]),
  resolution = cellsize_m,
  crs = moll_crs
)

forest_fraction <- terra::project(
  forest_2000, target_grid,
  method = "average",
  use_gdal = TRUE,
  filename = file.path(out_dir, "hansen_forest_proportion_500m.tif"),
  overwrite = FALSE,
  wopt = list(gdal = "COMPRESS=LZW")
)

names(forest_fraction) <- "treecover_2000"
cov_data$treecover_2000 <- extract_one(forest_fraction, method = "simple")

rm(tree_raw, forest_2000, forest_fraction, target_grid)
gc()

# Cropland cover in 2000 (GACED30)

crop_raw <- rast(gaced_file)[["crop_2000"]]
crop_baseline <- terra::ifel(crop_raw == 10, 1, terra::ifel(crop_raw == 0, 0, NA))

cov_data$cropland_baseline <- extract_one(crop_baseline)

# Built-up surface in 2000 (GHSL)

built_files <- sort(list.files(
  built_dir,
  pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
  full.names = TRUE
))

built_baseline <- rast(built_files[[1]])[[1]]
cov_data$baseline_built <- extract_one(built_baseline)

# Matching input

matching_input <- left_join(assignment, cov_data, by = "global_id")

saveRDS(cov_data, file.path(out_dir, "covariates_unique_cells.rds"))
saveRDS(matching_input, file.path(out_dir, "matching_input_revised_500m.rds"))

# Covariate completeness by protected area

covars <- c(
  "elevation", "slope", "dist_roads", "dist_settlements",
  "treecover_2000", "cropland_baseline", "baseline_built"
)

report <- bind_rows(lapply(split(matching_input, matching_input$NAME_ENG), function(x) {

  complete <- complete.cases(x[, covars, drop = FALSE])

  data.frame(
    NAME_ENG = x$NAME_ENG[[1]],
    treated_total = sum(x$treated),
    treated_complete = sum(x$treated & complete),
    controls_20_50_total = sum(!x$treated & x$pool_20_50),
    controls_20_50_complete = sum(!x$treated & x$pool_20_50 & complete),
    controls_20_100_total = sum(!x$treated & x$pool_20_100),
    controls_20_100_complete = sum(!x$treated & x$pool_20_100 & complete)
  )
})) |>
  arrange(NAME_ENG)

write.csv(
  report,
  file.path(out_dir, "covariate_completeness_by_park.csv"),
  row.names = FALSE
)
