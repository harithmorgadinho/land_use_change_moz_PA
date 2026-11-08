
# 05. Temporal land-use change within and around protected areas

library(sf)
library(terra)
library(dplyr)
library(tidyr)

# Inputs and settings

parks_file <- "data/processed/parks_core.Rdata"
rings_file <- "Mozambique_PA_equal_area_surroundings_nationally_clipped.shp"
treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
lossyear_file <- "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"
gaced_file <- "data/raw/GACED30/Mozambique_GACED30_cropland_2000_2024.tif"
built_dir <- "data/raw/ghsl"

output_dir <- "outputs/temporal"
checkpoint_dir <- file.path(output_dir, "checkpoints_binary30_tiled_v1")
tmp_dir <- file.path(output_dir, "terra_tmp_binary30_tiled")

for (p in c(output_dir, checkpoint_dir, tmp_dir)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

terra::terraOptions(
  memfrac = 0.15,
  memmax = 4,
  tempdir = normalizePath(tmp_dir),
  progress = 0
)

TILE_PIXELS <- 256L
TREE_COVER_THRESHOLD <- 30

forest_years <- 2001:2025
crop_years <- 2000:2024
built_years <- c(2000, 2005, 2010, 2015, 2020, 2025)

RUN_FOREST <- TRUE
RUN_CROP <- TRUE
RUN_BUILT <- TRUE

# Protected areas and surrounding regions

park_env <- new.env(parent = globalenv())
load(parks_file, envir = park_env)

stopifnot(exists("parks_core", envir = park_env))

parks_core <- get("parks_core", envir = park_env) |>
  dplyr::filter(
    DESIG %in% c(
      "Parque Nacional",
      "National Park",
      "Reserva Nacional",
      "Reserva Especial"
    )
  ) |>
  sf::st_make_valid()

stopifnot(
  nrow(parks_core) == 13L,
  !anyDuplicated(parks_core$NAME_ENG)
)

rings_all <- sf::st_read(rings_file, quiet = TRUE) |>
  dplyr::rename(NAME_ENG = park_name) |>
  sf::st_make_valid() |>
  dplyr::filter(NAME_ENG %in% parks_core$NAME_ENG)

stopifnot(
  !anyDuplicated(rings_all$NAME_ENG),
  setequal(rings_all$NAME_ENG, parks_core$NAME_ENG)
)

moll_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

region_areas <- dplyr::bind_rows(
  parks_core |>
    sf::st_transform(moll_crs) |>
    dplyr::mutate(
      region_area_ha = as.numeric(sf::st_area(geometry)) / 10000,
      location = "Inside"
    ) |>
    sf::st_drop_geometry() |>
    dplyr::select(NAME_ENG, location, region_area_ha),

  rings_all |>
    sf::st_transform(moll_crs) |>
    dplyr::mutate(
      region_area_ha = as.numeric(sf::st_area(geometry)) / 10000,
      location = "Outside"
    ) |>
    sf::st_drop_geometry() |>
    dplyr::select(NAME_ENG, location, region_area_ha)
)

get_region_area <- function(park_name, location_name) {

  x <- region_areas |>
    dplyr::filter(
      NAME_ENG == park_name,
      location == location_name
    ) |>
    dplyr::pull(region_area_ha)

  stopifnot(length(x) == 1L, is.finite(x), x > 0)

  x
}

polygon_for_raster <- function(polygon, raster) {
  terra::vect(
    sf::st_transform(
      polygon,
      sf::st_crs(terra::crs(raster))
    )
  )
}

# Checkpoint functions

checkpoint_path <- function(metric, location, park, year = NULL) {

  key <- paste(
    metric,
    location,
    park,
    if (is.null(year)) "all" else year,
    sep = "__"
  )

  key <- gsub(
    "[^[:alnum:]_-]",
    "_",
    iconv(key, to = "ASCII//TRANSLIT", sub = "_")
  )

  file.path(checkpoint_dir, paste0(key, ".rds"))
}

cached <- function(path, compute) {

  if (file.exists(path)) {
    return(readRDS(path))
  }

  value <- compute()
  pending <- paste0(path, ".pending")

  saveRDS(value, pending)

  if (!file.rename(pending, path)) {
    stop("Could not save checkpoint: ", path)
  }

  value
}

# Source-aligned tiled processing

visit_tiles <- function(reference_raster, polygon_v, tile_fun) {

  window <- terra::crop(
    reference_raster,
    polygon_v,
    snap = "out"
  )

  if (is.null(window)) {
    stop("Region does not overlap source raster")
  }

  nr <- terra::nrow(window)
  nc <- terra::ncol(window)

  row_starts <- seq.int(1L, nr, by = TILE_PIXELS)
  col_starts <- seq.int(1L, nc, by = TILE_PIXELS)

  accumulator <- NULL

  for (r1 in row_starts) {

    r2 <- min(nr, r1 + TILE_PIXELS - 1L)

    for (c1 in col_starts) {

      c2 <- min(nc, c1 + TILE_PIXELS - 1L)

      tile_ext <- terra::ext(
        terra::xmin(window) + (c1 - 1L) * terra::res(window)[1],
        terra::xmin(window) + c2 * terra::res(window)[1],
        terra::ymax(window) - r2 * terra::res(window)[2],
        terra::ymax(window) - (r1 - 1L) * terra::res(window)[2]
      )

      tile <- terra::crop(
        reference_raster,
        tile_ext,
        snap = "near"
      )

      piece <- tile_fun(tile, tile_ext, polygon_v)

      if (is.null(accumulator)) {
        accumulator <- piece
      } else {
        accumulator <- accumulator + piece
      }

      rm(tile, piece)
    }
  }

  if (is.null(accumulator)) {
    stop("No tiles were processed")
  }

  accumulator
}

weighted_total <- function(area_raster, polygon_v, value_name) {

  names(area_raster) <- value_name

  x <- terra::extract(
    area_raster,
    polygon_v,
    weights = TRUE
  )

  if (nrow(x) == 0L) {
    return(0)
  }

  sum(
    x[[value_name]] * x$weight,
    na.rm = TRUE
  )
}

# Tree-cover loss, 2001-2025

forest_series_one <- function(
    polygon, park, location, region_area_ha,
    treecover, lossyear
) {

  path <- checkpoint_path("forest", location, park)

  cached(path, function() {

    poly <- polygon_for_raster(polygon, treecover)

    totals <- visit_tiles(
      treecover,
      poly,
      function(tile, tile_ext, polygon_v) {

        loss_tile <- terra::crop(
          lossyear,
          tile_ext,
          snap = "near"
        )

        if (!terra::compareGeom(
          tile, loss_tile,
          stopOnError = FALSE
        )) {
          stop("Forest and loss tiles are not aligned")
        }

        baseline <- terra::ifel(
          is.na(tile),
          NA,
          terra::ifel(
            tile >= TREE_COVER_THRESHOLD,
            1,
            0
          )
        )

        area <- terra::cellSize(tile, unit = "ha")
        baseline_area <- area * baseline

        annual <- terra::ifel(
          baseline == 1 &
            loss_tile >= 2001 &
            loss_tile <= 2025,
          loss_tile,
          NA
        )

        names(baseline_area) <- "baseline_area"
        names(annual) <- "loss_year"

        extracted <- terra::extract(
          c(baseline_area, annual),
          polygon_v,
          weights = TRUE
        )

        result <- numeric(26L)
        names(result) <- c(
          "baseline",
          as.character(forest_years)
        )

        if (nrow(extracted)) {

          result[1] <- sum(
            extracted$baseline_area * extracted$weight,
            na.rm = TRUE
          )

          valid <- !is.na(extracted$loss_year) &
            extracted$loss_year %in% forest_years &
            !is.na(extracted$baseline_area)

          if (any(valid)) {

            year_sums <- rowsum(
              extracted$baseline_area[valid] *
                extracted$weight[valid],
              as.character(extracted$loss_year[valid]),
              reorder = FALSE
            )

            result[rownames(year_sums)] <- as.numeric(
              year_sums[, 1]
            )
          }
        }

        result
      }
    )

    baseline_ha <- unname(totals["baseline"])
    annual <- unname(totals[as.character(forest_years)])
    cumulative <- c(0, cumsum(annual))

    tibble::tibble(
      NAME_ENG = park,
      location = location,
      year = 2000:2025,
      region_area_ha = region_area_ha,
      baseline_forest_ha = baseline_ha,
      annual_loss_ha = c(0, annual),
      cumulative_loss_ha = cumulative,
      change_pp = if (baseline_ha > 0) {
        100 * cumulative / baseline_ha
      } else {
        rep(NA_real_, 26L)
      },
      threat = "Tree-cover loss"
    )
  })
}

# Cropland, 2000-2024

crop_series_one <- function(
    polygon, park, location, region_area_ha, gaced
) {

  annual <- lapply(crop_years, function(year) {

    path <- checkpoint_path(
      "crop", location, park, year
    )

    cached(path, function() {

      source <- gaced[[paste0("crop_", year)]]
      poly <- polygon_for_raster(polygon, source)

      area_ha <- visit_tiles(
        source,
        poly,
        function(tile, tile_ext, polygon_v) {

          crop_binary <- terra::ifel(
            tile == 10,
            1,
            terra::ifel(tile == 0, 0, NA)
          )

          weighted_total(
            terra::cellSize(tile, unit = "ha") *
              crop_binary,
            polygon_v,
            "crop_area_ha"
          )
        }
      )

      tibble::tibble(
        year = year,
        crop_area_ha = area_ha
      )
    })
  }) |>
    dplyr::bind_rows()

  baseline <- annual$crop_area_ha[
    annual$year == 2000
  ]

  annual |>
    dplyr::mutate(
      crop_change_ha = crop_area_ha - baseline,
      change_pp = crop_change_ha / region_area_ha * 100,
      NAME_ENG = park,
      location = location,
      region_area_ha = region_area_ha,
      threat = "Cropland"
    ) |>
    dplyr::select(
      NAME_ENG, location, year, region_area_ha,
      crop_area_ha, crop_change_ha, change_pp, threat
    )
}

# Built-up surface, 2000-2025

find_year_file <- function(files, year) {

  hits <- files[
    grepl(
      paste0("(^|[^0-9])", year, "([^0-9]|$)"),
      basename(files)
    )
  ]

  if (length(hits) != 1L) {
    stop(
      "Expected exactly one GHSL raster for ",
      year, "; found ", length(hits)
    )
  }

  hits[[1]]
}

built_series_one <- function(
    polygon, park, location, region_area_ha, year_files
) {

  annual <- lapply(built_years, function(year) {

    path <- checkpoint_path(
      "built", location, park, year
    )

    cached(path, function() {

      source <- terra::rast(
        year_files[[as.character(year)]]
      )[["built_surface"]]

      poly <- polygon_for_raster(polygon, source)

      area_ha <- visit_tiles(
        source,
        poly,
        function(tile, tile_ext, polygon_v) {

          valid <- terra::ifel(
            tile <= 10000,
            tile,
            NA
          )

          weighted_total(
            valid,
            polygon_v,
            "built_surface"
          ) / 10000
        }
      )

      tibble::tibble(
        year = year,
        built_area_ha = area_ha
      )
    })
  }) |>
    dplyr::bind_rows()

  baseline <- annual$built_area_ha[
    annual$year == 2000
  ]

  annual |>
    dplyr::mutate(
      built_change_ha = built_area_ha - baseline,
      change_pp = built_change_ha / region_area_ha * 100,
      NAME_ENG = park,
      location = location,
      region_area_ha = region_area_ha,
      threat = "Built-up"
    ) |>
    dplyr::select(
      NAME_ENG, location, year, region_area_ha,
      built_area_ha, built_change_ha, change_pp, threat
    )
}

# Run analyses and export results

regions <- list(
  Inside = parks_core,
  Outside = rings_all
)

run_regions <- function(fun) {

  result <- list()

  for (location in names(regions)) {

    polys <- regions[[location]]

    for (i in seq_len(nrow(polys))) {

      park <- as.character(
        polys$NAME_ENG[i]
      )

      result[[paste(location, park, sep = "__")]] <- fun(
        polys[i, ],
        park,
        location,
        get_region_area(park, location)
      )

      gc(verbose = FALSE)
    }
  }

  dplyr::bind_rows(result) |>
    dplyr::arrange(NAME_ENG, location, year)
}

if (RUN_FOREST) {

  treecover <- terra::rast(treecover_file)
  lossyear <- terra::rast(lossyear_file)

  if (!terra::compareGeom(
    treecover, lossyear,
    stopOnError = FALSE
  )) {
    stop("Tree-cover and loss-year rasters are not aligned")
  }

  tree_ts <- run_regions(function(poly, park, loc, area) {
    forest_series_one(
      poly, park, loc, area,
      treecover, lossyear
    )
  })

  saveRDS(
    tree_ts,
    file.path(
      output_dir,
      "PA_treecover_temporal_inside_outside.rds"
    )
  )

  write.csv(
    tree_ts,
    file.path(
      output_dir,
      "PA_treecover_temporal_inside_outside.csv"
    ),
    row.names = FALSE
  )

  rm(treecover, lossyear, tree_ts)
  gc()
}

if (RUN_CROP) {

  gaced <- terra::rast(gaced_file)

  stopifnot(
    all(paste0("crop_", crop_years) %in% names(gaced))
  )

  crop_ts <- run_regions(function(poly, park, loc, area) {
    crop_series_one(
      poly, park, loc, area, gaced
    )
  })

  saveRDS(
    crop_ts,
    file.path(
      output_dir,
      "PA_cropland_temporal_inside_outside_gaced.rds"
    )
  )

  write.csv(
    crop_ts,
    file.path(
      output_dir,
      "PA_cropland_temporal_inside_outside_gaced.csv"
    ),
    row.names = FALSE
  )

  rm(gaced, crop_ts)
  gc()
}

if (RUN_BUILT) {

  built_files <- list.files(
    built_dir,
    pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
    full.names = TRUE
  )

  year_files <- setNames(
    vapply(
      built_years,
      function(y) find_year_file(built_files, y),
      character(1)
    ),
    as.character(built_years)
  )

  built_ts <- run_regions(function(poly, park, loc, area) {
    built_series_one(
      poly, park, loc, area, year_files
    )
  })

  saveRDS(
    built_ts,
    file.path(
      output_dir,
      "PA_builtup_temporal_inside_outside.rds"
    )
  )

  write.csv(
    built_ts,
    file.path(
      output_dir,
      "PA_builtup_temporal_inside_outside.csv"
    ),
    row.names = FALSE
  )

  rm(built_ts)
  gc()
}
