
# 06. National land-use change

library(sf)
library(terra)
library(dplyr)

RUN_FOREST <- TRUE
RUN_CROP <- TRUE
RUN_BUILT <- TRUE

TILE_SIZE <- 512L

terraOptions(
  progress = 0,
  memfrac = 0.15,
  memmax = 4
)

out <- "outputs/national"
ckroot <- file.path(out, "checkpoints_national_binary30_tiles512_v1")
tmp <- file.path(out, "terra_tmp")

dir.create(ckroot, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = normalizePath(tmp))

# National boundary

boundary <- "data/raw/moz_admin_boundaries.shp/moz_admin0.shp"

moz0 <- st_make_valid(st_read(boundary, quiet = TRUE))

moz <- st_sf(
  geometry = st_sfc(
    st_union(st_geometry(moz0)),
    crs = st_crs(moz0)
  )
)

land_ha <- as.numeric(
  st_area(st_transform(moz, 6933))
) / 10000

forest_years <- 2001:2025
crop_years <- 2000:2024
built_years <- c(2000, 2005, 2010, 2015, 2020, 2025)

# Tiled processing with checkpoints

tile_run <- function(r, tag, keys, calc) {

  poly <- vect(st_transform(moz, st_crs(crs(r))))

  rows <- seq(1, nrow(r), by = TILE_SIZE)
  cols <- seq(1, ncol(r), by = TILE_SIZE)

  ckdir <- file.path(ckroot, tag)
  dir.create(ckdir, recursive = TRUE, showWarnings = FALSE)

  total <- setNames(numeric(length(keys)), keys)

  for (ri in seq_along(rows)) {
    for (ci in seq_along(cols)) {

      file <- file.path(
        ckdir,
        sprintf("tile_%05d_%05d.rds", ri, ci)
      )

      if (file.exists(file)) {

        z <- readRDS(file)

      } else {

        rr <- rows[ri]
        cc <- cols[ci]

        re <- min(rr + TILE_SIZE - 1L, nrow(r))
        ce <- min(cc + TILE_SIZE - 1L, ncol(r))

        e <- ext(
          xmin(r) + (cc - 1) * res(r)[1],
          xmin(r) + ce * res(r)[1],
          ymax(r) - re * res(r)[2],
          ymax(r) - (rr - 1) * res(r)[2]
        )

        tilepoly <- as.polygons(e, crs = crs(r))

        z <- setNames(numeric(length(keys)), keys)

        if (any(relate(tilepoly, poly, "intersects"))) {

          block <- crop(r, e, snap = "near")

          vals <- extract(
            block,
            poly,
            weights = TRUE,
            cells = TRUE
          )

          if (nrow(vals)) {
            z <- calc(block, vals)
          }

          rm(block, vals)
        }

        if (!identical(names(z), keys) || any(!is.finite(z))) {
          stop("Invalid tile ", tag, " ", ri, "/", ci)
        }

        partial <- paste0(file, ".partial")
        saveRDS(z, partial)

        if (!file.rename(partial, file)) {
          stop("Checkpoint rename failed: ", file)
        }
      }

      if (!identical(names(z), keys)) {
        stop("Checkpoint schema mismatch: ", file)
      }

      total <- total + z
    }
  }

  total
}

# Tree-cover loss, 2001-2025

if (RUN_FOREST) {

  tc <- rast(
    "data/raw/hansen/Mozambique_treecover2000_percent.tif"
  )

  ly <- rast(
    "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"
  )

  if (!compareGeom(tc, ly, stopOnError = FALSE)) {
    stop("Hansen grids differ")
  }

  h <- c(tc, ly)
  names(h) <- c("tree", "loss")

  keys <- c("baseline", paste0("loss_", forest_years))

  ft <- tile_run(
    h,
    "forest",
    keys,
    function(block, v) {

      area <- cellSize(block[[1]], unit = "ha")
      names(area) <- "area"

      a <- extract(area, v$cell)

      z <- setNames(numeric(length(keys)), keys)

      forest <- !is.na(v$tree) & v$tree >= 30
      w <- a$area * v$weight

      z[1] <- sum(w[forest], na.rm = TRUE)

      for (i in seq_along(forest_years)) {

        z[i + 1L] <- sum(
          w[
            forest &
              !is.na(v$loss) &
              v$loss == forest_years[i]
          ],
          na.rm = TRUE
        )
      }

      z
    }
  )

  base_ha <- unname(ft["baseline"])

  forest_annual <- tibble(
    year = forest_years,
    annual_loss_ha = unname(
      ft[paste0("loss_", forest_years)]
    )
  ) |>
    mutate(
      cumulative_loss_ha = cumsum(annual_loss_ha)
    )

  loss_ha <- sum(forest_annual$annual_loss_ha)
  final_forest_ha <- base_ha - loss_ha

  if (final_forest_ha < -1e-6) {
    stop("Forest loss exceeds baseline")
  }

  loss_pct <- if (base_ha > 0) {
    100 * loss_ha / base_ha
  } else {
    NA_real_
  }

  write.csv(
    forest_annual,
    file.path(
      out,
      "Mozambique_treecover_annual_loss_2001_2025.csv"
    ),
    row.names = FALSE
  )

  rm(tc, ly, h)
  gc()
}

# Cropland, 2000-2024

if (RUN_CROP) {

  g <- rast(
    "data/raw/GACED30/Mozambique_GACED30_cropland_2000_2024.tif"
  )

  keys <- paste0("crop_", crop_years)

  if (!all(keys %in% names(g))) {
    stop("GACED layers missing")
  }

  g <- g[[keys]]

  ct <- tile_run(
    g,
    "crop",
    keys,
    function(block, v) {

      a <- cellSize(block[[1]], unit = "ha")
      names(a) <- "area"

      av <- extract(a, v$cell)$area
      w <- av * v$weight

      z <- setNames(numeric(length(keys)), keys)

      for (i in seq_along(keys)) {

        z[i] <- sum(
          w[
            !is.na(v[[keys[i]]]) &
              v[[keys[i]]] == 10
          ],
          na.rm = TRUE
        )
      }

      z
    }
  )

  crop_annual <- tibble(
    year = crop_years,
    cropland_ha = unname(ct[keys])
  ) |>
    mutate(
      change_from_2000_ha = cropland_ha - cropland_ha[1],
      change_from_2000_percent = if (cropland_ha[1] > 0) {
        100 * change_from_2000_ha / cropland_ha[1]
      } else {
        NA_real_
      },
      percent_land = 100 * cropland_ha / land_ha,
      change_from_2000_pp = percent_land - percent_land[1]
    )

  crop_2000_ha <- crop_annual$cropland_ha[1]
  crop_2024_ha <- tail(crop_annual$cropland_ha, 1)

  crop_delta <- crop_2024_ha - crop_2000_ha

  crop_pct <- if (crop_2000_ha > 0) {
    100 * crop_delta / crop_2000_ha
  } else {
    NA_real_
  }

  write.csv(
    crop_annual,
    file.path(
      out,
      "Mozambique_cropland_annual_2000_2024.csv"
    ),
    row.names = FALSE
  )

  rm(g)
  gc()
}

# Built-up surface, 2000-2025

find_year <- function(files, y) {

  f <- files[
    grepl(
      paste0("(^|[^0-9])", y, "([^0-9]|$)"),
      basename(files)
    )
  ]

  if (length(f) != 1L) {
    stop("GHSL file count for ", y, ": ", length(f))
  }

  f
}

if (RUN_BUILT) {

  fs <- list.files(
    "data/raw/ghsl",
    pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
    full.names = TRUE
  )

  rasters <- lapply(
    built_years,
    function(y) {
      rast(find_year(fs, y))[["built_surface"]]
    }
  )

  for (i in seq_along(rasters)[-1L]) {

    if (!compareGeom(
      rasters[[1]],
      rasters[[i]],
      stopOnError = FALSE
    )) {
      stop("GHSL year rasters have differing grids")
    }
  }

  b <- do.call(c, rasters)

  keys <- paste0("built_", built_years)
  names(b) <- keys

  bt <- tile_run(
    b,
    "built",
    keys,
    function(block, v) {

      z <- setNames(numeric(length(keys)), keys)

      for (i in seq_along(keys)) {

        x <- v[[keys[i]]]

        ok <- !is.na(x) &
          x >= 0 &
          x <= 10000

        z[i] <- sum(
          x[ok] * v$weight[ok],
          na.rm = TRUE
        )
      }

      z
    }
  )

  built_temporal <- tibble(
    year = built_years,
    built_surface_m2 = unname(bt[keys])
  ) |>
    mutate(
      built_surface_ha = built_surface_m2 / 10000,
      percent_land = 100 * built_surface_ha / land_ha,
      change_from_2000_ha = built_surface_ha - built_surface_ha[1],
      change_from_2000_percent = if (built_surface_ha[1] > 0) {
        100 * change_from_2000_ha / built_surface_ha[1]
      } else {
        NA_real_
      },
      change_from_2000_pp = percent_land - percent_land[1]
    )

  built_2000_ha <- built_temporal$built_surface_ha[1]
  built_2025_ha <- tail(built_temporal$built_surface_ha, 1)

  built_delta <- built_2025_ha - built_2000_ha

  built_pct <- if (built_2000_ha > 0) {
    100 * built_delta / built_2000_ha
  } else {
    NA_real_
  }

  write.csv(
    built_temporal,
    file.path(
      out,
      "Mozambique_builtup_temporal_2000_2025.csv"
    ),
    row.names = FALSE
  )

  rm(b, rasters)
  gc()
}

# National summary table

if (all(c(RUN_FOREST, RUN_CROP, RUN_BUILT))) {

  summary <- tibble(
    Threat = c("Forest", "Cropland", "Built-up surface"),
    `Baseline year` = c(2000, 2000, 2000),
    `Final year` = c(2025, 2024, 2025),
    `Baseline extent (ha)` = c(
      base_ha, crop_2000_ha, built_2000_ha
    ),
    `Final extent (ha)` = c(
      final_forest_ha, crop_2024_ha, built_2025_ha
    ),
    `Loss / expansion (ha)` = c(
      loss_ha, crop_delta, built_delta
    ),
    `Loss / expansion (%)` = c(
      loss_pct, crop_pct, built_pct
    )
  ) |>
    mutate(
      `Baseline extent (Mha)` = `Baseline extent (ha)` / 1e6,
      `Final extent (Mha)` = `Final extent (ha)` / 1e6,
      `Loss / expansion (Mha)` = `Loss / expansion (ha)` / 1e6,
      `Baseline land share (%)` = 100 * `Baseline extent (ha)` / land_ha,
      `Final land share (%)` = 100 * `Final extent (ha)` / land_ha,
      `Change in land share (pp)` =
        `Final land share (%)` - `Baseline land share (%)`
    ) |>
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
      `Change in land share (pp)` = round(
        `Change in land share (pp)`, 3
      )
    )

  write.csv(
    summary,
    file.path(
      out,
      "Table_S1_national_land_use_change_summary_fixed.csv"
    ),
    row.names = FALSE
  )
}
