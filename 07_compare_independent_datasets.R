# 07. Comparison with six independent datasets

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)

root <- getwd()

primary_dir <- file.path(root, "data/raw")
independent_dir <- file.path(primary_dir, "independent")
output_dir <- file.path(root, "outputs/dataset_comparison")
cache_dir <- file.path(output_dir, "cache_5km_six_datasets_v1")
tmp_dir <- file.path(output_dir, "terra_tmp")

for (p in c(output_dir, cache_dir, tmp_dir)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}



# Common 5-km equal-area grid

moz <- st_read(
  file.path(primary_dir, "moz_admin_boundaries.shp/moz_admin0.shp"),
  quiet = TRUE
) |>
  st_make_valid() |>
  st_transform(6933)

bbox <- st_bbox(moz)
cell_m <- 5000

common_grid <- rast(
  ext(
    floor(as.numeric(bbox["xmin"]) / cell_m) * cell_m,
    ceiling(as.numeric(bbox["xmax"]) / cell_m) * cell_m,
    floor(as.numeric(bbox["ymin"]) / cell_m) * cell_m,
    ceiling(as.numeric(bbox["ymax"]) / cell_m) * cell_m
  ),
  resolution = cell_m,
  crs = "EPSG:6933"
)

# Raster processing

project_cached <- function(r, key) {

  filename <- file.path(cache_dir, paste0(key, ".tif"))

  if (file.exists(filename)) {
    return(rast(filename))
  }


  project(
    r,
    common_grid,
    method = "average",
    filename = filename,
    overwrite = FALSE,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
}

rasterize_cached <- function(v, key) {

  filename <- file.path(cache_dir, paste0(key, ".tif"))

  if (file.exists(filename)) {
    return(rast(filename))
  }


  rasterize(
    v,
    common_grid,
    cover = TRUE,
    background = 0,
    filename = filename,
    overwrite = FALSE,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
}

make_comparison <- function(x, y, x_name, y_name, zero_x = FALSE) {

  if (zero_x) {
    x <- ifel(!is.na(y) & is.na(x), 0, x)
  }

  names(x) <- x_name
  names(y) <- y_name

  d <- as.data.frame(c(x, y), na.rm = TRUE)
  names(d) <- c("x", "y")

  d[
    is.finite(d$x) &
      is.finite(d$y) &
      (d$x > 0 | d$y > 0),
    ,
    drop = FALSE
  ]
}

# FNDS vs Hansen, 2024


fnds_dir <- file.path(root, "FNDS_deforestation_2024")

fnds_files <- list.files(
  fnds_dir,
  pattern = "_2024\\.gpkg$",
  full.names = TRUE
)

national_gpkg <- file.path(
  cache_dir,
  "FNDS_deforestation_2024_Mozambique.gpkg"
)

if (!file.exists(national_gpkg)) {

  provinces <- lapply(fnds_files, function(f) {


    z <- st_read(f, quiet = TRUE)

    if ("areaHa" %in% names(z)) {
      names(z)[names(z) == "areaHa"] <- "area_ha"
    }

    z$province <- sub("_2024\\.gpkg$", "", basename(f))

    st_transform(z, 6933)
  })

  fnds_sf <- bind_rows(provinces)

  st_write(fnds_sf, national_gpkg, quiet = TRUE)

  rm(provinces, fnds_sf)
  gc()
}

fnds <- vect(national_gpkg)

fnds_5km <- rasterize_cached(
  fnds, "fnds_deforestation_2024_fraction"
)

rm(fnds)
gc()

treecover <- rast(
  file.path(primary_dir, "hansen/Mozambique_treecover2000_percent.tif")
)

lossyear <- rast(
  file.path(primary_dir, "hansen/Mozambique_forest_loss_year_2001_2025.tif")
)

hansen24 <- ifel(
  is.na(treecover) | is.na(lossyear),
  NA,
  ifel(treecover >= 30 & lossyear == 2024, 1, 0)
)

hansen24_5km <- project_cached(
  hansen24, "hansen_binary30_loss_2024"
)

national_forest <- make_comparison(
  fnds_5km, hansen24_5km, "FNDS", "Hansen",
  zero_x = TRUE
)

rm(hansen24, hansen24_5km, fnds_5km)
gc()

# MozFields vs GACED30, 2023


mozfields <- rast(
  file.path(primary_dir, "mozambique_crop_fraction_p065_2023.tif")
)

mozfields_5km <- project_cached(
  mozfields, "mozfields_crop_fraction_2023"
)

rm(mozfields)
gc()

gaced <- rast(
  file.path(primary_dir, "GACED30/Mozambique_GACED30_cropland_2000_2024.tif")
)

gaced23 <- gaced[["crop_2023"]]

gaced23_binary <- ifel(
  gaced23 == 10, 1,
  ifel(gaced23 == 0, 0, NA)
)

gaced23_5km <- project_cached(
  gaced23_binary, "gaced30_crop_fraction_2023"
)

national_crop <- make_comparison(
  mozfields_5km, gaced23_5km,
  "MozFields", "GACED30"
)

rm(mozfields_5km, gaced23_5km, gaced23, gaced23_binary)
gc()

# GRID3 vs GHSL, 2025


grid3 <- vect(
  file.path(primary_dir, "GRID3_MOZ_settlement_extents_v3_0.gpkg")
)

grid3 <- project(grid3, "EPSG:6933")
grid3$presence <- 1

grid3_5km <- rasterize_cached(
  grid3, "grid3_settlement_fraction_5km"
)

rm(grid3)
gc()

ghsl25 <- rast(
  file.path(primary_dir, "ghsl/Mozambique_GHSL_built_surface_2025.tif")
)[["built_surface"]]

ghsl25_fraction <- ifel(
  ghsl25 >= 0 & ghsl25 <= 10000,
  ghsl25 / 10000,
  NA
)

ghsl25_5km <- project_cached(
  ghsl25_fraction, "ghsl_built_fraction_2025"
)

national_built <- make_comparison(
  grid3_5km, ghsl25_5km,
  "GRID3", "GHSL",
  zero_x = TRUE
)

rm(grid3_5km, ghsl25, ghsl25_fraction, ghsl25_5km)
gc()

# GLAD vs Hansen, 2001-2020


hansen_binary <- ifel(
  is.na(treecover) | is.na(lossyear),
  NA,
  ifel(
    treecover >= 30 & lossyear >= 2001 & lossyear <= 2020,
    1, 0
  )
)

hansen_5km <- project_cached(
  hansen_binary, "hansen_loss_2001_2020"
) * 100

rm(treecover, lossyear, hansen_binary)
gc()

glad <- rast(
  file.path(
    independent_dir,
    "Mozambique_GLAD_forest_extent_loss_2000_2020.tif"
  )
)

glad_5km <- project_cached(
  glad, "glad_forest_loss_2000_2020"
) * 100

rm(glad)
gc()

international_forest <- make_comparison(
  glad_5km, hansen_5km,
  "GLAD", "Hansen"
)

rm(glad_5km, hansen_5km)
gc()

# WorldCereal vs GACED30, 2021


gaced21 <- gaced[["crop_2021"]]

gaced21_binary <- ifel(
  gaced21 == 10, 100,
  ifel(gaced21 == 0, 0, NA)
)

gaced21_5km <- project_cached(
  gaced21_binary, "gaced30_crop_2021"
)

rm(gaced, gaced21, gaced21_binary)
gc()

wc <- rast(
  file.path(
    independent_dir,
    "Mozambique_WorldCereal_cropland_2021.tif"
  )
)

wc_5km <- project_cached(
  wc, "worldcereal_crop_2021"
)

rm(wc)
gc()

international_crop <- make_comparison(
  wc_5km, gaced21_5km,
  "WorldCereal", "GACED30"
)

rm(wc_5km, gaced21_5km)
gc()

# WSF vs GHSL, 2015


ghsl15 <- rast(
  file.path(primary_dir, "ghsl/Mozambique_GHSL_built_surface_2015.tif")
)[["built_surface"]]

ghsl15_pct <- ifel(
  ghsl15 >= 0 & ghsl15 <= 10000,
  ghsl15 / 100,
  NA
)

ghsl15_5km <- project_cached(
  ghsl15_pct, "ghsl_built_2015"
)

rm(ghsl15, ghsl15_pct)
gc()

wsf <- rast(
  file.path(
    independent_dir,
    "Mozambique_WSF2015_settlement_100m.tif"
  )
)

wsf_pct <- ifel(is.na(wsf), 0, wsf)

wsf_5km <- project_cached(
  wsf_pct, "wsf_settlement_2015"
)

rm(wsf, wsf_pct)
gc()

international_built <- make_comparison(
  wsf_5km, ghsl15_5km,
  "WSF", "GHSL"
)

rm(wsf_5km, ghsl15_5km)
gc()

# Convert national comparisons to percentages

national_forest <- national_forest * 100
national_crop <- national_crop * 100
national_built <- national_built * 100

# Statistical comparisons

fit_stats <- function(df, name) {

  fit <- lm(y ~ x, data = df)
  model <- summary(fit)

  data.frame(
    comparison = name,
    n_cells = nrow(df),
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2]),
    pearson_r = cor(df$x, df$y),
    spearman_rho = cor(df$x, df$y, method = "spearman"),
    R2 = model$r.squared,
    p_value = coef(model)[2, 4],
    MAE = mean(abs(df$y - df$x)),
    bias = mean(df$y - df$x),
    RMSE = sqrt(mean((df$y - df$x)^2))
  )
}

datasets <- list(
  national_forest,
  national_crop,
  national_built,
  international_forest,
  international_crop,
  international_built
)

comparison_names <- c(
  "FNDS vs Hansen 2024",
  "MozFields vs GACED30 2023",
  "GRID3 vs GHSL 2025",
  "GLAD vs Hansen 2001-2020",
  "WorldCereal vs GACED30 2021",
  "WSF vs GHSL 2015"
)

comparison_summary <- bind_rows(
  Map(fit_stats, datasets, comparison_names)
)

print(comparison_summary)

write.csv(
  comparison_summary,
  file.path(output_dir, "independent_dataset_comparison_six_datasets.csv"),
  row.names = FALSE
)

comparison_data <- bind_rows(
  Map(
    function(d, nm) {
      data.frame(comparison = nm, x = d$x, y = d$y)
    },
    datasets,
    comparison_names
  )
)

write.csv(
  comparison_data,
  file.path(output_dir, "independent_dataset_comparison_six_datasets_data.csv"),
  row.names = FALSE
)

# Figure 4

panel <- function(df, x_label, y_label, title, colour, idx) {

  ggplot(df, aes(x = x, y = y)) +
    geom_point(
      alpha = 0.15,
      size = 0.8,
      colour = colour
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = "black"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      hjust = 1.1,
      vjust = -0.5,
      size = 3.4,
      label = sprintf(
        "r = %.3f\nR² = %.3f",
        comparison_summary$pearson_r[idx],
        comparison_summary$R2[idx]
      )
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

p_a <- panel(
  national_forest,
  "FNDS deforestation (%)",
  "Hansen forest loss (%)",
  "a. Forest loss (2024)",
  "#238B45", 1
)

p_b <- panel(
  national_crop,
  "MozFields cropland (%)",
  "GACED30 cropland (%)",
  "b. Cropland (2023)",
  "#E68A00", 2
)

p_c <- panel(
  national_built,
  "GRID3 settlement extent (%)",
  "GHSL built-up surface (%)",
  "c. Built-up (2025)",
  "#7A3DB8", 3
)

p_d <- panel(
  international_forest,
  "GLAD forest loss (%)",
  "Hansen forest loss (%)",
  "d. Forest loss (2001-2020)",
  "#238B45", 4
)

p_e <- panel(
  international_crop,
  "WorldCereal cropland (%)",
  "GACED30 cropland (%)",
  "e. Cropland (2021)",
  "#E68A00", 5
)

p_f <- panel(
  international_built,
  "WSF settlement (%)",
  "GHSL built-up surface (%)",
  "f. Built-up (2015)",
  "#7A3DB8", 6
)

figure_4 <- (p_a + p_b + p_c) / (p_d + p_e + p_f)


ggsave(
  file.path(output_dir, "Figure4_independent_dataset_comparison_2x3.png"),
  figure_4,
  width = 12,
  height = 8,
  dpi = 600
)

