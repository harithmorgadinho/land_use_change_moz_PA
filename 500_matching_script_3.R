# 04. Matched land-use effects - 500 m sensitivity analysis

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(MatchIt)

matching_data_dir <- "data/processed/matching_GACED30/revised_500m_forest_proportion_exclusion_exception"
matching_output_dir <- "outputs/matching/revised_no_replacement_500m_forest_proportion_exclusion_exception"
output_dir <- "outputs/results/revised_no_replacement_500m_forest_proportion_cell_area_exclusion_exception"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(memfrac = 0.5, tempdir = tempdir(), progress = 1)

treecover_file <- "data/raw/hansen/Mozambique_treecover2000_percent.tif"
lossyear_file <- "data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif"
gaced_file <- "data/raw/GACED30/Mozambique_GACED30_cropland_2000_2024.tif"
built_dir <- "data/raw/ghsl"

main_balance_threshold <- 0.25
strict_balance_threshold <- 0.10

all_parks <- c(
  "Banhine", "Chimanimani", "Gilé", "Gorongosa", "Limpopo",
  "Magoe", "Maputo", "Marromeu", "Niassa", "Pomene",
  "Quirimbas", "Serra da Gorongosa", "Zinave"
)

# Matching results

usage <- read.csv(file.path(matching_output_dir, "matching_usage_all.csv")) |>
  filter(pool == "20_100")

m_tree <- readRDS(file.path(matching_output_dir, "matching_tree_20_100.rds"))
m_crop <- readRDS(file.path(matching_output_dir, "matching_crop_20_100.rds"))
m_built <- readRDS(file.path(matching_output_dir, "matching_built_20_100.rds"))

# Grid and cell extraction

grid_500 <- readRDS(file.path(matching_data_dir, "grid_500m_mollweide.rds"))
grid_extent <- ext(vect(grid_500))

r_500 <- rast(
  xmin = xmin(grid_extent),
  xmax = xmax(grid_extent),
  ymin = ymin(grid_extent),
  ymax = ymax(grid_extent),
  resolution = 500,
  crs = st_crs(grid_500)$wkt
)

grid_points <- vect(st_centroid(grid_500))
grid_ids <- st_drop_geometry(grid_500) |> select(global_id)

extract_outcome <- function(r, name) {
  values <- terra::extract(r, grid_points)
  result <- grid_ids
  result[[name]] <- values[[2]]
  result
}

get_matched_pairs <- function(match_object, outcome_table, outcome_col) {

  md <- MatchIt::match.data(match_object, drop.unmatched = TRUE) |>
    left_join(outcome_table, by = "global_id")

  md |>
    select(subclass, treated, global_id, all_of(outcome_col)) |>
    mutate(group = if_else(treated == 1, "inside", "control")) |>
    select(subclass, group, all_of(outcome_col)) |>
    pivot_wider(
      names_from = group,
      values_from = all_of(outcome_col)
    ) |>
    mutate(difference = inside - control)
}

# Matched differences and confidence intervals

calculate_matched_effects <- function(matches, outcome_table, outcome_col, analysis) {

  results <- list()
  paired_data <- list()

  for (park in all_parks) {

    m <- matches[[park]]
    if (is.null(m)) next

    pairs <- get_matched_pairs(m, outcome_table, outcome_col)
    n_pairs_total <- nrow(pairs)

    pairs_valid <- pairs[
      is.finite(pairs$inside) & is.finite(pairs$control),
      , drop = FALSE
    ]

    n_pairs <- nrow(pairs_valid)
    if (n_pairs == 0) next

    differences <- pairs_valid$difference
    effect <- mean(differences)

    if (n_pairs >= 2) {
      se <- sd(differences) / sqrt(n_pairs)
      ci_halfwidth <- qt(0.975, df = n_pairs - 1) * se
      ci_lower <- effect - ci_halfwidth
      ci_upper <- effect + ci_halfwidth
    } else {
      se <- NA_real_
      ci_lower <- NA_real_
      ci_upper <- NA_real_
    }

    results[[park]] <- data.frame(
      NAME_ENG = park,
      analysis = analysis,
      inside_pp = mean(pairs_valid$inside),
      control_pp = mean(pairs_valid$control),
      effect_pp = effect,
      se_pp = se,
      ci_lower_pp = ci_lower,
      ci_upper_pp = ci_upper,
      matched_pairs_total = n_pairs_total,
      matched_pairs_with_outcome = n_pairs,
      matched_pairs_missing_outcome = n_pairs_total - n_pairs
    )

    pairs_valid$NAME_ENG <- park
    pairs_valid$analysis <- analysis
    paired_data[[park]] <- pairs_valid

    message(
      analysis, " | ", park,
      " | Pairs: ", n_pairs,
      " | Effect: ", round(effect, 4)
    )
  }

  list(
    effects = bind_rows(results),
    pairs = bind_rows(paired_data)
  )
}

# Tree-cover loss, 2001-2025

treecover <- rast(treecover_file)
lossyear <- rast(lossyear_file)

forest_2000 <- ifel(
  is.na(treecover), NA,
  ifel(treecover >= 30, 1, 0)
)

forest_lost <- ifel(
  forest_2000 == 1 & lossyear >= 2001 & lossyear <= 2025,
  1, 0
)

forest_lost <- ifel(
  is.na(forest_2000) | is.na(lossyear),
  NA, forest_lost
)

# Forest loss as percentage of total cell area

forest_lost_500 <- project(
  forest_lost,
  r_500,
  method = "average",
  filename = file.path(output_dir, "forest_lost_fraction_of_cell_500m_mollweide.tif"),
  overwrite = TRUE
)

forest_loss_cell_pp <- forest_lost_500 * 100
names(forest_loss_cell_pp) <- "forest_loss_cell_pp"

writeRaster(
  forest_loss_cell_pp,
  file.path(output_dir, "forest_loss_pct_of_cell_500m_mollweide.tif"),
  overwrite = TRUE
)

tree_cell_values <- extract_outcome(forest_loss_cell_pp, "tree_loss_pp")

tree_results <- calculate_matched_effects(
  m_tree, tree_cell_values, "tree_loss_pp", "tree"
)

tree_effects <- tree_results$effects

# Forest loss relative to baseline forest cover

forest_baseline_500 <- project(
  forest_2000,
  r_500,
  method = "average",
  filename = file.path(output_dir, "forest_baseline_fraction_500m_mollweide.tif"),
  overwrite = TRUE
)

forest_loss_baseline_pct <- ifel(
  forest_baseline_500 > 0,
  100 * forest_lost_500 / forest_baseline_500,
  NA
)

names(forest_loss_baseline_pct) <- "forest_loss_baseline_pct"

writeRaster(
  forest_loss_baseline_pct,
  file.path(output_dir, "forest_loss_pct_of_baseline_500m_mollweide.tif"),
  overwrite = TRUE
)

tree_sensitivity_values <- extract_outcome(
  forest_loss_baseline_pct, "tree_loss_pp"
)

tree_sensitivity <- calculate_matched_effects(
  m_tree, tree_sensitivity_values, "tree_loss_pp", "tree"
)

write.csv(
  tree_sensitivity$effects,
  file.path(output_dir, "SENSITIVITY_PA_forest_loss_pct_of_baseline.csv"),
  row.names = FALSE
)

saveRDS(
  tree_sensitivity$pairs,
  file.path(output_dir, "SENSITIVITY_matched_pairs_tree_baseline.rds")
)

saveRDS(
  tree_sensitivity_values,
  file.path(output_dir, "SENSITIVITY_forest_baseline_cell_values.rds")
)

rm(
  treecover, lossyear, forest_2000, forest_lost,
  forest_baseline_500, forest_lost_500,
  forest_loss_cell_pp, forest_loss_baseline_pct
)
gc()

# Cropland change, 2000-2024

gaced <- rast(gaced_file)

crop_2000 <- ifel(
  gaced[["crop_2000"]] == 10, 1,
  ifel(gaced[["crop_2000"]] == 0, 0, NA)
)

crop_2024 <- ifel(
  gaced[["crop_2024"]] == 10, 1,
  ifel(gaced[["crop_2024"]] == 0, 0, NA)
)

crop_change <- crop_2024 - crop_2000

crop_change_500 <- project(
  crop_change,
  r_500,
  method = "average",
  filename = file.path(output_dir, "crop_change_500m_mollweide.tif"),
  overwrite = TRUE
)

crop_cell_values <- extract_outcome(crop_change_500, "crop_change_pp") |>
  mutate(crop_change_pp = crop_change_pp * 100)

crop_results <- calculate_matched_effects(
  m_crop, crop_cell_values, "crop_change_pp", "crop"
)

crop_effects <- crop_results$effects

rm(gaced, crop_2000, crop_2024, crop_change, crop_change_500)
gc()

# Built-up surface change, 2000-2025

built_files <- list.files(
  built_dir,
  pattern = "^Mozambique_GHSL_built_surface_.*\\.tif$",
  full.names = TRUE
)

built_years <- as.integer(
  sub(".*?(20[0-9]{2}).*", "\\1", basename(built_files))
)

last_built_year <- max(built_years)

built_2000 <- rast(built_files[built_years == 2000])[[1]]
built_last <- rast(built_files[built_years == last_built_year])[[1]]

built_2000 <- ifel(built_2000 > 10000, NA, built_2000)
built_last <- ifel(built_last > 10000, NA, built_last)

built_2000_pct <- built_2000 / 10000 * 100
built_last_pct <- built_last / 10000 * 100

built_change_100m <- built_last_pct - built_2000_pct

built_change_500m_native <- aggregate(
  built_change_100m,
  fact = 5,
  fun = "mean",
  na.rm = TRUE
)

built_change_500 <- project(
  built_change_500m_native,
  r_500,
  method = "bilinear",
  filename = file.path(
    output_dir,
    paste0("built_change_2000_", last_built_year, "_500m_mollweide.tif")
  ),
  overwrite = TRUE
)

built_cell_values <- extract_outcome(
  built_change_500, "built_change_pp"
)

built_results <- calculate_matched_effects(
  m_built, built_cell_values, "built_change_pp", "built"
)

built_effects <- built_results$effects

rm(
  built_2000, built_last, built_2000_pct, built_last_pct,
  built_change_100m, built_change_500m_native, built_change_500
)
gc()

# Combined results

effects_all <- bind_rows(
  tree_effects, crop_effects, built_effects
) |>
  complete(
    NAME_ENG = all_parks,
    analysis = c("tree", "crop", "built")
  )

usage_summary <- usage |>
  select(
    NAME_ENG, analysis, treated_total, treated_matched,
    treated_retained_percent, unique_controls,
    control_effective_sample_size, mean_abs_SMD, max_abs_SMD
  ) |>
  distinct(NAME_ENG, analysis, .keep_all = TRUE)

effects_all <- effects_all |>
  left_join(usage_summary, by = c("NAME_ENG", "analysis")) |>
  mutate(
    balance_025 = if_else(
      is.na(max_abs_SMD), NA, max_abs_SMD <= main_balance_threshold
    ),
    balance_010 = if_else(
      is.na(max_abs_SMD), NA, max_abs_SMD <= strict_balance_threshold
    ),
    pressure = case_when(
      analysis == "tree" ~ "a. Baseline forest loss",
      analysis == "crop" ~ "b. Cropland expansion",
      analysis == "built" ~ "c. Built-up expansion"
    ),
    status = case_when(
      is.na(effect_pp) ~ "Effect unavailable",
      is.na(max_abs_SMD) ~ "Balance unavailable",
      balance_025 ~ "Adequate balance",
      TRUE ~ "Inadequate balance"
    )
  ) |>
  arrange(pressure, match(NAME_ENG, all_parks))

# Export results

write.csv(
  tree_effects,
  file.path(output_dir, "PA_matched_effects_treecover_2000_2025.csv"),
  row.names = FALSE
)

write.csv(
  crop_effects,
  file.path(output_dir, "PA_matched_effects_cropland_2000_2024.csv"),
  row.names = FALSE
)

write.csv(
  built_effects,
  file.path(output_dir, paste0("PA_matched_effects_built_2000_", last_built_year, ".csv")),
  row.names = FALSE
)

write.csv(
  effects_all,
  file.path(output_dir, "PA_matching_effects_all_pressures.csv"),
  row.names = FALSE
)

write.csv(
  effects_all,
  file.path(output_dir, "effects_all.csv"),
  row.names = FALSE
)

saveRDS(
  effects_all,
  file.path(output_dir, "PA_matching_effects_all_pressures.rds")
)

saveRDS(
  tree_cell_values,
  file.path(output_dir, "Hansen_tree_loss_matching_cell_values.rds")
)

saveRDS(
  crop_cell_values,
  file.path(output_dir, "GACED_crop_change_matching_cell_values.rds")
)

saveRDS(
  built_cell_values,
  file.path(output_dir, "GHSL_built_change_matching_cell_values.rds")
)

saveRDS(tree_results$pairs, file.path(output_dir, "matched_pairs_tree.rds"))
saveRDS(crop_results$pairs, file.path(output_dir, "matched_pairs_crop.rds"))
saveRDS(built_results$pairs, file.path(output_dir, "matched_pairs_built.rds"))
