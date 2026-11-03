# 09. SUPPLEMENTARY TABLES


library(sf)
library(dplyr)
library(tidyr)

output_dir <- "outputs/supplementary_tables/wdpa_exclusion_exception"
matching_dir <- "outputs/matching/revised_no_replacement_forest_proportion_exclusion_exception"
results_dir <- "outputs/results/revised_no_replacement_forest_proportion_cell_area_exclusion_exception"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Protected areas and surrounding landscapes
load("data/processed/parks_core.Rdata")

parks_core <- parks_core |>
  filter(DESIG %in% c("Parque Nacional", "National Park",
                      "Reserva Nacional", "Reserva Especial")) |>
  st_make_valid()

rings <- st_read(
  "data/processed/Mozambique_PA_equal_area_surroundings_nationally_clipped.shp",
  quiet = TRUE
) |>
  rename(NAME_ENG = park_name) |>
  filter(NAME_ENG %in% parks_core$NAME_ENG) |>
  st_make_valid()

moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

pa_area <- parks_core |>
  st_transform(moll) |>
  mutate(protected_area_ha = as.numeric(st_area(geometry)) / 10000) |>
  st_drop_geometry() |>
  select(NAME_ENG, DESIG, protected_area_ha)

ring_area <- rings |>
  st_transform(moll) |>
  mutate(surrounding_area_ha = as.numeric(st_area(geometry)) / 10000) |>
  st_drop_geometry() |>
  select(NAME_ENG, surrounding_area_ha)

table_pa_areas <- pa_area |>
  left_join(ring_area, by = "NAME_ENG") |>
  transmute(
    `Protected area` = NAME_ENG,
    Designation = DESIG,
    `Protected area (ha)` = round(protected_area_ha),
    `Surrounding landscape (ha)` = round(surrounding_area_ha)
  ) |>
  arrange(`Protected area`)

write.csv(table_pa_areas,
          file.path(output_dir, "protected_areas_and_surrounding_areas.csv"),
          row.names = FALSE)


# Temporal land-use change
temporal_dir <- "outputs/temporal"

get_final_change <- function(file, year, prefix) {
  read.csv(file.path(temporal_dir, file)) |>
    filter(.data$year == .env$year) |>
    select(NAME_ENG, location, change_pp) |>
    pivot_wider(names_from = location, values_from = change_pp,
                names_prefix = prefix)
}

tree_final <- get_final_change(
  "PA_treecover_temporal_inside_outside.csv", 2025, "tree_"
)

crop_final <- get_final_change(
  "PA_cropland_temporal_inside_outside_gaced.csv", 2024, "crop_"
)

built_final <- get_final_change(
  "PA_builtup_temporal_inside_outside.csv", 2025, "built_"
)

table_landuse <- tree_final |>
  full_join(crop_final, by = "NAME_ENG") |>
  full_join(built_final, by = "NAME_ENG") |>
  transmute(
    `Protected area` = NAME_ENG,
    `Tree-cover loss within PA (pp)` = round(tree_Inside, 4),
    `Tree-cover loss surrounding PA (pp)` = round(tree_Outside, 4),
    `Cropland change within PA (pp)` = round(crop_Inside, 4),
    `Cropland change surrounding PA (pp)` = round(crop_Outside, 4),
    `Built-up change within PA (pp)` = round(built_Inside, 4),
    `Built-up change surrounding PA (pp)` = round(built_Outside, 4)
  ) |>
  arrange(`Protected area`)

write.csv(table_landuse,
          file.path(output_dir, "land_use_change_within_and_around_protected_areas_fixed.csv"),
          row.names = FALSE)


# Matching balance
usage <- read.csv(file.path(matching_dir, "matching_usage_all.csv"))

table_balance <- usage |>
  mutate(
    Pressure = recode(analysis,
                      tree = "Tree-cover loss",
                      crop = "Cropland expansion",
                      built = "Built-up expansion"),
    `Main balance criterion` = case_when(
      is.na(max_abs_SMD) ~ NA_character_,
      max_abs_SMD <= 0.25 ~ "Yes",
      TRUE ~ "No"
    ),
    `Strict balance criterion` = case_when(
      is.na(max_abs_SMD) ~ NA_character_,
      max_abs_SMD <= 0.10 ~ "Yes",
      TRUE ~ "No"
    )
  ) |>
  transmute(
    `Protected area` = NAME_ENG,
    Pressure,
    `Mean absolute SMD` = round(mean_abs_SMD, 4),
    `Maximum absolute SMD` = round(max_abs_SMD, 4),
    `Main balance criterion`,
    `Strict balance criterion`
  ) |>
  arrange(Pressure, `Protected area`)

write.csv(table_balance,
          file.path(output_dir, "matching_covariate_balance.csv"),
          row.names = FALSE)


# Matched protected-area effects
effects <- read.csv(
  file.path(results_dir, "PA_matching_effects_all_pressures.csv")
)

table_effects <- effects |>
  mutate(
    Pressure = recode(analysis,
                      tree = "Tree-cover loss",
                      crop = "Cropland expansion",
                      built = "Built-up expansion")
  ) |>
  transmute(
    `Protected area` = NAME_ENG,
    Pressure,
    `Difference (PA - control, pp)` = round(effect_pp, 4),
    `95% CI lower (pp)` = round(ci_lower_pp, 4),
    `95% CI upper (pp)` = round(ci_upper_pp, 4),
    `Matched pairs` = matched_pairs_total,
    `Mean absolute SMD` = round(mean_abs_SMD, 4),
    `Maximum absolute SMD` = round(max_abs_SMD, 4),
    `Adequate balance` = case_when(
      is.na(balance_025) ~ NA_character_,
      balance_025 ~ "Yes",
      TRUE ~ "No"
    ),
    Status = status
  ) |>
  arrange(`Protected area`, Pressure)

write.csv(table_effects,
          file.path(output_dir, "matched_protected_area_effects.csv"),
          row.names = FALSE)


# Independent dataset comparison
comparison <- read.csv(
  "outputs/dataset_comparison/independent_dataset_comparison_summary_fixed.csv"
)

table_comparison <- comparison |>
  transmute(
    Pressure = comparison,
    `Primary product` = y_product,
    `Independent product` = x_product,
    `Number of cells` = n_cells,
    `Pearson's r` = round(pearson_r, 3),
    `R-squared` = round(R2, 3),
    Intercept = round(intercept, 3),
    Slope = round(slope, 3)
  )

write.csv(table_comparison,
          file.path(output_dir, "independent_dataset_comparison.csv"),
          row.names = FALSE)

