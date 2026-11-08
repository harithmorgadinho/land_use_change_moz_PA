# 05. Spatial resolution sensitivity - 1 km vs 500 m

library(dplyr)
library(tidyr)

results_1km <- "outputs/results/revised_no_replacement_forest_proportion_cell_area_exclusion_exception/PA_matching_effects_all_pressures.csv"
results_500m <- "outputs/results/revised_no_replacement_500m_forest_proportion_cell_area_exclusion_exception/PA_matching_effects_all_pressures.csv"

out_dir <- "outputs/results/sensitivity_500m_vs_1km"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c(
  "NAME_ENG", "analysis", "effect_pp", "ci_lower_pp", "ci_upper_pp",
  "matched_pairs_with_outcome", "treated_matched", "max_abs_SMD", "status"
)

read_effects <- function(path, suffix) {

  x <- read.csv(path, check.names = FALSE)
  x <- x[, cols, drop = FALSE]

  names(x)[!names(x) %in% c("NAME_ENG", "analysis")] <-
    paste0(names(x)[!names(x) %in% c("NAME_ENG", "analysis")], "_", suffix)

  x
}

one_km <- read_effects(results_1km, "1km")
five_hundred <- read_effects(results_500m, "500m")

# Compare matched effects

comparison <- full_join(
  one_km, five_hundred,
  by = c("NAME_ENG", "analysis")
) |>
  mutate(
    difference_500m_minus_1km_pp = effect_pp_500m - effect_pp_1km,

    same_effect_direction = case_when(
      is.na(effect_pp_1km) | is.na(effect_pp_500m) ~ NA,
      effect_pp_1km == 0 | effect_pp_500m == 0 ~
        effect_pp_1km == effect_pp_500m,
      TRUE ~ sign(effect_pp_1km) == sign(effect_pp_500m)
    ),

    adequate_balance_both_025 =
      !is.na(max_abs_SMD_1km) &
      !is.na(max_abs_SMD_500m) &
      max_abs_SMD_1km <= 0.25 &
      max_abs_SMD_500m <= 0.25,

    adequate_balance_both_010 =
      !is.na(max_abs_SMD_1km) &
      !is.na(max_abs_SMD_500m) &
      max_abs_SMD_1km <= 0.10 &
      max_abs_SMD_500m <= 0.10
  ) |>
  arrange(match(analysis, c("tree", "crop", "built")), NAME_ENG)

write.csv(
  comparison,
  file.path(out_dir, "supplementary_matching_resolution_comparison.csv"),
  row.names = FALSE,
  na = ""
)

# Summary by land-use pressure

summary <- comparison |>
  group_by(analysis) |>
  summarise(
    n_with_both_effects = sum(is.finite(effect_pp_1km) & is.finite(effect_pp_500m)),
    n_same_direction = sum(same_effect_direction, na.rm = TRUE),
    n_adequate_balance_both_025 = sum(adequate_balance_both_025),
    median_abs_effect_difference_pp = median(
      abs(difference_500m_minus_1km_pp), na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  summary,
  file.path(out_dir, "matching_resolution_summary.csv"),
  row.names = FALSE
)

