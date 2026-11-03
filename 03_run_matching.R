# 03. MATCHING WITHOUT REPLACEMENT

library(MatchIt)
library(cobalt)
library(dplyr)

# Directories
input_dir <- "data/processed/matching_GACED30/revised_20_50_100_forest_proportion_exclusion_exception"
output_dir <- "outputs/matching/revised_no_replacement_forest_proportion_exclusion_exception"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(file.path(input_dir, "matching_input_revised_1km.rds"))

# Matching covariates
analysis_spec <- list(
  tree = list(
    label = "Tree-cover loss",
    covars = c(
      "elevation", "slope", "dist_roads",
      "dist_settlements", "treecover_2000"
    )
  ),
  crop = list(
    label = "Cropland expansion",
    covars = c(
      "elevation", "slope", "dist_roads",
      "dist_settlements", "cropland_baseline"
    )
  ),
  built = list(
    label = "Built-up expansion",
    covars = c(
      "elevation", "slope", "dist_roads",
      "dist_settlements", "baseline_built"
    )
  )
)

park_names <- sort(unique(as.character(dat$NAME_ENG)))

# Matching function
run_matching <- function(d0, covars, park, analysis) {

  d <- d0[
    d0$treated | (!d0$treated & d0$pool_20_100),
    , drop = FALSE
  ]

  d <- d[order(d$global_id, d$treated), , drop = FALSE]

  n_treat_total <- sum(d$treated)
  n_control_total <- sum(!d$treated)

  complete <- complete.cases(d[, covars, drop = FALSE])

  n_treat_complete <- sum(d$treated & complete)
  n_control_complete <- sum(!d$treated & complete)

  d <- d[complete, , drop = FALSE]
  d$treated <- as.integer(d$treated)
  rownames(d) <- paste0("cell_", d$global_id)

  # Matching diagnostics
  usage <- data.frame(
    NAME_ENG = park,
    analysis = analysis,
    pool = "20_100",
    treated_total = n_treat_total,
    treated_complete = n_treat_complete,
    treated_missing_covariates = n_treat_total - n_treat_complete,
    controls_eligible = n_control_total,
    controls_complete = n_control_complete,
    controls_missing_covariates = n_control_total - n_control_complete,
    treated_matched = NA_integer_,
    treated_unmatched = NA_integer_,
    treated_retained_percent = NA_real_,
    unique_controls = NA_integer_,
    control_effective_sample_size = NA_real_,
    max_control_reuse = NA_real_,
    mean_abs_SMD = NA_real_,
    max_abs_SMD = NA_real_,
    n_over_0_10 = NA_integer_,
    n_over_0_25 = NA_integer_,
    status = "not_run",
    detail = "",
    stringsAsFactors = FALSE
  )

  # Matching
  form <- reformulate(covars, response = "treated")

  m <- MatchIt::matchit(
    formula = form,
    data = d,
    method = "nearest",
    distance = "glm",
    mahvars = covars,
    caliper = 0.2,
    std.caliper = TRUE,
    replace = FALSE,
    ratio = 1,
    m.order = "closest"
  )

  # Matched pairs and control usage
  matched_controls <- as.character(m$match.matrix[, 1])
  matched_controls <- matched_controls[!is.na(matched_controls)]

  reuse <- table(matched_controls)

  usage$treated_matched <- length(matched_controls)
  usage$treated_unmatched <- n_treat_complete - usage$treated_matched
  usage$treated_retained_percent <- 100 * usage$treated_matched / n_treat_complete
  usage$unique_controls <- length(reuse)

  usage$control_effective_sample_size <- if (length(reuse) > 0) {
    sum(as.numeric(reuse))^2 / sum(as.numeric(reuse)^2)
  } else {
    0
  }

  usage$max_control_reuse <- if (length(reuse) > 0) {
    max(as.numeric(reuse))
  } else {
    0
  }

  usage$status <- if (usage$treated_matched > 0) "matched" else "no_matches"

  # Covariate balance
  b <- cobalt::bal.tab(
    m,
    un = TRUE,
    stats = "mean.diffs",
    binary = "std"
  )

  bal <- b$Balance
  bal <- bal[rownames(bal) != "distance", , drop = FALSE]

  cov_balance <- data.frame(
    NAME_ENG = park,
    analysis = analysis,
    pool = "20_100",
    covariate = rownames(bal),
    SMD_before = as.numeric(bal$Diff.Un),
    SMD_after = as.numeric(bal$Diff.Adj),
    abs_SMD_after = abs(as.numeric(bal$Diff.Adj)),
    row.names = NULL
  )

  smd <- cov_balance$abs_SMD_after

  usage$mean_abs_SMD <- mean(smd, na.rm = TRUE)
  usage$max_abs_SMD <- max(smd, na.rm = TRUE)
  usage$n_over_0_10 <- sum(smd > 0.10, na.rm = TRUE)
  usage$n_over_0_25 <- sum(smd > 0.25, na.rm = TRUE)

  list(
    match = m,
    usage = usage,
    balance = cov_balance
  )
}

# Run matching
usage_list <- list()
balance_list <- list()

for (analysis in names(analysis_spec)) {

  covars <- analysis_spec[[analysis]]$covars
  matches <- setNames(vector("list", length(park_names)), park_names)

  message("Matching: ", analysis)

  for (park in park_names) {

    d0 <- dat[dat$NAME_ENG == park, , drop = FALSE]

    res <- run_matching(d0, covars, park, analysis)

    matches[[park]] <- res$match
    usage_list[[length(usage_list) + 1L]] <- res$usage
    balance_list[[length(balance_list) + 1L]] <- res$balance

    message(
      park,
      " | Matched: ", res$usage$treated_matched,
      " | Max SMD: ", round(res$usage$max_abs_SMD, 3)
    )
  }

  saveRDS(
    matches,
    file.path(output_dir, paste0("matching_", analysis, "_20_100.rds"))
  )

  write.csv(
    bind_rows(usage_list),
    file.path(output_dir, "matching_usage_all.csv"),
    row.names = FALSE
  )

  write.csv(
    bind_rows(balance_list),
    file.path(output_dir, "covariate_balance_all.csv"),
    row.names = FALSE
  )

  gc()
}

# Final diagnostics
usage_all <- bind_rows(usage_list)
balance_all <- bind_rows(balance_list)

saveRDS(
  list(
    usage = usage_all,
    balance = balance_all,
    analysis_spec = analysis_spec,
    method = "nearest_without_replacement",
    pool = "20_100"
  ),
  file.path(output_dir, "matching_diagnostics_all.rds")
)

# Comparison with previous WDPA exclusion
original_file <- file.path(
  "outputs/matching/revised_no_replacement_forest_proportion",
  "matching_usage_all.csv"
)

original <- read.csv(original_file) |>
  filter(pool == "20_100") |>
  mutate(covariate_version = "Previous WDPA exclusion")

alternative <- usage_all |>
  mutate(covariate_version = "WDPA exclusion exception")

comparison <- bind_rows(original, alternative) |>
  mutate(
    retention_pct = if_else(
      treated_total > 0,
      100 * treated_matched / treated_total,
      NA_real_
    )
  ) |>
  select(
    NAME_ENG, analysis, covariate_version, treated_matched,
    retention_pct, unique_controls,
    control_effective_sample_size, max_abs_SMD, status
  ) |>
  arrange(analysis, NAME_ENG, covariate_version)

write.csv(
  comparison,
  file.path(output_dir, "matching_wdpa_exclusion_comparison.csv"),
  row.names = FALSE
)
