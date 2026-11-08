# 03. Matching without replacement - 500 m sensitivity analysis

library(MatchIt)
library(cobalt)
library(dplyr)

input_dir <- "data/processed/matching_GACED30/revised_500m_forest_proportion_exclusion_exception"
output_dir <- "outputs/matching/revised_no_replacement_500m_forest_proportion_exclusion_exception"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(file.path(input_dir, "matching_input_revised_500m.rds"))

analysis_spec <- list(
  tree = c("elevation", "slope", "dist_roads", "dist_settlements", "treecover_2000"),
  crop = c("elevation", "slope", "dist_roads", "dist_settlements", "cropland_baseline"),
  built = c("elevation", "slope", "dist_roads", "dist_settlements", "baseline_built")
)

park_names <- sort(unique(as.character(dat$NAME_ENG)))

run_matching <- function(d0, covars, park, analysis) {

  d <- d0[d0$treated | (!d0$treated & d0$pool_20_100), , drop = FALSE]
  d <- d[order(d$global_id, d$treated), , drop = FALSE]

  n_treat_total <- sum(d$treated)
  n_control_total <- sum(!d$treated)

  complete <- complete.cases(d[, covars, drop = FALSE])
  n_treat_complete <- sum(d$treated & complete)
  n_control_complete <- sum(!d$treated & complete)

  d <- d[complete, , drop = FALSE]
  d$treated <- as.integer(d$treated)
  rownames(d) <- paste0("cell_", d$global_id)

  m <- MatchIt::matchit(
    formula = reformulate(covars, response = "treated"),
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

  matched_controls <- as.character(m$match.matrix[, 1])
  matched_controls <- matched_controls[!is.na(matched_controls)]

  n_matched <- length(matched_controls)

  b <- cobalt::bal.tab(m, un = TRUE, stats = "mean.diffs", binary = "std")
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
    treated_matched = n_matched,
    treated_unmatched = n_treat_complete - n_matched,
    treated_retained_percent = 100 * n_matched / n_treat_complete,
    unique_controls = length(unique(matched_controls)),
    control_effective_sample_size = n_matched,
    max_control_reuse = as.numeric(n_matched > 0),
    mean_abs_SMD = mean(smd, na.rm = TRUE),
    max_abs_SMD = max(smd, na.rm = TRUE),
    n_over_0_10 = sum(smd > 0.10, na.rm = TRUE),
    n_over_0_25 = sum(smd > 0.25, na.rm = TRUE),
    status = if (n_matched > 0) "matched" else "no_matches",
    detail = ""
  )

  list(match = m, usage = usage, balance = cov_balance)
}

usage_list <- list()
balance_list <- list()

for (analysis in names(analysis_spec)) {

  covars <- analysis_spec[[analysis]]
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
}

usage_all <- bind_rows(usage_list)
balance_all <- bind_rows(balance_list)

write.csv(
  usage_all,
  file.path(output_dir, "matching_usage_all.csv"),
  row.names = FALSE
)

write.csv(
  balance_all,
  file.path(output_dir, "covariate_balance_all.csv"),
  row.names = FALSE
)

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
