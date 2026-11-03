#08. ROSENBAUM SENSITIVITY ANALYSIS

library(rbounds)
library(dplyr)

input_dir <- "outputs/results/revised_no_replacement_forest_proportion_cell_area_exclusion_exception"
matching_dir <- "outputs/matching/revised_no_replacement_forest_proportion_exclusion_exception"
output_dir <- file.path(input_dir, "rosenbaum_sensitivity", "final")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

max_gamma <- 3
step_gamma <- 0.1
alpha <- 0.05
balance_main <- 0.25
balance_strict <- 0.10

all_parks <- c(
  "Banhine", "Chimanimani", "Gilé", "Gorongosa", "Limpopo",
  "Magoe", "Maputo", "Marromeu", "Niassa", "Pomene",
  "Quirimbas", "Serra da Gorongosa", "Zinave"
)
pressures <- c("tree", "crop", "built")

usage <- read.csv(file.path(matching_dir, "matching_usage_all.csv"))
if ("pool" %in% names(usage)) usage <- usage[usage$pool == "20_100", ]

# Rosenbaum bounds for one park and pressure
run_rosenbaum <- function(d, park, analysis, smd) {

  n_total <- nrow(d)
  d <- d[is.finite(d$inside) & is.finite(d$control), ]
  diffs <- d$inside - d$control

  n_valid <- length(diffs)
  n_nonzero <- sum(diffs != 0)
  mean_diff <- if (n_valid) mean(diffs) else NA_real_

  direction <- if (is.na(mean_diff) || mean_diff == 0) NA_character_ else
    if (mean_diff > 0) "protected_higher" else "protected_lower"

  status <- "ok"
  crossing <- p_gamma1 <- p_gamma3 <- p_wilcox <- gamma1_difference <- NA_real_
  bounds <- NULL

  if (n_valid < 2) {
    status <- "fewer_than_two_pairs"
  } else if (n_nonzero < 2) {
    status <- "fewer_than_two_nonzero_differences"
  } else {

    result <- tryCatch({

      rb_plus <- rbounds::psens(d$inside, d$control,
                                Gamma = max_gamma, GammaInc = step_gamma)
      rb_minus <- rbounds::psens(d$control, d$inside,
                                 Gamma = max_gamma, GammaInc = step_gamma)

      b1 <- as.data.frame(rb_plus$bounds)
      b2 <- as.data.frame(rb_minus$bounds)

      names(b1)[1:3] <- c("Gamma", "p_lower_plus", "p_upper_plus")
      names(b2)[1:3] <- c("Gamma", "p_lower_minus", "p_upper_minus")

      b <- merge(b1[, 1:3], b2[, 1:3], by = "Gamma")
      b[] <- lapply(b, as.numeric)

      b$p_upper_two_sided <- pmin(
        1, 2 * pmin(b$p_upper_plus, b$p_upper_minus)
      )

      b$NAME_ENG <- park
      b$analysis <- analysis

      i1 <- which.min(abs(b$Gamma - 1))
      i3 <- which.min(abs(b$Gamma - max_gamma))

      p_w <- wilcox.test(d$inside, d$control, paired = TRUE,
                         alternative = "two.sided", exact = FALSE,
                         correct = FALSE)$p.value

      hit <- which(b$p_upper_two_sided >= alpha)

      list(
        bounds = b,
        p1 = b$p_upper_two_sided[i1],
        p3 = b$p_upper_two_sided[i3],
        wilcox = p_w,
        difference = abs(p_w - b$p_upper_two_sided[i1]),
        crossing = if (length(hit)) b$Gamma[hit[1]] else NA_real_
      )

    }, error = function(e) e)

    if (inherits(result, "error")) {
      status <- paste0("rbounds_error: ", conditionMessage(result))
    } else {
      bounds <- result$bounds
      p_gamma1 <- result$p1
      p_gamma3 <- result$p3
      p_wilcox <- result$wilcox
      gamma1_difference <- result$difference
      crossing <- result$crossing
    }
  }

  threshold_note <- if (status != "ok") status else
    if (is.na(crossing)) paste0("No crossing through Gamma=", max_gamma) else
      if (crossing == 1) "Already p >= 0.05 at Gamma=1" else
        "First crossing on tested Gamma grid"

  summary <- data.frame(
    NAME_ENG = park,
    analysis = analysis,
    matched_pairs_in_file = n_total,
    matched_pairs_with_finite_outcomes = n_valid,
    nonzero_differences = n_nonzero,
    mean_effect_pp = mean_diff,
    effect_direction = direction,
    max_abs_SMD = smd,
    balance_025 = if (is.na(smd)) NA else smd <= balance_main,
    balance_010 = if (is.na(smd)) NA else smd <= balance_strict,
    p_wilcox_Gamma_1 = p_wilcox,
    p_upper_two_sided_Gamma_1 = p_gamma1,
    gamma1_absolute_difference = gamma1_difference,
    p_upper_two_sided_Gamma_3 = p_gamma3,
    first_Gamma_two_sided_p_upper_ge_0_05 = crossing,
    gamma_threshold_note = threshold_note,
    status = status
  )

  list(summary = summary, bounds = bounds)
}

# Run analysis
summaries <- list()
curves <- list()

for (analysis in pressures) {

  dat <- readRDS(file.path(input_dir, paste0("matched_pairs_", analysis, ".rds")))
  message("\nRosenbaum analysis: ", analysis)

  for (park in all_parks) {

    d <- dat[as.character(dat$NAME_ENG) == park, ]

    u <- usage[
      as.character(usage$NAME_ENG) == park &
        as.character(usage$analysis) == analysis, ]

    smd <- if (nrow(u)) as.numeric(u$max_abs_SMD[1]) else NA_real_
    res <- run_rosenbaum(d, park, analysis, smd)

    summaries[[length(summaries) + 1L]] <- res$summary
    if (!is.null(res$bounds)) curves[[length(curves) + 1L]] <- res$bounds

    message(park, " | Pairs: ", res$summary$matched_pairs_with_finite_outcomes,
            " | Gamma: ", res$summary$first_Gamma_two_sided_p_upper_ge_0_05,
            " | ", res$summary$status)
  }
}

# Save results
summary_all <- bind_rows(summaries)

curve_all <- if (length(curves)) bind_rows(curves) else data.frame(
  Gamma = numeric(), p_lower_plus = numeric(), p_upper_plus = numeric(),
  p_lower_minus = numeric(), p_upper_minus = numeric(),
  p_upper_two_sided = numeric(), NAME_ENG = character(), analysis = character()
)

write.csv(summary_all, file.path(output_dir, "rosenbaum_summary_final.csv"),
          row.names = FALSE)
write.csv(curve_all, file.path(output_dir, "rosenbaum_pvalue_bounds_final.csv"),
          row.names = FALSE)

saveRDS(
  list(summary = summary_all, bounds = curve_all,
       max_gamma = max_gamma, step_gamma = step_gamma, alpha = alpha),
  file.path(output_dir, "rosenbaum_results_final.rds")
)




