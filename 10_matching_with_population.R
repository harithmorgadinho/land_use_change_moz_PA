
# Population-density sensitivity analysis

library(MatchIt)
library(cobalt)
library(dplyr)
library(tidyr)

main_input <- "data/processed/matching_GACED30/revised_20_50_100_forest_proportion_exclusion_exception/matching_input_revised_1km.rds"
old_pop_input <- "data/processed/matching_GACED30/revised_20_50_100_population/matching_input_revised_1km.rds"

out <- "outputs/matching/population_sensitivity_exclusion_exception"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(main_input)
popdat <- readRDS(old_pop_input)

key <- c("NAME_ENG", "global_id")
popvar <- "population_density_2000"

needed <- c(
  key, "treated", "pool_20_100",
  "elevation", "slope", "dist_roads", "dist_settlements",
  "treecover_2000", "cropland_baseline", "baseline_built"
)

if (length(setdiff(needed, names(dat)))) {
  stop("Current grid lacks: ", paste(setdiff(needed, names(dat)), collapse = ", "))
}

if (length(setdiff(c(key, popvar), names(popdat)))) {
  stop("Population grid lacks required key/population fields")
}

if (anyDuplicated(dat[key]) || anyDuplicated(popdat[key])) {
  stop("Duplicate NAME_ENG/global_id assignments")
}

if (anyNA(dat$treated) || anyNA(dat$pool_20_100)) {
  stop("Missing treatment/pool flags")
}

if (any(dat$treated & dat$pool_20_100)) {
  stop("Treated cells marked as controls")
}

# Transfer population density without changing current grid eligibility

if (!(popvar %in% names(dat))) {

  shared <- intersect(c("elevation", "slope"), names(popdat))

  if (length(shared)) {

    verify <- inner_join(
      dat[, c(key, shared), drop = FALSE],
      popdat[, c(key, shared), drop = FALSE],
      by = key,
      suffix = c(".current", ".old")
    )

    for (v in shared) {

      a <- verify[[paste0(v, ".current")]]
      z <- verify[[paste0(v, ".old")]]

      different <- is.finite(a) & is.finite(z) & abs(a - z) > 1e-6

      if (any(different)) {
        stop(
          "Cell ID / grid covariate mismatch for ", v,
          ". Verify spatial alignment before transferring population."
        )
      }
    }
  }

  pop_lookup <- popdat[, c(key, popvar)]

  dat <- left_join(dat, pop_lookup, by = key)
}

if (!is.numeric(dat[[popvar]]) ||
    any(dat[[popvar]] < 0, na.rm = TRUE)) {
  stop("Invalid population density values")
}

# Matching specifications

spec <- list(
  tree = c(
    "elevation", "slope", "dist_roads",
    "dist_settlements", "treecover_2000"
  ),
  crop = c(
    "elevation", "slope", "dist_roads",
    "dist_settlements", "cropland_baseline"
  ),
  built = c(
    "elevation", "slope", "dist_roads",
    "dist_settlements", "baseline_built"
  )
)

parks <- sort(unique(as.character(dat$NAME_ENG)))

# Matching function

run_one <- function(d0, covars, park, analysis, variant) {

  d <- d0[
    d0$treated | (!d0$treated & d0$pool_20_100),
    ,
    drop = FALSE
  ]

  d <- d[order(d$global_id, d$treated), , drop = FALSE]

  total_t <- sum(d$treated)

  complete <- complete.cases(d[, covars, drop = FALSE])

  n_t <- sum(d$treated & complete)
  n_c <- sum(!d$treated & complete)

  missing_pop <- sum(
    !is.finite(d[[popvar]]) | is.na(d[[popvar]])
  )

  d <- d[complete, , drop = FALSE]

  usage <- data.frame(
    NAME_ENG = park,
    analysis = analysis,
    variant = variant,
    treated_total = total_t,
    treated_complete = n_t,
    controls_complete = n_c,
    population_missing_in_eligible = missing_pop,
    treated_matched = NA_integer_,
    max_abs_SMD = NA_real_,
    status = "failed",
    detail = ""
  )

  if (!n_t || !n_c) {
    usage$detail <- "No complete treated or control observations"
    return(list(match = NULL, usage = usage, balance = data.frame()))
  }

  d$treated <- as.integer(d$treated)
  rownames(d) <- paste0("cell_", d$global_id)

  m <- tryCatch(
    MatchIt::matchit(
      reformulate(covars, response = "treated"),
      data = d,
      method = "nearest",
      distance = "glm",
      mahvars = covars,
      caliper = 0.2,
      std.caliper = TRUE,
      replace = FALSE,
      ratio = 1,
      m.order = "closest"
    ),
    error = function(e) e
  )

  if (inherits(m, "error")) {
    usage$detail <- conditionMessage(m)
    return(list(match = NULL, usage = usage, balance = data.frame()))
  }

  mm <- m$match.matrix

  usage$treated_matched <- if (is.null(mm)) {
    0L
  } else {
    sum(!is.na(mm[, 1]))
  }

  usage$status <- if (usage$treated_matched) {
    "matched"
  } else {
    "no_matches"
  }

  b <- tryCatch(
    cobalt::bal.tab(
      m,
      un = TRUE,
      stats = "mean.diffs",
      binary = "std",
      addl = if (popvar %in% covars) NULL else popvar
    ),
    error = function(e) e
  )

  if (inherits(b, "error")) {
    usage$detail <- paste("Balance error:", conditionMessage(b))
    return(list(match = m, usage = usage, balance = data.frame()))
  }

  bal <- b$Balance
  bal <- bal[rownames(bal) != "distance", , drop = FALSE]

  if (!all(c("Diff.Un", "Diff.Adj") %in% names(bal))) {
    stop("Unexpected cobalt balance columns")
  }

  balance <- data.frame(
    NAME_ENG = park,
    analysis = analysis,
    variant = variant,
    covariate = rownames(bal),
    SMD_before = as.numeric(bal$Diff.Un),
    SMD_after = as.numeric(bal$Diff.Adj),
    row.names = NULL
  )

  balance$abs_SMD_after <- abs(balance$SMD_after)

  orig <- balance$abs_SMD_after[
    balance$covariate %in% covars
  ]

  usage$max_abs_SMD <- if (
    length(orig) && any(is.finite(orig))
  ) {
    max(orig, na.rm = TRUE)
  } else {
    NA_real_
  }

  list(match = m, usage = usage, balance = balance)
}

# Run five- and six-covariate matching

usage <- list()
balances <- list()

for (analysis in names(spec)) {

  for (variant in c(
    "five_original",
    "six_population",
    "five_common_sample",
    "six_common_sample"
  )) {

    matches <- setNames(
      vector("list", length(parks)),
      parks
    )

    for (park in parks) {

      d0 <- dat[
        dat$NAME_ENG == park,
        ,
        drop = FALSE
      ]

      covars <- if (grepl("^six", variant)) {
        c(spec[[analysis]], popvar)
      } else {
        spec[[analysis]]
      }

      if (grepl("common_sample", variant)) {

        eligible <- d0$treated |
          (!d0$treated & d0$pool_20_100)

        keep <- complete.cases(
          d0[, c(spec[[analysis]], popvar), drop = FALSE]
        )

        d0 <- d0[!eligible | keep, , drop = FALSE]
      }

      r <- run_one(
        d0, covars, park, analysis, variant
      )

      matches[[park]] <- r$match

      usage[[length(usage) + 1L]] <- r$usage

      if (nrow(r$balance)) {
        balances[[length(balances) + 1L]] <- r$balance
      }
    }

    saveRDS(
      matches,
      file.path(
        out,
        paste0("matching_", analysis, "_", variant, ".rds")
      )
    )
  }
}

# Export matching diagnostics

u <- bind_rows(usage)
b <- bind_rows(balances)

write.csv(
  u,
  file.path(out, "matching_usage_all.csv"),
  row.names = FALSE
)

write.csv(
  b,
  file.path(out, "covariate_balance_all.csv"),
  row.names = FALSE
)

pop_balance <- b |>
  filter(covariate == popvar) |>
  select(
    NAME_ENG, analysis, variant,
    SMD_before, SMD_after, abs_SMD_after
  ) |>
  pivot_wider(
    names_from = variant,
    values_from = c(
      SMD_before, SMD_after, abs_SMD_after
    )
  )

write.csv(
  pop_balance,
  file.path(out, "population_balance_comparison.csv"),
  row.names = FALSE
)

matched <- u |>
  select(
    NAME_ENG, analysis, variant,
    treated_complete, controls_complete,
    treated_matched, status, max_abs_SMD
  ) |>
  pivot_wider(
    names_from = variant,
    values_from = c(
      treated_complete, controls_complete,
      treated_matched, status, max_abs_SMD
    )
  )

write.csv(
  matched,
  file.path(out, "matching_comparison.csv"),
  row.names = FALSE
)
