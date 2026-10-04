# ============================================================
# 03. RUN MATCHING
# Mozambique protected areas
# ============================================================

library(MatchIt)
library(cobalt)
library(dplyr)
library(sf)

matching_dir <- "data/processed/matching"
wdpa_dir <- "data/raw/wdpa"
output_dir <- "outputs/matching"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load prepared 1-km data
# ------------------------------------------------------------

dat_1000 <- readRDS(file.path(matching_dir, "matching_input_1000m.rds"))
grid_1000 <- readRDS(file.path(matching_dir, "matching_grid_1000m.rds"))

# ------------------------------------------------------------
# Exclude candidate controls that fall inside ANY terrestrial PA
# ------------------------------------------------------------

wdpa_files <- list.files(
  wdpa_dir,
  pattern = ".shp$",
  recursive = TRUE,
  full.names = TRUE
)


parks_all <- do.call(rbind, lapply(wdpa_files, st_read, quiet = TRUE))
parks_all_terrestrial <- parks_all %>%
  filter(REALM == "Terrestrial") %>%
  st_make_valid()

grid_centres <- st_centroid(grid_1000)
parks_all_terrestrial <- st_transform(parks_all_terrestrial, st_crs(grid_centres))

inside_any_pa <- lengths(st_intersects(grid_centres, parks_all_terrestrial)) > 0

grid_1000$inside_any_pa <- inside_any_pa

dat_1000$inside_any_pa <- inside_any_pa

dat_1000_clean <- dat_1000 %>%
  filter(treated | !inside_any_pa)

saveRDS(
  dat_1000_clean,
  file.path(matching_dir, "matching_input_1000m_unprotected_controls.rds")
)

# ------------------------------------------------------------
# Matching specification
# 1:1 nearest-neighbour matching with replacement.
# Mahalanobis distance on the matching covariates within a
# propensity-score caliper of 0.2 standard deviations.
# ------------------------------------------------------------

analysis_spec <- list(
  tree = list(
    label = "Tree-cover loss",
    covars = c("elevation", "slope", "dist_roads", "dist_settlements", "treecover_2000")
  ),
  crop = list(
    label = "Cropland expansion",
    covars = c("elevation", "slope", "dist_roads", "dist_settlements", "cropland_baseline")
  ),
  built = list(
    label = "Built-up expansion",
    covars = c("elevation", "slope", "dist_roads", "dist_settlements", "baseline_built")
  )
)

run_matching <- function(dat, covars) {

  form <- reformulate(covars, response = "treated")
  vars_needed <- all.vars(form)
  parks <- unique(dat$NAME_ENG)
  results <- setNames(vector("list", length(parks)), parks)

  for (park in parks) {

    message("Matching: ", park)

    d <- dat[dat$NAME_ENG == park, , drop = FALSE]
    d <- d[complete.cases(d[, vars_needed, drop = FALSE]), , drop = FALSE]
    d$treated <- as.integer(d$treated)

    if (length(unique(d$treated)) < 2) {
      message("  MATCHING FAILED: treated and control cells are not both available")
      next
    }

    results[[park]] <- tryCatch(
      MatchIt::matchit(
        formula = form,
        data = d,
        method = "nearest",
        distance = "glm",
        mahvars = covars,
        caliper = 0.2,
        std.caliper = TRUE,
        replace = TRUE,
        ratio = 1
      ),
      error = function(e) {
        message("  MATCHING FAILED: ", conditionMessage(e))
        NULL
      }
    )
  }

  results[!vapply(results, is.null, logical(1))]
}

get_balance <- function(matches, analysis_label) {
  bind_rows(lapply(names(matches), function(park) {
    b <- cobalt::bal.tab(matches[[park]], un = TRUE)
    x <- b$Balance[rownames(b$Balance) != "distance", , drop = FALSE]
    smd <- abs(x$Diff.Adj)

    data.frame(
      NAME_ENG = park,
      analysis = analysis_label,
      resolution_m = 1000,
      mean_abs_SMD = mean(smd, na.rm = TRUE),
      max_abs_SMD = max(smd, na.rm = TRUE),
      n_over_0.10 = sum(smd > 0.10, na.rm = TRUE)
    )
  }))
}

get_covariate_balance <- function(matches) {
  bind_rows(lapply(names(matches), function(park) {
    b <- cobalt::bal.tab(matches[[park]], un = TRUE)
    x <- b$Balance[
      rownames(b$Balance) != "distance",
      c("Diff.Un", "Diff.Adj"),
      drop = FALSE
    ]

    data.frame(
      NAME_ENG = park,
      covariate = rownames(x),
      SMD_before = x$Diff.Un,
      SMD_after = x$Diff.Adj,
      abs_SMD_after = abs(x$Diff.Adj),
      row.names = NULL
    )
  }))
}

get_matching_usage <- function(matches) {
  bind_rows(lapply(names(matches), function(park) {
    m <- matches[[park]]
    md <- MatchIt::match.data(m)
    treated <- md[md$treated == 1, , drop = FALSE]
    controls <- md[md$treated == 0, , drop = FALSE]

    data.frame(
      NAME_ENG = park,
      treated_original = sum(m$treat == 1),
      treated_matched = nrow(treated),
      unique_controls = nrow(controls),
      treated_retained_percent = 100 * nrow(treated) / sum(m$treat == 1),
      treated_per_unique_control = nrow(treated) / nrow(controls)
    )
  }))
}

# ------------------------------------------------------------
# Run and save each pressure separately.
# ------------------------------------------------------------

for (analysis in names(analysis_spec)) {

  spec <- analysis_spec[[analysis]]
  matches <- run_matching(dat_1000_clean, spec$covars)
  balance <- get_balance(matches, spec$label)
  covariate_balance <- get_covariate_balance(matches)
  usage <- get_matching_usage(matches)

  saveRDS(matches, file.path(output_dir, paste0("matching_", analysis, "_1000m_final.rds")))
  write.csv(balance, file.path(output_dir, paste0("matching_", analysis, "_balance_final.csv")), row.names = FALSE)
  write.csv(covariate_balance, file.path(output_dir, paste0("matching_", analysis, "_covariate_balance_final.csv")), row.names = FALSE)
  write.csv(usage, file.path(output_dir, paste0("matching_", analysis, "_control_usage_final.csv")), row.names = FALSE)
}
