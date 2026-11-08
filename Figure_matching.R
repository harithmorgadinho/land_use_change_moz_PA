# FIGURE BALANCED MATCHED COMPARISONS ----

results_file <- paste0(
  "outputs/results/",
  "revised_no_replacement_forest_proportion_cell_area_exclusion_exception/",
  "PA_matching_effects_all_pressures.csv"
)

#from script 04
effects_all <- read.csv(
  results_file,
  stringsAsFactors = FALSE
)
effects_all$pressure <- dplyr::recode(
  effects_all$analysis,
  tree  = "a. Baseline forest loss",
  crop  = "b. Cropland expansion",
  built = "c. Built-up expansion"
)

all_parks <- c(
  "Banhine", "Chimanimani", "Gilé", "Gorongosa", "Limpopo",
  "Magoe", "Maputo", "Marromeu", "Niassa", "Pomene",
  "Quirimbas", "Serra da Gorongosa", "Zinave"
)

effects_plot <- effects_all |>
  filter(
    !is.na(effect_pp),
    !is.na(ci_lower_pp), !is.na(ci_upper_pp),
    !is.na(balance_025),
    balance_025
  ) |>
  mutate(
    NAME_ENG = factor(NAME_ENG, levels = rev(all_parks)),
    pressure = factor(
      pressure,
      levels = c(
        "a. Tree-cover loss",
        "b. Cropland expansion",
        "c. Built-up expansion"
      )
    ),
    effect_label = sprintf("%+.2f", effect_pp),
    ci_lower_plot = ci_lower_pp,
    ci_upper_plot = ci_upper_pp
  )

ci_candidates <- list(
  c("ci_low_pp", "ci_high_pp"),
  c("ci_lower_pp", "ci_upper_pp"),
  c("lower_ci", "upper_ci"),
  c("CI_lower", "CI_upper")
)

ci_pair <- NULL

for (pair in ci_candidates) {
  if (all(pair %in% names(effects_plot))) {
    ci_pair <- pair
    break
  }
}

effects_plot <- effects_plot |>
  group_by(pressure) |>
  mutate(
    label_offset = 0.03 * diff(
      range(
        c(ci_lower_plot, ci_upper_plot),
        na.rm = TRUE
      )
    )
  ) |>
  ungroup()

label_offset <- 0.03 * diff(
  range(
    c(effects_plot$ci_lower_plot,
      effects_plot$ci_upper_plot),
    na.rm = TRUE
  )
)

shade_rows <- data.frame(
  ymin = seq(0.5, length(park_order) - 0.5, by = 2),
  ymax = seq(1.5, length(park_order) + 0.5, by = 2)
)

library(ggh4x)


#Figure 3 bottom ----
figure_3 <- ggplot(
  effects_plot,
  aes(x = effect_pp, y = NAME_ENG, fill = pressure)
) +
  geom_rect(
    data = shade_rows,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = "grey95", colour = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey40"
  ) +
  geom_col(
    width = 0.65,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower_plot,
      xmax = ci_upper_plot
    ),
    orientation = "y",
    width = 0.22,
    linewidth = 0.55,
    colour = "black"
  ) +
  geom_text(
    aes(
      x = ifelse(
        effect_pp >= 0,
        pmax(effect_pp, ci_upper_plot) + label_offset,
        pmin(effect_pp, ci_lower_plot) - label_offset
      ),
      label = effect_label,
      hjust = ifelse(effect_pp >= 0, 0, 1)
    ),
    size = 2.7,
    colour = "black"
  ) +
  facet_wrap(
    ~ pressure,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(c(
      "a. Tree-cover loss" = "b.",
      "b. Cropland expansion" = "d.",
      "c. Built-up expansion" = "f."
    ))
  ) +
  scale_fill_manual(
    values = c(
      # "a. Tree-cover loss" = "forestgreen",
      # "b. Cropland expansion" = "red3",
      # "c. Built-up expansion" = "grey40"
      "a. Tree-cover loss" = "#238B45",
      "b. Cropland expansion" = "#E68A00",
      "c. Built-up expansion" = "#7A3DB8")
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Protection effect (percentage points)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    #panel.background = element_rect(colour = 'grey90', fill = NA),
    strip.text = element_text(
      size = 17,
      face = "bold",hjust = 0
    ),
    axis.text.y = element_text(size = 10),
    #axis.text.x = element_text(size = 8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.spacing = grid::unit(0.1, "cm"),
    plot.margin = margin(0, 18, 0, 10)
  )

#Figure 3 top

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MatchIt)

matching_data_dir <- "data/processed/matching_GACED30/revised_20_50_100_forest_proportion_exclusion_exception"

matching_output_dir <- "outputs/matching/revised_no_replacement_forest_proportion_exclusion_exception"

results_file <- paste0(
  "outputs/results/",
  "revised_no_replacement_forest_proportion_cell_area_exclusion_exception/",
  "PA_matching_effects_all_pressures.csv"
)

wdpa_dir <- "data/raw/wdpa"

moz_file <- paste0(
  "data/raw/moz_admin_boundaries.shp/",
  "moz_admin0.shp"
)

figure_dir <- "outputs/figures/forest_proportion_wdpa_exception"

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

matching_data_dir <- "data/processed/matching_GACED30/revised_20_50_100_forest_proportion_exclusion_exception"
matching_output_dir <- "outputs/matching/revised_no_replacement_forest_proportion_exclusion_exception"

grid_1000 <- readRDS(file.path(
  matching_data_dir,
  "grid_1km_mollweide.rds"
))

m_tree <- readRDS(file.path(
  matching_output_dir,
  "matching_tree_20_100.rds"
))

m_crop <- readRDS(file.path(
  matching_output_dir,
  "matching_crop_20_100.rds"
))

m_built <- readRDS(file.path(
  matching_output_dir,
  "matching_built_20_100.rds"
))

effects_final <- read.csv(
  results_file,
  stringsAsFactors = FALSE
)
effects_final$pressure <- dplyr::recode(
  effects_final$analysis,
  tree  = "a. Tree-cover loss",
  crop  = "b. Cropland expansion",
  built = "c. Built-up expansion"
)

pressure_levels <- c(
  "a. Tree-cover loss",
  "b. Cropland expansion",
  "c. Built-up expansion"
)

pressure_colours <- c(
  "a. Tree-cover loss" = "#238B45",
  "b. Cropland expansion" = "#E68A00",
  "c. Built-up expansion" = "#7A3DB8"
)


extract_controls <- function(matches, pressure_label) {

  bind_rows(lapply(names(matches), function(park) {

    obj <- matches[[park]]

    if (is.null(obj)) return(NULL)

    md <- MatchIt::match.data(obj)

    md$global_id <- as.integer(
      sub("^cell_", "", rownames(md))
    )

    md |>
      filter(
        treated == 0,
        weights > 0
      ) |>
      transmute(
        global_id,
        NAME_ENG = park,
        pressure = pressure_label
      )
  }))
}

controls <- bind_rows(
  extract_controls(m_tree, "a. Tree-cover loss"),
  extract_controls(m_crop, "b. Cropland expansion"),
  extract_controls(m_built, "c. Built-up expansion")
)


#Retain controls from adequately balanced comparisons

balanced_parks <- effects_final |>
  filter(
    !is.na(effect_pp),
    !is.na(balance_025),
    balance_025
  ) |>
  select(
    NAME_ENG,
    pressure
  )

controls <- controls |>
  inner_join(
    balanced_parks,
    by = c("NAME_ENG", "pressure")
  )


#Coordinates of matched cells

grid_points <- st_centroid(grid_1000)

controls_sf <- grid_points |>
  inner_join(
    controls,
    by = "global_id"
  )


controls_sf <- st_transform(
  controls_sf,
  4326
)

controls_sf$pressure <- factor(
  controls_sf$pressure,
  levels = pressure_levels
)


moz <- st_read(
  moz_file,
  quiet = TRUE
)

moz_map <- st_transform(
  moz,
  4326
)


# protected-area analysis.


parks_map <- st_transform(
  parks_core,
  4326
)

# All other protected areas

wdpa_files <- list.files(
  wdpa_dir,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)


all_pas <- lapply(wdpa_files, function(f) {
  st_read(f, quiet = TRUE) |>
    st_transform(4326)
}) |>
  bind_rows()

# Remove study PA footprints from the background
# layer, preserving the remaining PA geometries.

all_pas[all_pas$NAME_ENG %in% 'Quirimbas',]
unique(all_pas$DESIG)

all_pas = all_pas[!all_pas$DESIG %in% "Zona Tampao" ,]

other_pas <- st_difference(
  st_union(st_geometry(all_pas)),
  st_union(st_geometry(parks_map))
)

other_pas <- st_sf(
  geometry = other_pas
)


# Create a copy of each study PA for each pressure.

park_effects <- bind_rows(
  lapply(pressure_levels, function(p) {

    parks_map |>
      mutate(pressure = p)

  })
)

park_effects <- park_effects |>
  left_join(
    effects_final |>
      select(
        NAME_ENG,
        pressure,
        effect_pp,
        balance_025
      ),
    by = c("NAME_ENG", "pressure")
  ) |>
  mutate(
    effect_map = if_else(
      !is.na(balance_025) & balance_025,
      effect_pp,
      NA_real_
    ),
    pressure = factor(
      pressure,
      levels = pressure_levels
    )
  )

effect_limit <- max(
  abs(park_effects$effect_map),
  na.rm = TRUE
)

# Plot

#theme_harith----
theme_harith <- theme(
  strip.text = element_text(
    size = 16,
    face = "bold"
  ),
  strip.background = element_blank(),
  axis.title.y = element_text(
    size = 13,
    face = "bold",
    margin = margin(r = 8)
  ),
  legend.position = c(-0.02, 0.12),
  legend.direction = "horizontal",
  legend.title.position = "top",
  legend.key.height = unit(
    0.75,
    "cm"),
  legend.key.width = unit(
    0.05,
    "cm"),
  legend.key.spacing = unit(
    0.15,
    "cm"),
  legend.title = element_text(
    size = 9,
    hjust = 0.5),
  legend.text = element_text(
    size = 10),
  panel.spacing = unit(
    0.25,
    "lines"),
  plot.margin = margin(0,0,0,0),
  plot.tag.position = c(0.03,0.95),
  plot.tag = element_text(size = 12,face = 'bold')
)
pal = colorRampPalette(hcl.colors(n = 5,palette = 'RdYlBu',rev = T))

library(png)
library(grid)
library(cowplot)

img <- readPNG('miombo_forest_silhouette_transparent.png', native = TRUE)
miombo_icon <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Three-Stalk Maize Silhouette.png', native = TRUE)
crop_icon <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Urban Settlement Silhouette Icon.png', native = TRUE)
built_up_icon <- rasterGrob(img, interpolate=TRUE)

library(ggspatial)
library(ggspatial)

map_annotations <- data.frame(
  pressure = "a. Tree-cover loss"
)

#plot 3a ----
counterfactual_map <- ggplot() +
  geom_sf(
    data = moz_map,
    fill = "white",
    colour = NA,
    linewidth = 0
  ) +
  geom_sf(
    data = other_pas,
    fill = "grey30",
    colour = NA,
    linewidth = 0
  ) +
  geom_sf(
    data = controls_sf,
    aes(colour = pressure),
    size = 0.12,
    alpha = 0.55
  ) +
  geom_sf(
    data = park_effects,
    aes(fill = effect_map),
    colour = NA,
    linewidth = 0
  ) +
  geom_sf(
    data = moz_map,
    fill = NA,
    colour = 'black',
    linewidth = 0.025
  ) +

  facet_wrap(
    ~pressure,
    nrow = 1,
    labeller = as_labeller(c(
      "a. Tree-cover loss" = "a. Tree-cover loss",
      "b. Cropland expansion" = "c. Cropland change",
      "c. Built-up expansion" = "e. Built-up expansion"
    ))
  ) +

  scale_fill_gradientn(
    colours = pal(10),
    limits = c(-effect_limit, effect_limit),
    breaks = c(-effect_limit, 0, effect_limit),
    labels = c(
      "Reduced\npressure",
      "0",
      "Increased\npressure"
    ),
    na.value = "grey65",
    name = "Protection effect"
  ) +

  scale_colour_manual(
    values = pressure_colours,
    name = "Matched control cells",
    labels = c('Tree cover loss','Cropland change', 'Built-up expansion')
  ) +

  coord_sf(
    datum = NA
  ) +

  annotation_scale(
    data = map_annotations,
    location = "br",
    width_hint = 0.25,
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm"),
    text_cex = 0.6,
    line_width = 0.4
  ) +

  annotation_north_arrow(
    data = map_annotations,
    location = "tl",
    which_north = "true",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.3, "cm"),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(
      text_size = 6,
      line_width = 0.5
    )
  ) +

  theme_void() +

  theme(
    strip.text = element_text(
      size = 8,
      face = "bold",
      hjust = 0,
      margin = margin(t = 0, b = 2, l = 0, r = 0)),
    plot.title = element_text(
      size = 8,
      face = "bold"
    ),
    panel.spacing = grid::unit(
      0.8,
      "cm"
    ),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  theme_harith+
  guides(
    fill = guide_colourbar(
      direction = "vertical",
      barheight = unit(3, "cm"),
      barwidth = unit(0.4, "cm"),
      title.position = "left",
      reverse = TRUE,
      title.theme = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      )
    ),
    colour = guide_legend(
      direction = "vertical",
      title.position = "top",
      title.theme = element_text(
        angle = 0,
        hjust = 0.5,
        size = 10
      ),
      nrow = 3,
      theme = theme(
        legend.key.height = unit(0.35, "cm"),
        legend.key.spacing.y = unit(0, "cm"),
        legend.text = element_text(
          margin = margin(0, 0, 0, 0)
        )
      ),
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  )


library(cowplot)

legend_gradient <- get_legend(
  counterfactual_map +
    guides(colour = "none") +
    theme(legend.position = "right")
)


map_no_legend <- counterfactual_map +
  guides(fill = "none")


library(cowplot)

figure_3_final <- ggdraw() +
  draw_plot(
    map_no_legend,
    x = 0.1,
    y = 0.30,
    width = 0.9,
    height = 0.70
  ) +
  draw_plot(
    figure_3,
    x = 0,
    y = 0.05,
    width = 1,
    height = 0.3
  )







final_plot <- ggdraw() +
  draw_plot(figure_3_final, 0, 0, 1, 1) + # Main plot
  draw_plot(miombo_icon, x = 0.12, y = 0.8, width = 0.075, height = 0.075)+
  draw_plot(crop_icon, x = 0.425, y = 0.8, width = 0.075, height = 0.075)+
  draw_plot(built_up_icon, x = 0.72, y = 0.79, width = 0.095, height = 0.095)+
  draw_grob(
    legend_gradient,
    x = 0.01,
    y = 0.50,
    width = 0.12,
    height = 0.25
  )


ggsave(
  file.path(
    figure_dir,
    "Figure3_final_nov2v2.png"
  ),
  plot = final_plot,
  width = 12,
  height = 10,
  dpi = 600
)
