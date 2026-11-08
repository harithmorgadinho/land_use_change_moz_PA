
# New Figure 1 script

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)
library(png)
library(cowplot)

output_dir <- "outputs/figures/dominant_pressure"
cache_dir <- file.path(output_dir, "cache")

pressures <- rast(file.path(output_dir, "pressure_percentages.tif"))
dominant <- rast(file.path(output_dir, "dominant_pressure.tif"))

map_df <- readRDS(file.path(cache_dir, "dominant_pressure_map_df.rds"))
moz <- readRDS(file.path(cache_dir, "moz_boundary.rds"))
parks <- readRDS(file.path(cache_dir, "parks.rds"))

bb <- st_bbox(moz)
x_limits <- c(bb["xmin"] - 0.05, bb["xmax"] + 0.05)
y_limits <- c(bb["ymin"] - 0.05, bb["ymax"] + 0.05)

# Satellite background

satellite_fixed <- satellite_aligned

black_pixels <- (
  satellite_fixed[[1]] == 0 &
    satellite_fixed[[2]] == 0 &
    satellite_fixed[[3]] == 0
)

satellite_fixed <- ifel(
  black_pixels,
  NA,
  satellite_fixed
)

satellite_fixed <- app(
  satellite_fixed,
  fun = function(x) {
    pmin(255, 255 * (x / 255)^0.75)
  }
)

names(satellite_fixed) <- c("red", "green", "blue")

satellite_vivid <- satellite_fixed / 255
satellite_vivid <- satellite_vivid ^ 0.85

satellite_vivid <- clamp(
  (satellite_vivid - 0.5) * 1.12 + 0.5,
  lower = 0,
  upper = 1
)

rgb_mean <- mean(satellite_vivid)

satellite_vivid <- clamp(
  rgb_mean + (satellite_vivid - rgb_mean) * 1.30,
  lower = 0,
  upper = 1
)

satellite_vivid[[2]] <- clamp(
  satellite_vivid[[2]] * 1.08,
  lower = 0,
  upper = 1
)

satellite_vivid <- round(satellite_vivid * 255)
names(satellite_vivid) <- c("red", "green", "blue")

# Protected area labels

park_labels <- parks |>
  st_transform(32737) |>
  st_point_on_surface() |>
  st_transform(4326)

coords <- st_coordinates(park_labels)

park_labels <- park_labels |>
  mutate(
    lon = coords[, 1],
    lat = coords[, 2]
  ) |>
  arrange(desc(lat)) |>
  mutate(park_id = row_number())

park_legend <- park_labels |>
  st_drop_geometry() |>
  transmute(
    park_id,
    park_name = NAME
  ) |>
  arrange(park_id)

park_legend_text <- paste(
  paste0(park_legend$park_id, ". ", park_legend$park_name),
  collapse = "\n"
)

levels(map_df$pressure)[
  levels(map_df$pressure) == "Deforestation"
] <- "Tree-cover loss"

# Dominant pressure map

p <- ggplot() +
  geom_sf(
    data = moz,
    fill = "grey95",
    colour = NA
  ) +
  geom_spatraster_rgb(
    data = satellite_vivid,
    r = 1,
    g = 2,
    b = 3,
    maxcell = 100000000
  ) +
  geom_raster(
    data = map_df,
    aes(x = x, y = y, fill = pressure),
    alpha = 0.65
  ) +
  geom_sf(
    data = parks,
    fill = NA,
    colour = "black",
    linewidth = 0.15
  ) +
  geom_sf(
    data = moz,
    fill = NA,
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Tree-cover loss" = "red2",
      "Cropland" = "#E68A00",
      "Built-up surface" = "#7A3DB8"
    ),
    name = "Dominant pressure",
    drop = FALSE
  ) +
  geom_text_repel(
    data = park_labels,
    aes(
      x = lon,
      y = lat,
      label = park_id
    ),
    size = 5,
    fontface = "bold",
    colour = "black",
    label.size = 0.15,
    label.padding = unit(0.15, "lines"),
    box.padding = 0.15,
    point.padding = 0,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  coord_sf(
    crs = st_crs(4326),
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE,
    datum = NA
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = c(0.7, 0.4),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(0, 0, 0, 0),
    plot.tag.position = c(0, 0.999),
    plot.tag = element_text(
      face = "bold",
      hjust = 0,
      size = 13
    )
  ) +
  labs(tag = "a. Dominant pressure") +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm"),
    pad_x = unit(0.45, "cm"),
    pad_y = unit(0.9, "cm"),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.35,
    text_cex = 0.55,
    line_width = 0.4,
    pad_x = unit(0.15, "cm"),
    pad_y = unit(0.15, "cm")
  ) +
  annotate(
    "text",
    x = x_limits[2] - 4.6,
    y = y_limits[1] + 0.45,
    label = park_legend_text,
    hjust = 0,
    vjust = 0,
    size = 3.5,
    lineheight = 1.1,
    colour = "black"
  )

# National land-use trends

root <- getwd()
national_dir <- file.path(root, "outputs/national")
out_dir <- file.path(national_dir, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

forest <- read.csv(
  file.path(national_dir, "Mozambique_treecover_annual_loss_2001_2025.csv")
)

crop <- read.csv(
  file.path(national_dir, "Mozambique_cropland_annual_2000_2024.csv")
)

built <- read.csv(
  file.path(national_dir, "Mozambique_builtup_temporal_2000_2025.csv")
)

national_summary <- read.csv(
  file.path(
    national_dir,
    "Table_S1_national_land_use_change_summary_fixed.csv"
  ),
  check.names = FALSE
)

baseline_forest_ha <- national_summary[
  national_summary$Threat == "Forest",
  "Baseline extent (ha)"
]

forest_df <- forest |>
  transmute(
    year = as.integer(year),
    change_ha = cumulative_loss_ha,
    relative_change = 100 * cumulative_loss_ha / baseline_forest_ha
  ) |>
  bind_rows(
    data.frame(
      year = 2000,
      change_ha = 0,
      relative_change = 0
    )
  ) |>
  arrange(year)

crop_baseline <- crop$cropland_ha[crop$year == 2000]

crop_df <- crop |>
  transmute(
    year = as.integer(year),
    change_ha = cropland_ha - crop_baseline,
    relative_change = 100 *
      (cropland_ha - crop_baseline) / crop_baseline
  ) |>
  arrange(year)

built_baseline <- built$built_surface_ha[built$year == 2000]

built_df <- built |>
  transmute(
    year = as.integer(year),
    change_ha = built_surface_ha - built_baseline,
    relative_change = 100 *
      (built_surface_ha - built_baseline) / built_baseline
  ) |>
  arrange(year)

# Endpoint labels

format_ha <- function(x) {
  paste0(
    ifelse(x >= 0, "+", "-"),
    format(
      round(abs(x)),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    " ha"
  )
}

format_pct <- function(x) {
  sprintf("%+.1f%%", x)
}

make_panel <- function(df, title, colour) {

  final <- df[nrow(df), , drop = FALSE]

  final$label <- paste0(
    format_pct(final$relative_change),
    "\n",
    format_ha(final$change_ha)
  )

  ggplot(df, aes(x = year, y = relative_change)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey65",
      linewidth = 0.4
    ) +
    geom_area(
      fill = colour,
      alpha = 0.25
    ) +
    geom_line(
      colour = colour,
      linewidth = 0.9
    ) +
    geom_point(
      data = final,
      colour = colour,
      size = 2.8
    ) +
    geom_text(
      data = final,
      aes(label = label),
      hjust = 1.05,
      vjust = -0.4,
      size = 3.5,
      fontface = "bold",
      colour = colour,
      lineheight = 1.2
    ) +
    scale_x_continuous(
      breaks = seq(2000, 2025, 5),
      limits = c(2000, 2025)
    ) +
    scale_y_continuous(
      labels = label_number(
        accuracy = 0.1,
        suffix = "%"
      ),
      expand = expansion(mult = c(0.06, 0.30))
    ) +
    labs(
      title = title,
      x = "Year",
      y = "Change relative to 2000 (%)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = 13
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

p_forest <- make_panel(
  forest_df,
  "b. Tree-cover loss",
  "#238B45"
)

p_crop <- make_panel(
  crop_df,
  "c. Cropland",
  "#E68A00"
)

p_built <- make_panel(
  built_df,
  "d. Built-up surface",
  "#7A3DB8"
)

# Export national trends data

all_data <- bind_rows(
  forest_df |> mutate(pressure = "Tree-cover loss"),
  crop_df |> mutate(pressure = "Cropland"),
  built_df |> mutate(pressure = "Built-up surface")
)

write.csv(
  all_data,
  file.path(out_dir, "Mozambique_national_relative_land_use_change.csv"),
  row.names = FALSE
)

# Combine map and national trends

figure_national <- (
  p_forest /
    p_crop /
    p_built
) +
  plot_layout(
    ncol = 1,
    heights = c(1, 1, 1)
  )

miombo_icon <- rasterGrob(
  readPNG("miombo_forest_silhouette_transparent.png", native = TRUE),
  interpolate = TRUE
)

crop_icon <- rasterGrob(
  readPNG("Three-Stalk Maize Silhouette.png", native = TRUE),
  interpolate = TRUE
)

built_up_icon <- rasterGrob(
  readPNG("Urban Settlement Silhouette Icon.png", native = TRUE),
  interpolate = TRUE
)

figure_combined <- (p | figure_national) +
  plot_layout(widths = c(2, 1))

final_plot <- ggdraw() +
  draw_plot(figure_combined, 0, 0, 1, 1) +
  draw_plot(
    miombo_icon,
    x = 0.7, y = 0.87,
    width = 0.075, height = 0.075
  ) +
  draw_plot(
    crop_icon,
    x = 0.7, y = 0.54,
    width = 0.075, height = 0.075
  ) +
  draw_plot(
    built_up_icon,
    x = 0.69, y = 0.2,
    width = 0.095, height = 0.095
  )

ggsave(
  "outputs/figures/Figure_dominant_pressure_national_trends.png",
  final_plot,
  width = 10,
  height = 9,
  dpi = 1000,
  bg = "white"
)

ggsave(
  "outputs/figures/Figure_dominant_pressure_national_trends_low.png",
  final_plot,
  width = 10,
  height = 9,
  dpi = 600,
  bg = "white"
)
