# 01 FIGURE 1: SPATIAL PATTERNS OF LAND-USE CHANGE ACROSS TIME
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggspatial)
library(png)
library(grid)
library(cowplot)

figure_dir <- "outputs/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Protected areas and Mozambique boundary

load("data/processed/parks_core.Rdata")

parks <- parks_core |>
  filter(DESIG %in% c(
    "Parque Nacional", "National Park",
    "Reserva Nacional", "Reserva Especial"
  )) |>
  st_make_valid() |>
  st_transform(4326)

moz <- st_read(
  "/Users/gdt366/Library/CloudStorage/Dropbox/gee_deforestation_mozambique/moz_admin_boundaries.shp/moz_admin0.shp",
  quiet = TRUE
) |>
  st_make_valid() |>
  st_transform(4326)

bb <- st_bbox(moz)
x_limits <- c(bb["xmin"] - 0.05, bb["xmax"] + 0.05)
y_limits <- c(bb["ymin"] - 0.05, bb["ymax"] + 0.05)

# theme

theme_harith <- theme(
  strip.text = element_text(size = 16, face = "bold"),
  strip.background = element_blank(),
  axis.title.y = element_text(
    size = 13, face = "bold", angle = 90,
    margin = margin(r = 8)
  ),
  legend.position = c(0.95, 0.25),
  legend.direction = "vertical",
  legend.title.position = "left",
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.title = element_text(size = 8.5, angle = 90, hjust = 0.5),
  legend.text = element_text(size = 7.5),
  panel.spacing = unit(0.25, "lines"),
  plot.margin = margin(0, 0, 0, 0),
  plot.tag.position = c(0.03, 0.95),
  plot.tag = element_text(size = 12, face = "bold")
)

# REMAINING BASELINE FOREST COVER (2000, 2005, 2015, 2025)


TREE_COVER_THRESHOLD <- 30
display_fact <- 5
tree_years <- c(2000, 2005, 2015, 2025)

tree <- rast("data/raw/hansen/Mozambique_treecover2000_percent.tif")
loss <- rast("data/raw/hansen/Mozambique_forest_loss_year_2001_2025.tif")


forest_2000 <- ifel(
  is.na(tree), NA,
  ifel(tree >= TREE_COVER_THRESHOLD, 1, 0)
)
names(forest_2000) <- "forest_2000"

baseline_forest_pct <- aggregate(
  forest_2000, fact = display_fact, fun = "mean", na.rm = TRUE
) * 100

make_tree_state <- function(year) {
  if (year == 2000) {
    remaining_forest <- forest_2000
  } else {
    lost_by_year <- forest_2000 == 1 & loss >= 2001 & loss <= year
    remaining_forest <- ifel(lost_by_year, 0, forest_2000)
    remaining_forest <- ifel(
      is.na(forest_2000) | is.na(loss), NA, remaining_forest
    )
  }

  state <- aggregate(
    remaining_forest, fact = display_fact, fun = "mean", na.rm = TRUE
  ) * 100

  state <- ifel(baseline_forest_pct > 0, state, NA)
  names(state) <- "state"
  state
}

tree_states <- lapply(tree_years, make_tree_state)
names(tree_states) <- as.character(tree_years)

tree_df <- bind_rows(lapply(seq_along(tree_states), function(i) {
  r <- tree_states[[i]]
  if (!is.lonlat(r)) {
    r <- project(r, "EPSG:4326", method = "near")
  }
  x <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(x)[3] <- "state"
  x$year <- factor(tree_years[i], levels = tree_years)
  x
}))

pal_tree <- colorRampPalette(c(
  "#f1f5e9", "#c8dfb4", "#85bb7b", "#3f8d58"
))

p_tree <- ggplot() +
  geom_sf(data = moz, fill = "white", colour = NA) +
  geom_raster(data = tree_df, aes(x, y, fill = state)) +
  geom_sf(data = parks, fill = NA, colour = "black", linewidth = 0.12) +
  geom_sf(data = moz, fill = NA, colour = 'black',linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  annotation_north_arrow(
    data = data.frame(year = factor(2000, levels = tree_years)),
    location = "br",
    which_north = "true",
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm"),
    pad_x = unit(0.45, "cm"),
    pad_y = unit(0.45, "cm"),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    data = data.frame(year = factor(2000, levels = tree_years)),
    location = "br",
    width_hint = 0.35,
    text_cex = 0.55,
    line_width = 0.4,
    pad_x = unit(0.15, "cm"),
    pad_y = unit(0.15, "cm")
  ) +
  scale_fill_gradientn(
    colours = pal_tree(100),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    name = "Remaining 2000 forest (%)"
  ) +
  coord_sf(
    crs = st_crs(4326), xlim = x_limits, ylim = y_limits,
    expand = FALSE, datum = NA
  ) +
  labs(x = NULL, y = "Forest cover") +
  theme_void() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold", angle = 90),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top",
    legend.key.height = grid::unit(0.42, "cm"),
    legend.key.width = grid::unit(0.35, "cm"),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7),
    panel.spacing = grid::unit(0.12, "lines")
  ) +
  theme_harith

# BUILT-UP SURFACE

built_years <- c(2000, 2005, 2015, 2025)

built_df <- bind_rows(lapply(built_years, function(y) {
  r <- rast(file.path(
    "data/raw/ghsl",
    paste0("Mozambique_GHSL_built_surface_", y, ".tif")
  ))[["built_surface"]]

  r <- aggregate(r, fact = 5, fun = "sum", na.rm = TRUE)
  if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "bilinear")
  r <- ifel(r > 0, r, NA)
  names(r) <- "built"

  x <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(x)[3] <- "built"
  x$year <- factor(y, levels = built_years)
  x
}))

built_limit <- quantile(built_df$built, 0.99, na.rm = TRUE)
pal_built_up <- colorRampPalette(hcl.colors(5, "Sunset", rev = TRUE))
pal_built_up <- colorRampPalette(hcl.colors(5, "purples", rev = TRUE))

p_built <- ggplot() +
  geom_sf(data = moz, fill = "white", colour = NA) +
  geom_raster(data = built_df, aes(x, y, fill = built)) +
  geom_sf(data = parks, fill = NA, colour = "black", linewidth = 0.12) +
  geom_sf(data = moz, fill = NA, colour = 'black',linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  scale_fill_gradientn(
    colours = pal_built_up(10),
    trans = scales::pseudo_log_trans(sigma = 500),
    breaks = c(500, 2500, 10000, 40000, 100000),
    labels = c("0.5k", "2.5k", "10k", "40k", "100k"),
    na.value = NA,
    name = "Built surface (m²)"
  ) +
  coord_sf(
    crs = st_crs(4326), xlim = x_limits, ylim = y_limits,
    expand = FALSE, datum = NA
  ) +
  labs(x = NULL, y = "Built-up area") +
  theme_void() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold", angle = 90),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top",
    legend.key.height = grid::unit(0.75, "cm"),
    legend.key.width = grid::unit(0.35, "cm"),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7),
    panel.spacing = grid::unit(0.12, "lines")
  ) +
  theme_harith

# CROPLAND: GACED30 (2000-2024; UNCHANGED CALCULATION)

crop_years <- c(2000, 2005, 2015, 2024)

gaced <- rast(
  "data/raw/GACED30/Mozambique_GACED30_cropland_2000_2024.tif"
)

read_crop <- function(y) {
  gaced[[paste0("crop_", y)]]
}

# GACED30: 10 = cropland, 0 = non-cropland.
crop_binary <- lapply(crop_years, function(y) {
  r <- read_crop(y)
  ifel(r == 10, 1, ifel(r == 0, 0, NA))
})

crop_pct <- lapply(crop_binary, function(r) {
  aggregate(r, fact = 5, fun = "mean", na.rm = TRUE) * 100
})

crop_baseline <- crop_pct[[1]]
crop_presence_2000 <- ifel(crop_baseline > 0, 1, NA)
names(crop_presence_2000) <- "state"

make_crop_change <- function(current) {
  state <- current - crop_baseline
  agriculture_present <- crop_baseline > 0 | current > 0
  state <- ifel(agriculture_present, state, NA)
  names(state) <- "state"
  state
}

# The original script referenced crop_changes but did not create it.
crop_changes <- lapply(crop_pct[-1], make_crop_change)

to_df <- function(r, year) {
  if (!is.lonlat(r)) {
    r <- project(r, "EPSG:4326", method = "bilinear")
  }
  x <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(x)[3] <- "crop"
  x$year <- factor(year, levels = crop_years)
  x
}

# Baseline cropland is shown at the midpoint of the change scale.
crop_presence_2000 <- ifel(!is.na(crop_presence_2000), 0, NA)

df_crop <- bind_rows(
  to_df(crop_presence_2000, 2000),
  to_df(crop_changes[[1]], 2005),
  to_df(crop_changes[[2]], 2015),
  to_df(crop_changes[[3]], 2024)
)

pal_crop <- colorRampPalette(rev(hcl.colors(5, "RdYlBu")))
pal_crop <- colorRampPalette(hcl.colors(5, "Zissou1"))
pal_crop <- colorRampPalette(rev(hcl.colors(5, "Oranges")))

p_crop <- ggplot() +
  geom_sf(data = moz, fill = "white", colour = NA) +
  geom_raster(data = df_crop, aes(x = x, y = y, fill = crop)) +
  geom_sf(data = parks, fill = NA, colour = "black", linewidth = 0.12) +
  geom_sf(data = moz, fill = NA, colour = 'black',linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  scale_fill_gradientn(
    colours = pal_crop(100),
    limits = c(-100, 100),
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100),
    labels = c(
      "-100%", "-75%", "-50%", "-25%", "0",
      "+25%", "+50%", "+75%", "+100%"
    ),
    name = "Cropland change"
  ) +
  coord_sf(
    crs = st_crs(4326), xlim = x_limits, ylim = y_limits,
    expand = FALSE, datum = NA
  ) +
  labs(x = NULL, y = "Cropland") +
  theme_void() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold", angle = 90),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top",
    legend.key.height = grid::unit(0.5, "cm"),
    legend.key.width = grid::unit(0.25, "cm"),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.5),
    panel.spacing = grid::unit(0.12, "lines")
  ) +
  theme_harith


figure_1 <- p_tree / p_crop / p_built +
  plot_layout(heights = c(1, 1, 1))

library(grid)
library(png)
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

final_plot <- ggdraw() +
  draw_plot(figure_1, 0, 0, 1, 1) +
  draw_plot(miombo_icon, x = 0.01, y = 0.9, width = 0.075, height = 0.075) +
  draw_plot(crop_icon, x = 0.01, y = 0.57, width = 0.075, height = 0.075) +
  draw_plot(built_up_icon, x = 0.01, y = 0.25, width = 0.095, height = 0.095)

ggsave(
  file.path(figure_dir, "Figure1_threshold_30_remaining_forest.png"),
  final_plot,
  width = 9,
  height = 10.5,
  dpi = 600,
  bg = "white"
)

#saving files for future runs
save(tree_df,file = 'tree_df.Rdata')
save(built_df,file = 'built_df.Rdata')
save(crop_df,file = 'crop_df.Rdata')

