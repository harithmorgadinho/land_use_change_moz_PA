# Figure 2: Temporal change within and around protected areas

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

figure_dir <- "outputs/figures/forest_proportion_update"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)



tree_ts <- read.csv("outputs/temporal/PA_treecover_temporal_inside_outside.csv")
crop_ts <- read.csv("outputs/temporal/PA_cropland_temporal_inside_outside_gaced.csv")
built_ts <- read.csv("outputs/temporal/PA_builtup_temporal_inside_outside.csv")

landuse_ts <- bind_rows(
  tree_ts |> select(NAME_ENG, location, year, change_pp, threat),
  crop_ts |> select(NAME_ENG, location, year, change_pp, threat),
  built_ts |> select(NAME_ENG, location, year, change_pp, threat)
) |>
  mutate(
    threat = factor(threat,
                    levels = c("Tree-cover loss", "Cropland", "Built-up")),
    location = factor(location, levels = c("Inside", "Outside"))
  )

p_time <- ggplot(
  landuse_ts,
  aes(
    x = year, y = change_pp,
    colour = threat, linetype = location,
    group = interaction(threat, location)
  )
) +
  geom_line(linewidth = 0.8) +
  #geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ NAME_ENG, scales = "free_y", ncol = 4) +
  scale_linetype_manual(values = c("Inside" = "solid", "Outside" = "dashed")) +
  #scale_colour_manual(values = c("forestgreen", "red3", "grey40")) +
  scale_colour_manual(values = c("#238B45", "#E68A00", "#7A3DB8")) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(
    x = "Protection effect (percentage points)",
    colour = NULL, linetype = NULL, tag = "a"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.01,1),
    legend.key.width = unit(1, "cm"),
    axis.title.y = element_text(margin = margin(r=0,l = 0,t = 0,b = 0))
  )

final_change_ha <- bind_rows(
  tree_ts |> filter(year == 2025) |>
    transmute(NAME_ENG, location, threat = "Tree-cover loss",
              change_ha = cumulative_loss_ha),
  crop_ts |> filter(year == 2024) |>
    transmute(NAME_ENG, location, threat = "Cropland",
              change_ha = crop_change_ha),
  built_ts |> filter(year == 2025) |>
    transmute(NAME_ENG, location, threat = "Built-up",
              change_ha = built_change_ha)
) |>
  pivot_wider(names_from = location, values_from = change_ha)

park_order <- sort(unique(final_change_ha$NAME_ENG))
final_change_ha$NAME_ENG <- factor(
  final_change_ha$NAME_ENG, levels = rev(park_order)
)
final_change_ha$threat <- factor(
  final_change_ha$threat,
  levels = c("Tree-cover loss", "Cropland", "Built-up")
)

shade_rows <- data.frame(
  ymin = seq(0.5, length(park_order) - 0.5, by = 2),
  ymax = seq(1.5, length(park_order) + 0.5, by = 2)
)

library(ggh4x)

p_total <- ggplot(final_change_ha, aes(y = NAME_ENG)) +
  geom_rect(
    data = shade_rows,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = "grey95", colour = NA
  ) +
  geom_segment(
    aes(x = Inside, xend = Outside, yend = NAME_ENG, colour = threat),
    linewidth = 0.7, show.legend = FALSE
  ) +
  geom_point(aes(x = Inside, shape = "Inside"), size = 3) +
  geom_point(aes(x = Outside, shape = "Outside"), size = 3) +
  facet_wrap(~ threat, scales = "free_x", nrow = 1,
             labeller = as_labeller(c(
               "Tree-cover loss" = "Tree-cover loss (ha)",
               "Cropland" = "Net cropland change (ha)",
               "Built-up" = "Built-up (ha)"
             ))) +
  #scale_colour_manual(values = c("forestgreen", "red2", "grey40")) +
  scale_colour_manual(values = c("#238B45", "#E68A00", "#7A3DB8")) +
  scale_shape_manual(values = c("Inside" = 16, "Outside" = 1)) +
  labs(x = NULL, y = NULL, shape = NULL, tag = "b") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.02,0.97),
    legend.text = element_text(size = 11),
    legend.key.spacing.y = unit(0, "cm"),
    legend.box.spacing = unit(0.0, "cm") ,
    legend.key.spacing.x = unit(0.0, "cm"),
    legend.key.height = unit(0,'cm'),
    legend.key.width = unit(1.5,'cm'),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(colour = "none")


library(ggh4x)

p_total <- p_total +
  ggh4x::facetted_pos_scales(
    x = list(
      threat == "Tree-cover loss" ~
        scale_x_continuous(
          limits = c(0, 100000),
          breaks = c(0, 50000, 100000),
          labels = scales::label_comma(),
          expand = expansion(mult = c(0.02, 0.02)),

        )
      ,threat == "Cropland" ~
        scale_x_continuous(
          breaks = scales::breaks_pretty(n = 3),
          labels = scales::label_comma()
        )
    )
  )

figure_2 <- p_time +
  inset_element(
    p_total,
    left = 0.28, bottom = 0.04,
    right = 1.00, top = 0.29,
    align_to = "full"
  )

ggsave(
  file.path(figure_dir, "Figure2__20oct.png"),
  figure_2, width = 9, height = 10, dpi = 600
)

library(grid)
library(png)

img <- readPNG('Banhine.png', native = TRUE)
park1 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Chimanimani.png', native = TRUE)
park2 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Gilé.png', native = TRUE)
park3 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Gorongosa.png', native = TRUE)
park4 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Limpopo.png', native = TRUE)
park5 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Magoe.png', native = TRUE)
park6 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Maputo.png', native = TRUE)
park7 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Marromeu.png', native = TRUE)
park8 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Niassa.png', native = TRUE)
park9 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Pomene.png', native = TRUE)
park10 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Quirimbas.png', native = TRUE)
park11 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Serra da Gorongosa.png', native = TRUE)
park12 <- rasterGrob(img, interpolate=TRUE)
img <- readPNG('Zinave.png', native = TRUE)
park13 <- rasterGrob(img, interpolate=TRUE)


library(grid)
library(png)

make_icon <- function(filename, opacity = 0.15) {
  img <- png::readPNG(filename)

  # Ensure the image has an alpha channel
  if (length(dim(img)) == 3 && dim(img)[3] == 3) {
    alpha <- matrix(1, nrow(img), ncol(img))
    img <- abind::abind(img, alpha, along = 3)
  }

  # Apply opacity to the existing alpha channel
  img[, , 4] <- img[, , 4] * opacity

  grid::rasterGrob(img, interpolate = TRUE)
}

park1  <- make_icon("Banhine.png")
park2  <- make_icon("Chimanimani.png")
park3  <- make_icon("Gilé.png")
park4  <- make_icon("Gorongosa.png")
park5  <- make_icon("Limpopo.png")
park6  <- make_icon("Magoe.png")
park7  <- make_icon("Maputo.png")
park8  <- make_icon("Marromeu.png")
park9  <- make_icon("Niassa.png")
park10 <- make_icon("Pomene.png")
park11 <- make_icon("Quirimbas.png")
park12 <- make_icon("Serra da Gorongosa.png")
park13 <- make_icon("Zinave.png")


library(cowplot)
final_plot <- ggdraw() +
  draw_plot(figure_2, 0, 0, 1, 1)+
  draw_plot(park1, x = 0.085, y = 0.75, width = 0.2, height = 0.2)+
  draw_plot(park2, x = 0.3, y = 0.75, width = 0.2, height = 0.2)+
  draw_plot(park3, x = 0.55, y = 0.75, width = 0.2, height = 0.2)+
  draw_plot(park4, x = 0.8, y = 0.75, width = 0.2, height = 0.2)+
  draw_plot(park5, x = 0.085, y = 0.53, width = 0.2, height = 0.2)+
  draw_plot(park6, x = 0.32, y = 0.53, width = 0.2, height = 0.2)+
  draw_plot(park7, x = 0.55, y = 0.53, width = 0.2, height = 0.2)+
  draw_plot(park8, x = 0.785, y = 0.53, width = 0.2, height = 0.2)+
  draw_plot(park9, x = 0.085, y = 0.30, width = 0.2, height = 0.2)+
  draw_plot(park10, x = 0.3, y = 0.30, width = 0.2, height = 0.2)+
  draw_plot(park11, x = 0.55, y = 0.30, width = 0.2, height = 0.2)+
  draw_plot(park12, x = 0.8, y = 0.25, width = 0.2, height = 0.2)+
  draw_plot(park13, x = 0.085, y = 0.1, width = 0.2, height = 0.2)




#final_plot

ggsave(plot = final_plot,filename = file.path(figure_dir, 'figure2_icons_20oct.png'),width = 9,height = 10,dpi = 600)


