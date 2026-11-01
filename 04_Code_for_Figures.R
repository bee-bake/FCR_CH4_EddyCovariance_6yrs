#This script contains code to produce all the manuscript figures including supplemental figures and tables.

# Create folder for figures to organize everything in one folder
figure_dir <- "Figures"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
}

#Figure 1
# Falling Creek Reservoir bathymetry with the eddy covariance flux footprint.

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2) # >= 3.5, for legend.position.inside
})

#paths
shape_dir <- "FCR bathymetry shapefile"
dir.create("figures", showWarnings = FALSE)
stem <- file.path("figures", "fcr_bathymetry_fcr_2d_manuscript_R")

#map extent
xlim <- c(602805, 603375)
ylim <- c(4129077, 4129839)
span_x <- diff(xlim)
span_y <- diff(ylim)

#colors
depth_bins <- c("0-2 m", "2-4 m", "4-6 m", "6-8 m", "8-9.3 m")
depth_cols <- setNames(c("#F7FBFF", "#D7E6F5", "#AED1E7", "#71B1D7", "#3B8BC2"), depth_bins)

fp_cols <- c(`10` = "#FFE87D", `20` = "#FDC502", `30` = "#F2A102", `40` = "#E47E01",
             `50` = "#D55901", `60` = "#C32C06", `70` = "#A30218", `80` = "#780220")

#sampling sites
sites <- data.frame(
  site = c("inflow weir", "eddy covariance", "met station"),
  x = c(603243.949, 603066.340, 603072.080),
  y = c(4129741.320, 4129121.070, 4129090.407),
  shape = c(22, 24, 21),
  size = c(5.0, 3.4, 4.1),
  fill = c("#CC79A7", "#000000", "#009E73"),
  edge = c("grey10", "white", "grey10")
)

# EC tower position
tower <- unlist(sites[sites$site == "eddy covariance", c("x", "y")])

#footprint isolines

iso70 <- matrix(c(
  0.0, 0.9, -22.1, 7.9, -48.7, 19.3, -75.2, 35.9, -95.6, 58.6,
  -108.0, 84.9, -109.7, 104.2, -99.1, 119.9, -77.9, 130.4, -48.7, 136.5,
  -19.5, 137.4, 7.1, 132.1, 28.3, 119.9, 42.5, 102.4, 48.7, 82.3,
  46.0, 60.4, 35.4, 38.5, 21.2, 19.3, 8.8, 6.1), ncol = 2, byrow = TRUE)
iso80 <- matrix(c(
  0.0, -2.6, -26.5, 3.5, -61.9, 16.6, -97.3, 34.1, -128.3, 60.4,
  -148.7, 91.0, -155.8, 121.7, -148.7, 147.9, -123.9, 165.4, -88.5, 174.2,
  -48.7, 178.6, -8.8, 175.9, 22.1, 167.2, 44.2, 150.5, 58.4, 128.7,
  63.7, 104.2, 60.2, 77.9, 48.7, 51.6, 33.6, 27.1, 17.7, 9.6,
  4.4, -0.9), ncol = 2, byrow = TRUE)

# closed periodic spline through the digitised points
smooth_ring <- function(v, n = 400) {
  v <- rbind(v, v[1, ])
  d <- c(0, cumsum(sqrt(diff(v[, 1])^2 + diff(v[, 2])^2)))
  u <- seq(0, 1, length.out = n)
  cbind(spline(d / max(d), v[, 1], method = "periodic", xout = u)$y,
        spline(d / max(d), v[, 2], method = "periodic", xout = u)$y)
}

ring70 <- smooth_ring(iso70)
shrink <- function(sx, sy) cbind(ring70[, 1] * sx, ring70[, 2] * sy)
rings <- list(
  `10` = shrink(0.112, 0.090), `20` = shrink(0.184, 0.153), `30` = shrink(0.268, 0.223),
  `40` = shrink(0.346, 0.299), `50` = shrink(0.436, 0.395), `60` = shrink(0.559, 0.535),
  `70` = ring70, `80` = smooth_ring(iso80)
)

fp <- do.call(rbind, Map(function(xy, lv) {
  data.frame(x = xy[, 1] + tower[1], y = xy[, 2] + tower[2], level = lv)
}, rings, names(rings)))
fp$level <- factor(fp$level, levels = rev(names(rings)))
fp$traced <- fp$level %in% c("70", "80")

#shapefiles

contours <- c("506_9", "504_9", "502_9", "500_9", "498_9")
bathy <- do.call(rbind, Map(function(f, bin) {
  g <- st_read(file.path(shape_dir, paste0(f, ".shp")), quiet = TRUE)
  st_sf(depth = bin, geometry = st_geometry(g))
}, contours, depth_bins))
inflow <- st_read(file.path(shape_dir, "inf2.shp"), quiet = TRUE)
crs_utm <- st_crs(bathy)

#scale bar, north arrow and footprint key
bar_x <- xlim[2] - 0.055 * span_x - 200
bar_y <- ylim[1] + span_y * c(0.036, 0.048)
arrow_x <- bar_x + 100
arrow_y <- bar_y[2] + span_y * c(0.055, 0.089)

key_x <- xlim[1] + span_x * c(0.62, 0.98)
key_y <- ylim[1] + span_y * c(0.228, 0.252)
cell <- diff(key_x) / length(fp_cols)
key_left <- key_x[1] + cell * (seq_along(fp_cols) - 1)
key_mid <- key_left + cell / 2

#plot
p <- ggplot() +
  
  # Bathymetry
  geom_sf(
    data = bathy,
    aes(fill = depth),
    colour = "black",
    linewidth = 0.13
  ) +
  
  # Primary inflow stream
  geom_sf(
    data = inflow,
    aes(linetype = "primary inflow stream"),
    colour = "darkgray",
    linewidth = 0.45
  ) +
  
  # Footprint isolines
  # traced 70/80% lines heavier than the approximate inner rings
  geom_path(
    data = subset(fp, traced),
    aes(x, y, group = level, colour = level),
    linewidth = 1.3
  ) +
  
  geom_path(
    data = subset(fp, !traced),
    aes(x, y, group = level, colour = level),
    linewidth = 0.95
  ) +
  
  # Sampling sites
  geom_point(
    data = sites,
    aes(x, y, shape = site),
    size = sites$size,
    fill = sites$fill,
    colour = sites$edge,
    stroke = 0.5
  ) +
  
  #scale bar

annotate(
  "rect",
  xmin = bar_x + 50 * 0:3,
  xmax = bar_x + 50 * 1:4,
  ymin = bar_y[1],
  ymax = bar_y[2],
  fill = c(
    "grey35",
    "white",
    "grey35",
    "white"
  ),
  colour = "grey25",
  linewidth = 0.15
) +
  
  annotate(
    "text",
    x = bar_x + 50 * 0:4,
    y = bar_y[1] - 0.008 * span_y,
    label = c(
      "0",
      "50",
      "100",
      "150",
      "200 m"
    ),
    vjust = 1,
    size = 13,
    size.unit = "pt",
    family = "serif",
    colour = "grey15"
  ) +
  
  #north arrow

annotate(
  "segment",
  x = arrow_x,
  xend = arrow_x,
  y = arrow_y[1],
  yend = arrow_y[2],
  colour = "grey35",
  linewidth = 1.25,
  arrow = arrow(
    length = unit(8, "pt"),
    type = "closed"
  )
) +
  
  annotate(
    "text",
    x = arrow_x,
    y = arrow_y[2] + 0.006 * span_y,
    label = "N",
    vjust = 0,
    size = 15,
    size.unit = "pt",
    family = "serif",
    fontface = "bold",
    colour = "grey15"
  ) +
  
  #footprint contribution key

annotate(
  "rect",
  xmin = key_left,
  xmax = key_left + cell,
  ymin = key_y[1],
  ymax = key_y[2],
  fill = fp_cols
) +
  
  annotate(
    "rect",
    xmin = key_x[1],
    xmax = key_x[2],
    ymin = key_y[1],
    ymax = key_y[2],
    fill = NA,
    colour = "grey35",
    linewidth = 0.2
  ) +
  
  annotate(
    "segment",
    x = key_mid,
    xend = key_mid,
    y = key_y[1],
    yend = key_y[1] - 0.003 * span_y,
    colour = "grey20",
    linewidth = 0.2
  ) +
  
  annotate(
    "text",
    x = key_mid,
    y = key_y[1] - 0.006 * span_y,
    label = names(fp_cols),
    vjust = 1,
    size = 13,
    size.unit = "pt",
    family = "serif",
    colour = "grey20"
  ) +
  
  annotate(
    "text",
    x = mean(key_x),
    y = key_y[2] + 0.005 * span_y,
    label = "flux contribution (%)",
    vjust = 0,
    size = 13,
    size.unit = "pt",
    family = "serif",
    colour = "grey20"
  ) +
  
  #scales

scale_fill_manual(
  values = depth_cols,
  breaks = rev(depth_bins),
  name = "Depth (m)"
) +
  
  scale_colour_manual(
    values = fp_cols,
    guide = "none"
  ) +
  
  scale_linetype_manual(
    values = "solid",
    name = NULL
  ) +
  
  scale_shape_manual(
    values = setNames(
      sites$shape,
      sites$site
    ),
    breaks = sites$site,
    name = NULL
  ) +
  
  #legends

guides(
  fill = guide_legend(
    order = 1,
    override.aes = list(
      colour = "black",
      linewidth = 0.5
    )
  ),
  
  linetype = guide_legend(
    order = 2
  ),
  
  shape = guide_legend(
    order = 3,
    override.aes = list(
      size = sites$size,
      fill = sites$fill,
      colour = sites$edge,
      stroke = 0.5
    )
  )
) +
  
#coordinate system
# lon/lat labels come from the coordinate scales

coord_sf(
  xlim = xlim,
  ylim = ylim,
  expand = FALSE,
  crs = crs_utm
) +
  
  scale_x_continuous(
    breaks = seq(
      -79.839,
      -79.834,
      by = 0.001
    ),
    labels = function(x) {
      sprintf(
        "%.3f\u00b0W",
        abs(x)
      )
    }
  ) +
  
  scale_y_continuous(
    breaks = seq(
      37.303,
      37.309,
      by = 0.001
    ),
    labels = function(y) {
      sprintf(
        "%.3f\u00b0N",
        y
      )
    }
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  #theme

theme_bw(
  base_size = 14,
  base_family = "serif"
) +
  
  theme(
    
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "grey35",
      linewidth = 0.8
    ),
    
    #axis text
    
    axis.text = element_text(
      size = 14,
      colour = "grey20"
    ),
    
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 14
    ),
    
    axis.text.y = element_text(
      size = 14
    ),
    
    axis.title = element_text(
      size = 15
    ),
    
    #legend
    
    legend.position = "inside",
    
    legend.position.inside = c(
      0.975,
      0.52
    ),
    
    legend.justification.inside = c(
      1,
      0.5
    ),
    
    legend.background = element_blank(),
    
    legend.margin = margin(
      1,
      0,
      1,
      0
    ),
    
    legend.box.background = element_rect(
      fill = "white",
      colour = "grey55",
      linewidth = 0.4
    ),
    
    legend.box.margin = margin(
      3,
      5,
      3,
      5
    ),
    
    legend.spacing.y = unit(
      2,
      "pt"
    ),
    
    legend.key = element_blank(),
    
    legend.key.width = unit(
      20,
      "pt"
    ),
    
    legend.key.height = unit(
      17,
      "pt"
    ),
    
    legend.key.spacing.y = unit(
      2,
      "pt"
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    legend.title = element_text(
      size = 15,
      hjust = 0.5
    ),
    
    # margins
    
    plot.margin = margin(
      5,
      6,
      3,
      3
    )
  )


# Save JPG
ggsave(
  filename = file.path(figure_dir, paste0(stem, ".jpg")),
  plot = p,
  width = 210,
  height = 270,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# PDF device
pdf_device <- if (capabilities("aqua")) {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}

# Save PDF
ggsave(
  filename = file.path(figure_dir, paste0(stem, ".pdf")),
  plot = p,
  device = pdf_device,
  width = 210,
  height = 270,
  units = "mm",
  bg = "white"
)

############
#Figure 2
# Prepare data
ch4_timeseries <- comp_data %>%
  mutate(date = as.Date(date)) %>%
  filter(
    !is.na(date),
    !is.na(daily_CH4)
  ) %>%
  arrange(date)

# Define x-axis limits and yearly breaks
x_limits <- as.Date(c("2020-01-01", "2026-12-31"))

year_breaks <- seq(
  as.Date("2020-01-01"),
  as.Date("2026-01-01"),
  by = "1 year"
)

# Plot methane flux time series
p_ch4_timeseries <- ggplot(
  ch4_timeseries,
  aes(x = date, y = daily_CH4)
) +
  geom_point(
    color = "black",
    alpha = 0.85,
    size = 1
  ) +
  scale_x_date(
    limits = x_limits,
    breaks = year_breaks,
    date_labels = "%Y",
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Year",
    y = expression(CH[4]~flux~(mg~m^{-2}~d^{-1})),
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.ticks.x = element_line(color = "black")
  )

# Display plot
print(p_ch4_timeseries)

# Save high-resolution figure
ggsave(
  filename = file.path(figure_dir,"FCR_CH4_timeseries.png"),
  plot = p_ch4_timeseries,
  width = 200,
  height = 110,
  units = "mm",
  dpi = 600,
  bg = "white"
)


############
#Figure 3:
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

#PREPARE DATA
comp_data <- comp_data %>%
  mutate(
    date = as.Date(date)
  )

# NOTE: this filter selects the DATA range and stays at 2020-01-01.
# The plotted x-axis is extended further left (see section 6b) so the
# panel labels sit in empty space instead of on top of the data.
plot_data <- comp_data %>%
  filter(
    date >= as.Date("2020-01-01"),
    date <= as.Date("2026-12-31")
  )


#Select DO data for different depths from plot data created above
do_data <- plot_data %>%
  select(
    date,
    daily_CH4,
    DOsat_surface,
    DOsat_5m,
    DOsat_9m
  ) %>%
  filter(
    !is.na(daily_CH4)
  ) %>%
  pivot_longer(
    cols = c(
      DOsat_surface,
      DOsat_5m,
      DOsat_9m
    ),
    names_to = "depth",
    values_to = "DO"
  ) %>%
  mutate(
    depth = case_when(
      depth == "DOsat_surface" ~ "1.6 m",
      depth == "DOsat_5m"      ~ "5 m",
      depth == "DOsat_9m"      ~ "9 m"
    )
  ) %>%
  filter(
    !is.na(DO)
  )


#Select temperature data for different depths from plot data created above
temp_data <- plot_data %>%
  select(
    date,
    daily_CH4,
    temp_1.6m,
    temp_5m,
    temp_9m
  ) %>%
  filter(
    !is.na(daily_CH4)
  ) %>%
  pivot_longer(
    cols = c(
      temp_1.6m,
      temp_5m,
      temp_9m
    ),
    names_to = "depth",
    values_to = "temperature"
  ) %>%
  mutate(
    depth = case_when(
      depth == "temp_1.6m" ~ "1.6 m",
      depth == "temp_5m"   ~ "5 m",
      depth == "temp_9m"   ~ "9 m"
    )
  ) %>%
  filter(
    !is.na(temperature)
  )

#Create color scale for the plot

ch4_lim <- quantile(
  plot_data$daily_CH4,
  probs = c(0.10, 0.90),
  na.rm = TRUE
)

ch4_colors <- c(
  "#542788",   # dark purple - low
  "#998EC3",   # light purple
  "#009E73",   # green
  "#F1A340",   # light orange
  "#E08214",   # orange
  "#B35806"    # dark orange - high
)

#Adjust text size for the plot

figure_text_size <- 17
panel_label_size <- figure_text_size / 3

#Create common theme for the plot
common_theme <- theme_classic(
  base_size = figure_text_size,
  base_family = "serif"
) +
  theme(
    # Axis scale numbers
    axis.text = element_text(
      size = figure_text_size,
      colour = "black"
    ),
    
    # Individual y-axis titles OFF by default
    axis.title.y = element_blank(),
    
    # Individual x-axis titles OFF by default
    axis.title.x = element_blank(),
    
    # Axis ticks
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.ticks.length = unit(
      0.18,
      "cm"
    ),
    
    # No gridlines
    panel.grid = element_blank(),
    
    # Legend
    legend.title = element_text(
      size = figure_text_size,
      colour = "black",
      face = "plain"
    ),
    
    legend.text = element_text(
      size = figure_text_size,
      colour = "black"
    ),
    
    legend.key.width = unit(
      1.6,
      "cm"
    ),
    
    legend.key.height = unit(
      0.5,
      "cm"
    ),
    
    
    # Margins
    plot.margin = margin(
      t = 10,
      r = 14,
      b = 7,
      l = 14
    )
  )


#X-axis settings (extra space on the left for panel labels)
#The axis starts a year before the data so the (a)-(f) labels sit in
#empty space rather than on top of the 2020 points.
x_start <- as.Date("2019-09-01")
x_end   <- as.Date("2026-12-31")
label_x <- as.Date("2019-10-01")

#Ticks only on the real data years (no empty "2019" tick)
year_breaks <- seq(
  as.Date("2020-01-01"),
  as.Date("2026-01-01"),
  by = "1 year"
)


#1.6 m DO
p_do_1.6m <- ggplot(
  do_data %>% filter(depth == "1.6 m"),
  aes(
    x = date,
    y = DO,
    colour = daily_CH4
  )
) +
  
  annotate(
    "rect",
    xmin = x_start,
    xmax = x_end,
    ymin = 0,
    ymax = 20,
    fill = "grey75",
    alpha = 0.35
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 180),
    breaks = seq(0, 180, 30),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  common_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(a) 1.6 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )

#5 m DO


p_do_5m <- ggplot(
  do_data %>% filter(depth == "5 m"),
  aes(
    x = date,
    y = DO,
    colour = daily_CH4
  )
) +
  
  annotate(
    "rect",
    xmin = x_start,
    xmax = x_end,
    ymin = 0,
    ymax = 20,
    fill = "grey75",
    alpha = 0.35
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 180),
    breaks = seq(0, 180, 30),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "DO saturation (%)"
  ) +
  
  common_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_text(
      angle = 90,
      size = figure_text_size,
      colour = "black",
      face = "plain",
      margin = margin(r = 18)
    )
  ) +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(b) 5 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )

#9 m DO
p_do_9m <- ggplot(
  do_data %>% filter(depth == "9 m"),
  aes(
    x = date,
    y = DO,
    colour = daily_CH4
  )
) +
  
  annotate(
    "rect",
    xmin = x_start,
    xmax = x_end,
    ymin = 0,
    ymax = 20,
    fill = "grey75",
    alpha = 0.35
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 180),
    breaks = seq(0, 180, 30),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  common_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(c) 9 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )


#1.6 m temperature
p_temp_1.6m <- ggplot(
  temp_data %>% filter(depth == "1.6 m"),
  aes(
    x = date,
    y = temperature,
    colour = daily_CH4
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  common_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(d) 1.6 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )

#5 m temperature
p_temp_5m <- ggplot(
  temp_data %>% filter(depth == "5 m"),
  aes(
    x = date,
    y = temperature,
    colour = daily_CH4
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Temperature (°C)"
  ) +
  
  common_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_text(
      angle = 90,
      size = figure_text_size,
      colour = "black",
      face = "plain",
      margin = margin(r = 18)
    )
  ) +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(e) 5 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )

#9 m temperature
p_temp_9m <- ggplot(
  temp_data %>% filter(depth == "9 m"),
  aes(
    x = date,
    y = temperature,
    colour = daily_CH4
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_colour_gradientn(
    colours = ch4_colors,
    limits = ch4_lim,
    oob = scales::squish,
    na.value = NA,
    name = expression(
      CH[4]~flux~"(mg"~m^{-2}~d^{-1}*")"
    )
  ) +
  
  scale_x_date(
    limits = c(x_start, x_end),
    breaks = year_breaks,
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +
  
  labs(
    x = "Date",
    y = NULL
  ) +
  
  common_theme +
  
  annotate(
    "text",
    x = label_x,
    y = Inf,
    label = "(f) 9 m",
    hjust = 0,
    vjust = 1.3,
    size = panel_label_size,
    fontface = "bold"
  )

#combine all the plots
hypolimnetic_figure <-
  p_do_1.6m /
  p_do_5m /
  p_do_9m /
  p_temp_1.6m /
  p_temp_5m /
  p_temp_9m +
  
  plot_layout(
    ncol = 1,
    guides = "collect"
  ) &
  
  theme(
    legend.position = "top"
  )

#display
hypolimnetic_figure

#save as png
ggsave(
  filename = file.path(figure_dir,"hypolimnetic_DO_temperature_CH4_6x1.png"),
  plot = hypolimnetic_figure,
  width = 275,
  height = 330,
  units = "mm",
  dpi = 600,
  bg = "white"
)

#save as pdf
ggsave(
  filename = file.path(figure_dir,"hypolimnetic_DO_temperature_CH4_6x1.pdf"),
  plot = hypolimnetic_figure,
  width = 275,
  height = 330,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

###############
#Figure 4: Hypolimnetic DO series for all years
# 1. Prepare data for 5 m and 9 m
phase_data <- bind_rows(
  
  comp_data %>%
    mutate(
      year = year(date),
      depth = "5 m",
      temp = temp_5m,
      DOsat = DOsat_5m
    ) %>%
    select(year, depth, temp, DOsat),
  
  comp_data %>%
    mutate(
      year = year(date),
      depth = "9 m",
      temp = temp_9m,
      DOsat = DOsat_9m
    ) %>%
    select(year, depth, temp, DOsat)
  
) %>%
  filter(
    !is.na(temp),
    !is.na(DOsat)
  )

# 2. Determine temperature range where DO < 20%
hypoxic_ranges <- phase_data %>%
  filter(DOsat < 20) %>%
  group_by(depth) %>%
  summarise(
    temp_min_hypoxic = min(temp, na.rm = TRUE),
    temp_max_hypoxic = max(temp, na.rm = TRUE),
    .groups = "drop"
  )

print(hypoxic_ranges)


# 3. Function to create each panel
make_phase_plot <- function(depth_label, panel_label,
                            show_y_title = TRUE) {
  
  depth_data <- phase_data %>%
    filter(depth == depth_label)
  
  depth_range <- hypoxic_ranges %>%
    filter(depth == depth_label)
  
  ggplot(
    depth_data,
    aes(
      x = temp,
      y = DOsat,
      color = factor(year)
    )
  ) +
    
    # Highlight temperature range in which hypoxia occurs
    geom_rect(
      data = depth_range,
      aes(
        xmin = temp_min_hypoxic,
        xmax = temp_max_hypoxic,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.35
    ) +
    
    # Observations
    geom_point(
      size = 2.5,
      alpha = 0.70
    ) +
    
    # Hypoxic threshold
    geom_hline(
      yintercept = 20,
      linetype = "dashed",
      linewidth = 0.9,
      color = "black"
    ) +
    
    # Panel label
    labs(
      title = panel_label,
      x = "Hypolimnetic temperature (°C)",
      y = if (show_y_title) {
        "Hypolimnetic oxygen (% saturation)"
      } else {
        NULL
      },
      color = "Year"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      # Place (a) and (b) at the upper-left of each panel
      plot.title = element_text(
        face = "plain",
        size = 13,
        hjust = 0,
        margin = margin(b = 8)
      ),
      
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11),
      
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      
      # Panel border
      panel.border = element_rect(
        colour = "grey40",
        fill = NA,
        linewidth = 0.6
      ),
      
      # Space around each panel
      plot.margin = margin(8, 8, 8, 8)
    )
}


# 4. Create the two panels
p5 <- make_phase_plot(
  depth_label = "5 m",
  panel_label = "(a) 5m",
  show_y_title = TRUE
)

p9 <- make_phase_plot(
  depth_label = "9 m",
  panel_label = "(b) 9m",
  show_y_title = FALSE
)


# 5. Combine panels and collect one shared legend
p_phase_year <- (p5 + p9) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )


# 6. Display Figure 4
print(p_phase_year)

ggsave(
  filename = file.path(figure_dir,"temp_oxygen.png"),  plot = p_phase_year, width = 11, height = 6, dpi = 300)

##############

###########
#Figure 5: SHAP dependency and beeswarm plots
#Note: To run this, you'll need to run '03_XGBoost_SHAP_model' first!!!
p_imp_A <- sv_importance(shp_A, show_numbers = TRUE) + scale_y_discrete(
  labels = shap_labels) +
  labs(title = "") +
  theme_classic()

p_bee_A <- sv_importance(shp_A, kind = "beeswarm") + scale_y_discrete(
  labels = shap_labels) +
  labs(title = "") +
  scale_color_gradient(low = "blue", high = "red",
                       na.value = NA, name = "Feature value") +
  labs(title = "") + theme_classic()

print(p_imp_A)
print(p_bee_A)

#Figure for paper
p1 <- sv_dependence(shp_A, v = "DOsat_5m", color_var = NULL) +
  labs(x = "% DO saturation (5 m)", y = "SHAP value") 
p2 <- sv_dependence(shp_A, v = "AirTemp_C_Average", color_var = NULL)
p3 <- sv_dependence(shp_A, v = "temp_5m", color_var = NULL) +
  labs(x = "Temperature (5 m)", y = "SHAP value")
p4 <- sv_dependence(shp_A, v = "DOsat_9m", color_var = NULL)
p5 <- sv_dependence(shp_A, v = "temp_9m", color_var = NULL)
p6 <- sv_dependence(shp_A, v = "WindSpeed_Average_m_s", color_var = NULL)

p1 + p2 + p3 + p6

# Figure 3
p1 <- p1 + theme_classic()
p3 <- p3 + theme_classic()

shap_plots <- 
  (p_imp_A + p_bee_A) /
  p1 /
  p3 +
  plot_layout(
    heights = c(1, 1.3, 1.3)
  ) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )

ggsave(
  filename = file.path(figure_dir,"Figure_3_SHAP.png"),
  shap_plots,
  width = 11,
  height = 12,
  units = "in",
  dpi = 600
)


# SENSITIVITY ANALYSIS: Hypoxic threshold at 15, 20, 25 % DO saturation
model_data_raw2 <- comp_data |>
  mutate(
    month = month(date),
    year  = year(date)
  ) |>
  dplyr::filter(
    month %in% 5:9,
    !is.na(daily_CH4)
  ) |>
  arrange(date) |>
  select(
    all_of(c(
      "daily_CH4",
      xvars_clean,
      "DOsat_5m",
      "month",
      "year"
    ))
  )

# The same observations are used for all three DO thresholds.
# This is important because it makes the sensitivity analysis
# test only the effect of changing the threshold.

set.seed(42)

sample_index <- sample(
  seq_len(nrow(model_data_raw2)),
  min(149, nrow(model_data_raw2))
)

X_explain_5m <- model_data_raw2[
  sample_index,
  ,
  drop = FALSE
]

#Create the SHAP prediction matrix
# Since, fit_raw was trained using xvars_clean only,
# DOsat_5m is NOT an additional model predictor.
# It is retained separately for DO classification.
X_pred_5m <- as.matrix(
  X_explain_5m[
    ,
    xvars_clean,
    drop = FALSE
  ]
)


# Verify predictor structure
cat(
  "Number of predictors supplied to XGBoost:",
  ncol(X_pred_5m),
  "\n"
)

cat(
  "Number of predictors expected by model:",
  fit_raw$nfeatures,
  "\n"
)

print(colnames(X_pred_5m))

print(
  identical(
    colnames(X_pred_5m),
    xvars_clean
  )
)

#Calculate SHAP values
shp_2020_2026 <- shapviz(
  fit_raw,
  X_pred = X_pred_5m,
  X      = X_explain_5m[, xvars_clean, drop = FALSE]
)

#Create common 5-m SHAP dataset
# DOsat_5m remains here specifically for the sensitivity analysis.
plot_5m_base <- data.frame(
  temperature = X_explain_5m$temp_5m,
  SHAP        = shp_2020_2026$S[, "temp_5m"],
  DOsat_5m    = X_explain_5m$DOsat_5m
) |>
  filter(
    !is.na(temperature),
    !is.na(SHAP),
    !is.na(DOsat_5m)
  )

# Check the data
summary(plot_5m_base)


#Sensitivity analysis
thresholds <- c(15, 20, 25)
set.seed(42)
n_perm <- 5000
sensitivity_results <- list()
for (thresh in thresholds) {
  
  dat <- plot_5m_base |>
    mutate(DOgroup = ifelse(DOsat_5m < thresh, "Hypoxic", "Oxic"),
           DOgroup = factor(DOgroup, levels = c("Oxic", "Hypoxic"))) |>
    filter(!is.na(SHAP), !is.na(temperature), !is.na(DOgroup))
  
  # sample sizes
  n_tab <- dat |> count(DOgroup)
  n_oxic    <- n_tab$n[n_tab$DOgroup == "Oxic"]
  n_hypoxic <- n_tab$n[n_tab$DOgroup == "Hypoxic"]
  
  # fit interaction model
  m <- lm(SHAP ~ temperature * DOgroup, data = dat)
  coefs <- summary(m)$coefficients
  
  slope_oxic     <- coefs["temperature", "Estimate"]
  interaction_b  <- coefs["temperature:DOgroupHypoxic", "Estimate"]
  interaction_p  <- coefs["temperature:DOgroupHypoxic", "Pr(>|t|)"]
  slope_hypoxic  <- slope_oxic + interaction_b
  
  # permutation test
  observed_diff <- interaction_b
  perm_diffs <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
      shuffled <- dat
        shuffled$DOgroup <- sample(
      shuffled$DOgroup
    )
    pm <- lm(
      SHAP ~ temperature * DOgroup,
      data = shuffled
    )
    perm_diffs[i] <- coef(pm)[
      grepl(
        "temperature:DOgroup",
        names(coef(pm))
      )
    ]
  }
  
  # Two-sided permutation p-value
  perm_p <- mean(
    abs(perm_diffs) >= abs(observed_diff)
  )
  
  # Format permutation p-value
  perm_p_display <- if (perm_p == 0) {
    
    paste0(
      "< ",
      format(
        1 / n_perm,
        scientific = FALSE
      )
    )
    
  } else {
    
    formatC(
      perm_p,
      digits = 4,
      format = "f"
    )
  }
  
  # report as "< 1/n_perm" when zero, since a permutation test can't
  # certify an exact zero probability; otherwise show 4 decimal places
  perm_p_display <- if (perm_p == 0) {
    paste0("< ", format(1 / n_perm, scientific = FALSE))
  } else {
    formatC(perm_p, digits = 4, format = "f")
  }
  
  sensitivity_results[[as.character(thresh)]] <- data.frame(
    threshold_mgL      = thresh,
    n_oxic              = n_oxic,
    n_hypoxic           = n_hypoxic,
    slope_oxic          = round(slope_oxic, 3),
    slope_hypoxic       = round(slope_hypoxic, 3),
    slope_diff          = round(interaction_b, 3),
    parametric_p        = signif(interaction_p, 3),
    permutation_p       = perm_p_display
  )
}
sensitivity_table <- bind_rows(sensitivity_results)
print(sensitivity_table)

library(knitr)
library(kableExtra)

sensitivity_table |>
  kable(
    caption = "Sensitivity of temperature–SHAP slope divergence to hypoxia threshold",
    col.names = c("Threshold (mg/L)", "n (Oxic)", "n (Hypoxic)",
                  "Slope (Oxic)", "Slope (Hypoxic)", "Slope Diff.",
                  "Parametric p", "Permutation p"),
    align = "c"
  ) |>
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)


#Across hypoxia thresholds of 15%, 20%, and 25% DO saturation, the temperature–SHAP
#relationship was consistently steeper under hypoxic than oxic conditions, although
#the strength of this interaction varied among thresholds. At the 15% threshold, the
#temperature–SHAP slopes were 1.307 and 1.732 SHAP units per °C under oxic and hypoxic
#conditions, respectively, corresponding to a slope difference of 0.426 SHAP units per 
#°C (parametric p = 0.044; permutation p = 0.0422). At the 20% threshold, the divergence
#was strongest, with slopes of 1.105 and 2.032 SHAP units per °C under oxic and hypoxic 
#conditions, respectively (slope difference = 0.927; parametric p = <0.0001; permutation p < 0.0002).
#At the 25% threshold, the hypoxic slope remained steeper than the oxic slope (1.845 vs. 1.130 SHAP units per °C;
#slope difference = 0.715; parametric p < 0.0001; permutation p = 0.0002). 
#Thus, the direction and statistical significance of the temperature–SHAP interaction were robust to the choice 
#of hypoxia threshold, although the magnitude of the slope difference was greatest at the 20% threshold. 
#The changing group sizes across thresholds (57–75 hypoxic observations and 71–89 oxic observations) should be 
#considered when comparing the magnitude of the interaction among thresholds.


############
#Figure 6: Temperature-DO interaction plot 
# Entire 5 m dataset, 2020–2026
# Summer CH4 fluxes only (May–September)

set.seed(42)

# Prepare data

model_data_raw2 <- comp_data |>
  mutate(
    month = month(date),
    year  = year(date)
  ) |>
  dplyr::filter(
    month %in% 5:9,
    !is.na(daily_CH4)
  ) |>
  arrange(date) |>
  select(
    all_of(c(
      "daily_CH4",
      xvars_clean,
      "DOsat_5m",
      "month",
      "year"
    ))
  )


# Sample the entire dataset
# Use the same maximum sample size as before

X_explain_5m <- model_data_raw2[
  sample(
    nrow(model_data_raw2),
    min(411, nrow(model_data_raw2))
  ),
  c(xvars_clean, "DOsat_5m")
]


# Calculate SHAP values for the full dataset
# Uses the existing model: fit_raw

X_pred_5m <- as.matrix(
  X_explain_5m[, xvars_clean]
)

shp_5m <- shapviz(
  fit_raw,
  X_pred = X_pred_5m
)

#Build plotting data

hypoxic_threshold <- 20  # % saturation

plot_5m <- data.frame(
  temperature = X_explain_5m$temp_5m,
  SHAP        = shp_5m$S[, "temp_5m"],
  DO_mgL      = X_explain_5m$DOsat_5m
) |>
  mutate(
    DOgroup = ifelse(
      DO_mgL < hypoxic_threshold,
      "Hypoxic",
      "Oxic"
    ),
    DOgroup = factor(
      DOgroup,
      levels = c("Oxic", "Hypoxic")
    )
  )


#Sample size per DO group

n_by_group <- plot_5m |>
  filter(
    !is.na(SHAP),
    !is.na(temperature)
  ) |>
  count(
    DOgroup,
    name = "n"
  )

print(n_by_group)


#Temperature range by DO group

range_by_group <- plot_5m |>
  filter(
    !is.na(SHAP),
    !is.na(temperature)
  ) |>
  group_by(DOgroup) |>
  summarise(
    n         = n(),
    temp_min  = round(min(temperature, na.rm = TRUE), 1),
    temp_max  = round(max(temperature, na.rm = TRUE), 1),
    temp_mean = round(mean(temperature, na.rm = TRUE), 1),
    .groups   = "drop"
  )

#Supplementary material Table 1
print(range_by_group)


#Temperature × DO interaction model

interaction_model_5m <- lm(
  SHAP ~ temperature * DOgroup,
  data = plot_5m
)

#Supplementary material Table 2
summary(interaction_model_5m)


#Plot

merged_plot <- ggplot(
  plot_5m,
  aes(
    x = temperature,
    y = SHAP,
    color = DOgroup
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.75
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  
  scale_color_manual(
    values = c(
      "Oxic" = "#E8A33D",
      "Hypoxic" = "#7B52AE"
    ),
    na.translate = FALSE
  ) +
  
  labs(
    title = "",
    x = "Temperature (5 m, °C)",
    y = "SHAP value",
    color = "DO group"
  ) +
  
  theme_classic()

ggsave(
  filename = file.path(figure_dir,"SHAP_temp_by_DOgroup.png"), plot = merged_plot, width = 6, height = 5, dpi = 300)

#The temperature–SHAP relationship was significantly steeper under hypoxic conditions than 
#oxic conditions, indicating that temperature sensitivity of predicted CH4 flux is substantially amplified when bottom-water oxygen is limiting — consistent with 
#reduced methanotrophic consumption under hypoxia allowing a more direct translation of thermally-stimulated
#methanogenesis into surface flux. Oxic (n = 75) and hypoxic (n = 71) observations spanned overlapping temperature
#ranges (oxic: 8.8–13.8°C, mean 11.5°C; hypoxic: 9.4–12.6°C, mean 10.7°C), indicating that the steeper temperature
#sensitivity observed under hypoxic conditions was not attributable to a systematic difference in the temperature
#range sampled between groups.

#Because SHAP values derived from a single fitted model are not independent observations, we complemented the parametric
#interaction-model test with a permutation test (5,000 iterations, shuffling DO group labels) to obtain a distribution-free
#significance estimate robust to this non-independence.
# Permutation test (same logic as before, shuffling DOgroup instead of Period)
set.seed(42)
n_perm <- 5000
observed_diff_do <- coef(interaction_model_5m)[grepl("temperature:DOgroup",
                                                     names(coef(interaction_model_5m)))]
perm_diffs_do <- numeric(n_perm)
for (i in 1:n_perm) {
  shuffled <- plot_5m
  shuffled$DOgroup <- sample(shuffled$DOgroup)
  m <- lm(SHAP ~ temperature * DOgroup, data = shuffled)
  perm_diffs_do[i] <- coef(m)[grepl("temperature:DOgroup", names(coef(m)))]
}
perm_p_do <- mean(abs(perm_diffs_do) >= abs(observed_diff_do))

# +1 correction so p is never reported as exactly zero
perm_p_do_corrected <- (
  sum(abs(perm_diffs_do) >= abs(observed_diff_do)) + 1
) / (n_perm + 1)

cat(
  "Observed slope difference (Oxic vs Hypoxic):",
  round(observed_diff_do, 3), "\n"
)

cat(
  "Permutation p-value:",
  round(perm_p_do, 4), "\n"
)

cat(
  "Corrected permutation p-value:",
  round(perm_p_do_corrected, 4), "\n"
)


# Visualize the null distribution vs observed statistic
ggplot(data.frame(perm_diffs), aes(x = perm_diffs_do)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = observed_diff_do, color = "firebrick", 
             linewidth = 1, linetype = "dashed") +
  labs(title = "Permutation null distribution of slope difference",
       subtitle = paste("Observed difference =", round(observed_diff_do, 3),
                        "| p =", round(perm_p, 4)),
       x = "Slope difference (temperature:Period)", y = "Count") +
  theme_classic()


############
#Supplemental Information Figures
#Figure S1
#Check correlation
#prepare data
cor_data <- comp_data %>%
  select(temp_1.6m, temp_5m, temp_9m, DOsat_surface, DOsat_9m, DOsat_5m, mean_fDOM, mean_Chla_lag1,
         WindSpeed_Average_m_s, AirTemp_C_Average, DOY, daily_CH4, sch_stab, thermo.depth, DC_mgL)

#correlation matrix
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

#cleaner plot
corrplot(
  cor_matrix,
  method   = "color",
  type     = "upper",
  order    = "hclust",
  
  # appearance
  col      = colorRampPalette(
    c("#3B4CC0", "white", "#B40426")
  )(80),
  
  # labels
  tl.col   = "black",
  tl.cex   = 0.6,
  tl.srt   = 35,
  
  # show correlation coefficients
  diag        = FALSE,
  addCoef.col = "black",
  number.cex  = 0.6,     # size of the coefficient text
  number.digits = 2,      # decimal places shown
  
  # remove visual clutter
  cl.cex   = 0.7,
  mar      = c(0, 0, 1, 0),
  
  title    = "Correlation Matrix"
)

############
#Figure S2
#Check autocorrelation and lags
AR_df <- comp_data |>
  dplyr::select(-date, -Date)

# remove NA values
chla_series <- na.omit(AR_df$mean_Chla)

# PACF plot
pacf(
  chla_series,
  lag.max = 30,   # adjust as needed
  main    = "",
  col     = "#2166AC",
  lwd     = 2
)

#Figure S3
sch_stab_ts <- comp_data %>%
  mutate(date = as.Date(date)) %>%
  filter(
    !is.na(date),
    !is.na(sch_stab)
  )

ggplot(sch_stab_ts, aes(x = date, y = sch_stab)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  labs(
    x = "Date",
    y = "Schmidt stability",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

###########
#Figure S4
set.seed(42)
X_explain_A <- model_data_raw[sample(nrow(model_data_raw),
                                     min(415, nrow(model_data_raw))),
                              xvars_clean]

shp_A <- shapviz(fit_raw,
                 X_pred = data.matrix(X_explain_A),
                 X      = X_explain_A)

#rename variables for figures
shap_labels <- c(
  "WindSpeed_Average_m_s"     = "Wind speed",
  "AirTemp_C_Average"         = "Air temperature",
  "thermo.depth"              = "Thermocline depth",
  "mean_Chla_lag1"            = "Chla lag1"
)

p1 <- sv_dependence(shp_A, v = "WindSpeed_Average_m_s", color_var = NULL) +
  labs(x = "Wind speed (m/s)", y = "SHAP value") 
p2 <- sv_dependence(shp_A, v = "AirTemp_C_Average", color_var = NULL) +
  labs(x = "Air temperature (C)", y = "SHAP value") 
p3 <- sv_dependence(shp_A, v = "mean_Chla_lag1", color_var = NULL) +
  labs(x = "Chla lag1", y = "SHAP value")
p4 <- sv_dependence(shp_A, v = "thermo.depth", color_var = NULL) +
  labs(x = "Thermocline depth", y = "SHAP value")

p1 + p2 + p3 + p4

p1 <- p1 + theme_classic()
p2 <- p2 + theme_classic()
p3 <- p3 + theme_classic()
p4 <- p4 + theme_classic()

shap_plots <-
  (p1 | p2) /
  (p3 | p4) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )


ggsave(
  filename = file.path(figure_dir,"Figure_SI5_SHAP.png"),shap_plots,width = 11,height = 12,units = "in",dpi = 600)


############
#Figure S5
plot_data <- comp_data %>%
  dplyr::select(date, mean_Chla, DOsat_surface, mean_fDOM) %>%
  mutate(
    date = as.Date(date)
  ) %>%
  pivot_longer(
    cols = c(mean_Chla, DOsat_surface, mean_fDOM),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("mean_Chla", "DOsat_surface", "mean_fDOM"),
      labels = c(
        "Chlorophyll-a (µg/L)",
        "DO surface saturation (%)",
        "Dissolved Organic Matter (RFU)"
      )
    )
  )

ggplot(plot_data, aes(x = date, y = value)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    x = "Date",
    y = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.8, "lines")
  )


###############
#Figure S6
#DOY distribution of fluxes
doy_data <- comp_data %>%
  mutate(
    date = as.Date(date),
    year = format(date, "%Y"),
    DOY = as.integer(format(date, "%j")),
    status = case_when(
      daily_CH4 > 0 ~ "Source",
      daily_CH4 < 0 ~ "Sink",
      daily_CH4 == 0 ~ "Zero",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(date),
    !is.na(daily_CH4)
  )


flux_distribution <- ggplot(doy_data, aes(y = year, x = DOY, color = status)) +
  
  geom_jitter(
    height = 0.15,
    width = 0,
    size = 2,
    alpha = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "Sink" = "#542788",
      "Zero" = "#F7F7F7",
      "Source" = "#E08214"
    ),
    name = "CH4 status"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 365, 30),
    limits = c(0, 365)
  ) +
  
  labs(
    x = "Day of year (DOY)",
    y = "Year"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "top",
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

ggsave(
  filename = file.path(figure_dir,"Flux distribution across years.png"), plot = flux_distribution, width = 6, height = 5, dpi = 300)

############
#Figure SI8: Time-series of CH4 fluxes, DO and temperature across different depths
x_limits <- as.Date(c("2020-01-01", "2026-12-31"))

#Plot appearance
pt_size   <- 0.8
flux_size <- 1.3
bg_alpha  <- 0.5


#Define the DO and temperature range
# Pool all water variables across depths.
# The primary axis is forced to start at 0.

water_vars <- c(
  comp_data$DOmgL_surface,
  comp_data$DOmgL_5m,
  comp_data$DOmgL_9m,
  comp_data$temp_1.6m,
  comp_data$temp_5m,
  comp_data$temp_9m
)

water_min <- 0
water_max <- max(water_vars, na.rm = TRUE)


#Define the CH4 flux range
# CH4 flux is allowed to contain both negative and positive values.

ch4_min <- min(comp_data$daily_CH4, na.rm = TRUE)
ch4_max <- max(comp_data$daily_CH4, na.rm = TRUE)

# Add a small amount of padding to the CH4 range
ch4_range <- ch4_max - ch4_min

ch4_min_pad <- ch4_min - 0.05 * ch4_range
ch4_max_pad <- ch4_max + 0.05 * ch4_range


#Scale and offset CH4 onto the primary axis
# This transformation allows the primary axis to start at 0
# while preserving negative CH4 values on the secondary axis.

scale_factor <- water_max / (ch4_max_pad - ch4_min_pad)

ch4_offset <- -ch4_min_pad * scale_factor


# Create scaled CH4 variable for plotting
comp_data <- comp_data |>
  mutate(
    CH4_scaled = daily_CH4 * scale_factor + ch4_offset
  )


#Define shared y-axis limits
# Start DO/temperature axis at exactly zero.
# Add headroom at the top for depth labels.

y_limits <- c(
  0,
  water_max * 1.15
)


#Common theme

base_theme <- theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


#Function to add depth labels

add_depth_label <- function(p, tag_label) {
  
  p +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = tag_label,
      hjust = -0.1,
      vjust = 1.5,
      size = 3.8,
      fontface = "italic"
    )
}


#Panel 1: 1.6 m

p1 <- ggplot(comp_data, aes(x = date)) +
  
  geom_point(
    aes(y = DOmgL_surface, color = "Dissolved Oxygen"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = temp_1.6m, color = "Temperature"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = CH4_scaled, color = "Atmospheric CH4 Flux"),
    size = flux_size,
    shape = 21,
    fill = "black",
    stroke = 0.4
  ) +
  
  scale_color_manual(
    values = c(
      "Atmospheric CH4 Flux" = "black",
      "Dissolved Oxygen" = "steelblue",
      "Temperature" = "red"
    )
  ) +
  scale_x_date(
    limits = x_limits,
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = NULL,
    limits = y_limits,
    
    sec.axis = sec_axis(
      ~ (. - ch4_offset) / scale_factor,
      name = NULL
    )
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "",
    color = ""
  ) +
  
  base_theme

p1 <- add_depth_label(p1, "(a) 1.6 m")


#Panel 2: 5 m

p2 <- ggplot(comp_data, aes(x = date)) +
  
  geom_point(
    aes(y = DOmgL_5m, color = "Dissolved Oxygen"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = temp_5m, color = "Temperature"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = CH4_scaled, color = "Atmospheric CH4 Flux"),
    size = flux_size,
    shape = 21,
    fill = "black",
    stroke = 0.4
  ) +
  
  scale_color_manual(
    values = c(
      "Atmospheric CH4 Flux" = "black",
      "Dissolved Oxygen" = "steelblue",
      "Temperature" = "red"
    )
  ) +
  scale_x_date(
    limits = x_limits,
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = "Dissolved Oxygen (mg/L) & Temperature (°C)",
    limits = y_limits,
    
    sec.axis = sec_axis(
      ~ (. - ch4_offset) / scale_factor,
      name = "Atmospheric CH4 Flux (mg/m²/d)"
    )
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "",
    color = ""
  ) +
  
  base_theme +
  
  theme(
    axis.title.y = element_text(
      margin = margin(r = 13)
    )
  )

p2 <- add_depth_label(p2, "(b) 5 m")


#Panel 3: 9 m

p3 <- ggplot(comp_data, aes(x = date)) +
  
  geom_point(
    aes(y = DOmgL_9m, color = "Dissolved Oxygen"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = temp_9m, color = "Temperature"),
    size = pt_size,
    alpha = bg_alpha
  ) +
  
  geom_point(
    aes(y = CH4_scaled, color = "Atmospheric CH4 Flux"),
    size = flux_size,
    shape = 21,
    fill = "black",
    stroke = 0.4
  ) +
  
  scale_color_manual(
    values = c(
      "Atmospheric CH4 Flux" = "black",
      "Dissolved Oxygen" = "steelblue",
      "Temperature" = "red"
    )
  ) +
  
  scale_x_date(
    limits = x_limits,
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_continuous(
    name = NULL,
    limits = y_limits,
    
    sec.axis = sec_axis(
      ~ (. - ch4_offset) / scale_factor,
      name = NULL
    )
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "",
    color = "",
    x = "Date"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "top"
  )

p3 <- add_depth_label(p3, "(c) 9 m")


#Combine panels

combined <- (p1 / p2 / p3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )


# Display the combined figure
combined

ggsave(
  filename = file.path(figure_dir,"Flux_timeseries.png"), plot = combined, width = 8, height = 6, dpi = 300)
###########










