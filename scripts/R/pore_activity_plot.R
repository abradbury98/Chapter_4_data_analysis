library(ggplot2)
library(dplyr)

# ── Data ──────────────────────────────────────────────────────────────────────
df <- read.csv("pore_activity_data.csv", stringsAsFactors = FALSE)

# Ordered factor for stacking (bottom to top)
state_levels <- c("Unclassified", "Inactive", "Unavailable", "Pore available", "Sequencing")
df$state <- factor(df$state, levels = state_levels)

# MinKNOW colours (from HTML source)
state_colours <- c(
  "Sequencing"     = "#00ff00",
  "Pore available" = "#00cc00",
  "Unavailable"    = "#0084a9",
  "Inactive"       = "#90c6e7",
  "Unclassified"   = "#315BAD"
)

# Display order: EVI, ASC, EVI_BFBB, ASC_BFBB, EVI_RVS, ASC_RVS (row-by-row, left=EVI, right=ASC)
run_levels <- c("EVI", "ASC", "EVI_BFBB", "ASC_BFBB", "EVI_RVS", "ASC_RVS")
run_labels <- c(
  "EVI"      = "Post-EVI",
  "ASC"      = "Post-ASC",
  "EVI_BFBB" = "Post-EVI BF Bolton base + CAT",
  "ASC_BFBB" = "Post-ASC BF Bolton base + CAT",
  "EVI_RVS"  = "Post-EVI BPW + RVS",
  "ASC_RVS"  = "Post-ASC BPW + RVS"
)
df$run <- factor(df$run, levels = run_levels)

# ── Compute bar centres and widths for geom_col ───────────────────────────────
bar_gap <- 0.15   # 15% gap between bars

df <- df %>%
  mutate(
    x_center  = (xmin_h + xmax_h) / 2,
    bar_width = (xmax_h - xmin_h) * (1 - bar_gap)
  )

# ── Theme colours ─────────────────────────────────────────────────────────────
bg_col        <- "white"
panel_col     <- "white"
text_col      <- "#333333"
grid_col      <- "#dddddd"
strip_bg      <- "#f0f0f0"
strip_text    <- "#111111"

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = x_center, y = pct, fill = state, width = bar_width)) +
  geom_col(position = "stack", colour = NA, linewidth = 0) +
  facet_wrap(~ run, ncol = 2, labeller = labeller(run = run_labels), axes = "all") +
  scale_fill_manual(
    values = state_colours,
    breaks = rev(state_levels),
    labels = rev(state_levels),
    name   = "Pore state"
  ) +
  scale_x_continuous(
    name   = "Run time (hours)",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "Percentage of pores (%)",
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal(base_size = 11) +
  theme(
    # Overall background
    plot.background   = element_rect(fill = bg_col, colour = NA),
    panel.background  = element_rect(fill = panel_col, colour = NA),
    # Grid
    panel.grid.major  = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    # Axes
    axis.text         = element_text(colour = text_col, size = 9),
    axis.title        = element_text(colour = text_col, size = 10),
    axis.ticks        = element_line(colour = grid_col),
    axis.line         = element_line(colour = grid_col, linewidth = 0.3),
    # Facet strips
    strip.background  = element_rect(fill = strip_bg, colour = NA),
    strip.text        = element_text(colour = strip_text, face = "bold", size = 10),
    # Legend
    legend.background = element_rect(fill = bg_col, colour = NA),
    legend.key        = element_rect(fill = bg_col, colour = NA),
    legend.text       = element_text(colour = text_col, size = 10),
    legend.title      = element_text(colour = text_col, face = "bold", size = 10),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.5, "cm"),
    # Title
    plot.title        = element_text(colour = strip_text, face = "bold", size = 13, hjust = 0.5),
    plot.subtitle     = element_text(colour = text_col, size = 10, hjust = 0.5),
    # Spacing
    panel.spacing     = unit(0.8, "lines"),
    plot.margin       = margin(12, 12, 12, 12)
  ) +
  labs(
    title    = "Nanopore pore activity over time",
    subtitle = "Percentage of pores in each state at each pore scan"
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(
  "pore_activity_all_runs.png",
  plot   = p,
  width  = 12,
  height = 10,
  dpi    = 300,
  bg     = bg_col
)

message("Saved: pore_activity_all_runs.png")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 – Pore scan (absolute pore counts, 1.5-hour scan interval)
# ══════════════════════════════════════════════════════════════════════════════

ps <- read.csv("pore_scan_data.csv", stringsAsFactors = FALSE)

# Stacking order bottom to top (matches MinKNOW display)
ps_state_levels <- c("Inactive", "Zero", "Saturated", "Unavailable", "Reserved pore", "Pore available")
ps$state <- factor(ps$state, levels = ps_state_levels)

ps_colours <- c(
  "Pore available" = "#00CC00",
  "Reserved pore"  = "#EDE797",
  "Unavailable"    = "#54B8B1",
  "Saturated"      = "#68767E",
  "Zero"           = "#90C6E7",
  "Inactive"       = "#0084A9"
)

ps$run <- factor(ps$run, levels = run_levels)

ps <- ps %>%
  mutate(
    x_center  = (xmin_h + xmax_h) / 2,
    bar_width = (xmax_h - xmin_h) * (1 - bar_gap)
  )

total_pores <- unique(ps$total)   # 2048 for all runs

p2 <- ggplot(ps, aes(x = x_center, y = count, fill = state, width = bar_width)) +
  geom_col(position = "stack", colour = NA, linewidth = 0) +
  facet_wrap(~ run, ncol = 2, labeller = labeller(run = run_labels), axes = "all") +
  scale_fill_manual(
    values = ps_colours,
    breaks = rev(ps_state_levels),
    labels = rev(ps_state_levels),
    name   = "Pore state"
  ) +
  scale_x_continuous(
    name   = "Run time (hours)",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "Number of pores",
    expand = c(0, 0),
    limits = c(0, total_pores)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = bg_col, colour = NA),
    panel.background  = element_rect(fill = panel_col, colour = NA),
    panel.grid.major  = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(colour = text_col, size = 9),
    axis.title        = element_text(colour = text_col, size = 10),
    axis.ticks        = element_line(colour = grid_col),
    axis.line         = element_line(colour = grid_col, linewidth = 0.3),
    strip.background  = element_rect(fill = strip_bg, colour = NA),
    strip.text        = element_text(colour = strip_text, face = "bold", size = 10),
    legend.background = element_rect(fill = bg_col, colour = NA),
    legend.key        = element_rect(fill = bg_col, colour = NA),
    legend.text       = element_text(colour = text_col, size = 10),
    legend.title      = element_text(colour = text_col, face = "bold", size = 10),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.5, "cm"),
    plot.title        = element_text(colour = strip_text, face = "bold", size = 13, hjust = 0.5),
    plot.subtitle     = element_text(colour = text_col, size = 10, hjust = 0.5),
    panel.spacing     = unit(0.8, "lines"),
    plot.margin       = margin(12, 12, 12, 12)
  ) +
  labs(
    title    = "Nanopore pore scan over time",
    subtitle = "Number of pores in each state at each 1.5-hour scan"
  )

ggsave(
  "pore_scan_all_runs.png",
  plot   = p2,
  width  = 12,
  height = 10,
  dpi    = 300,
  bg     = bg_col
)

message("Saved: pore_scan_all_runs.png")
