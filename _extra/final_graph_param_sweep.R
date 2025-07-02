# ────────────────────────────────────────────────
# Load libraries
library(tidyverse)
library(here)
library(patchwork)
library(viridis)    # for turbo colors
library(ggrepel)    # for repel annotations

# ────────────────────────────────────────────────
# Load hyperparameter sweep results
results_tbl <- read_csv(here("_extra", "hyper_param_sweep.csv"), show_col_types = FALSE)

# ────────────────────────────────────────────────
# Shared plot theme
base_theme <- theme_minimal(base_size = 13)

# ────────────────────────────────────────────────
# Set cutoff for significance
z_cutoff <- 1.96

# ────────────────────────────────────────────────
# Define shading colors from viridis turbo palette
shade_hot_color <- viridis(10, option = "turbo")[7]  # bright yellow-ish
shade_cold_color <- viridis(10, option = "turbo")[4] # teal-ish

# ────────────────────────────────────────────────
# Hotspot plot (positive z)
p_hotspot <- ggplot(
  results_tbl %>% filter(!is.na(`Hotspot (Worsening)`)),
  aes(x = z, y = `Hotspot (Worsening)`, color = factor(snap))
) +
  # Shade right of cutoff line (zone of significance)
  annotate(
    "rect",
    xmin = z_cutoff, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = shade_hot_color,
    alpha = 0.7
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(~ queen, labeller = label_both) +
  geom_vline(xintercept = z_cutoff, linetype = "dashed", color = "red") +
  # Repel annotation for vertical line
  geom_text_repel(
    data = tibble(
      x = z_cutoff,
      y = max(results_tbl$`Hotspot (Worsening)`, na.rm = TRUE) * 0.3,
      label = paste0("z = ", z_cutoff)
    ),
    aes(x = x, y = y, label = label),
    color = "black",
    nudge_x = -0.32,
    nudge_y = -0.4 * max(results_tbl$`Hotspot (Worsening)`, na.rm = TRUE),
    direction = "y",
    hjust = 1,
    segment.size = 0.3,
    segment.color = "grey50",
    size = 4
  ) +
  annotate(
    "text",
    x = z_cutoff + 0.2,
    y = max(results_tbl$`Hotspot (Worsening)`, na.rm = TRUE) * 0.8,
    label = "Statistical Significance\nHotspots (p ≤ 0.05)",
    color = "black",
    hjust = 0.2,
    size = 4
  ) +
  labs(
    title = "Hotspot Counts vs Positive Z-Threshold",
    x = "Z-Threshold",
    y = "Hotspot Count",
    color = "Snap Distance"
  ) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "lines"),  # increased horizontal spacing
    strip.text = element_text(size = 12, face = "bold"),  # Larger & bold facet titles
    plot.title = element_text(hjust = 0.5)
    )

# ────────────────────────────────────────────────
# Coldspot plot (negative z)
coldspot_data <- results_tbl %>%
  filter(!is.na(`Coldspot (Improving)`)) %>%
  mutate(z_neg = -z)

p_coldspot <- ggplot(
  coldspot_data,
  aes(x = z_neg, y = `Coldspot (Improving)`, color = factor(snap))
) +
  # Shade left of cutoff line (zone of significance)
  annotate(
    "rect",
    xmin = -Inf, xmax = -z_cutoff,
    ymin = -Inf, ymax = Inf,
    fill = shade_cold_color,
    alpha = 0.7
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(~ queen, labeller = label_both) +
  geom_vline(xintercept = -z_cutoff, linetype = "dashed", color = "red") +
  # Repel annotation for vertical line
  geom_text_repel(
    data = tibble(
      x = -z_cutoff,
      y = max(coldspot_data$`Coldspot (Improving)`, na.rm = TRUE) * 0.3,
      label = paste0("z = -", z_cutoff)
    ),
    aes(x = x, y = y, label = label),
    color = "black",
    nudge_x = 0.3,
    nudge_y = -0.5 * max(coldspot_data$`Coldspot (Improving)`, na.rm = TRUE),
    direction = "y",
    hjust = 0,
    segment.size = 0.3,
    segment.color = "grey50",
    size = 4
  ) +
  annotate(
    "text",
    x = -z_cutoff - 0.2,
    y = max(coldspot_data$`Coldspot (Improving)`, na.rm = TRUE) * 0.85,
    label = "Statistical Significance\nColdspots (p ≤ 0.05)",
    color = "black",
    hjust = 0.9,
    size = 4
  ) +
  labs(
    title = "Coldspot Counts vs Negative Z-Threshold",
    x = "Z-Threshold (Negative for Coldspots)",
    y = "Coldspot Count",
    color = "Snap Distance"
  ) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "lines"),  # increased horizontal spacing
    strip.text = element_text(size = 12, face = "bold"),  # Larger & bold facet titles
    plot.title = element_text(hjust = 0.5)
  )
# ────────────────────────────────────────────────
# Combine plots vertically with equal heights
combined_plot <- p_hotspot / p_coldspot + plot_layout(heights = c(1, 1))

# ────────────────────────────────────────────────
# Save and print
ggsave(
  filename = here("_extra", "topo_lines_hotspot_coldspot_shaded_annotated.png"),
  plot = combined_plot,
  width = 8,
  height = 10
)

print(combined_plot)

cat("✅ Hotspot and Coldspot count visualizations with shading and repel annotations saved and rendered.\n")
