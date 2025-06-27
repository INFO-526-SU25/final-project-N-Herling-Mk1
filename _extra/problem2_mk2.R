library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(tmap)

options(tigris_use_cache = TRUE)

# Load and clean water insecurity data (2023 only)
water_2023_clean <- read.csv("Water_insecurity/water_insecurity_2023.csv") %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

# Load counties shapefile
counties_sf <- counties(cb = TRUE, year = 2020, class = "sf")

# Join 2023 data to counties
spatial_data <- counties_sf %>%
  left_join(water_2023_clean %>% select(geoid, percent_lacking_plumbing),
            by = c("GEOID" = "geoid")) %>%
  filter(!is.na(percent_lacking_plumbing))  # Remove counties with missing plumbing data

# Project to Albers Equal Area for spatial analysis
spatial_data <- st_transform(spatial_data, 5070)

# Create neighbors list (Queen contiguity with snap to reduce isolation)
nb <- poly2nb(spatial_data, queen = TRUE, snap = 5000)

# Remove polygons with no neighbors
no_neighbors <- which(card(nb) == 0)
if (length(no_neighbors) > 0) {
  spatial_data <- spatial_data[-no_neighbors, ]
  nb <- poly2nb(spatial_data, queen = TRUE, snap = 5000)
}

# Create spatial weights
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Compute Local G* statistic on plumbing levels
local_gstar <- localG(spatial_data$percent_lacking_plumbing, listw = lw, zero.policy = TRUE)
spatial_data$gstar_z <- as.numeric(local_gstar)

# Classify clusters
spatial_data <- spatial_data %>%
  mutate(cluster_type = case_when(
    gstar_z >= 1.96  ~ "Hotspot (Low Plumbing Access)",
    gstar_z <= -1.96 ~ "Coldspot (High Plumbing Access)",
    TRUE             ~ "Not Significant"
  ))

# Plot
print(
tmap_mode("plot") +
tm_shape(spatial_data) +
  tm_polygons("cluster_type",
              palette = c("#FDE725", "#21908C", "#BFBFBF"),  # Hotspot, Coldspot, Not Significant,
              fill.legend = tm_legend(title = "Plumbing Access Clusters (2023)")) +
  tm_layout(legend.outside = TRUE)
)