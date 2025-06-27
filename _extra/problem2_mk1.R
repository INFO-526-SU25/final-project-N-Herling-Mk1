# Load libraries
library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(tmap)

# Cache tigris shapefiles
options(tigris_use_cache = TRUE)

# ────────────────────────────────────────────────
# Load and clean water insecurity data
water_2022_clean <- read.csv("Water_insecurity/water_insecurity_2022.csv") %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

water_2023_clean <- read.csv("Water_insecurity/water_insecurity_2023.csv") %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

# Filter to common GEOIDs
common_geoids <- intersect(water_2022_clean$geoid, water_2023_clean$geoid)
water_2022_filtered <- water_2022_clean %>% filter(geoid %in% common_geoids)
water_2023_filtered <- water_2023_clean %>% filter(geoid %in% common_geoids)

# ────────────────────────────────────────────────
# Load and join shapefile
counties_sf <- counties(cb = TRUE, year = 2020, class = "sf")

counties_2022 <- counties_sf %>%
  left_join(water_2022_filtered %>%
              select(geoid, percent_lacking_plumbing_2022 = percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

counties_2023 <- counties_sf %>%
  left_join(water_2023_filtered %>%
              select(geoid, percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

# ────────────────────────────────────────────────
# Merge and calculate change
spatial_data <- counties_2023 %>%
  left_join(counties_2022 %>% st_drop_geometry() %>%
              select(GEOID, percent_lacking_plumbing_2022),
            by = "GEOID") %>%
  mutate(change = percent_lacking_plumbing - percent_lacking_plumbing_2022) %>%
  filter(!is.na(change))  # Remove rows with NA in change

# ────────────────────────────────────────────────
# Project to Albers Equal Area CRS
spatial_data <- st_transform(spatial_data, crs = 5070)

# ────────────────────────────────────────────────
# Create neighbors list (with Queen contiguity and snap distance)
nb <- poly2nb(spatial_data, queen = TRUE, snap = 5000)

# Identify and exclude observations with no neighbors
no_neighbors <- which(card(nb) == 0)
cat("Number of observations with no neighbors (excluded):", length(no_neighbors), "\n")

# Remove polygons with no neighbors from spatial data and neighbor list
if(length(no_neighbors) > 0) {
  spatial_data <- spatial_data[-no_neighbors, ]
  nb <- poly2nb(spatial_data, queen = TRUE, snap = 5000)  # Recalculate after filtering
}

# Now continue safely with weights and analysis
lw <- nb2listw(nb, style = "W", zero.policy = FALSE)

# Compute Local G* (all units now have neighbors)
local_gstar <- localG(spatial_data$change, listw = lw, zero.policy = FALSE)
spatial_data$gstar_z <- as.numeric(local_gstar)

# ────────────────────────────────────────────────
# Create spatial weights list
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Calculate Local G* statistic
local_gstar <- localG(spatial_data$change, listw = lw, zero.policy = TRUE)
spatial_data$gstar_z <- as.numeric(local_gstar)

# ────────────────────────────────────────────────
# Classify spatial clusters
spatial_data <- spatial_data %>%
  mutate(cluster_type = case_when(
    gstar_z >= 1.96  ~ "Hotspot (Worsening)",
    gstar_z <= -1.96 ~ "Coldspot (Improving)",
    TRUE             ~ "Not Significant"
  ))

# ────────────────────────────────────────────────
# Visualization
print(
  tm_shape(spatial_data) +
    tm_polygons("cluster_type",
                palette = c("red", "blue", "grey80"),
                fill.legend = tm_legend(title = "Plumbing Access Change Clusters")) +
    tm_layout(legend.outside = TRUE)
)
