library(tidyverse)
library(sf)
library(tigris)        # For county shapefiles
library(spdep)
library(tmap)

options(tigris_use_cache = TRUE)

# Load and clean water insecurity data
water_2022_clean <- read.csv("Water_insecurity/water_insecurity_2022.csv") %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

water_2023_clean <- read.csv("Water_insecurity/water_insecurity_2023.csv") %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

# Find common GEOIDs
common_geoids <- intersect(water_2022_clean$geoid, water_2023_clean$geoid)

water_2022_filtered <- water_2022_clean %>% filter(geoid %in% common_geoids)
water_2023_filtered <- water_2023_clean %>% filter(geoid %in% common_geoids)

# Load US counties shapefile (2020 Census counties)
counties_sf <- counties(cb = TRUE, year = 2020, class = "sf")

# Join water data to counties (2022 and 2023)
counties_2022 <- counties_sf %>%
  left_join(water_2022_filtered %>% select(geoid, percent_lacking_plumbing_2022 = percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

counties_2023 <- counties_sf %>%
  left_join(water_2023_filtered %>% select(geoid, percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

# Combine into one spatial dataframe
spatial_data <- counties_2023 %>%
  left_join(counties_2022 %>% st_drop_geometry() %>% select(GEOID, percent_lacking_plumbing_2022),
            by = "GEOID") %>%
  mutate(change = percent_lacking_plumbing - percent_lacking_plumbing_2022)

# Reproject to Albers Equal Area (meters) for neighbor calculations
spatial_data <- st_transform(spatial_data, 5070)

# Remove rows with NA in 'change' before neighbors calculation
spatial_data <- spatial_data %>% filter(!is.na(change))

# Buffer polygons slightly to fix minor topology issues (0.1 meter buffer)
spatial_data_buff <- st_buffer(spatial_data, dist = 0.1)

# Create neighbors list with larger snap distance (3000 meters)
nb <- poly2nb(spatial_data_buff, queen = TRUE, snap = 3000)

# Identify isolated polygons with no neighbors
isolated <- which(card(nb) == 0)
cat("Number of isolated polygons with no neighbors:", length(isolated), "\n")

# Remove isolated polygons to avoid issues in spatial stats
if(length(isolated) > 0){
  spatial_data <- spatial_data[-isolated, ]
  spatial_data_buff <- spatial_data_buff[-isolated, ]
  nb <- poly2nb(spatial_data_buff, queen = TRUE, snap = 3000)
}

# Create spatial weights list
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Calculate Local G* statistic (Getis-Ord)
local_gstar <- localG(spatial_data$change, listw = lw, zero.policy = TRUE)
spatial_data$gstar_z <- as.numeric(local_gstar)

# Classify cluster types based on significance thresholds
spatial_data <- spatial_data %>%
  mutate(cluster_type = case_when(
    gstar_z >= 1.96  ~ "Hotspot (Worsening)",
    gstar_z <= -1.96 ~ "Coldspot (Improving)",
    TRUE             ~ "Not Significant"
  ))

# Visualize results with tmap
tmap_mode("plot")
tm_shape(spatial_data) +
  tm_polygons("cluster_type",
              palette = c("red", "blue", "grey80"),
              title = "Plumbing Access Change Clusters") +
  tm_layout(legend.outside = TRUE)
