#
# - gets geoid information from .csv files.
#
library(readr)
library(here)
library(dplyr)
library(beepr)  # Load the beepr package

# Read the CSV files
water_insecurity_2022 <- read_csv(here("water_insecurity", "water_insecurity_2022.csv"))
water_insecurity_2023 <- read_csv(here("water_insecurity", "water_insecurity_2023.csv"))

# Extract unique geoids
unique_geoids_2022 <- unique(water_insecurity_2022$geoid)
unique_geoids_2023 <- unique(water_insecurity_2023$geoid)

# Write to .txt files
writeLines(unique_geoids_2022, here("_extra", "unique_geoids_2022.txt"))
writeLines(unique_geoids_2023, here("_extra", "unique_geoids_2023.txt"))

# Find matching geoids
matching_geoids <- intersect(unique_geoids_2022, unique_geoids_2023)

# Write matching geoids to file
writeLines(matching_geoids, here("_extra", "geo_ids_present_in_both_22_23.txt"))

# User alert
message("✅ Unique geoids have been written to text files.")
beep(sound = 1)  # Plays default sound (you can change the number for other sounds)
