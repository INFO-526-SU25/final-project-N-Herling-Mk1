library(tidyverse)
library(usmap)
library(here)
library(stringr)

# Read your data (adjust paths as needed)
water_insecurity_2022 <- read_csv(here("water_insecurity", "water_insecurity_2022.csv"))
water_insecurity_2023 <- read_csv(here("water_insecurity", "water_insecurity_2023.csv"))

# Helper function to process a dataset for states/counties stats
process_year_data <- function(df) {
  df %>%
    distinct(name) %>%
    separate(name, into = c("county", "state"), sep = ",") %>%
    mutate(
      county = str_trim(str_replace(county, " County$", "")),
      state = str_trim(state)
    ) %>%
    distinct(state, county) %>%
    count(state, name = "represented_counties") %>%
    mutate(
      total_counties = sum(represented_counties),
      `% of total` = round(100 * represented_counties / total_counties, 3)
    ) %>%
    select(state, represented_counties, `% of total`)
}

# Process 2022 and 2023 separately
stats_2022 <- process_year_data(water_insecurity_2022)
stats_2023 <- process_year_data(water_insecurity_2023)

# Add total rows for each year
total_2022 <- tibble(
  state = "TOTAL",
  represented_counties = sum(stats_2022$represented_counties),
  `% of total` = sum(stats_2022$`% of total`)
)
stats_2022 <- bind_rows(stats_2022, total_2022)

total_2023 <- tibble(
  state = "TOTAL",
  represented_counties = sum(stats_2023$represented_counties),
  `% of total` = sum(stats_2023$`% of total`)
)
stats_2023 <- bind_rows(stats_2023, total_2023)

# Rename columns except 'state' for joining side-by-side
stats_2022_renamed <- stats_2022 %>%
  rename_with(~ paste0(.x, "_2022"), -state)

stats_2023_renamed <- stats_2023 %>%
  rename_with(~ paste0(.x, "_2023"), -state)

# Join on state (which remains the same column name in both)
combined_stats <- full_join(stats_2022_renamed, stats_2023_renamed, by = "state")

# Add difference flag column for represented_counties comparison, ignoring last row
combined_stats <- combined_stats %>%
  mutate(
    difference_flag = ""
  )

# For all but the last row, set '*' if represented counties differ
if (nrow(combined_stats) > 1) {
  combined_stats$difference_flag[1:(nrow(combined_stats) - 1)] <- ifelse(
    combined_stats$represented_counties_2022[1:(nrow(combined_stats) - 1)] != 
      combined_stats$represented_counties_2023[1:(nrow(combined_stats) - 1)],
    "*", ""
  )
}

# Print combined stats to console
cat("State and County Representation Stats by Year\n")
cat("--------------------------------------------\n")
print(combined_stats)

# Ensure _extra directory exists
if (!dir.exists("_extra")) dir.create("_extra")

# Write to CSV
write_csv(combined_stats, file = "_extra/state_count_stats_by_year.csv")

cat("\n✅ State and county representation stats saved to '_extra/state_count_stats_by_year.csv'\n")
