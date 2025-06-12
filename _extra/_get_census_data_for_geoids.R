#
# - get census data for geoids
# https://api.census.gov/data/create_success.html

library(tidycensus)
library(dplyr)
library(readr)
library(stringr)
library(here)
library(purrr)

#---------------------------
# Step 1: Set Census API Key
#---------------------------
tryCatch({
  census_api_key("81107bf14154f7b68ed09fe6c5289cb54be4c5f5", install = TRUE, overwrite = TRUE)
  readRenviron("~/.Renviron")
}, error = function(e) {
  stop("❌ [Line ~10] Error setting up Census API key: ", e$message)
})

#---------------------------
# Step 2: Read GEOID list
#---------------------------
geoids <- tryCatch({
  lines <- read_lines(here("_extra", "geo_ids_present_in_both_22_23.txt"))
  trimmed <- str_trim(lines)
  geoids_clean <- trimmed[trimmed != ""]  # remove empty strings
  geoids_clean
}, error = function(e) {
  stop("❌ [Line ~18] Failed to read GEOID list from file: ", e$message)
})

if (length(geoids) == 0) stop("⚠️ No valid GEOIDs found in input file.")

#---------------------------
# Step 3: Define ACS Variables
#---------------------------
acs_vars <- c(
  total_pop       = "B01003_001",
  median_income   = "B19013_001",
  white           = "B02001_002",
  black           = "B02001_003",
  native_american = "B02001_004",
  asian           = "B02001_005",
  pacific_islander = "B02001_006",
  other_race      = "B02001_007",
  two_or_more     = "B02001_008",
  hispanic        = "B03003_003",
  age_median      = "B01002_001"
)
#---------------------------
# Step 4: Fetch Function
#---------------------------
fetch_acs_data <- function(year) {
  tryCatch({
    message(glue::glue("🔄 Fetching ACS data for {year}..."))
    get_acs(
      geography = "county",
      variables = acs_vars,
      year = year,
      survey = "acs5",
      output = "wide"
    ) %>% filter(GEOID %in% geoids)
  }, error = function(e) {
    stop(glue::glue("❌ [Line ~44] Failed to fetch ACS data for {year}: {e$message}"))
  })
}

#---------------------------
# Step 5: Fetch Data
#---------------------------
acs_2022 <- fetch_acs_data(2022)
acs_2023 <- fetch_acs_data(2023)

#---------------------------
# Step 5.5: Rename Columns to Human-Readable Names
#---------------------------
rename_acs_columns <- function(df) {
  # Human-readable labels for estimate columns (E)
  est_names <- names(acs_vars)
  est_codes <- paste0(acs_vars, "E")
  est_map <- setNames(est_names, est_codes)
  
  # Human-readable labels for margin of error columns (M)
  moe_map <- c(
    B01003_001M = "Margin of error for total population",
    B19013_001M = "Margin of error for median income",
    B02001_002M = "Margin of error for white population",
    B02001_003M = "Margin of error for black population",
    B02001_004M = "Margin of error for native population",
    B02001_005M = "Margin of error for asian population",
    B02001_006M = "Margin of error for pacific islander population",
    B02001_007M = "Margin of error for other race",
    B02001_008M = "Margin of error for two or more races",
    B01002_001M = "Margin of error for median age"
  )
  
  # Always include GEOID and NAME labels too
  geo_map <- c(
    GEOID = "Geographic identifier (e.g., FIPS code for county)",
    NAME = "Name of the county/region"
  )
  
  rename_vec <- c(geo_map, est_map, moe_map)
  df %>% rename(any_of(rename_vec))
}

acs_2022 <- rename_acs_columns(acs_2022)
acs_2023 <- rename_acs_columns(acs_2023)
#---------------------------
# Step 6: Save Output
#---------------------------
tryCatch({
  if (!dir.exists("_extra")) dir.create("_extra")
  
  write_csv(acs_2022, "_extra/acs_demographics_2022.csv")
  write_csv(acs_2023, "_extra/acs_demographics_2023.csv")
  
  cat("✅ ACS data saved:\n")
  cat("   - _extra/acs_demographics_2022.csv\n")
  cat("   - _extra/acs_demographics_2023.csv\n")
}, error = function(e) {
  stop("❌ [Line ~65] Failed to write CSV files: ", e$message)
})

# User alert
message("✅ Unique geoids have been written to text files.")
beep(sound = 1)  # Plays default sound (you can change the number for other sounds)
