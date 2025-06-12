#---------------------------
# Save Column Descriptions to CSV
#---------------------------
column_descriptions <- tibble::tibble(
  column_name = c(
    "GEOID", "NAME",
    "total_popE", "total_popM",
    "median_incomeE", "median_incomeM",
    "whiteE", "whiteM",
    "blackE", "blackM",
    "native_americanE", "native_americanM",
    "asianE", "asianM",
    "pacific_islanderE", "pacific_islanderM",
    "other_raceE", "other_raceM",
    "two_or_moreE", "two_or_moreM",
    "hispanicE", "hispanicM",
    "age_medianE", "age_medianM"
  ),
  description = c(
    "Geographic identifier (e.g., FIPS code for county)",
    "Name of the county/region",
    "Estimated total population",
    "Margin of error for total population",
    "Median household income",
    "Margin of error for median income",
    "White alone population",
    "Margin of error for white population",
    "Black or African American alone population",
    "Margin of error for black population",
    "American Indian and Alaska Native alone population",
    "Margin of error for native population",
    "Asian alone population",
    "Margin of error for asian population",
    "Native Hawaiian or Pacific Islander alone population",
    "Margin of error for pacific islander population",
    "Some other race alone population",
    "Margin of error for other race",
    "Two or more races population",
    "Margin of error for two or more races",
    "Hispanic or Latino (any race) population",
    "Margin of error for Hispanic population",
    "Median age",
    "Margin of error for median age"
  )
)

# Save to CSV
readr::write_csv(column_descriptions, "_extra/acs_column_descriptions.csv")
