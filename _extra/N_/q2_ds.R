# Define U.S. states grouped by Census regions and divisions
us_regions <- list(
  Northeast = list(
    New_England = c("Maine", "New Hampshire", "Vermont", 
                    "Massachusetts", "Rhode Island", "Connecticut"),
    Mid_Atlantic = c("New York", "New Jersey", "Pennsylvania")
  ),
  
  Midwest = list(
    East_North_Central = c("Ohio", "Indiana", "Illinois", "Michigan", "Wisconsin"),
    West_North_Central = c("Minnesota", "Iowa", "Missouri", 
                           "North Dakota", "South Dakota", "Nebraska", "Kansas")
  ),
  
  South = list(
    South_Atlantic = c("Delaware", "Maryland", "District of Columbia", 
                       "Virginia", "West Virginia", "North Carolina", 
                       "South Carolina", "Georgia", "Florida"),
    East_South_Central = c("Kentucky", "Tennessee", "Mississippi", "Alabama"),
    West_South_Central = c("Arkansas", "Louisiana", "Oklahoma", "Texas")
  ),
  
  West = list(
    Mountain = c("Montana", "Idaho", "Wyoming", "Colorado", 
                 "Utah", "Nevada", "Arizona", "New Mexico"),
    Pacific = c("Washington", "Oregon", "California", "Alaska", "Hawaii")
  )
)

# Flatten the nested list to get all state-like entries
all_states <- unlist(us_regions)

# Count all entries (includes District of Columbia)
total_with_dc <- length(all_states)

# Exclude DC to get only the 50 official states
all_states_no_dc <- all_states[all_states != "District of Columbia"]
total_states_only <- length(all_states_no_dc)

# Identify non-state entities (territories or federal districts)
non_states <- unique(all_states[!(all_states %in% state.name)])

# Print summary
cat("Total entries (states + DC):", total_with_dc, "\n")
cat("Total U.S. states only:", total_states_only, "\n")
cat("Non-state entities present:", paste(non_states, collapse = ", "), "\n")
