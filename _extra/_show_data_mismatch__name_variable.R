## 0 - Setup

#################
# Package Setup #
#################
#Check if pacman [package manager] is installed, if not install it.
#throw [FYI] alert either way.
if (!requireNamespace("pacman", quietly = TRUE)) {
  message("Installing 'pacman' (not found locally)...")
  install.packages("pacman")
} else {
  message("[FYI]\n'pacman' already installed — skipping install.")
}
# use this line for installing/loading
# pacman::p_load()
# - packages to load stored in a variable (vector)
pkgs <- c(
  "tidyverse",
  "glue",
  "scales",
  "lubridate",
  "patchwork",
  "ggh4x",
  "ggrepel",
  "openintro",
  "here"
)
# - load from the character array/vector
pacman::p_load(char = pkgs)

# - install tidyverse/dsbox directly from Git Hub
# - this allows for the possible need to install on a repo. pull.
# - and, if it's already installed just thorw an alert.
if (!requireNamespace("dsbox", quietly = TRUE)) {
  message("Installing 'dsbox' from GitHub (not found locally)...")
  suppressMessages(devtools::install_github("tidyverse/dsbox"))
} else {
  message("[FYI]\n'dsbox' already installed — skipping GitHub install.")
}

# - alert to user packages loaded.
cat(paste(
  "The packages loaded:",
  paste("-", pkgs, collapse = "\n"),
  sep = "\n"
))

#-------------------------->
######################
# Basic set Theme up #
######################
# set theme for ggplot2
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 14))

# set width of code output
options(width = 65)

# set figure parameters for knitr
knitr::opts_chunk$set(
  fig.width = 7,        # 7" width
  fig.asp = 0.618,      # the golden ratio
  fig.retina = 3,       # dpi multiplier for displaying HTML output on retina
  fig.align = "center", # center align figures
  dpi = 300             # higher dpi, sharper image
)
#######
water_insecurity_2022 <- read_csv(here("water_insecurity",
                                       "water_insecurity_2022.csv"))
water_insecurity_2023 <- read_csv(here("water_insecurity",
                                       "water_insecurity_2023.csv"))



# Extract unique names
unique_names_22 <- unique(water_insecurity_2022$name)
unique_names_23 <- unique(water_insecurity_2023$name)

# Find all unique names from both years combined
all_names <- union(unique_names_22, unique_names_23)

# Create a data frame to store results
compendium <- data.frame(
  `22_names` = character(length(all_names)),
  `23_names` = character(length(all_names)),
  stringsAsFactors = FALSE
)

for (i in seq_along(all_names)) {
  name <- all_names[i]
  in_22 <- name %in% unique_names_22
  in_23 <- name %in% unique_names_23
  
  compendium$`22_names`[i] <- if (in_22) name else "X"
  compendium$`23_names`[i] <- if (in_23) name else "X"
}

# Optional: sort so matching pairs come first and mismatches at bottom
compendium$match <- compendium$`22_names` == compendium$`23_names`
compendium <- compendium[order(-compendium$match), ]
compendium$match <- NULL

# Make sure the folder exists
if (!dir.exists("_extra")) {
  dir.create("_extra")
}

# Prepare file path
file_path <- file.path("_extra", "county_name_compendium.csv")

# Write header with two empty leading columns
header_line <- paste(c(",","22_names", "23_names"), collapse = ",")
writeLines(header_line, con = file_path)

# Prepare data with two empty columns at start
data_to_write <- cbind( compendium)

# Append data to file, no row names or column names, quotes around entries
write.table(data_to_write, 
            file = file_path, 
            sep = ",", 
            append = TRUE, 
            col.names = FALSE, 
            row.names = FALSE,
            quote = TRUE)

# Alert user
message("✅ CSV file '_extra/county_name_compendium.csv' has been successfully created.")
