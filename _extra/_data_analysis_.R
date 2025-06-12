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


# Define a function to detect missing-like values in a data frame
detect_missing_values <- function(df) {
  # Define what counts as missing (besides NA)
  missing_patterns <- c(NA, "", " ", ".", "..", "...", "NA", "na", "Na", "\t", "\n")
  
  # Create a logical matrix TRUE if missing-like value, FALSE otherwise
  is_missing <- sapply(df, function(col) {
    # If factor, convert to character
    col_chr <- if (is.factor(col)) as.character(col) else col
    # Check if each element matches any missing pattern or is NA
    sapply(col_chr, function(x) {
      is.na(x) || trimws(x) %in% missing_patterns
    })
  })
  
  # Total cells
  total_cells <- nrow(df) * ncol(df)
  
  # Count missing
  missing_count <- sum(is_missing)
  
  # Percentage missing
  missing_pct <- 100 * missing_count / total_cells
  
  # Find positions of missing data (row and column)
  missing_positions <- which(is_missing, arr.ind = TRUE)
  missing_info <- data.frame(
    row = missing_positions[,1],
    col = colnames(df)[missing_positions[,2]],
    stringsAsFactors = FALSE
  )
  
  # Add name field if available
  if ("name" %in% colnames(df)) {
    missing_info$name <- df[missing_info$row, "name", drop = TRUE]
  } else {
    missing_info$name <- NA
  }
  
  list(
    missing_percentage = missing_pct,
    missing_count = missing_count,
    total_cells = total_cells,
    missing_positions = missing_info
  )
}

# Open a connection to the output file
output_file <- file("_extra/missing_data_info.txt", open = "wt")

# Redirect output to the file
sink(output_file)

# Run for 2022
missing_22 <- detect_missing_values(water_insecurity_2022)
cat("2022 Data Missing Percentage:", round(missing_22$missing_percentage, 4), "%\n")
print(missing_22$missing_positions)  # Show all missing locations

# Run for 2023
missing_23 <- detect_missing_values(water_insecurity_2023)
cat("2023 Data Missing Percentage:", round(missing_23$missing_percentage, 4), "%\n")
print(missing_23$missing_positions)  # Show all missing locations

# Stop redirecting output
sink()
close(output_file)

# Also print to screen
cat("2022 Data Missing Percentage:", round(missing_22$missing_percentage, 4), "%\n")
print(missing_22$missing_positions)

cat("2023 Data Missing Percentage:", round(missing_23$missing_percentage, 4), "%\n")
print(missing_23$missing_positions)

# Confirm to user
message("✅ Missing data info saved to '_extra/missing_data_info.txt'")
