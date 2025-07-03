# Data
- **Dataset**: This project uses two datasets from the **American Community Survey (ACS)** 5-Year Estimates, accessed via the **TidyTuesday** project dated **2025-01-28**. These datasets contain county-level data on indoor plumbing access across the U.S. for the years **2022** and **2023**. The data allow for analysis of water insecurity through the lens of plumbing access, population demographics, and geographic patterns. The data were originally compiled by the U.S. Census Bureau and visualized by the USGS Vizlab team.

# Codebook for Water Insecurity Dataset (`water_insecurity_2022` and `water_insecurity_2023`)

## Variable Names and Descriptions:

- **`geoid`**: Unique identifier for the county (FIPS code)
- **`name`**: County name
- **`year`**: The year of the ACS estimate (2022 or 2023)
- **`geometry`**: County boundary geometry as an `sf` object (used for mapping)
- **`total_pop`**: Total population of the county
- **`plumbing`**: Number of households without complete indoor plumbing
- **`percent_lacking_plumbing`**: Percentage of the population lacking complete plumbing facilities

## Data Types:

- **`geoid`**: `character`
- **`name`**: `character`
- **`year`**: `numeric`
- **`geometry`**: `sf / geometry`
- **`total_pop`**: `numeric`
- **`plumbing`**: `numeric`
- **`percent_lacking_plumbing`**: `numeric`
