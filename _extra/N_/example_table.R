# Load the gt package
library(gt)

# Create base table data
table_data <- data.frame(
  Variables = c("geoid", "name", "year", "geometry", 
                "total_pop", "plumbing", "percent_lacking_plumbing"),
  Q1 = c("✅", "", "✅", "", "✅", "", ""),
  Q2 = c("", "✅", "", "✅", "", "✅", "✅")
)

# Build and style the table
table_data |>
  gt() |>
  tab_header(
    title = "Variable Availability by Quarter"
  ) |>
  cols_label(
    Variables = "Variables",
    Q1 = "Q1",
    Q2 = "Q2"
  ) |>
  text_transform(
    locations = cells_body(columns = c(Q1, Q2)),
    fn = function(x) {
      dplyr::case_when(
        x == "✅" ~ "<span style='color: green; font-size: 20px;'>✅</span>",
        TRUE ~ ""
      )
    }
  ) |>
  fmt_markdown(columns = c(Q1, Q2))  # Render HTML inside cells

print(table_data)