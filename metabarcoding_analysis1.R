# Load packages
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)

# Path to Excel file
file_path <- "all_reports.xlsx"

# Function to extract species names from one sheet
extract_species <- function(df) {
  # Drop the "root" and "unclassified" rows
  df <- df %>% filter(!(rank %in% c("R", "U")))
  
  # Identify numeric taxonomy columns
  tax_cols <- grep("^[0-9]+$", names(df), value = TRUE)
  
  # Extract the last non-NA taxon name (the species)
  df %>%
    rowwise() %>%
    mutate(
      species = {
        lineage <- na.omit(c_across(all_of(tax_cols)))
        if (length(lineage) == 0) NA else lineage[length(lineage)]
      }
    ) %>%
    ungroup() %>%
    select(species)
}
# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Read and extract species from all sheets
all_species <- map_dfr(sheet_names, ~ {
  df <- read_excel(file_path, sheet = .x)
  extract_species(df) %>% mutate(sheet = .x)
})

# Remove empty species names
all_species <- all_species %>% filter(!is.na(species))

# Summarize: count how many sheets each species appears in
species_summary <- all_species %>%
  distinct(sheet, species) %>%    # each species counted once per sheet
  count(species, name = "n_sheets_present") %>%
  arrange(desc(n_sheets_present))

# View result
print(species_summary)

View(species_summary)

# Optionally, save to CSV
write.csv(species_summary, "sp_occ_summary.csv", row.names = FALSE)
