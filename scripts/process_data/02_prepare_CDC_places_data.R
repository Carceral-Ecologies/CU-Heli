# /////////////////////////////////////////////////////////////////////////////
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
#
# Prepare a data set for modeling the relationship between aircraft noise and
# insufficient sleep.
# This script performs two critical merges of data. First, CDC places data set
# to noise data set. Second, Census ACS data for covariates to the output of 
# the first merge. Saves as merged places data as a .RDS files to maintain sf
# object features.
# /////////////////////////////////////////////////////////////////////////////

# Load Packages ---------------------------------------------------------------
library(tidyverse)
library(sf)
library(units)
library(tidycensus)

# Set options for tigris
options(tigris_progress_bar = FALSE)
options(tigris_use_cache = TRUE)


# Read noise data -------------------------------------------------------------
noise <- readRDS("D:\\rwjf_noise\\data\\noise\\lac_pw_db.rds")

# Create a TractFIPS column for merging the cdc places data set to the noise 
# data set
noise <- noise %>%
  rename(TractFIPS = GEOID)


# Load the places data --------------------------------------------------------
# The 2022 Places Data Set corresponds to BRFSS data from 2020, the year that 
# the current average annual noise come from.
fpath <- "D:\\rwjf_noise\\data\\places\\PLACES__Census_Tract_Data_(GIS_Friendly_Format),_2022_release_20260605.csv"

# shp <- st_read("fpath")
places <- read_csv(fpath, show_col_types = FALSE)

# Places data contain data at the country level. Filter data to the county of
# interest. Results in 2,324 tracts within LA County
places <- places %>% 
  filter(CountyName == "Los Angeles", StateAbbr == "CA")

# Places data has a TractFIPS at 11 digits
# Census data has a GEOID at 12 digits.
# Those with the extra digit in the GEOID is for the block identifier
# We may be able to create TractFIPS in Census to merge in to Places.


# Merge 1 of 2. Places to noise -----------------------------------------------
data <- left_join(
  noise,
  places,
  by = "TractFIPS"
)



# Download ACS data code book for verifying variables -------------------------
# Use acs5 for B variables
# Use acs5/subject for S variables
# Use acs5/profile for D variables
code_book <- load_variables(2022, "acs5/profile", cache = TRUE)

# use name or label to look for variables or labels
code_book %>% 
  filter(grepl("DP05_0077", name))

acs_vars <- c(
  total_pop = "B01003_001E", # OK
  median_value = "B25077_001E", # OK
  median_year_built = "B25037_001E", # OK
  median_income = "S1901_C01_012E", # OK
  pct_homeowners = "DP04_0046PE", # OK
  pct_foreign_born = "DP02_0094PE", # OK
  pct_college = "DP02_0068PE", # OK, subsets to 25+ yrs age 
  latino_count = "DP05_0071E", # OK hispanic or latino of any race
  white_count = "DP05_0077E", # OK, non-hispanic white
  black_count = "DP05_0078E", # OK, non-hispanic black
  native_count = "DP05_0079E", # OK, non-hispanic native
  asian_count = "DP05_0080E", # OK, non-hispanic native
  median_age = "B01002_001"

)

# Load the census data --------------------------------------------------------
census <- get_acs(
  geography = "tract",
  variables = acs_vars,
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  output = "wide",
  year = 2020
)

# Set the land area and the population density in per km squared units
census <- census %>%
  mutate(land_area_km2 = set_units(st_area(census), km^2)) %>%
  mutate(pop_density = as.numeric(total_pop / land_area_km2))

# Create TractFIPS. This is the variable name used in the CDC places data set
census <- census %>%
  mutate(TractFIPS = GEOID) %>%
  select(-NAME)

# Drop the geometry information and convert to data frame for merging
census <- st_drop_geometry(census)


# Merge 2 of 2 census to data -------------------------------------------------
data <- data %>%
  left_join(
    census,
    by = "TractFIPS"
  )


# Write out the merged data ---------------------------------------------------
saveRDS(data, "D:\\rwjf_noise\\data\\prepped_data\\merged_places_2026-06-22.rds")



# Plot the tracts that don't have places data
# noise %>%
#   drop_na(pw_exposure) %>%
#   mutate(has_places = ifelse(TractFIPS %in% lac$TractFIPS, "Yes", "No")) %>%
#   ggplot() +
#     geom_sf(aes(fill = has_places)) +
#     labs(title = "Availability of Places Data by Tract (dB)") +
#     theme_minimal()
# ggsave("places_data_availability_x_tract.pdf")


# Check the number of tracts that have Places Data values
# Only 2000 tracts have places data, whereas the exposure data contains 2,493 tracts
# noise %>%
#   drop_na(pw_exposure) %>%
#   mutate(has_places = ifelse(TractFIPS %in% lac$TractFIPS, "Yes", "No")) %>%
#   pull(has_places) %>%
#   table()