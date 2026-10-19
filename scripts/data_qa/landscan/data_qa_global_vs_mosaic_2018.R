# /////////////////////////////////////////////////////////////////////////////
# Check the landscan mosaic against the landscan global data sets 
# Check both data sets to the ACS5.
# Global is much closer to the ACS5 estimates. The global could be off by 
# 250,000 when compared to mosaic. The global was off by as much as 40,000 vs
# the ACS5 and the mosaic was off by about 200,000 compared to the ACS5.
# /////////////////////////////////////////////////////////////////////////////
library(tidyverse)
library(tidycensus)
library(sf)
library(tidyterra)
library(terra)
library(exactextractr)

year <- 2018
noise_path_prefix <- "D:\\rwjf_noise\\data\\noise\\annual\\Noise_"
noise_path_suffix <- "_Noise.csv"
noise_file <- str_c(noise_path_prefix, year, noise_path_suffix)

exposure_df <-read_csv(noise_file, show_col_types = FALSE) %>%
    rename(lon = Longitude, lat = Latitude) %>%
    mutate(L_night_ge32 = ifelse(L_night >= 32, 1, 0))

ls_mosaic <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-mosaic_timeseries-2018-assets\\landscan_mosaic_timeseries_2018.tif")

ls_global <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-global-2018-assets\\landscan-global-2018.tif")


#%% Prepare LA County tracts object
options(tigris_progress_bar = FALSE)
options(tigris_use_cache = TRUE)

suppressMessages({
# Load LA County tracts as simple features (sf) from tidycensus
# crs loads as 4269
tracts <- get_acs(
  geography = "tract",

  # Corresponds to total population variable, not used but argument needs to be
  # set for the function to work
  variables = "B01003_001",   
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  year = 2018,
  survey = "acs5"
)
})

exposure_sf <- st_as_sf(
  exposure_df %>% 
    sample_frac(.1), 
  coords = c("lon", "lat"), 
  crs = 4326)

tracts <- st_transform(tracts, 4326)

# Get the minimum latitude value from the exposure data to use as 
# a cutoff for removing the islands
lat_cutoff <- min(exposure_df$lat)

# Create a bb object of LA County Tracts
bb <- st_bbox(tracts)

# Modify the bb object's y minimum value
# uses the previously defined cut off value above to raise the souther edge
bb["ymin"] <- lat_cutoff   

# Create an sfc object from the modified bb object that will be used to clip
# away the islands
clip_poly <- st_as_sfc(bb)

# and the clipping polygon.
la_tracts_cropped <- st_intersection(tracts, clip_poly)

# Convert to spatial vector to be used for cropping
la_vect <- vect(la_tracts_cropped)

la_global <- mask(crop(ls_global, la_vect), la_vect)

la_mosaic <- mask(crop(ls_mosaic, la_vect), la_vect)



global(la_global, "isNA", na.rm = FALSE)

global(la_mosaic, "isNA", na.rm = FALSE)

global_sums <- exact_extract(la_global, la_tracts_cropped, 'sum', progress = FALSE)

mosaic_sums <- exact_extract(la_mosaic, la_tracts_cropped, 'sum', progress = FALSE)

la_tracts_cropped <- la_tracts_cropped %>%
  mutate(
    global_sum = global_sums,
    mosaic_sum = mosaic_sums,
    diff =  mosaic_sum - global_sum
  )

la_tracts_cropped %>%
  ggplot() +
    geom_sf(aes(fill = diff)) +
    scale_fill_viridis_c() +
    labs(title = "Mosaic - Global Cell Values") +
    theme_minimal()


la_tracts_cropped %>%
  mutate(diff = mosaic_sum - estimate) %>%
  ggplot() +
    geom_sf(aes(fill = diff)) +
    scale_fill_viridis_c() +
    labs(title = "Mosaic Cell Values - ACS5 Estimate ") +
    theme_minimal()

la_tracts_cropped %>%
  mutate(diff = global_sum - estimate) %>%
  ggplot() +
    geom_sf(aes(fill = diff)) +
    scale_fill_viridis_c() +
    labs(title = "Global Cell Values - ACS5 Estimate ") +
    theme_minimal()