# -----------------------------------------------------------------------------
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
# RWJF HelicopterSleep

# Data processing
#  Tract & county maps are from year = 2020
#  Noise is from 2020 annual
#  Population data is from

# N.B. When data workflow is finalized, comment out code that randomly samples
# a subset of the noise data for ease of processing. 
# # 
# I think the hang up there for me is the approach. We currently start off with
# noise levels in latitude longitude coordinates, we then transform those levels
# to match up with the land scan grid. For those grid cells that have more than
# one noise level, we average the noise levels per cell. Then we use the
# population value to weight the noise level. 

# From there, we calculate a population weighted value per tract by taking the
# sum of noise levels multiplied by the population in a cell and divide by the
# sum of population in all cells within a tract

# This version gets counts
# Considering that 
# The average of 2020
# count returns the number of times where the yearly average estimated noise level was estimated to be above 50 dB

# Using monthly data
# count returns the number of times over a year where the monthly average noise level was estimated to be above 50 dB.
# -----------------------------------------------------------------------------

#%% Libraries
library(tidyverse)
library(tidycensus)
library(sf)
library(terra)
library(exactextractr)

#%% Load the landscan population raster stored locally
pop <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-mosaic-unitedstates-v1-assets\\landscan-mosaic-unitedstates-v1.tif")

# Set File
noise_file <- "D:\\rwjf_noise\\data\\noise\\annual\\Noise_Annual_2020.csv"

#%% Load and prepare noise data 
# File is .csv format with a noise level for each lat and long coordinates.
# Rename coordinate columns for convenience.
exposure_df <-read_csv(noise_file, show_col_types = FALSE) %>%
    rename(lon = Longitude, lat = Latitude) %>%
    # slice_sample(.2) %>% # comment/uncomment as needed to speed up development
    mutate(
      L_night_ge32 = ifelse(L_night >= 32, 1, 0),
      L_night_gt50 = ifelse(L_night > 50, 1, 0),
      L_night_gt40 = ifelse(L_night > 40, 1,0),
      L_night_gt60 = ifelse(L_night > 60, 1,0)
    )

# Convert noise to sf and set coordinate system
exposure_sf <- st_as_sf(exposure_df,
                       coords = c("lon", "lat"),
                       crs = 4326)

# Convert the sf object to terra SpatVector
# needed for vector projection to landscan pop dataset
exposure_vect <- vect(exposure_sf)

# Vector projection to pop crs
exposure_vect <- project(exposure_vect, pop)


# Rasterize the exposure vector using mean L_night
# names(exposure_vect)
exposure_raster <- rasterize(
  exposure_vect,
  pop,
  field = "L_night_gt50",
  fun = sum,
  background = NA
)


#%% Prepare LA County tract objects
# Load LA County tracts as sf, from tidycensus 
tracts <- get_acs(
  geography = "tract",
  variables = "B01001_001E",   # total population (any variable works)
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  year = 2020
)

# transform the tracts to the common coordinate system
tracts <- st_transform(tracts, crs(pop))

# Convert to spatial vector
tracts_vect <- vect(tracts)

# Crop to the LA County tracts
pop_crop <- mask(crop(pop, tracts_vect), tracts_vect)

# exposure cropped, number of measurements above 50 per cell
exp_crop <- mask(crop(exposure_raster, tracts_vect), tracts_vect)
compareGeom(pop_crop, exp_crop)

# is pixel valid, is not missing 
valid <- !is.na(exp_crop)

# Adjust to valid non-missing pixels
pop_adj <- pop_crop * valid
weighted <- exp_crop * pop_adj

# remove empty tracts
tracts <- tracts[!st_is_empty(tracts), ] 

#%% Calculated a weighted exposure variable & plot
# num <- exact_extract(weighted, tracts, 'sum') # Sum of counts multiplied by population where pixels are valid within a tract
# den <- exact_extract(pop_adj, tracts, 'sum') # Sum of population where pixels are valid

tracts$n_gt_50 <- exact_extract(exp_crop, tracts, "sum")
tracts$population <- exact_extract(pop_adj, tracts, 'sum')

# tracts$pw_exposure <- num / den

compareGeom(pop_crop, exp_crop)

# Plot population weighted exposure by tract
tracts %>%
  drop_na(pw_exposure) %>%
  ggplot() +
    geom_sf(aes(fill = pw_exposure)) +
    scale_fill_viridis_c() +
    labs(title = "Population-Weighted Exposure by Tract (dB)") +
    theme_minimal()
  ggsave("pw_exposure_x_tract.pdf")

# diagnose bad tract, could be due to the sampling, and one time it didn't capture all of the noise
# problem_tract <- tracts[is.na(tracts$pw_exposure), ]

# nrow(problem_tract)

# plot(exp_crop)
# plot(vect(problem_tract), add = TRUE, col = "red", lwd = 2)

# plot(pop_crop)
# plot(vect(problem_tract), add = TRUE, col = "blue", lwd = 2)

# https://www.esri.com/about/newsroom/blog/la-county-maps-equitable-nature-access


#%% Diagnostic plot
# % of pixels in a tract that have a noise value AND have people living in it.

# Use max, because it will either be one or 0
exposure_presence_raster <- rasterize(
  exposure_vect,
  pop,
  field = "L_night_ge32",
  fun = max,
  background = 0
)

# Mask & crop exposure presence to tract
exp_present_crop <- mask(crop(exposure_presence_raster, tracts_vect), tracts_vect)

# Binarize the population values
pop_presence <- pop > 0

# Mask & crop population presence to tract
pop_present_crop <- mask(crop(pop_presence, tracts_vect), tracts_vect)

compareGeom(pop_present_crop, exp_present_crop)


num <- exact_extract(pop_present_crop, tracts, 'sum')
den <- exact_extract(pop_crop, tracts, function(values, coverage_fraction) {
  sum(coverage_fraction)
})

tracts$prop_pop_pixels <- num / den

# Output file to .rds
saveRDS(tracts, "D:\\rwjf_noise\\data\\noise\\lac_pw_noise.rds")

# Plot the 
tracts %>%
  drop_na(pw_exposure) %>% #used for filtering only
  ggplot() +
    geom_sf(aes(fill = prop_pop_pixels)) +
    scale_fill_viridis_c() +
    labs(title = "Proportion of Populated Pixels") +
    theme_minimal()
ggsave("prop_populated_pixels.pdf")

mean(tracts$prop_pixels_w_value)

# /////////////////////////////////////////////////////////////////////////////
# Both noise greater than 32 and population present
both <- pop_present_crop * exp_present_crop

num <- exact_extract(both, tracts, 'sum')
den <- exact_extract(pop_present_crop, tracts, 'sum')

tracts$prop_exposed_pixels <- num / den

tracts %>%
  drop_na(pw_exposure) %>% #used for filtering only
  ggplot() +
    geom_sf(aes(fill = prop_exposed_pixels)) +
    scale_fill_viridis_c() +
    labs(title = "Proportion of Populated Pixels with Exposure >= 32dB") +
    theme_minimal()
ggsave("prop_exposed_pixels_ge32.pdf")


ggplot() +
  geom_sf(
    data = exposure_sf %>% slice_sample(prop = 0.1), 
    aes(color = L_night), 
    size = 0.5) +
  geom_sf(
    data = tracts, 
    fill = NA, 
    color = "black") +
  scale_color_viridis_c() +
  theme_minimal()