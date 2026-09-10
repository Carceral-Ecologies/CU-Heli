# -----------------------------------------------------------------------------
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
# RWJF HelicopterSleep

# Noise data processing
#  Tract & county maps are from year = 2023 
#  Noise is from 2020
#  Population data is from

# N.B. When data workflow is finalized, comment out code that randomly samples
# a subset of the noise data for ease of processing. 
# -----------------------------------------------------------------------------

#%% Libraries
library(tidyverse)
library(tidycensus)
library(sf)
library(terra)
library(exactextractr)
# library(tigris)
# options(tigris_use_cache = TRUE)

#%%
# Load the landscan population raster
pop <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-mosaic-unitedstates-v1-assets\\landscan-mosaic-unitedstates-v1.tif")


# Noise data processing
# Load the noise dataset (.csv with lon/lat + values) and rename coordinate
# columns for convenience.
# *** Currently set to sample 20% of noise data to speed up processing while
# prototyping, since the noise data set is very large
exposure_df <-read_csv("D:\\rwjf_noise\\data\\noise\\Noise_Apr_2020.csv", show_col_types = FALSE) %>%
    rename(lon = Longitude, lat = Latitude) %>%
    # slice_sample(.2) %>% # comment uncomment as needed to speed up processing while prototyping
    mutate(L_night_ge32 = ifelse(L_night >= 32, 1, 0))

# Convert noise to sf and set coordinate system
exposure_sf <- st_as_sf(exposure_df,
                       coords = c("lon", "lat"),
                       crs = 4326)

# Convert to terra SpatVector
exposure_vect <- vect(exposure_sf)

# Vector projection to pop crs
exposure_vect <- project(exposure_vect, pop)

# show variable names
names(exposure_vect)

# Rasterize the exposure vector using mean L_night
# This SpatRaster represents the mean L_night value per pixel in the geometry of
# the pop SpatRaster.
exposure_raster <- rasterize(
  exposure_vect,
  pop,
  field = "L_night",
  fun = mean,
  background = NA
)

# Get a count of noise measurements per pixel
exp_count <- rasterize(
  exposure_vect,
  pop,
  fun = length,
  background = 0
)

# Load county map as sf with common coordinate system in order to crop the
# population and exposure rasters
# county <- tigris::counties(
#     state = "CA", 
#     cb = TRUE, 
#     year = 2023) %>%
#   filter(NAME == "Los Angeles") %>%
#   st_transform(4326)

# Convert to terra SpatVector
# county_vect <- vect(county)

# Vector projection to pop crs
# county_vect <- project(county_vect, pop)

# Crop & mask the population and exposure rasters to the county level
# *** I don't think these are necessary, they seem to veer off the goal
# pop_county <- mask(crop(pop, county_vect), county_vect)
# exp_county <- mask(crop(exposure_raster, county_vect), county_vect)


# Load LA County tracts as sf, from tidycensus 
tracts <- get_acs(
  geography = "tract",
  variables = "B01001_001E",   # total population (any variable works)
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  year = 2023
)

# transform the tracts to the common coordinate system
tracts <- st_transform(tracts, crs(pop))


# # *** This seems redundant, why project exp_county to pop, it should have already been in the 
# # same space
# exp_county <- project(exp_county, pop)
# exp_county <- resample(exp_county, pop, method = "bilinear") # This takes a really long time

# Convert to spatial vector
tracts_vect <- vect(tracts)


# # ----------- CHECK HERE
# # New version didnt work
# # Vector projection to pop crs
# # tracts_vect <- project(vect(tracts), pop)
# # pop_crop <- mask(crop(pop, tracts_vect), tracts_vect)
# # exp_crop <- mask(crop(exposure_raster, tracts_vect), tracts_vect) # This works for creating exp_crop

# Original version, crop to the LA County tracts
pop_crop <- mask(crop(pop, tracts_vect), tracts_vect)

# exposure cropped, average noise per pixel
# Instead of operating at the county level, it should be at the tract
exp_crop <- mask(crop(exposure_raster, tracts_vect), tracts_vect)
compareGeom(pop_crop, exp_crop)

# is pixel valid, is not missing 
valid <- !is.na(exp_crop)

# Adjust to valid non-missing pixels
pop_adj <- pop_crop * valid
weighted <- exp_crop * pop_adj

# *** Some thing breaks here with old and new version Error getting geometry extent
# loading library(raster) in workspace didn't do anything
tracts <- tracts[!st_is_empty(tracts), ] # error was happening because of this was not performed

num <- exact_extract(weighted, tracts, 'sum')
den <- exact_extract(pop_adj, tracts, 'sum')

tracts$pw_exposure <- num / den

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