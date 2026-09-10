# -----------------------------------------------------------------------------
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
# RWJF HelicopterSleep

# Noise data processing
#  Year = 2023 

# N.B. When data workflow is finalized, comment out code that randomly samples
# a subset of the noise data for ease of processing. 
# -----------------------------------------------------------------------------

#%% Libraries
library(tidyverse)
library(tidycensus)
library(sf)
library(terra)
library(exactextractr)
library(tigris)
library(raster) # Will over ride dplyr::select()
options(tigris_use_cache = TRUE)

#%% Load noise data
# Noise data delivered in .csv file(s) from UofI team
noise_full <- read_csv("D:\\rwjf_noise\\data\\noise\\Noise_Apr_2020.csv", show_col_types = FALSE)

# Rename columns, for friendlier coding
noise_full <- noise_full %>%
  rename("lat" = "Latitude", "long" = "Longitude")


#%% Noise data preparation
# Create a smaller data set for prototyping, (use less RAM, speed up code
# development). Take 20% sample
noise <- sample_frac(noise_full, .2)


#%% NOISE PLOT (used in development, not needed to prepare/process data)
# Generate a plot of the lat and long noise data
# ggplot(noise, aes(x = long, y = lat, color = L_night)) +
#   geom_point(size = 1) +
#   scale_color_viridis_c(option = "plasma") +
#   coord_equal() +
#   labs(title = "Noise Levels (L_night)",
#        x = "Longitude", y = "Latitude",
#        color = "Ldn (dB)") +
#   theme_minimal()


# COUNTY PLOT -----------------------------------------------------------------
# Overlay the LA County map: Performed in a series of steps that include:
# 1) Create an sf object from the noise. sf objects are created
# via the simple features package for spatial data. crs 4326 is a geographic 
# coordinate system.
# 2) Download the LA County shapefile, to overlay on top of the noise data. County
# data includes Santa Catalina and San Clemente islands. These are removed in the
# subsequent step.
# 3) In this step, the islands are removed by creating an overlay of the desired
# area that does not include them. This is accomplished by manipulating the
# lattitude value to use as a y limit in cartesian plotting. The intersection 
# of the county map and the desired area can then be used in downtream 
# processing.


# 1) Create Noise Level sf
noise_sf <- st_as_sf(noise, coords = c("long", "lat"), crs = 4326)

# 2) Download the LA County shapefile, and transform it the coordinate system
# uses tigris::counties()
la_county <- counties(state = "CA", cb = TRUE, year = 2023) %>%
  filter(NAME == "Los Angeles") %>%
  st_transform(4326)

# 3) Remove San Clemente and Santa Catalina Islands
# Create a lattitude cutoff to exclude San Clemente and Santa Catalina Islands
lat_cutoff <- min(noise$lat)

# Create a bb object of LA County
bb <- st_bbox(la_county)

# Modify the bb object's y minimum value
bb["ymin"] <- lat_cutoff   # raise the southern edge

# Create an sfc object from the modified bb object that will be used to clip
# away the islands
clip_poly <- st_as_sfc(bb)

# Crop the la county map by taking the overlap between the county map shape file
# and the clipping polygon.
la_county_cropped <- st_intersection(la_county, clip_poly)

# Plot noise level with LA County overlay
# Used in development (used in development, not needed to prepare/process data)
ggplot() +
  geom_sf(data = noise_sf, aes(color = L_night), size = 0.5, alpha = 0.8) +
  geom_sf(data = la_county_cropped, fill = NA, color = "green", size = 1) +
  scale_color_viridis_c(option = "plasma") +
  coord_sf() +
  labs(
    title = "Noise Levels (L_night) Over Los Angeles County",
    color = "Ldn (dB)"
  ) +
  theme_minimal()
ggsave("sample_l_night.pdf")

# -----------------------------------------------------------------------------
# Overlay the census tracts
# Link a census tract given a lattitude and longitude
# used tidycensus::get_acs()
# 1. Download LA County census tracts
la_tracts <- get_acs(
  geography = "tract",
  variables = "B01001_001E",   # total population (any variable works)
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  year = 2023
)

# Check the la tracts object
# Needs to have san clemente and santa catalina islands removed
# ggplot() +
#   geom_sf(data = la_tracts, fill = NA, color = "gray60")

# Create a bb object of LA County Tracts
bb <- st_bbox(la_tracts)

# Modify the bb object's y minimum value
# uses the previously defined cut off value with the county map above.
bb["ymin"] <- lat_cutoff   # raise the southern edge

# Create an sfc object from the modified bb object that will be used to clip
# away the islands
clip_poly <- st_as_sfc(bb)

# Crop the la county map by taking the overlap between the county map shape file
# and the clipping polygon.
la_tracts_cropped <- st_intersection(la_tracts, clip_poly)

# Check tracts visually
ggplot() +
  geom_sf(data = la_tracts_cropped, fill = NA, color = "gray60")


# Convert the noise_sf coordinate system to match that of what's listed in the
# LA tracts object.
noise_sf <- st_transform(noise_sf, st_crs(la_tracts))

# Joined data where a GEOID has NAs in some columns because some of the 
# locations where aircraft noise is measured do not fall within the LA
# tracts, after inspecting with subsequent plot.
noise_sf_tract <- st_join(noise_sf, la_tracts_cropped)

#%% Check to see where the noise values that have an NA in GEOID get plotted to
# see what didn't get linked. Appears as though they did not get assigned a
# GEOID because they fall outside the county area.
# Performed as a data qa step. Commented out to speed up development upon 
# reopening an R Session
# ggplot() +
#   geom_sf(
#     data = (noise_sf_tract %>% filter(is.na(GEOID))), # only where GEOID is NA
#     aes(color = L_night), 
#     size = 0.5, 
#     alpha = 0.8) +
#   geom_sf(
#     data = la_county_cropped, 
#     fill = NA, 
#     color = "green", 
#     size = 1) +
#   scale_color_viridis_c(option = "plasma") +
#   coord_sf() +
#   labs(
#     title = "Noise Levels (L_night) Over Los Angeles County",
#     color = "Ldn (dB)"
#   ) +
#   theme_minimal()


#%% Resume processing
# Filter out the values where GEOID is missing
# It's safe to drop those that did not get linked, after reviewing plot
noise_sf_tract <- noise_sf_tract %>%
  filter(!is.na(GEOID))

# Plott the subset of noise levels, tracts, and county boundary
# * Results in some jagged edges. For prototyping, we might be able to do
# something like expanding the line width of the county map and turning it to
# white. For publication, we may want to make crop the noise data by the LA
# count map object.
ggplot() +
  # L_night noise levels
  geom_sf(
    data = (noise_sf_tract),
    aes(color = L_night), 
    size = 0.5, 
    alpha = 0.8) +
  # Tract boundaries
  geom_sf(data = la_tracts_cropped, fill = NA, color = "gray60") +
  # # County boundary
  # geom_sf(
  #   data = la_county_cropped, 
  #   fill = NA, 
  #   color = "white", 
  #   size = 1) +
  scale_color_viridis_c(option = "plasma") +
  coord_sf() +
  labs(
    title = "Noise Levels (L_night) Over Los Angeles County",
    color = "Ldn (dB)"
  ) +
  theme_minimal()


# %% /////////////////////////////////////////////////////////////////////////////
# For dasymetric weighting, need additional information such as land use or 
# population to develop a scale that can then be used to weight each noise
# level. Varun mentioned using additional methods like a population weighted 
# centroid. Data can be found at landscan by OakRidge Natl. Lab.


# Import LandScan Data --------------------------------------------------------

landscan_us <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-mosaic-unitedstates-v1-assets\\landscan-mosaic-unitedstates-v1.tif")

# transform la cropped to landscan coordinates
la_tracts_ls <- st_transform(la_tracts_cropped, crs(landscan_us))

# This is the raw raster sum
pop_sum <- exact_extract(landscan_us, la_tracts_ls, 'sum')

# need to merge pop_sum to la_tracts_cropped
la_tracts_ls$raw_raster_sum <- pop_sum

# THIS WORK FLOW WORKS UP TO ABOUT HERE NEED NOISE IN LANDSCAN GRID, TO AVERAGE OUT NOISE PER CELL
# This could be replaced with a control level population
# la_tracts_ls$landscan_pop <- pop_sum

# # Create a scale factor
# la_tracts_ls$scale_factor <- as.numeric(la_tracts_ls$landscan_pop) / as.numeric(la_tracts_ls$raw_raster_sum)

# scale_rast <- rasterize(vect(la_tracts_ls), landscan_us, field = "scale_factor")

# landscan_scaled <- landscan_us * scale_rast

# /////////////////////////////////////////////////////////////

# LandScan raster (population)
pop <- rast("D:\\rwjf_noise\\data\\landscan\\landscan-mosaic-unitedstates-v1-assets\\landscan-mosaic-unitedstates-v1.tif")

# Exposure dataset (CSV with lon/lat + value)
exposure_df <-read_csv("D:\\rwjf_noise\\data\\noise\\Noise_Apr_2020.csv", show_col_types = FALSE)

exposure_df <- exposure_df %>%
  rename(lon = Longitude, lat = Latitude)

exposure_df <- sample_frac(exposure_df, .2)

# Convert to sf
exposure_sf <- st_as_sf(exposure_df,
                       coords = c("lon", "lat"),
                       crs = 4326)

# Convert to terra SpatVector
exposure_vect <- vect(exposure_sf)

exposure_vect <- project(exposure_vect, pop)

# Add exposure field name (assume column is "value")
names(exposure_vect)

# Rasterize: mean exposure per pixel
exposure_raster <- rasterize(
  exposure_vect,
  pop,
  field = "L_night",
  fun = mean,
  background = NA
)


# Step 4
county <- counties(
    state = "CA", 
    cb = TRUE, 
    year = 2023) %>%
  filter(NAME == "Los Angeles") %>%
  st_transform(4326)

county_vect <- vect(county)

# Ensure same CRS
county_vect <- project(county_vect, pop)

# Crop + mask
pop_county <- mask(crop(pop, county_vect), county_vect)
exp_county <- mask(crop(exposure_raster, county_vect), county_vect)

# Step 5 
# tracts <- st_read("tracts.shp")
tracts <- get_acs(
  geography = "tract",
  variables = "B01001_001E",   # total population (any variable works)
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE,
  year = 2023
)

tracts <- st_transform(tracts, crs(pop))
exp_county <- project(exp_county, pop)

exp_county <- resample(exp_county, pop, method = "bilinear")

tracts_vect <- vect(tracts)

pop_crop <- mask(crop(pop, tracts_vect), tracts_vect)
exp_crop <- mask(crop(exp_county, tracts_vect), tracts_vect)

valid <- !is.na(exp_crop)

pop_adj <- pop_crop * valid
weighted <- exp_crop * pop_adj


weighted_r <- raster(weighted)
pop_r <- raster(pop_adj)

# After running the lines below something worked here,
# could be manipulations of tracts
tracts <- st_make_valid(tracts)
tracts <- tracts[!st_is_empty(tracts), ] # error was happening because of this
num <- exact_extract(weighted_r, tracts, 'sum')
den <- exact_extract(pop_r, tracts, 'sum')


# this worked
tracts <- st_make_valid(tracts)
tracts <- tracts[!st_is_empty(tracts), ]
st_crs(tracts)
crs(weighted)
tracts <- st_transform(tracts, crs(weighted))
ext(weighted)
st_bbox(tracts)
num <- exact_extract(weighted, tracts, 'sum')
den <- exact_extract(pop_adj, tracts, 'sum')


tracts$pw_exposure <- num / den

compareGeom(pop_crop, exp_crop)

ext(pop_crop)
st_bbox(tracts)

plot(pop_crop)
plot(vect(tracts), add = TRUE)

summary(den)

plot(tracts["pw_exposure"])

ggplot(tracts) +
  geom_sf(aes(fill = pw_exposure)) +
  scale_fill_viridis_c() +
  labs(title = "Population-Weighted Exposure by Tract") +
  theme_minimal()

# diagnose bad tract
problem_tract <- tracts[is.na(tracts$pw_exposure), ]

nrow(problem_tract)

plot(exp_crop)
plot(vect(problem_tract), add = TRUE, col = "red", lwd = 2)

plot(pop_crop)
plot(vect(problem_tract), add = TRUE, col = "blue", lwd = 2)

# https://www.esri.com/about/newsroom/blog/la-county-maps-equitable-nature-access