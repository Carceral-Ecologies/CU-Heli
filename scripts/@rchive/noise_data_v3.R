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
    mutate(L_night_ge32 = ifelse(L_night >= 32, 1, 0))

# Convert exposure_df to terra spatial vector (SpatVector) with point geometry
exposure_vect <- vect(
  exposure_df,
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

# Check that the exposure vect is in the same coordinate reference system
same.crs(pop, exposure_vect)

#%% Prepare LA County tracts object
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
la_tracts <- st_transform(tracts, crs(pop))

# remove rows where geometry is empty
# Tracts 9901:9903 corresponding to water tracts with no geometry are removed
# Required to perform exact_extract() otherwise it will fail to get geometry
# extent
la_tracts <- la_tracts[!st_is_empty(la_tracts), ]


# Convert to spatial vector
la_vect <- vect(la_tracts)

# Convert tracts to the same pop crs
la_vect <- project(la_vect, pop)

# Check if the tracts and pop have the same crs
same.crs(la_vect, pop)


#%% Crop the landscan population raster
# crop and mask have different functions, crop is to establish the bounding box
# of pop to that of la_vect, and mask is to set any value outside the county, 
# but within the bbox to NA
la_pop <- mask(crop(pop, la_vect), la_vect)


#%% Get a count of noise measurements per pixel
exp_count <- rasterize(
  exposure_vect,
  la_pop,
  fun = length,
  background = 0
)

# Check the max number of points assigned to a cell.
# The max value is 2
global(exp_count, "max", na.rm = TRUE)

# Check the number of pixels where the number of points per pixel is greater
# than 1.
# Results in 383,180 cells with exactly 2 points per pixel.
global(exp_count > 1, "sum", na.rm = TRUE)

# Total number of pixels
# 1,599,360 in LA county
sum(values(exp_count) > 0, na.rm = TRUE)

# Proportion of cells that have 2 points (noise values) per pixel
# 383180/1599360

# Tabulation of the frequencies of the values in exp_count pixels
# 2180348 pixels have values of 0, 1216180 have values of 1, 383180 have values of 2
freq(exp_count)

# Visualize
dup_cells <- exp_count == 2

plot(dup_cells)
plot(la_vect, add = TRUE)

#%% Custom mean_db function
# names(exposure_vect) # To display the available variable names during
# development
# L_* variables are in dB scale.
# L_bar = L - 10 * log10(n)
# L = 10 * log10(sum(10^(x/10)))
# eq.s refer to Onosokki document from Yair
mean_db <- function(x) {
  # Calculate the average power given a vector of values in decibel(dB) scale 
  # Only average non-missing values
  x <- x[!is.na(x)]

  # if there are no values within a pixel return NA
  if (length(x) == 0) return(NA_real_)
  
  # eq. 12-5
  # 10*log10(sum(10^(x/10))) - (10 * log10(length(x))) 

  #Alternative form from eq. 12-6, bc each p_n^2 / p_0^2 = 10^(L_n/10)
  10 * log10(mean(10^(x/10))) 
}

# Check the log_mean_db function
# This should return 77.4, per example equation 12-6 from Onosokki document
# x  <- c(80, 70)
# mean_db(x)

#%%
# Rasterize() essentially assigns points to landscan cells/pixels
exposure_raster <- rasterize(
  exposure_vect,
  la_pop,
  field = "L_night",
  fun = mean_db,
  background = NA
)




#%% Crop the exposure
# exposure cropped, average dB per pixel
exp_crop <- mask(crop(exposure_raster, la_vect), la_vect)

compareGeom(la_pop, exp_crop)

# is pixel valid, is not missing 
# Creates another raster of logicals where a dB value is not missing
valid <- !is.na(exp_crop)

# Adjust to valid non-missing pixels
pop_adj <- la_pop * valid # removes cells where noise is missing
weighted <- exp_crop * pop_adj



#%% Calculated a weighted exposure variable & plot
num <- exact_extract(weighted, la_tracts, 'sum')
den <- exact_extract(pop_adj, la_tracts, 'sum')

# Set the population weighted exposure variable by dividing num by den
la_tracts$pw_exposure <- num / den


# Plot population weighted exposure by tract
la_tracts %>%
  drop_na(pw_exposure) %>%
  ggplot() +
    geom_sf(aes(fill = pw_exposure)) +
    scale_fill_viridis_c() +
    labs(title = "Population-Weighted Exposure by Tract (dB)") +
    theme_minimal()
  ggsave("pw_exposure_x_tract.pdf")


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

# Plot the proportion of populated pixels
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