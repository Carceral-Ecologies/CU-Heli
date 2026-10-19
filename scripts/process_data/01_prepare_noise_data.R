# /////////////////////////////////////////////////////////////////////////////
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
#
# This script contains a function that is designed to loope through different
# noise input data sets, and different years of the landscan global dataset
# to create the L_night, L_day, L_den, and other variables that are population
# weighted. For each noise variable, two additional variables are created, a
# population weighted variable in db and a linear form of the population 
# weighted variable. The outputs are to crete multiple data sets that can be
# ready for modeling.
# /////////////////////////////////////////////////////////////////////////////


library(tidyverse)
library(tidycensus)
library(sf)
library(tidyterra)
library(terra)
library(exactextractr)

# Specify the noise file path
# noise_file <- "D:\\rwjf_noise\\data\\noise\\annual\\2018_Noise.csv"

# Specify the year to use for the landscan
# year <- 2018

# Specify the out file path
# out_file <- "D:\\rwjf_noise\\data\\noise\\lac_pw_db_2018.rds"

process_noise <- function(noise_file, year, out_file) {

  # Specify the landscan file path
  landscan_path <- str_c("D:\\rwjf_noise\\data\\landscan\\landscan-global-", year, "-assets\\landscan-global-", year, ".tif")


  #%% Load the landscan population raster stored locally ----------------------
  # has resolution of 1/1200 degree
  # coordinate refeference system (crs) is 4326
  pop <- rast(landscan_path)

  #%% Load noise data --------------------------------------------------------- 
  # File is .csv format with a noise level for each lat and long coordinates.
  # Rename coordinate columns for convenience.
  exposure_df <- read_csv(noise_file, show_col_types = FALSE) %>%
      rename(lon = Longitude, lat = Latitude)
  
  #%% Prepare LA County tracts object -----------------------------------------
  options(tigris_progress_bar = FALSE)
  options(tigris_use_cache = TRUE)

  suppressMessages({
    # Load LA County tracts as simple features (sf) from tidycensus
    # crs loads as 4269
    # May need tidycensus API to be set
    tracts <- get_acs(
      geography = "tract",

    # Corresponds to total population variable, not used but argument needs to be
    # set for the function to work
    variables = "B01001_001E",   
    state = "CA",
    county = "Los Angeles",
    geometry = TRUE,
    year = year
  )
  })

  # Samples a fraction of the exposure data frame for plotting, otherwise it's too
  # computationally intensive to plot. Complete data set is over 2million rows.
  exposure_sf <- st_as_sf(
    exposure_df %>% 
      sample_frac(.1), 
    coords = c("lon", "lat"), 
    crs = 4326)

  # Transform tracts crs from 4269 to 4326
  tracts <- st_transform(tracts, 4326)

  # Check if both exposure and tracts have the same crs
  # If they don't
  if (st_crs(exposure_sf) != st_crs(tracts)){
    stop("Mismatch in the CRS of the exposure and census tracts.")
  }

  # Create a clipping polygon -------------------------------------------------
  # Clippling polygon to crop maps to excude San Celemente and Santa Catalina
  # Islands
  
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

  # Crop the tracts by taking the overlap between the tracts sf object
  # and the clipping polygon.
  la_tracts_cropped <- st_intersection(tracts, clip_poly)


  # remove rows where geometry is empty
  # Tracts where GEOID end in 9901:9903 with no geometry are removed
  # Empty tracts must be removed to successfully perform exact_extract() in a
  # downstream step, otherwise "fail to get geometry extent" errors arise.
  # la_tracts_cropped <- la_tracts_cropped %>%
  #   filter(!st_is_empty(la_tracts_cropped))

  # Convert to spatial vector to be used for cropping
  la_vect <- vect(la_tracts_cropped)

  # Check if the tracts and pop have the same crs
  same.crs(la_vect, pop)

  #%% Crop the landscan population raster
  # crop and mask have different functions, crop is to establish the bounding box
  # of pop to that of la_vect, and mask is to set any value outside the tracts, 
  # but within the bbox to NA
  la_pop <- mask(crop(pop, la_vect), la_vect)

  cells_w_missing_values <- global(la_pop, "isNA", na.rm = FALSE)


  variables <- c( "L_day", "L_evening", "L_night", "L_den", "L_dn")

  for (v in variables) {
    # Linearize the dB values
    exposure_df <- exposure_df %>%
      mutate(var := 10 ^ (!!rlang::sym(v)/10))

    # Create raster directly from points, instead of points to spatial vector to raster
    exposure_rast <- rast(
      exposure_df[, c("lon", "lat", "var")],
      type = "xyz",
      crs = "EPSG:4326"
    )

    # Resample the exposure raster to landscan grid
    exp_landscan_lin <- resample(
      exposure_rast,
      la_pop,
      method = "bilinear"
    )

    ## Transform the resampled exposure raster back to dB
    # Convert exposure raster back to dB
    exp_landscan_db <- 10 * log10(exp_landscan_lin)

    # Check the resolutions of exp_landscan and la_pop
    resolutions <- c(
      res(exp_landscan_lin) == res(exp_landscan_db),
      res(exp_landscan_lin) == res(la_pop),
      res(exp_landscan_db) == res(la_pop)
    )
    
    if (all(resolutions) != TRUE) {
      stop("Not all resolutions are the same.")
    }
      
    # Compare la_pop and exp_landscan geometries
    # At this point la_pop is cropped to exclude Santa Catalina and San Celemente islands
    compareGeom(la_pop, exp_landscan_lin)
    
    # Create a weighted value for each pixel by element-wise multiplication of the
    # population and exposure rasters
    weighted <- exp_landscan_lin * la_pop
    
    # Check the range of values
    minmax(weighted)

    #%% Calculated a weighted exposure variable
    num <- exact_extract(weighted, la_tracts_cropped, 'sum', progress = FALSE)

    den <- exact_extract(la_pop, la_tracts_cropped, 'sum', progress = FALSE)

    # Set the population weighted exposure variable by dividing num by den
    la_tracts_cropped <- la_tracts_cropped %>%
      mutate(pw_exp_lin = num / den)

    # Convert the linearized noise levels back to dB
    la_tracts_cropped <- la_tracts_cropped %>%
      mutate(pw_exp_db = 10 * log10(pw_exp_lin))
    
    # Rename the variable
    la_tracts_cropped <- la_tracts_cropped %>%
      rename(
        !!rlang::sym(str_c("pw_", v, "_db")) := pw_exp_db,
        !!rlang::sym(str_c("pw_", v, "_lin")) := pw_exp_lin)
    
  }


  # Output file to .rds
  saveRDS(la_tracts_cropped, out_file)
}


# Create 2018 for LACHS
process_noise(
  "D:\\rwjf_noise\\data\\noise\\annual\\2018_Noise.csv",
  2018,
  "D:\\rwjf_noise\\data\\noise\\lac_pw_db_2018.rds"
)


# Create 2023 for LACHS
process_noise(
  "D:\\rwjf_noise\\data\\noise\\annual\\2023_Noise.csv",
  2023,
  "D:\\rwjf_noise\\data\\noise\\lac_pw_db_2023.rds"
)
