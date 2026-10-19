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
options(tigris_use_cache = TRUE)

# Prepare LA tracts -----------------------------------------------------------
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

# Remove San Clemente and Santa Catalina Islands by setting a cutoff value
# Create a lattitude cutoff to exclude San Clemente and Santa Catalina Islands
lat_cutoff <- 33.5

# Create a bb object of LA County Tracts
bb <- st_bbox(la_tracts)

# Modify the bb object's y minimum value, raise the southern edge
bb["ymin"] <- lat_cutoff

# Create an sfc object from the modified bb object that will be used to clip
# away the islands
clip_poly <- st_as_sfc(bb)

# Crop the la county map by taking the overlap between the tracts and the
# clipping polygon.
la_tracts_cropped <- st_intersection(la_tracts, clip_poly)

# Check tracts visually
ggplot() +
  geom_sf(data = la_tracts_cropped, fill = NA, color = "gray60")


la_tracts <- la_tracts_cropped

rm(la_tracts_cropped)

# Import LandScan Data --------------------------------------------------------
# What is the units of the land scan raster? Should be population count per 
# pixel
ls_path <- "D:\\rwjf_noise\\data\\landscan\\landscan-mosaic-unitedstates-v1-assets\\landscan-mosaic-unitedstates-v1.tif"
landscan_us <- rast(ls_path)

# transform la tracts to landscan coordinates
la_tracts_ls <- st_transform(la_tracts, crs(landscan_us))

landscan_la <- crop(landscan_us, vect(la_tracts_ls))
landscan_la <- mask(landscan_la, vect(la_tracts_ls))



# ///////////////////////////////////////////////////////////////////////////
# # Get bounding box of the tracts
# ext_la_tracts_ls <- ext(la_tracts_ls)

# # Crop raster to that extent
# la_pop <- mask(
#   crop(landscan_us, la_tracts_ls), 
#   la_tracts_ls
# )

# r_df <- as.data.frame(la_pop, xy = TRUE, na.rm = TRUE)

# # Name the raster column (terra uses layer names)
# names(r_df)[3] <- "pop"


# # Landscan population overlaid with LA county census tracts
# ggplot() +
#   geom_raster(
#     data = r_df,
#     aes(x = x, y = y, fill = pop)) + 
#   geom_sf(data = la_tracts_ls, fill = NA, color = "gray60") +
#   labs(
#     title = "Landscan Population overlaid with Los Angeles County Census Tracts") +
#   theme_minimal()


# //////////////////////////////// Noise data /////////////////////////////////
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

# Create Noise Level sf and convert to coordinate system
noise_sf <- st_as_sf(noise, coords = c("long", "lat"), crs = 4326)

noise_ls <- st_transform(noise_sf, crs(landscan_la))

coords <- st_coordinates(noise_sf)

noise_sf$cell_id <- cellFromXY(landscan_la, coords)

noise_sf

noise_sf$ls_pop <- terra::extract(landscan_la, noise_sf)[,1]

# Get the average noise per cell
noise_by_cell <- noise_sf %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    noise = mean(L_night, na.rm = TRUE),
    ls_pop = first(ls_pop))

noise_by_cell <- noise_by_cell %>%
  mutate(weighted_noise = noise * ls_pop) %>%
  ungroup()


# vector of cell IDs we care about
cell_ids <- noise_by_cell$cell_id

xy <- terra::xyFromCell(landscan_la, cell_ids)

sum(is.na(xy[,1]) | is.na(xy[,2]))  

valid <- complete.cases(xy)

noise_by_cell_valid <- noise_by_cell[valid, ]
xy_valid <- xy[valid, ]

cells_sf <- st_as_sf(
  data.frame(
    cell_id = noise_by_cell_valid$cell_id,
    x = xy_valid[, "x"],
    y = xy_valid[, "y"]
  ),
  coords = c("x", "y"),
  crs = crs(landscan_la)
)


cells_with_tracts <- st_join(
  cells_sf,
  la_tracts_ls |> select(GEOID),
  join = st_intersects
)


cells_full <- cells_with_tracts |>
  st_drop_geometry() |>
  left_join(noise_by_cell_valid, by = "cell_id")
