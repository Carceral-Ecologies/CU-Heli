# Seems like the LACHS will have a census tract 

# Use the census tract to merge in census data

# Helicopter data should either have a census tract id or something else to link
# to the LACHS data

# Likely a series of models for each outcome variable.


# Inspect census files
library(tidyverse)
library(sf)
library(tigris)



noise_full <- read_csv("D:\\rwjf_noise\\data\\noise\\Noise_Apr_2020.csv", show_col_types = FALSE)
noise_full <- noise_full %>%
  rename("lat" = "Latitude", "long" = "Longitude")

# To make things easier for troublshooting/prototyping, sample a percentage of
# the data
noise <- sample_frac(noise_full, .2)

# Generate a plot of the lat and long data
ggplot(noise, aes(x = long, y = lat, color = L_dn)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  labs(title = "Noise Levels (Ldn)",
       x = "Longitude", y = "Latitude",
       color = "Ldn (dB)") +
  theme_minimal()


# Overlay the LA County map
options(tigris_use_cache = TRUE)

grid_sf <- st_as_sf(noise, coords = c("long", "lat"), crs = 4326)

la_county <- counties(state = "CA", cb = TRUE, year = 2023) %>%
  filter(NAME == "Los Angeles") %>%
  st_transform(4326)

lat_cutoff <- min(noise$lat)

bb <- st_bbox(la_county)
bb["ymin"] <- lat_cutoff   # raise the southern edge

clip_poly <- st_as_sfc(bb)

la_county_cropped <- st_intersection(la_county, clip_poly)


ggplot() +
  geom_sf(data = grid_sf, aes(color = L_dn), size = 0.5, alpha = 0.8) +
  geom_sf(data = la_county_cropped, fill = NA, color = "green", size = 1) +
  scale_color_viridis_c(option = "plasma") +
  coord_sf() +
  labs(
    title = "Noise Levels (Ldn) Over Los Angeles County",
    color = "Ldn (dB)"
  ) +
  theme_minimal()
