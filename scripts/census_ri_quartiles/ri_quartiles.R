# /////////////////////////////////////////////////////////////////////////////
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
#
# Calculate measures of racial ethnic concentration
# For each tract, and for each racial group, calculated the difference between
# the proportion of the racial group in the tract and the proportion of the 
# racial group in the county. Using census population counts.
 
# Approach from Kitchesn et. al, https://www.mdpi.com/1660-4601/21/5/613
# For each group, calculate C = A - B where A and B are the proportions of a 
# specific group at the tract and county levels. The authors stated:
# "Positive values suggest that a community is relatively overrepresented in a
# tract, while negative values indicate underrepresentation. Values close to
# zero imply that a community’s representation in the tract is relatively
# proportional to its representation in the city’s overall population."



# Currently set to use the ACS 2020 data.

# /////////////////////////////////////////////////////////////////////////////

# Load Packages ---------------------------------------------------------------
library(tidyverse)
library(sf)
library(units)
library(tidycensus)
library(broom)

# Set options for tigris
options(tigris_progress_bar = FALSE)
options(tigris_use_cache = TRUE)

# Download ACS data code book for verifying variables -------------------------
# Use acs5 for B variables
# Use acs5/subject for S variables
# Use acs5/profile for D variables
code_book <- load_variables(2020, "acs5/profile", cache = TRUE)

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

# Year 2025 does not work
years <- c(2018, 2023)

for (year in years) {
  # Load the census data --------------------------------------------------------
  census <- get_acs(
    geography = "tract",
    variables = acs_vars,
    state = "CA",
    county = "Los Angeles",
    geometry = TRUE,
    output = "wide",
    survey = "acs5",
    year = year
  )

  # Define groups to use for calculating representation index quartiles
  groups <- c("latino", "white", "black", "native", "asian")

  # County-level counts and proportions (B)
  county_props <- census %>%
    sf::st_drop_geometry() %>%
    summarise(across(all_of(paste0(groups, "_count")), sum, na.rm = TRUE),
              total_pop = sum(total_pop, na.rm = TRUE)) %>%
    mutate(across(all_of(paste0(groups, "_count")),
                  ~ . / total_pop,
                  .names = "{gsub('_count', '_county_prop', .col)}"))

  # Compute representation index (A - B) and quartiles for each group
  # There are 20 tracts with missing values for all of the variables pulled.
  census_ri <- census %>%
    sf::st_drop_geometry() %>%
    select(GEOID, total_pop, all_of(paste0(groups, "_count"))) %>%
    mutate(across(
      all_of(paste0(groups, "_count")),
      ~ . / total_pop,
      .names = "{gsub('_count', '_tract_prop', .col)}"
    )) %>%
    mutate(
      latino_ri = latino_tract_prop - county_props$latino_county_prop,
      white_ri  = white_tract_prop  - county_props$white_county_prop,
      black_ri  = black_tract_prop  - county_props$black_county_prop,
      native_ri = native_tract_prop - county_props$native_county_prop,
      asian_ri  = asian_tract_prop  - county_props$asian_county_prop
    ) %>%
    mutate(across(
      ends_with("_ri"),
      ~ ntile(., 4),
      .names = "{gsub('_ri', '_ri_quartile', .col)}"
    )) %>%
    select(GEOID, ends_with("_ri_quartile"))

  # Join the segregation values
  census <- census %>% 
    left_join(census_ri, by = "GEOID")

  # Load Noise data
  noise <- readRDS(str_c("D:\\rwjf_noise\\data\\noise\\lac_pw_db_", year, ".rds"))

  # Create a data set
  data <- 
    census %>%
    sf::st_drop_geometry() %>%
    select(-NAME) %>%
    left_join(., noise, by = "GEOID")

  # Set quartile vars
  quartile_vars <- paste0(groups, "_ri_quartile")
  l_vars <- data %>% select(ends_with("_db")) %>% names()

  tab_out <- data.frame()

  # For each group and for each noise variable
  for (qvar in quartile_vars){

    for (l in l_vars) {
      
     sub_tab <- data %>%
       sf::st_drop_geometry() %>%
       select(!!rlang::sym(l), !!rlang::sym(qvar)) %>%
       mutate(!!rlang::sym(qvar) := factor(!!rlang::sym(qvar))) %>%
       group_by(!!rlang::sym(qvar)) %>%
       summarise(
         mean_L = mean(!!rlang::sym(l), na.rm = TRUE),
         sd_L = sd(!!rlang::sym(l), na.rm = TRUE),
         n = n()
       ) %>%
       ungroup() %>%
       drop_na() %>%
       mutate(
         se = sd_L / sqrt(n),
         lwr = mean_L - (1.96 * se),
         upr = mean_L + (1.96 * se)) %>%
       mutate(variable = l, group = str_c(qvar, "_", !!rlang::sym(qvar))) %>%
       select(variable, group, everything()) %>%
       select(-!!rlang::sym(qvar))
      
      
      tab_out <- bind_rows(tab_out, sub_tab)
      
    }
  }


  tab_out %>% 
    arrange(variable, group) %>% 
    write_csv(., str_c("D:\\rwjf_noise\\data\\representation_index\\all_values_", year, ".csv"))


  p <- tab_out %>%
    filter(variable == "pw_L_night_db") %>%
    mutate(
      quartile = str_sub(group, -1),
      group = str_remove(group, "_ri_quartile") %>% str_to_title()) %>%
    mutate(group = str_sub(group, 1, nchar(group) - 2)) %>%
    select(group, quartile, mean_L, lwr, upr) %>%
    ggplot(aes(x = quartile, y = mean_L, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    facet_wrap(~ group, nrow = 1) +
    labs(
      x = "Representation Quartile",
      y = "Mean L_night (dB)",
      title = str_c("Noise exposure by racial/ethnic representation quartile ", year)
    )
  
  ggsave(
    str_c("D:\\rwjf_noise\\data\\representation_index\\L_night_ri_", year, ".png"), 
    plot = p, 
    width = 8, 
    height = 6, dpi = 600)

}


#  ////////////////////////////////////////////////////////////////////////////
# Calculate the census over-representation quartiles
# census_seg <- census %>%
#   sf::st_drop_geometry() %>%
#   select(GEOID, total_pop, ends_with("_count")) %>%
#   mutate(across(
#     all_of(paste0(groups, "_count")),
#     ~ . / total_pop,
#     .names = "{.col}_pct"
#   )) %>%
#   mutate(across(
#     all_of(paste0(groups, "_count_pct")),
#     ~ ntile(., 4),
#     .names = "{gsub('_count_pct', '_seg_quartile', .col)}"
#   )) %>%
#   select(GEOID, ends_with("_seg_quartile"))

# # Summary
# census %>% 
#   sf::st_drop_geometry() %>% 
#   select(GEOID, ends_with("_seg_quartile")) %>% 
#   summary()


# quartile_vars <- paste0(groups, "_seg_quartile")

# models <- quartile_vars %>%
#   set_names() %>%
#   map(\(qvar) {
#     df <- data[, c("pw_exp_db", qvar)] %>% sf::st_drop_geometry()
#     df[[qvar]] <- factor(df[[qvar]])
#     lm(pw_exp_db ~ ., data = df)
#   })

# models %>% 
#   map(summary)


# # For each mod in model, set a group variable that is equal to the name of the
# # element in each tidy output. Once the tidy output has a group variable, the
# # data can be stacked, the intercept removed, and the quartile is extracted
# # from the term variable. Then the group variable is modified to strip the 
# # suffix and convert it to title case. Then plot.
# models %>%
#   imap(\(mod, nm) tidy(mod, conf.int = TRUE) %>% mutate(group = nm)) %>%
#   bind_rows() %>%
#   filter(term != "(Intercept)") %>%
#   mutate(
#     quartile = str_extract(term, "[0-9]+$"),
#     group = str_remove(group, "_seg_quartile") %>% str_to_title()
#   ) %>%
#   ggplot(aes(x = quartile, y = estimate, ymin = conf.low, ymax = conf.high)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_pointrange() +
#   facet_wrap(~ group, nrow = 1) +
#   labs(
#     x = "Concentration Quartile (ref = Q1)",
#     y = "Estimated difference in noise exposure (dB)",
#     title = "Noise exposure by racial/ethnic concentration quartile",
#     subtitle = "Coefficients from 5 separate OLS models; reference category is quartile 1"
#   )
