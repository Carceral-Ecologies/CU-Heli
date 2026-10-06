# /////////////////////////////////////////////////////////////////////////////
# Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine
#
# Prepare a data sets for modeling the relationship between aircraft noise and
# lachs outcomes.

# Outputs two separate data sets one for 2018 and one for 2013
# /////////////////////////////////////////////////////////////////////////////

# Load Packages ---------------------------------------------------------------
library(tidyverse)


# For two years, 2018 and 2023

years <- c(2018, 2023)

prep_data <- function(year) {
  # Read noise data -------------------------------------------------------------
  noise <- readRDS(
    str_c(
      "D:\\rwjf_noise\\data\\noise\\lac_pw_db_", 
      year,
      ".rds"
    )
  )

  # Create a TractFIPS column for merging the cdc places data set to the noise 
  # data set
  noise <- noise %>%
    rename(TractFIPS = GEOID)

  # Load the LACHS data -------------------------------------------------------
  # The 2022 Places Data Set corresponds to BRFSS data from 2020, the year that 
  # the current average annual noise come from.
  fpath <- str_c("D:\\rwjf_noise\\data\\LACPH\\adult", year, "_ucamc.csv")

  if (year == 2023) {
    lachs <- read_csv(
    fpath, 
    show_col_types = FALSE,
    col_types = cols(
      AGE = col_character(),
      AGEGRP = col_character(),
      AGEUND = col_character()
      )
    ) 
  } else {
    lachs <- read_csv(
    fpath, 
    show_col_types = FALSE
  )
  }
   

  # Geocodes are 2010 tracts for 2018 survey and 2020 tracts for 2023 survey

  # 2018 data has 1,989 unique tracts, 2023 data has 2,163 tracts whereas noise
  # data is available for 2,342 and 2,493 for 2018 and 2023 data respectively
  lachs %>%
    select(GEO_CT) %>%
    summarise(n = n_distinct(GEO_CT))

  # Create GEOID by prepending 06 = CA, 037 = LA County to the GEO_CT variable
  # paste0 will prepend the GEOID to rows with missing values. These rows will
  # not get a noise variable merged in because there will not be a match for it
  # but explicitly setting as NA.
  lachs <- lachs %>%
  mutate(
    GEOID = paste0("06037", GEO_CT)
  ) %>%
  mutate(GEOID = ifelse(is.na(GEO_CT), NA, GEOID))

  # Merge the tract level population weight noise levels to the lachs
  lachs <- lachs %>%
    left_join(
      noise %>% select(TractFIPS, pw_exp_lin, pw_exp_db, pw_exp_db2),
      by = c("GEOID" = "TractFIPS")
    )

  
  # 2018 and 2023 have different column names. Harmonize to the 2023 names
  if (year == 2018) {
    lachs <- lachs %>%
      rename(
        HSPH = HS2PH, 
        HSMH = HS3MH, 
        HSALD = HS4ALD)
  } 

  # Set the Don't Know and Refused to Missing.
  lachs <- lachs %>%
      mutate(across(c(AGE, HSTOTDAY, PHQ, ASTHMA, DEPNOW, HSPH, HSMH, HSALD, HSGEN), ~ ifelse(.x %in% c("Do not know", "Refused"), NA, .x)))
  
  # Collapse Age group
  lachs <- lachs %>%
    mutate(AGEGROUP_c = 
      case_match(
        AGEGROUP,
        c("18-24", "25-29", "30-39") ~ "18-39",
        c("40-49", "50-59") ~ "40-59",
        .default = "60+")
      )

  # Set Refused to NA for GENDERO
  # 2018 and 2023 data have different columns for capturing gender 
  # (see code book)
  # 2018
  #   GENDERO - Original Birth Certificate
  #   GENDERC - Current
  #   GENGERN - DPH SOP
  # 2023 
  #   GENDERB - Assigned at birth
  #   GENDERC - Current
  #   GENGERN - DPH SOP
  
  # This creates a harmonized version of GENDERO and GENDERB and renames it to
  # GENDER, so that it can be programmed in SAS macros/loops.
  if (year == 2018) {
    lachs <- lachs %>%
      mutate(GENDER = GENDERO) %>%
      mutate(GENDER = ifelse(!GENDER %in% c("Male", "Female"), NA, GENDER))
  } else {
    lachs <- lachs %>%
      mutate(GENDER = GENDERB) %>%
      mutate(GENDER = ifelse(!GENDER %in% c("Male", "Female"), NA, GENDER))
  }


  # Collapse race
  # NHOPI is native hawaiian or pacific islander
  if (year == 2018) {
  lachs <- lachs %>%
    mutate(RACESOP_c = 
      case_match(
        RACESOP,
        c("AI/AN", "Do not know", "Multiracial/Other", "NHOPI", "Refused") ~ "Other/Unknown",
        .default = RACESOP)
      ) %>%
    mutate(RACESOP_c = factor(RACESOP_c, levels = c("White", "African American", "Asian", "Latino", "Other/Unknown")))
  } else {
    lachs <- lachs %>%
      mutate(RACESOP_c = 
        case_match(
          RACESOP,
          c("Non-Hispanic AI/AN", "Non-Hispanic Multi-Racial or Other", "Non-Hispanic NHPI", "Refused") ~ "Other/Unknown",
          "Non-Hispanic Black or African American" ~ "African American",
          "Non-Hispanic Asian" ~ "Asian",
          "Non-Hispanic White" ~ "White",
          "Latinx" ~ "Latino",
          .default = RACESOP)
        ) %>% 
      mutate(RACESOP_c = factor(RACESOP_c, levels = c("White", "African American", "Asian", "Latino", "Other/Unknown")))
    
  }
  # Convert to numeric
  # lachs <- lachs %>%
  #   mutate(across(HSPH:HSALD, ~ as.numeric(.x)))
    
  lachs <- lachs %>%
    mutate(PHQ_num = case_match(
      PHQ,
      "Yes" ~ 1,
      "No" ~ 0,
      .default = NA))
    
  lachs <- lachs %>%
    mutate(DEPNOW_num = case_match(
      DEPNOW,
      "Yes" ~ 1,
      "No" ~ 0,
      .default = NA))

  lachs <- lachs %>%
    mutate(ASTHMA_num = case_match(
      ASTHMA,
      "Yes" ~ 1,
      "No" ~ 0,
      .default = NA))
    
  lachs <- lachs %>%
    mutate(pw_exp_db_cen = pw_exp_db - mean(pw_exp_db, na.rm = TRUE))

  # Write out to a file for R
  saveRDS(lachs, str_c("D:\\rwjf_noise\\data\\prepped_data\\merged_lachs_", year, ".rds"))

  # Write out to a file for SAS
  # Remove geometry otherwise it will break the .csv file when importing to SAS
  lachs %>%
    write_csv(
      .,
      str_c("D:\\rwjf_noise\\data\\prepped_data\\merged_lachs_", year, ".csv"),
      na = ""
    )

}

# run the loop
walk(years, prep_data)