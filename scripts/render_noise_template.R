# /////////////////////////////////////////////////////////////////////////////
# Driver program to render 3 reports using 3 different year's worth of data
# Settings can be changed in the .qmd document to specify which data source
# to use
# 
# /////////////////////////////////////////////////////////////////////////////

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(furrr)

# Set the input clinics, output format parameters, and number of cores ---------
# by_vars <- c("by_timepoint")

by_vars <- c(2018, 2022, 2023)

format <- "html"

cores <- 3

# Create reports function ------------------------------------------------------
create_reports <- function(by_var) {

  # Relative to the root project directory, set the path to the master layout
  # .qmd file
  layout <-  "./scripts/noise_template.qmd"

  # Create a variable-specific .qmd file. Copy and rename the master layout
  # because the same file can't be read in multiple future sessions in parallel.

  modified_layout_file <- here("scripts", str_c(sub(" ", "", by_var), "_temp.qmd"))

  file.copy(
    from = layout,
    to = modified_layout_file,
    overwrite = TRUE
  )

  # Set file_in to the copied and renamed clinic-specific .qmd file. Serves as
  # an input to the quarto_render() function
  file_in <- modified_layout_file

  # Set file_out to the file name of the rendered report. Serves as a parameter
  # in the output-file option of the clinic-specific .qmd file.
  file_out <- str_c(by_var, ".", format)

  # Render the report, this will output the report to the root directory.
  quarto::quarto_render(
    input = file_in,
    execute_params = list(by_variable = by_var),
    output_format = "html",
    output_file = file_out
  )

  # The reports are initially placed in the root project directory. Copy the
  # output .html file to the reports directory
  file.copy(
    from = here("scripts", file_out),
    # to = here("deliverables/reports", str_c(by_var, ".", format)),
    to = here("scripts/noise_data_html_reports"),
    overwrite = TRUE
  )

  # Remove the output .html file from the root project directory. Uses a
  # relative path
  file.remove(here("scripts", file_out))

  # Remove the copy of the master layout .qmd file. Uses an absolute Path
  file.remove(file_in)

}

# Set furrr options ------------------------------------------------------------
options(future.rng.onMisuse = "ignore")
plan(multisession, workers = cores)

# Render reports in parallel ---------------------------------------------------
system.time(by_vars %>% future_walk(~ create_reports(.x)))

# Render reports in serial -----------------------------------------------------
# system.time(by_vars %>% walk(~ create_reports(.x)))