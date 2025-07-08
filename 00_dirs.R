# This sets up all the directories for the project so that all scripts use the same directory names
library(here)
library(stringr)
library(magrittr)

# Directory structure
output_path <- here() %>% file.path("outputs")

# Input paths
farms_locs_path <- file.path("data", "farm_locations")
farm_sst_path <- file.path("data", "SST")
input_species_param_path <- file.path("data", "atlantic_salmon_params")
feed_profile_path <- file.path("data", "diets")

# Output paths
output_farm_data_path <- file.path(output_path, "farm_data")
output_species_data_path <- file.path(output_path, "species_data")
output_sens_data_path <- file.path(output_path, "sensitivity_data")
output_growth_data_path <- file.path(output_path, "farm_growth_data")
output_cohorts_data_path <- file.path(output_path, "cohort_growth_data")
output_model_farm_path <- file.path(output_path, "all_outputs_farm")
output_model_cohort_path <- file.path(output_path, "all_outputs_cohort")
total_uneaten_path <- file.path(output_path, "total_uneaten_cohort")
total_excreted_path <- file.path(output_path, "total_excreted_cohort")

# Create output directories
dir.create(output_farm_data_path, recursive = T, showWarnings = F)
dir.create(output_species_data_path, showWarnings = F)
dir.create(output_sens_data_path, showWarnings = F)
dir.create(output_growth_data_path, showWarnings = F)
dir.create(output_cohorts_data_path, showWarnings = F)
dir.create(output_model_farm_path, showWarnings = F)
dir.create(output_model_cohort_path, showWarnings = F)
dir.create(total_uneaten_path, showWarnings = F)
dir.create(total_excreted_path, showWarnings = F)
