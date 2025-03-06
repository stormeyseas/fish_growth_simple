library(targets, warn.conflicts = F)
library(crew, warn.conflicts = F)
library(mirai, warn.conflicts = F)
library(here, warn.conflicts = F)
library(arrow, warn.conflicts = F)
library(sf, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(tidyr, warn.conflicts = F)
library(qs, warn.conflicts = F)
library(terra, warn.conflicts = F) # do not remove even though it looks unecessary
library(magrittr, warn.conflicts = F)
library(conflicted, warn.conflicts = F)
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "qs2", "here", "future", "furrr", "qs", "terra", "readxl"), 
  format = "qs", 
  controller = crew_controller_local(workers = 4, seconds_idle = 60),
  workspace_on_error = TRUE
)

tar_source(here("src/model_functions.R"))

# Paths ----------------------------------------------------------------------------------------------
this_species <- "atlantic_salmon"
this_path <- sprintf(here("data/%s"), this_species)
google_path <- "C:/Users/treimer/OneDrive - University of Tasmania/Sustainable Aquafeeds/Data from Google Drive"

# Begin targets pipeline -----------------------------------------------------------------------------
list(
  tar_target(prod_cycle_file, "data/_general_data/production_cycles/production_cycle.csv", format = "file"),
  # tar_target(day_num, seq(1:production_cycle)),
  # tar_target(yday_num, unique(ifelse(day_num > 365, day_num-365, day_num))),
  tar_target(farm_data_file, here("data/_general_data/SST/farm_SST_extracted.parquet"), format = "file"),
  tar_target(
    farm_IDs, 
    farm_data_file %>% 
      read_parquet() %>% 
      distinct(farm_id) %>% 
      slice_head(n = 250) %>% 
      as.vector() %>% unlist() %>% unname()
  ),
  tar_target(
    farm_data, 
    read_parquet(farm_data_file) %>% 
      mutate(day = str_split_i(day, "day_", 2) %>% as.numeric()) %>% 
      filter(farm_id == farm_IDs),
    pattern = farm_IDs,
    deployment = "main"
  ),
  tar_target(
    water_temp, 
    farm_data %>% select(day, temp_c),
    pattern = farm_IDs
  ),
  
  # Species parameters
  tar_target(
    production_cycle, 
    read.csv(prod_cycle_file) %>% 
      filter(species == this_species) %>% 
      rename(production_cycle_length = days) %>% 
      pull(production_cycle_length)
  ),
  tar_target(stocking_N, unique(farm_data$stocking_n), pattern = farm_data),
  
  tar_target(species_param_file, file.path(google_path,"Species Parameters.xlsx"), format = "file"),
  tar_target(species_params, get_species_params(species_param_file)),
  
  # Feed parameters
  tar_target(feed_types, c("reference", "past", "future")[1]),
  tar_target(feed_profile_file, file.path(google_path, "feed_profiles_and_composition.xlsx"), format = "file"),
  tar_target(
    feed_profile,
    read_excel(feed_profile_file)
  )
  
  
  

)
