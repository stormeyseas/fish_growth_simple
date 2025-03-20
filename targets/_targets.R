suppressMessages(suppressWarnings(library(targets)))
suppressMessages(suppressWarnings(library(crew)))
suppressMessages(suppressWarnings(library(mirai)))
suppressMessages(suppressWarnings(library(arrow)))
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(terra))) # do not remove
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(conflicted)))
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "future", "furrr", "ggplot2", "matrixStats"), 
  format = "qs", 
  controller = crew_controller_local(workers = 8),
  workspace_on_error = TRUE
)

tar_source("src/model_functions.R")

# Paths ----------------------------------------------------------------------------------------------
this_species <- "atlantic_salmon"
this_path <- file.path("data", this_species)
google_path <- "C:/Users/treimer/OneDrive - University of Tasmania/Sustainable Aquafeeds/Data from Google Drive"

# Begin targets pipeline -----------------------------------------------------------------------------
list(
  tar_target(prod_cycle_file, "data/_general_data/production_cycles/production_cycle.csv", format = "file"),
  # tar_target(day_num, seq(1:production_cycle)),
  # tar_target(yday_num, unique(ifelse(day_num > 365, day_num-365, day_num))),
  tar_target(farm_data_file, "data/_general_data/SST/farm_SST_extracted.parquet", format = "file"),
  tar_target(
    farm_IDs, 
    farm_data_file %>% 
      read_parquet() %>% 
      distinct(farm_id) %>% 
      # slice_sample(n = 250) %>%
      as.vector() %>% unlist() %>% unname()
  ),
  tar_target(
    farm_static_data, 
    read_parquet(farm_data_file) %>% 
      select(-c(day, temp_c)) %>% 
      filter(farm_id == farm_IDs) %>% 
      slice_head(n = 1),
    pattern = farm_IDs
  ),
  tar_target(
    farm_ts_data, 
    read_parquet(farm_data_file) %>% 
      select(c(farm_id, day, temp_c)) %>% 
      mutate(day = str_split_i(day, "day_", 2) %>% as.numeric()) %>% 
      filter(farm_id == farm_IDs),
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
  tar_target(times, c("t_start" = 1, "t_end" = 455, "dt" = 1)),
  tar_target(
    N_population, 
    generate_pop(
      N_seed = unique(farm_static_data$stocking_n), 
      mort = pop_params['mortmyt'],
      times = times
    ), 
    pattern = farm_static_data
  ),
  
  tar_target(species_param_file, file.path(google_path,"Species Parameters.xlsx"), format = "file"),
  tar_target(species_params, get_spec_params(species_param_file, "Atlantic salmon")),
  
  # Population parameters
  tar_target(pop_params_file, file.path(google_path,"Population parameters.xlsx"), format = "file"),
  # tar_target(pop_params, get_pop_params(pop_params_file, "Atlantic salmon")),
  tar_target(pop_params, c('meanW' = 1250, 'deltaW' = 100, 'Wlb' = 1.00e-04, 'meanImax' = 3.50e-02, 'deltaImax' = 3.50e-03, 'mortmyt' = 4.10e-04, 'nruns' = 5.00e+03)),
  
  # Feed parameters
  tar_target(feed_types, c("reference", "past", "future")),
  tar_target(feed_profile_file, "data/_general_data/diets/feed_profiles_and_composition_Atlantic_salmon.xlsx", format = "file"),
  tar_target(
    feed_params_protein, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("protein")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_protein, digest = ing_protein_digestibility), 
    pattern = feed_types
  ),
  tar_target(
    feed_params_carbs, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("carb")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_carb, digest = ing_carb_digestibility), 
    pattern = feed_types
  ),
  tar_target(
    feed_params_lipids, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("lipid")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_lipid, digest = ing_lipid_digestibility), 
    pattern = feed_types
  ),
  
  tar_target(
    example_individual,
    fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_ts_data$temp_c,
      feed_params = list(
        Proteins = feed_params_protein,
        Carbohydrates = feed_params_carbs,
        Lipids = feed_params_lipids
      ),
      times = times,
      init_weight = pop_params["meanW"],
      ingmax = pop_params["meanImax"]
    ),
    pattern = cross(map(feed_params_protein, feed_params_carbs, feed_params_lipids), map(farm_IDs, farm_ts_data)),
  ),
  
  tar_target(exind_weight, cbind(example_individual[, 1], example_individual[, 2], 
                                 matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                 matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_dw, cbind(example_individual[, 1], example_individual[, 3], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_water_temp, cbind(example_individual[, 1], example_individual[, 4], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_T_response, cbind(example_individual[, 1], example_individual[, 5], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_P_excr, cbind(example_individual[, 1], example_individual[, 6], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_L_excr, cbind(example_individual[, 1], example_individual[, 7], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_C_excr, cbind(example_individual[, 1], example_individual[, 8], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_P_uneat, cbind(example_individual[, 1], example_individual[, 9], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_L_uneat, cbind(example_individual[, 1], example_individual[, 10], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_C_uneat, cbind(example_individual[, 1], example_individual[, 11], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_food_prov, cbind(example_individual[, 1], example_individual[, 12], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_food_enc, cbind(example_individual[, 1], example_individual[, 13], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_rel_feeding, cbind(example_individual[, 1], example_individual[, 14], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_ing_pot, cbind(example_individual[, 1], example_individual[, 15], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_ing_act, cbind(example_individual[, 1], example_individual[, 16], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_E_assim, cbind(example_individual[, 1], example_individual[, 17], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_E_somat, cbind(example_individual[, 1], example_individual[, 18], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_anab, cbind(example_individual[, 1], example_individual[, 19], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_catab, cbind(example_individual[, 1], example_individual[, 20], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_O2, cbind(example_individual[, 1], example_individual[, 21], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  tar_target(exind_NH4, cbind(example_individual[, 1], example_individual[, 22], 
                                          matrix(farm_IDs, times['t_end'], 1, dimnames = list(NULL, 'farm_ID')), 
                                          matrix(feed_types, times['t_end'], 1, dimnames = list(NULL, 'feed'))), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),

  # This is the main farm growth target
  tar_target(
    main_farm_growth,
    farm_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_ts_data$temp_c,
      feed_params = list(
        Proteins = feed_params_protein,
        Carbohydrates = feed_params_carbs,
        Lipids = feed_params_lipids
      ),
      times = times,
      N_pop = N_population/2,
      nruns = 5000
    ),
    pattern = cross(
      map(feed_types, feed_params_protein, feed_params_carbs, feed_params_lipids),
      map(farm_ts_data, N_population)
    )
  ),

  tar_target(CommSize_file, "data/_general_data/harvest_sizes/all_harvest_sizes.csv", format = "file"),
  tar_target(CommSize, 4500),
             # read.csv(CommSize_file) %>% filter(species == this_species) %>% select(harvest_size_g) %>% as.numeric()),
  
  tar_target(
    days_to_CommSize,
    data.frame(
      farm_ID = farm_IDs,
      feed = feed_types,
      Ub = days_to_CS(weight_stat = (main_farm_growth[[1]][,1] + main_farm_growth[[1]][,2]), days = 1:455, CS = CommSize),
      Mean = days_to_CS(weight_stat = main_farm_growth[[1]][,1], days = 1:455, CS = CommSize),
      Lb = days_to_CS(weight_stat = (main_farm_growth[[1]][,1] - main_farm_growth[[1]][,2]), days = 1:455, CS = CommSize)
    ),
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs))
  )
  
)


