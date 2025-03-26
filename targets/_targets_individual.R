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
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "future", "furrr", "ggplot2", "matrixStats", "tibble"), 
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
  tar_target(farm_data_file, "data/_general_data/SST/farm_SST_extracted.parquet", format = "file"),
  tar_target(farm_coord_file, "data/_general_data/farm_locations/farm_coords.parquet", format = "file"),
  tar_target(farms_to_omit_file, sprintf("data/_general_data/farm_locations/%s_farms_to_omit.qs", this_species), format = "file"),
  tar_target(farms_to_omit, qs::qread(farms_to_omit_file)),
  tar_target(
    farm_IDs, 
    farm_data_file %>% 
      read_parquet() %>% 
      filter(!farm_id %in% farms_to_omit) %>% 
      distinct(farm_id) %>% 
      # slice_sample(n = 650) %>% # Only doing a subset of the farms because this is the individuals stuff, not necessary to do them all
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
  
  # starts 1 May in northern hemisphere, 1 October in southern hemisphere
  tar_target(times_N, c("t_start" = 121, "t_end" = 121+547, "dt" = 1)),
  tar_target(times_S, c("t_start" = 274, "t_end" = 274+547, "dt" = 1)),
  tar_target(
    farm_coords, 
    read_parquet(farm_coord_file) %>% 
      mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']),
             t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']))
  ),
  tar_target(farm_times, 
             c(farm_coords$t_start[farm_coords$farm_ID == farm_IDs], farm_coords$t_end[farm_coords$farm_ID == farm_IDs], dt = 1),
             pattern = farm_IDs),
  tar_target(farm_temp, farm_ts_data$temp_c[farm_times['t_start']:farm_times['t_end']], pattern = map(farm_ts_data, farm_times)),
  
  tar_target(species_param_file, file.path(google_path,"Species Parameters.xlsx"), format = "file"),
  tar_target(species_params, get_spec_params(species_param_file, "Atlantic salmon")),
  
  # Population parameters
  # tar_target(pop_params_file, file.path(google_path,"Population parameters.xlsx"), format = "file"),
  tar_target(pop_params, c('meanW' = 125, 'deltaW' = 10, 'Wlb' = 0.0001, 'meanImax' = 0.035, 'deltaImax' = 0.0035, 'mortmyt' = 0.00041, 'nruns' = 5000)),

  # Feed parameters
  tar_target(feed_types, c("reference", "past", "future")),
  tar_target(feed_profile_file, "data/_general_data/diets/feed_profiles_and_composition_Atlantic_salmon.xlsx", format = "file"),
  tar_target(
    feed_params_protein, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("protein")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_protein, digest = ing_protein_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  tar_target(
    feed_params_carbs, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("carb")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_carb, digest = ing_carb_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  tar_target(
    feed_params_lipids, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("lipid")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_lipid, digest = ing_lipid_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  
  tar_target(
    example_individual,
    fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_temp,
      feed_params = list(
        Proteins = feed_params_protein,
        Carbohydrates = feed_params_carbs,
        Lipids = feed_params_lipids
      ),
      times = farm_times,
      init_weight = pop_params["meanW"],
      ingmax = pop_params["meanImax"]
    ),
    pattern = cross(map(feed_params_protein, feed_params_carbs, feed_params_lipids), map(farm_IDs, farm_temp, farm_times)),
  ),
  
  tar_target(exind_weight, cbind(example_individual[, 1], example_individual[, 2]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, weight = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_SGR, 
             exind_weight %>% 
               mutate(SGR = 100 * (exp((log(weight)-log(exind_weight$weight[1]))/(days-exind_weight$days[1])) - 1)) %>% 
               relocate(SGR, .after = days) %>% select(-weight) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = exind_weight),
  
  tar_target(exind_dw, cbind(example_individual[, 1], example_individual[, 3]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, dw = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),

  tar_target(exind_TGC, 
             1000*((exind_weight$weight[nrow(exind_weight)])^(1/3) - (exind_weight$weight[1])^(1/3))/sum(farm_ts_data$temp_c[farm_times['t_start']:farm_times['t_end']]), 
             pattern = map(exind_weight, map(example_individual, cross(feed_types, map(farm_IDs, farm_times))))
             ),
  
  tar_target(exind_water_temp, cbind(example_individual[, 1], example_individual[, 4]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, water_temp = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_T_response, cbind(example_individual[, 1], example_individual[, 5]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, T_response = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_P_excr, cbind(example_individual[, 1], example_individual[, 6]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, P_excr = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_L_excr, cbind(example_individual[, 1], example_individual[, 7]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, L_excr = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_C_excr, cbind(example_individual[, 1], example_individual[, 8]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, C_excr = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_P_uneat, cbind(example_individual[, 1], example_individual[, 9]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, P_uneat = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_L_uneat, cbind(example_individual[, 1], example_individual[, 10]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, L_uneat = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_C_uneat, cbind(example_individual[, 1], example_individual[, 11]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, C_uneat = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_food_prov, cbind(example_individual[, 1], example_individual[, 12]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, food_prov = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_FCR, 
             cbind(example_individual[, 1], example_individual[, 12], example_individual[, 3]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types), FCR = V2/V3) %>% 
               relocate(FCR, .after = days) %>% select(-c(V2, V3)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(),
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_food_enc, cbind(example_individual[, 1], example_individual[, 13]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, food_enc = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_rel_feeding, cbind(example_individual[, 1], example_individual[, 14]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, rel_feeding = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_ing_pot, cbind(example_individual[, 1], example_individual[, 15]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, ing_pot = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_ing_act, cbind(example_individual[, 1], example_individual[, 16]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, ing_act = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_E_assim, cbind(example_individual[, 1], example_individual[, 17]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, E_assim = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_E_somat, cbind(example_individual[, 1], example_individual[, 18]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, E_somat = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_anab, cbind(example_individual[, 1], example_individual[, 19]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, anab = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_catab, cbind(example_individual[, 1], example_individual[, 20]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, catab = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_O2, cbind(example_individual[, 1], example_individual[, 21]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, O2 = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs))),
  
  tar_target(exind_NH4, cbind(example_individual[, 1], example_individual[, 22]) %>% 
               as.data.frame() %>% remove_rownames() %>% rename(days = V1, NH4 = V2) %>% 
               mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
               group_by(farm_ID, feed) %>% mutate(prod_days = days - min(days)+1) %>% ungroup(), 
             pattern = map(example_individual, cross(feed_types, farm_IDs)))
  
  ## Senstivities -------------------------------------------------------------------------------------------------
)


