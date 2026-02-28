# This sets up all the directories for the project so that all scripts use the same directory names
library(tidyverse)
library(here)

output_path <- here() %>% file.path("outputs")
data_path <- here() %>% file.path("data")

# Create directories
dir.create(output_path, recursive = T, showWarnings = F)
dir.create(data_path, recursive = T, showWarnings = F)
