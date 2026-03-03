# This sets up all the directories for the project so that all scripts use the same directory names
library(tidyverse)
library(here)

output_path <- here() %>% file.path("outputs")
data_path <- here() %>% file.path("data")

# Create directories
dir.create(output_path, recursive = T, showWarnings = F)
dir.create(data_path, recursive = T, showWarnings = F)

fixnum <- function(n, digits = 4) {
  vapply(n, function(x) {
    str_flatten(c(rep("0", digits-nchar(as.character(x))), as.character(x)))
  }, character(1))
}
meanna <- function(x, ...) mean(x, na.rm = TRUE, ...)
minna <- function(x, ...) min(x, na.rm = TRUE, ...)
maxna <- function(x, ...) max(x, na.rm = TRUE, ...)
sdna <- function(x, ...) sd(x, na.rm = TRUE, ...)
sumna <- function(x, ...) sum(x, na.rm = TRUE, ...)
medianna <- function(x, ...) median(x, na.rm = TRUE, ...)
