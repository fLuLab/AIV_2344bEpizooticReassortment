################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2024-XX-XX
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(cmdstanr)
library(posterior)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)
library(ggmcmc)


################################### DATA #######################################
# Read and inspect data
data <- read_csv('./data/countmodeldata_2025Apr23.csv')

################################### MAIN #######################################
# Main analysis or transformation steps

# Data pre-processing
data_processed_2 <- data %>%
  
  # remove reassortant classes
  dplyr::select(-c(ends_with('class'), time_since_last_dominant)) %>%
  
  # scaling
  mutate(across(c(collection_dec, collection_year),
                .fns = ~ subtract(.x, 2015))) %>%
  
  # set default NA values
  replace_na(list(woah_cases = 0, 
                  woah_susceptibles = 0,
                  woah_deaths = 0,
                  n_sequences = 0,
                  n_reassortants = 0,
                  time_since_last_dominant = 0, # this may be quite a strong assumption 
                  minor = 0,
                  major = 0,
                  dominant = 0)) %>%
  mutate(across(starts_with('woah'), ~.x/6, .names = "{.col}_monthly")) %>%
  
  mutate(across(ends_with('_monthly'), ~log1p(.x), .names = "{.col}_log1p")) %>%
  mutate(n_sequences_log1p = log1p(n_sequences)) %>%
  #rename_with(~gsub('_', '-' ,.x)) %>%
  
  filter(collection_regionname != 'south america')


# Prepate Data
numbers_data <- list(N = nrow(data_processed_2),
                     y = data_processed_2 %>% pull(n_reassortants),
                     K = data_processed_2 %>% pull(n_reassortants) %>% max(),
                     C = data_processed_2 %>% pull(collection_regionname) %>% n_distinct(),
                     Y = data_processed_2 %>% pull(collection_year) %>% n_distinct(),
                     continent_index = data_processed_2 %>% pull(collection_regionname) %>% as.factor() %>% as.numeric(),
                     year_index = data_processed_2 %>% pull(collection_year) %>% as.factor() %>% as.numeric(),
                     cases =  data_processed_2 %>% pull(woah_susceptibles_monthly_log1p),
                     sequences =  data_processed_2 %>% pull(n_sequences_log1p))


# Compile Stan Model
numbers_mod <- cmdstan_model('./scripts/numberofreassortants_model/n_reassortants_final.stan')


# Set Priors - within n_reassortants_final.stan


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472

# Fit model to data
numbers_model_2 <- numbers_mod$sample(
  data = numbers_data,
  seed = SEED,
  chains = CHAINS,
  parallel_chains = CORES,
  iter_warmup = BURNIN,
  iter_sampling = ITER)


numbers_parms <- numbers_model_2$summary()


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################