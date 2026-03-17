################################################################################
## Script Name:         Class Model (using cumulative distribution)
## Purpose:            to model the probability that any given reassortant belongs 
##                     to class X, stratified by continent.
## Author:             James Baxter
## Date Created:       2025-05-13
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
library(broom)
library(broom.mixed)
library(brms)
library(cmdstanr)
library(igraph)

################################### DATA #######################################
# Read and inspect data
reassortant_ancestral_changes <- read_csv('./beast_outputs_summary/reassortant_ancestral_changes.csv')%>%
  mutate(across(ends_with('class'), .fns = ~ case_when(.x == 1 ~ 'moderate',
                                                       .x == 2 ~ 'minor',
                                                       .x == 3 ~ 'major')))

combined_data <- read_csv('./beast_outputs_summary/2024-09-20_combined_data.csv')
summary_data <- read_csv('./beast_outputs_summary/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
                   clade)) 
meta <- read_csv('./data/2024-09-09_meta.csv')


################################### MAIN #######################################
# Main analysis or transformation steps

my_edgelist<- reassortant_ancestral_changes %>% 
  dplyr::select(ends_with('label'), cluster_class) %>%
  drop_na() %>%
  relocate(parent_label, cluster_label) %>% 
  graph_from_data_frame(.)
plot(my_edgelist, vertex.size = 7)


most_recent_major <- distances(my_edgelist,
                               weights = reassortant_ancestral_changes %>%  
                                 dplyr::select(ends_with('label'), segments_changed) %>% 
                                 drop_na() %>% 
                                 pull(segments_changed),
                               mode = 'in',
                               to =  c( "H5N8/2019/R7_Africa",
                                        "H5N1/2020/R1_Europe",
                                        "H5N1/2021/R1_Europe",
                                        "H5N1/2021/R3_Europe" , 
                                        "H5N1/2022/R7_NAmerica" ,
                                        "H5N1/2022/R12_Europe")) %>%
  as_tibble(rownames = 'cluster_label') %>%
  rowwise() %>%
  filter(cluster_label != "H5N8/2019/R7_Africa") %>%
  mutate(across(starts_with('H5'), .fns = ~ case_when(.x==0 ~ NaN, 
                                                      abs(.x) == Inf ~ NaN,
                                                      .default = .x))) %>%
  mutate(most_recent_major = list(names(.)[2:7][which(c_across(starts_with('H5')) == min(c_across(starts_with('H5')), na.rm = T))])) %>%
  mutate(most_recent_major = paste(most_recent_major, collapse = '')) %>%
  mutate(last_major_label = na_if(most_recent_major, '')) %>%
  as_tibble() %>%
  dplyr::select(cluster_label, last_major_label) %>% 
  left_join(meta %>% dplyr::select(last_major_profile = cluster_profile, last_major_label = cluster_label) %>% distinct())


updated <- reassortant_ancestral_changes %>%
  
  # Add  origin continent of current reassortant
  left_join(combined_data %>% filter(segment == 'ha') %>% 
              dplyr::select(cluster_profile,
                            cluster_region = collection_regionname) %>% 
              group_by(cluster_profile) %>% 
              #slice_min(cluster_tmrca, n = 1) %>% 
              distinct()) %>%
  
  # Add tmrca and origin continent of 'parent' reassortant (with respect to HA)
  left_join(combined_data %>% 
              filter(segment == 'ha') %>% 
              dplyr::select(cluster_profile, 
                            parent_region = collection_regionname) %>% 
              group_by(cluster_profile) %>% 
              distinct(),
            by = join_by(parent_profile == cluster_profile), 
            relationship = 'many-to-many') %>%
  
  rename(cluster_tmrca = height_median,
         parent_tmrca = parent_height_median) %>%
  
  # Add 'last' (with respect to HA) major reassortant
  left_join(most_recent_major) %>%
  left_join(., dplyr::select(., last_major_profile = cluster_profile, 
                             last_major_region = cluster_region, 
                             last_major_tmrca = cluster_tmrca)) %>%
  
  mutate(across(ends_with('tmrca'), .fns = ~subtract(2024.21, .x))) %>%
  
  # calculate time to most recent major reassortant
  rowwise() %>%
  mutate(time_since_last_major =  cluster_tmrca-last_major_tmrca) %>%
  as_tibble() %>% 
  
  # Count the number of segments changed since last major reassortant
  rowwise() %>%
  mutate(segments_changed_last_major = sum(as.integer(str_split_fixed(cluster_profile, '_', 8)) != as.integer(str_split_fixed(last_major_profile, '_', 8)), na.rm = TRUE)) %>%
  as_tibble() %>%
  
  # Has RNP changed since last major?
  rowwise() %>%
  mutate(rnp_changed = if_else(all(!str_split_fixed(cluster_profile, '_', 8)[c(1,2,3)] %in% str_split_fixed(last_major_profile, '_', 8)[c(1,2,3)]), '1', '0')) %>%
  as_tibble() %>%
  
  
  # group continents
  mutate(across(ends_with('region'),
                .fns = ~case_when(grepl('europe', .x) ~ 'europe',
                                  grepl('africa', .x) ~ 'africa',
                                  grepl('asia', .x) ~ 'asia',
                                  grepl('(central|northern) america', .x) ~ 'central & northern america',
                                  .default = .x))) %>%
  
  # Does reassortant coincide with region change
  mutate(region_changed_from_previous = if_else(parent_region != cluster_region, '1', '0')) 

updated %<>%
  mutate(cluster_region = case_when(cluster_label == 'H5N1/2020/R1_Europe' ~ 'europe', .default = cluster_region))

class_data <- updated %>%
  dplyr::select(cluster_class,
                parent_class,
                cluster_region, 
                segments_changed,
                region_changed_from_previous,
                time_since_last_major,
                time_since_parent,
                rnp_changed) %>%
  
  mutate(cluster_class = ordered(cluster_class, levels = c('minor', 'moderate', 'major')))


# Formula
ordinal_formula_priors <- get_prior(cluster_class~1 + parent_class + cluster_region + s(time_since_last_major) + segments_changed,
                                    data = class_data,
                                    family = cumulative("probit")) 



# Set Priors
ordinal_formula_priors <- c(set_prior("normal(0, 2)", class = 'b'),
                            set_prior('exponential(0.5)', class = c('sds')))


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 2000
BURNIN <- ITER/10 # Discard0% burn in from each chain
SEED <- 4472


# Fit model to data
ordinal_model <- brm(cluster_class~1 + parent_class + cluster_region + s(time_since_last_major) + segments_changed,
                     data=class_data,
                     prior = ordinal_formula_priors,
                     family =cumulative("probit"),
                     chains = CHAINS,
                     cores = CORES, 
                     threads = 2, 
                     backend = "cmdstanr",
                     iter = ITER,
                     warmup = BURNIN,
                     seed = SEED)


# Post-fitting checks (including inspection of ESS, Rhat and posterior predictive)
performance(ordinal_model)


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################