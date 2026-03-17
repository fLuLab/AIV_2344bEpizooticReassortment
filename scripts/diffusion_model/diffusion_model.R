################################################################################
## Script Name:        Diffusion Coefficient Model
## Purpose:           A Bayesian mixed model to determine associations between 
##                    the diffusion coefficient obtained from our phylogenetic 
##                    analysis and for each reassortant and 1. region of origin, 
##                    2. proportion of evolutionary time in wild birds, 
##                    3. persistence time and 4. total number of  species jumps
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
library(brms)
library(rstan)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(ggmcmc)
library(performance)
library(compositions)

################################### DATA #######################################
# Read and inspect data
combined_data <- read_csv('./beast_outputs_summary/2024-09-20_combined_data.csv')

reassortant_ancestral_changes <- read_csv('./beast_outputs_summary/reassortant_ancestral_changes.csv')

reassortant_stratifiedpersistence <- read_csv('./beast_outputs_summary/reassortant_stratifiedpersistence.csv')

meta <- read_csv('./data/2024-09-09_meta.csv') 

summary_data <- read_csv('./beast_outputs_summary/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
                   clade)) 

################################### MAIN #######################################
# Main analysis or transformation steps
# Data preprocessing

# Number and proportion of domestic and wild sequences per reassortant
wild_seqs <- meta %>% 
  drop_na(cluster_profile) %>%
  count(cluster_profile, host_isdomestic) %>% 
  rename(n_seq = n) %>%
  pivot_wider(names_from = host_isdomestic,
              values_from = n_seq,
              names_glue = 'n_{host_isdomestic}') %>%
  replace_na(list('n_domestic' = 0,
                  'n_wild' = 0)) %>%
  mutate(n_seq = n_domestic + n_wild,
         prop_wild = n_wild/n_seq)


# Summarise data across reassortants
diffusion_data <- combined_data %>%
  
  # select variables of interes
  dplyr::select(
    segment,
    cluster_profile,
    TMRCA,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    #persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  drop_na(collection_regionname) %>%
  
  rename_with(~gsub('-', '_', .x)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  
  mutate(collection_regionname = factor(collection_regionname, levels = c('asia', 'africa', 'europe', 'central & northern america'))) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  
  # join host richness
  left_join(summary_data %>% dplyr::select(c(cluster_profile, 
                                             host_richness)),
            by = join_by(cluster_profile)) %>%
  
  # join persistence proportions
  left_join(reassortant_stratifiedpersistence,
            by = c('segment', 'cluster_profile')) %>%
  rename_with(~gsub('-', '_', .x))%>%
  
  # season (breeding, migrating_spring, migrating_autumn, overwintering)
  mutate(collection_month = date_decimal(TMRCA) %>% format(., "%m") %>% as.integer(),
         collection_year = date_decimal(TMRCA) %>% format(., "%Y") %>% as.integer() %>% factor(levels = c(2017,2018,2019,2020,2021,2022,2023)),
         collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                                       collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                                       collection_month %in% c(6,7,8)  ~ 'breeding', 
                                       collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
         )) %>%
  
  # persistence time proportions
  #mutate(across(c('median_anseriformes_wild', 
  # 'median_charadriiformes_wild',
  # 'median_galliformes_domestic'), .fns= ~ ifelse(.x/persist.time > 1, 1, .x/persist.time ), .names = '{.col}_prop')) %>%
  
  
  # binary
  mutate(across(c('count_cross_species'), .fns= ~ log1p(.x), .names = '{.col}_log1p'))  %>%
  mutate(persistence = if_else(persistence < 5.5, persistence, 5.5)) %>%
  
  # Ancestral (WRT HA) Changes
  left_join(reassortant_ancestral_changes %>% dplyr::select(cluster_profile, segments_changed,cluster_class)) %>%
  replace_na(list(segments_changed = 0)) %>%
  
  dplyr::select(cluster_profile,
                weighted_diff_coeff, 
                persistence,
                collection_year,
                count_cross_species_log1p,
                path_prop_anseriformes_wild,
                path_prop_charadriiformes_wild,
                path_prop_galliformes_domestic,
                host_simplifiedhost,
                segments_changed,
                segment,
                cluster_class,
                collection_regionname) %>%
  mutate(path_prop_remainder = 1 - (path_prop_anseriformes_wild + path_prop_charadriiformes_wild + path_prop_galliformes_domestic )) %>%
  mutate(across(starts_with('path_prop'), .fns = ~ if_else(.x ==0, 1e-6, .x))) %>%
  mutate(across(starts_with('path_prop'), .fns = ~ .x/(path_prop_anseriformes_wild + path_prop_charadriiformes_wild + path_prop_galliformes_domestic + path_prop_remainder))) %>%
  filter(weighted_diff_coeff > 0) %>%
  filter(segment == 'ha') %>%
  drop_na()


# Define  contrast matrix defining the sequential binary partitions: 
# (i) (Anseriformes + Galliformes + Charadriiformes) / Other host
# (ii) (Anseriformes + Charadriiformes) / Galliformes
# (iii) Charadriiformes / Anseriformes

# row order = Anser, Charad, Gall, Other
m <- matrix(data=c(1, 1, -1,
                   1, 1, 1,
                   1, -1, 0,
                   -1, 0, 0),
            ncol=3,
            nrow=4, 
            byrow = T)


transformed_props <- diffusion_data %>%
  dplyr::select(starts_with('path_prop')) %>%
  ilr(V = gsi.buildilrBase(m)) %>%
  as_tibble()

diffusion_data %<>%
  mutate(collection_regionname = case_when(cluster_profile == '2_1_1_1_1_1_1_1' ~ 'europe', .default = collection_regionname)) %>%
  bind_cols(transformed_props)

diffusion_data%<>%
  left_join(wild_seqs) %>%
  mutate(n_seq = log10(n_seq))


###############################################################################
# Model Formula
diffusion_formula <- bf(weighted_diff_coeff ~ 1 +
                          n_seq * count_cross_species_log1p +
                          V1 +
                          V2  +
                          V3 + persistence + 
             
                          (1 + V2 + V3 || collection_regionname) +
                          (1 | collection_year),
                        shape ~  1 + (1 | collection_regionname) )



# Set Priors
diffusionmodel1_priors <- c(set_prior("normal(0, 2)", class = 'b'),
                            set_prior('normal(5,2)', class = 'Intercept', lb = 0),
                            set_prior('normal(0,2)', class = 'Intercept',  dpar = 'shape'),
                            set_prior('exponential(0.5)', class = 'sd'))


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472


# Fit model to data
diffusion_model<- brm(
  diffusion_formula,  
  data = diffusion_data,
  prior = diffusionmodel1_priors,
  family = Gamma(link = "log"),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  backend = "cmdstanr",
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95))


# Post-fitting checks (including inspection of ESS, Rhat and posterior predictive)
performance(diffusion_model)


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################