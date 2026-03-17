####################################################################################################
####################################################################################################
## Script name: Diffusion Model Interpretation and Figures
##
## Purpose of script: To extract, interpret and display the parameter values and predictions of the
## hurdle model fitted in ./scripts/diffusion_model.R
##
## Date created: 2025-01-22
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 
options('marginaleffects_posterior_interval' = 'hdi')
  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(tidybayes) 
library(marginaleffects) # missing
library(magrittr)
library(compositions)

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
  }

############################################## DATA ################################################
#diffusionmodel_fit <- readRDS('./saved_models/diffusion_model_2.rds')sw 
#diffusion_data <- read_csv('./saved_models/diffusion_model.csv')


############################################## MAIN ################################################
options("marginaleffects_posterior_interval" = "hdi")
options("marginaleffects_posterior_center" = "median")

# Posterior predictions marginalised over a balanced dataset
# Balanced = one row for each combination of unique values for the categorical 
# (or binary) predictors, holding numeric variables at their means
avg_predictions(diffusion_model, newdata = 'balanced') %>%
  as_tibble() %>%
  mutate(across(everything(), .fns = ~.x/365.25)) %>% as.numeric()


# Posterior predictions marginalised over the empirical data distribution, 
# stratified by Continent
avg_predictions(diffusion_model, by = 'collection_regionname' )%>%
as_tibble() %>%
  mutate(across(where(is.numeric), .fns = ~.x/365.25)) %>%
  view()


# Posterior predictions marginalised over a balanced dataset, stratified by 
# the number of host jumps
avg_predictions(diffusion_model, 
                variables = list('count_cross_species_log1p' = log1p(0:10)),
                newdata = 'balanced') %>%
  as_tibble() %>%
  mutate(across(where(is.numeric), .fns = ~.x/365.25))


# Posterior predictions marginalised over a balanced dataset, stratified by 
# the year of tMRCA
avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_year')


# Posterior predictions marginalised over a continent, stratified by proportions
# in host orders (using isomentric log ratio transform)

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

v3_preds <- expand_grid( 'path_prop_anseriformes_wild' = seq(0.1, 0.9, by = 0.1),
                         'path_prop_galliformes_domestic' = 1e-6,
                         'path_prop_remainder' = 1e-6) %>%
  
  mutate(path_prop_charadriiformes_wild = 1- path_prop_anseriformes_wild) %>%
  mutate(across(starts_with('path_prop'), 
                .fns = ~ .x/(path_prop_anseriformes_wild + 
                               path_prop_charadriiformes_wild + 
                               path_prop_galliformes_domestic + 
                               path_prop_remainder))) %>%
  
  dplyr::select(path_prop_anseriformes_wild,
                path_prop_charadriiformes_wild,
                path_prop_galliformes_domestic,
                path_prop_remainder)

v3_preds_transformed = v3_preds %>%
  ilr(V = gsi.buildilrBase(m)) %>%
  as_tibble()

v3_preds %<>% bind_cols(v3_preds_transformed) 

v3_predicted <- avg_predictions(diffusion_model, 
                                datagrid(V1 = v3_preds$V1,
                                         V2 = v3_preds$V2,
                                         V3 = v3_preds$V3,
                                         persistence = 1,
                                         collection_regionname = unique(diffusion_data$collection_regionname)),
                                by = c( 'V3')) %>%
  get_draws(shape = 'rvar') %>%
  as_tibble() %>%
  left_join(v3_preds, by = 'V3')  %>%
  mutate(across(any_of(c('estimate', 'conf.low', 'conf.high')), .fns = ~.x/365.25)) %>%
  view()


############################################## WRITE ###############################################

############################################## END #################################################
####################################################################################################
####################################################################################################