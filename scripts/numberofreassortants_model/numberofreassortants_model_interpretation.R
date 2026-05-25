####################################################################################################
####################################################################################################
## Script name: Count and Class Model Interpretation and Figures
##
## Purpose of script: To extract, interpret and display the parameter values and predictions of the
## hurdle model fitted in ./scripts/numberofreassortants_model.R
##
## Date created: 2025-01-22
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)
library(ggmcmc)

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}



############################################## DATA ################################################
# reads numbers_model_2 from ./numberofreassortants_model.R
                                             

############################################## MAIN ################################################
# interpretation of stan model


#1) continent stratified estimates of the observed number of reassortants per month
test_pred <- as_draws_df(numbers_model_2$draws('y_rep') ) %>%
  pivot_longer(cols = starts_with("y_rep["), names_to = "row", values_to = ".epred") %>%
  mutate(row = as.integer(str_extract(row, "\\d+"))) %>%
  left_join(data_processed_2 %>% rowid_to_column('row'),
            by = 'row')

test_pred %>%
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)


#1) continent stratified estimates of the true (latent) number of reassortants per month
test_pred_2 <- as_draws_df(numbers_model_2$draws('N_rep') ) %>%
  pivot_longer(cols = starts_with("N_rep["), names_to = "row", values_to = ".epred") %>%
  mutate(row = as.integer(str_extract(row, "\\d+"))) %>%
  left_join(data_processed_2 %>% rowid_to_column('row'),
            by = 'row')

test_pred_2 %>%
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)


rpois(100, exp(1.77 + -0.06697910 * 0.000000 + -0.03163115 ))

#2) ZI 
t <- get_variables(numbers_model_2)
# Trace plot
numbers_model_2 %>%
  gather_draws(., !!!syms(t)) %>%
  filter(grepl('theta', .variable)) %>%
  
  group_by(.draw, .variable) %>% 
  summarise(avg_value = mean(.value), .groups = "drop") %>% 
  group_by(.variable) %>%
  median_hdci(avg_value)


#2) Probability of detection


#3) Effect of number of sequences on p(detect)
inv_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

continent_specific_detection_samples <- numbers_model_2 %>%
  gather_draws(., continent_specific_detection[i]) 

beta_sequences_samples <-  numbers_model_2 %>%
  gather_draws(., beta_sequences) %>%
  ungroup() %>% 
  select(.draw,
         .chain, 
         .iteration, 
         beta_sequences = .value)

year_detection_samples <- numbers_model_2 %>%
  gather_draws(., year_detection[i]) %>% 
  select(.draw,
         .chain,
         .iteration,
         year_detection = .value)


detection_probabilities <- continent_specific_detection_samples %>%
  rename(continent_specific_detection = .value) %>%
  inner_join(beta_sequences_samples ,
             by = c(".draw", ".chain", ".iteration")) %>%
  inner_join(year_detection_samples, 
             by = c(".draw", ".chain", ".iteration"), 
             relationship = "many-to-many") %>%
  mutate(sequences = list(log1p(c(3, 12, 21, 30)))) %>%
  unnest(sequences) %>%
  mutate(
    probability = inv_logit(continent_specific_detection + beta_sequences * sequences + year_detection)
  )

detection_probabilities %>%
  group_by(sequences) %>% 
  median_hdci(probability, .width = 0.95)

#4) Effect of number of cases on N


#5) Grouped average across years



draws <- as_draws_df(numbers_model_2$draws())
newdata <- tibble(x = seq(-2, 2, length.out = 50))

# Add fitted draws
preds <- draws %>%
  select(.draw, alpha, beta) %>%
  crossing(newdata) %>%
  mutate(
    y_hat = alpha + beta * x
  )

avg_preds <- preds %>%
  group_by(x) %>%
  summarise(
    estimate = mean(y_hat),
    lower = quantile(y_hat, 0.025),
    upper = quantile(y_hat, 0.975)
  )

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################