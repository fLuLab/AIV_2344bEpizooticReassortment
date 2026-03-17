################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-05-13
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)
options('marginaleffects_posterior_interval' = 'hdi')

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(ggraph)
library(tidygraph)
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes) # missing
library(bayesplot)
library(emmeans) # missing
library(marginaleffects) # missing
library(magrittr)
library(ggmcmc) # missing
library(bayestestR)
library(modelbased)


################################### DATA #######################################
# Read and inspect data

################################### MAIN #######################################
# Main analysis or transformation steps
# average posterior predictions 
avg_predictions(ordinal_model)
avg_predictions(ordinal_model, by = 'cluster_region') # empirical distribution (ie mean of predictions)


#percentage point change of reassortant class ~ previous class 
avg_comparisons(ordinal_model, variables = list("parent_class" = 'revpairwise'), 
                newdata = 'balanced')


# average prediction of time since last major (1, 2, 3yr interals)
avg_predictions(ordinal_model, variables = list('time_since_last_major' = c(0.5,1,3,5)))

# slops of time since last major (1, 2, 3yr interals)
ordinal_model %>%
  avg_comparisons(variables = list('time_since_last_major' = 1), 
                  by = "time_since_last_major",
                  newdata = datagrid(time_since_last_major = c(0.5, 1, 3),
                                     grid_type = 'counterfactual'), 
                  type = 'response')



# average prediction for 1,3 segmetns changing
avg_predictions(ordinal_model, variables = list('segments_changed' = c(1,2,4)))


# effect of eac hadditional segment
ordinal_model %>%
  avg_comparisons(variables = list('segments_changed' = 1), 
                  #by = "segments_changed",
                  type = 'response')


# contrasts for continent, all else equal
avg_comparisons(ordinal_model, variables = list("cluster_region" = 'pairwise'),
                newdata = 'balanced')

avg_predictions(ordinal_model, by = 'parent_class') 


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################