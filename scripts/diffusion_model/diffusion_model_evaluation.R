####################################################################################################
####################################################################################################
## Script name: Diffusion Model Evaluation
##
## Purpose of script:
##
## Date created: 2025-04-11
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(scales)
library(DHARMa)
library(ggdist)
library(ggmcmc)


# User functions

######################################### DATA & MODEL #############################################
fitted_model <- diffusionmodel1_fit_gamma_19

diffusion_data <- diffusion_data
###################################### MCMC Diagnostics #############################################

#MCMC chain convergence (Visual Inspection)
# Plot Chains
ggs(diffusionmodel1_fit_gamma_19) %>% 
  from_ggmcmc_names() %>%
  filter(.iteration > 400) %>%
  mutate(label = case_when(.variable == 'b_collection_regionnameasia'~ "alpha['asia']",
                           .variable == 'b_collection_regionnameafrica'~ "alpha['africa']",
                           .variable == 'b_collection_regionnameeurope'~ "alpha['europe']",
                           .variable == "b_collection_regionnamecentral&northernamerica"~ "alpha['americas']",
                           
                           .variable == 'b_median_anseriformes_wild_prop'~ "beta['anseriformes']",
                           .variable == 'b_Imedian_anseriformes_wild_propEQEQ0TRUE'~ "beta['step-anseriformes']",
                           
                           .variable == 'b_median_charadriiformes_wild_prop'~ "beta['charadriiformes']",
                           .variable == 'b_Imedian_charadriiformes_wild_propEQEQ0TRUE'~ "beta['step-charadriiformes']",
                           
                           #.variable == 'b_int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                           .variable == 'b_count_cross_species_log1p'~ "beta['hostjump']",
                           .variable == 'b_persist.time_log1p'~ "beta['persistence']",
                           
                          # .variable == 'b_collection_regionnameafrica:median_anseriformes_wild_prop'~ "gamma['africa-anser']",
                           #.variable == 'b_collection_regionnameeurope:median_anseriformes_wild_prop'~ "gamma['europe-anser']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:median_anseriformes_wild_prop'~ "gamma['americas-anser']",
                           
                           #.variable == 'b_collection_regionnameafrica:median_charadriiformes_wild_prop'~ "gamma['africa-charad']",
                          # .variable == 'b_collection_regionnameeurope:median_charadriiformes_wild_prop'~ "gamma['europe-charad']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:median_charadriiformes_wild_prop'~ "gamma['americas-charad']",
                           
                           #.variable == 'b_collection_regionnameafrica:persist.time_log1p'~ "gamma['africa-persist']",
                          # .variable == 'b_collection_regionnameasia:persist.time_log1p'~ "gamma['asia-persist']",
                           #.variable == 'b_collection_regionnameeurope:persist.time_log1p'~ "gamma['europe-persist']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:persist.time_log1p'~ "gamma['americas-persist']",
                          
                          .variable == "sd_collection_regionname__Intercept" ~ "sigma[continent]",                                                                           
                          .variable == "sd_collection_year__Intercept" ~ "sigma[year]",  
                          
                          #.variable ==  "r_collection_regionname[asia,Intercept]"~ "gamma[asia]",                                                                         
                          #.variable ==  "r_collection_regionname[africa,Intercept]"~ "gamma[africa]",                                                                       
                          #.variable ==  "r_collection_regionname[europe,Intercept]"~ "gamma[europe]",                                                                       
                          #.variable ==  "r_collection_regionname[central.&.northern.america,Intercept]"~ "gamma[americas]",                                                   
                          #.variable ==  "r_collection_regionname[asia,median_anseriformes_wild_prop]"~ "gamma[asia-anseriformes]",                                                     
                          #.variable ==  "r_collection_regionname[africa,median_anseriformes_wild_prop]"~ "gamma[africa-anseriformes]",                                                   
                          #.variable ==  "r_collection_regionname[europe,median_anseriformes_wild_prop]"~ "gamma[europe-anseriformes]",                                                   
                          #.variable ==  "r_collection_regionname[central.&.northern.america,median_anseriformes_wild_prop]"~ "gamma[americas-anseriformes]",                               
                          #.variable ==  "r_collection_regionname[asia,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[asia-step-anseriformes]",                                           
                          #.variable ==  "r_collection_regionname[africa,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-anseriformes]",                                         
                          #.variable ==  "r_collection_regionname[europe,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-anseriformes]",                                         
                          #.variable ==  "r_collection_regionname[central.&.northern.america,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[americas-step-anseriformes]",                     
                          #.variable ==  "r_collection_regionname[asia,median_charadriiformes_wild_prop]"~ "gamma[asia-charadriiformes]",                                                  
                          #.variable ==  "r_collection_regionname[africa,median_charadriiformes_wild_prop]"~ "gamma[africa-charadriiformes]",                                                
                          #.variable ==  "r_collection_regionname[europe,median_charadriiformes_wild_prop]"~ "gamma[europe-charadriiformes]",                                                
                          #.variable ==  "r_collection_regionname[central.&.northern.america,median_charadriiformes_wild_prop]"~ "gamma[americas-charadriiformes]",                            
                          #.variable ==  "r_collection_regionname[asia,Imedian_charadriiformes_wild_propEQEQ0TRUE]" ~ "gamma[asia-step-charadriiformes]",                                       
                          #.variable ==  "r_collection_regionname[africa,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-charadriiformes]",                                      
                          #.variable ==  "r_collection_regionname[europe,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-charadriiformes]",                                      
                          #.variable ==  "r_collection_regionname[central.&.northern.america,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~  "gamma[americas-step-charadriiformes]",                  
                          #.variable ==  "r_collection_regionname[asia,persist.time_log1p]" ~ "gamma[asia-persist]",                                                               
                          #.variable ==  "r_collection_regionname[africa,persist.time_log1p]"~ "gamma[africa-persist]",                                                              
                          #.variable ==  "r_collection_regionname[europe,persist.time_log1p]"~ "gamma[europe-persiste]",                                                              
                          #.variable ==  "r_collection_regionname[central.&.northern.america,persist.time_log1p]" ~ "gamma[americas-persist]",                                         
                          #.variable ==  "r_collection_year[2017,Intercept]" ~ "delta[2017]",                                                                              
                          #.variable ==  "r_collection_year[2018,Intercept]"~ "delta[2018]",                                                                               
                          #.variable ==  "r_collection_year[2019,Intercept]"~ "delta[2019]",                                                                               
                          #.variable ==  "r_collection_year[2020,Intercept]"~ "delta[2020]",                                                                               
                          #.variable ==  "r_collection_year[2021,Intercept]"~ "delta[2021]",                                                                               
                          #.variable ==  "r_collection_year[2022,Intercept]"~ "delta[2022]",                                                                               
                          #.variable ==  "r_collection_year[2023,Intercept]"~ "delta[2023]",      
                           
                           
                           .variable == 'b_shape_collection_regionnameasia'~ "rho['asia']",
                           .variable == 'b_shape_collection_regionnameafrica'~ "rho['africa']",
                           .variable == 'b_shape_collection_regionnameeurope'~ "rho['europe']",
                           .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "rho['americas']",
                          .variable == 'lprior' ~ 'prior',
                          .variable == 'lp__' ~ 'log~probability')) %>%
  drop_na(label) %>%
  ggplot(aes(x = .iteration,
             y = .value, 
             col = as.factor(.chain)))+
  geom_line(alpha = 0.5)+
  facet_wrap(~ label,
             ncol = 2,
             scale  = 'free_y',
             strip.position = 'left',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  scale_x_continuous('Iteration') + 
  scale_y_continuous('Parameter Value') + 
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7))

ggsave('~/Downloads/flu_plots/diffusion_trace.jpeg',
       dpi = 360,
       height = 25,
       width = 17,
       units = 'cm')

# Plot ranked traces
# If chains are exploring the same space efficiently, the traces should be similar to one another 
# and largely overlapping.
as_draws_df(diffusionmodel1_fit_gamma_19) %>% 
  
  # Selecting only beta coefficients
  mcmc_rank_overlay(regex_pars = '^b_') %>%
  
  # Extract data
  .$data %>%
  rename(.variable = parameter) %>%
  mutate(label = case_when(.variable == 'b_collection_regionnameasia'~ "alpha['asia']",
                           .variable == 'b_collection_regionnameafrica'~ "alpha['africa']",
                           .variable == 'b_collection_regionnameeurope'~ "alpha['europe']",
                           .variable == "b_collection_regionnamecentral&northernamerica"~ "alpha['americas']",
                           
                           .variable == 'b_median_anseriformes_wild_prop'~ "beta['anseriformes']",
                           .variable == 'b_Imedian_anseriformes_wild_propEQEQ0TRUE'~ "beta['step-anseriformes']",
                           
                           .variable == 'b_median_charadriiformes_wild_prop'~ "beta['charadriiformes']",
                           .variable == 'b_Imedian_charadriiformes_wild_propEQEQ0TRUE'~ "beta['step-charadriiformes']",
                           
                           #.variable == 'b_int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                           .variable == 'b_count_cross_species_log1p'~ "beta['hostjump']",
                           .variable == 'b_persist.time_log1p'~ "beta['persistence']",
                           
                           # .variable == 'b_collection_regionnameafrica:median_anseriformes_wild_prop'~ "gamma['africa-anser']",
                           #.variable == 'b_collection_regionnameeurope:median_anseriformes_wild_prop'~ "gamma['europe-anser']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:median_anseriformes_wild_prop'~ "gamma['americas-anser']",
                           
                           #.variable == 'b_collection_regionnameafrica:median_charadriiformes_wild_prop'~ "gamma['africa-charad']",
                           # .variable == 'b_collection_regionnameeurope:median_charadriiformes_wild_prop'~ "gamma['europe-charad']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:median_charadriiformes_wild_prop'~ "gamma['americas-charad']",
                           
                           #.variable == 'b_collection_regionnameafrica:persist.time_log1p'~ "gamma['africa-persist']",
                           # .variable == 'b_collection_regionnameasia:persist.time_log1p'~ "gamma['asia-persist']",
                           #.variable == 'b_collection_regionnameeurope:persist.time_log1p'~ "gamma['europe-persist']",
                           #.variable == 'b_collection_regionnamecentral&northernamerica:persist.time_log1p'~ "gamma['americas-persist']",
                           
                           .variable == "sd_collection_regionname__Intercept" ~ "sigma[continent]",                                                                           
                           .variable == "sd_collection_year__Intercept" ~ "sigma[year]",  
                           
                           #.variable ==  "r_collection_regionname[asia,Intercept]"~ "gamma[asia]",                                                                         
                           #.variable ==  "r_collection_regionname[africa,Intercept]"~ "gamma[africa]",                                                                       
                           #.variable ==  "r_collection_regionname[europe,Intercept]"~ "gamma[europe]",                                                                       
                           #.variable ==  "r_collection_regionname[central.&.northern.america,Intercept]"~ "gamma[americas]",                                                   
                           #.variable ==  "r_collection_regionname[asia,median_anseriformes_wild_prop]"~ "gamma[asia-anseriformes]",                                                     
                           #.variable ==  "r_collection_regionname[africa,median_anseriformes_wild_prop]"~ "gamma[africa-anseriformes]",                                                   
                           #.variable ==  "r_collection_regionname[europe,median_anseriformes_wild_prop]"~ "gamma[europe-anseriformes]",                                                   
                           #.variable ==  "r_collection_regionname[central.&.northern.america,median_anseriformes_wild_prop]"~ "gamma[americas-anseriformes]",                               
                           #.variable ==  "r_collection_regionname[asia,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[asia-step-anseriformes]",                                           
                           #.variable ==  "r_collection_regionname[africa,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-anseriformes]",                                         
                           #.variable ==  "r_collection_regionname[europe,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-anseriformes]",                                         
                           #.variable ==  "r_collection_regionname[central.&.northern.america,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[americas-step-anseriformes]",                     
                           #.variable ==  "r_collection_regionname[asia,median_charadriiformes_wild_prop]"~ "gamma[asia-charadriiformes]",                                                  
                           #.variable ==  "r_collection_regionname[africa,median_charadriiformes_wild_prop]"~ "gamma[africa-charadriiformes]",                                                
                           #.variable ==  "r_collection_regionname[europe,median_charadriiformes_wild_prop]"~ "gamma[europe-charadriiformes]",                                                
                           #.variable ==  "r_collection_regionname[central.&.northern.america,median_charadriiformes_wild_prop]"~ "gamma[americas-charadriiformes]",                            
                           #.variable ==  "r_collection_regionname[asia,Imedian_charadriiformes_wild_propEQEQ0TRUE]" ~ "gamma[asia-step-charadriiformes]",                                       
                           #.variable ==  "r_collection_regionname[africa,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-charadriiformes]",                                      
                           #.variable ==  "r_collection_regionname[europe,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-charadriiformes]",                                      
                           #.variable ==  "r_collection_regionname[central.&.northern.america,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~  "gamma[americas-step-charadriiformes]",                  
                           #.variable ==  "r_collection_regionname[asia,persist.time_log1p]" ~ "gamma[asia-persist]",                                                               
                           #.variable ==  "r_collection_regionname[africa,persist.time_log1p]"~ "gamma[africa-persist]",                                                              
                           #.variable ==  "r_collection_regionname[europe,persist.time_log1p]"~ "gamma[europe-persiste]",                                                              
                           #.variable ==  "r_collection_regionname[central.&.northern.america,persist.time_log1p]" ~ "gamma[americas-persist]",                                         
                           #.variable ==  "r_collection_year[2017,Intercept]" ~ "delta[2017]",                                                                              
                           #.variable ==  "r_collection_year[2018,Intercept]"~ "delta[2018]",                                                                               
                           #.variable ==  "r_collection_year[2019,Intercept]"~ "delta[2019]",                                                                               
                           #.variable ==  "r_collection_year[2020,Intercept]"~ "delta[2020]",                                                                               
                           #.variable ==  "r_collection_year[2021,Intercept]"~ "delta[2021]",                                                                               
                           #.variable ==  "r_collection_year[2022,Intercept]"~ "delta[2022]",                                                                               
                           #.variable ==  "r_collection_year[2023,Intercept]"~ "delta[2023]",      
                           
                           
                           .variable == 'b_shape_collection_regionnameasia'~ "rho['asia']",
                           .variable == 'b_shape_collection_regionnameafrica'~ "rho['africa']",
                           .variable == 'b_shape_collection_regionnameeurope'~ "rho['europe']",
                           .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "rho['americas']",
                           .variable == 'lprior' ~ 'prior',
                           .variable == 'lp__' ~ 'log~probability')) %>%
  
  # Plot 
  ggplot() + 
  geom_step(aes(x = bin_start, 
                y = n,
                colour = chain),
            linewidth = 0.8) + 
  scale_y_continuous(expand = c(0,0), limits = c(150, 225)) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~label,  ncol = 3, labeller = label_parsed) +
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.title = element_blank())


#MCMC chain resolution (ESS)
model_resolution <- list(tidy(diffusionmodel1_fit_gamma_19, ess = T, rhat = T, effects = 'ran_vals') %>%
                           dplyr::select( ess, rhat, group, level) %>%
                           unite(term, group, level, sep = '_'),
                         tidy(diffusionmodel1_fit_gamma_19, ess = T, rhat = T) %>%  dplyr::select(term, ess, rhat,term) ) %>%
  bind_rows() %>%
  
  mutate(term = case_when(term == 'collection_regionnameasia'~ "beta['asia']",
                          term == 'collection_regionnameafrica'~ "beta['africa']",
                          term == 'collection_regionnameeurope'~ "beta['europe']",
                          term == "collection_regionnamecentral&northernamerica"~ "beta['americas']",
                           
                          term == 'int_stepmedian_anseriformes_wild_log1p'~ "beta['step_anseriformes']",
                          term == 'median_anseriformes_wild_log1p'~ "beta['persist_anseriformes']",
                          term == 'int_stepmedian_charadriiformes_wild_log1p'~ "beta['step_charadriiformes']",
                          term == 'median_charadriiformes_wild_log1p'~ "beta['persist_charadriiformes']",
                          term == 'int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                          term == 'count_cross_species_log1p'~ "beta['num_hostjump']",
                           
                          term == 'collection_regionnameafrica:median_anseriformes_wild_log1p'~ "beta['africa-anser']",
                          term == 'collection_regionnameeurope:median_anseriformes_wild_log1p'~ "beta['europe-anser']",
                          term == 'collection_regionnamecentral&northernamerica:median_anseriformes_wild_log1p'~ "beta['americas-anser']",
                           
                          term == 'collection_regionnameafrica:median_charadriiformes_wild_log1p'~ "beta['africa-charad']",
                          term == 'collection_regionnameeurope:median_charadriiformes_wild_log1p'~ "beta['europe-charad']",
                          term == 'collection_regionnamecentral&northernamerica:median_charadriiformes_wild_log1p'~ "beta['americas-charad']",
                           
                          term == 'shape_collection_regionnameasia'~ "alpha['asia']",
                          term == 'shape_collection_regionnameafrica'~ "alpha['africa']",
                          term == 'shape_collection_regionnameeurope'~ "alpha['europe']",
                          term == 'shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                          
                          term == 'sd__(Intercept)'~ "sigma[gamma]",
                          term == 'sd__shape_(Intercept)'~ "sigma[zeta]",
                          
                          term == 'segment_ha'~ "gamma['ha']",
                          term == 'segment_mp'~ "gamma['mp']",
                          term == 'segment_np'~ "gamma['np']",
                          term == 'segment_ns'~ "gamma['ns']",
                          term == 'segment_nx'~ "gamma['nx']",
                          term == 'segment_pa'~ "gamma['pa']",
                          term == 'segment_pb1'~ "gamma['pb1']",
                          term == 'segment_pb2'~ "gamma['pb2']",
                          term == 'collection_regionname__shape_asia'~ "zeta['asia']",
                          term == 'collection_regionname__shape_africa'~ "zeta['africa']",
                          term == 'collection_regionname__shape_europe'~ "zeta['europe']",
                          term == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']")) 
  
 

model_resolution %>% 
  ggplot(aes(x = term, y = ess)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels= label_parse(), 'Parameter', expand = c(0.05, 0)) + 
  scale_y_continuous(expand = c(0,0), 'Effective Sample Size')


# NB - technically a convergence check
model_resolution %>% 
  ggplot(aes(x = term, y = rhat)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels= label_parse(), 'Parameter', expand = c(0.05, 0)) + 
  scale_y_continuous(expand = c(0,0), expression(paste("Potential Reduction in Scale Factor (", hat(R), ')')))


#MCMC Autocorrelation
diffusionmodel1_fit_gamma_19 %>% 
  mcmc_acf() %>% 
  .$data %>%
  as_tibble() %>%
  mutate(Parameter = case_when(Parameter == 'b_collection_regionnameasia'~ "alpha['asia']",
                               Parameter == 'b_collection_regionnameafrica'~ "alpha['africa']",
                               Parameter == 'b_collection_regionnameeurope'~ "alpha['europe']",
                               Parameter == "b_collection_regionnamecentral&northernamerica"~ "alpha['americas']",
                               
                               Parameter == 'b_median_anseriformes_wild_prop'~ "beta['anseriformes']",
                               Parameter == 'b_Imedian_anseriformes_wild_propEQEQ0TRUE'~ "beta['step-anseriformes']",
                               
                               Parameter == 'b_median_charadriiformes_wild_prop'~ "beta['charadriiformes']",
                               Parameter == 'b_Imedian_charadriiformes_wild_propEQEQ0TRUE'~ "beta['step-charadriiformes']",
                               
                               #Parameter == 'b_int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                               Parameter == 'b_count_cross_species_log1p'~ "beta['hostjump']",
                               Parameter == 'b_persist.time_log1p'~ "beta['persistence']",
                               
                               # Parameter == 'b_collection_regionnameafrica:median_anseriformes_wild_prop'~ "gamma['africa-anser']",
                               #Parameter == 'b_collection_regionnameeurope:median_anseriformes_wild_prop'~ "gamma['europe-anser']",
                               #Parameter == 'b_collection_regionnamecentral&northernamerica:median_anseriformes_wild_prop'~ "gamma['americas-anser']",
                               
                               #Parameter == 'b_collection_regionnameafrica:median_charadriiformes_wild_prop'~ "gamma['africa-charad']",
                               # Parameter == 'b_collection_regionnameeurope:median_charadriiformes_wild_prop'~ "gamma['europe-charad']",
                               #Parameter == 'b_collection_regionnamecentral&northernamerica:median_charadriiformes_wild_prop'~ "gamma['americas-charad']",
                               
                               #Parameter == 'b_collection_regionnameafrica:persist.time_log1p'~ "gamma['africa-persist']",
                               # Parameter == 'b_collection_regionnameasia:persist.time_log1p'~ "gamma['asia-persist']",
                               #Parameter == 'b_collection_regionnameeurope:persist.time_log1p'~ "gamma['europe-persist']",
                               #Parameter == 'b_collection_regionnamecentral&northernamerica:persist.time_log1p'~ "gamma['americas-persist']",
                               
                               Parameter == "sd_collection_regionname__Intercept" ~ "sigma[continent]",                                                                           
                               Parameter == "sd_collection_year__Intercept" ~ "sigma[year]",  
                               
                               #Parameter ==  "r_collection_regionname[asia,Intercept]"~ "gamma[asia]",                                                                         
                               #Parameter ==  "r_collection_regionname[africa,Intercept]"~ "gamma[africa]",                                                                       
                               #Parameter ==  "r_collection_regionname[europe,Intercept]"~ "gamma[europe]",                                                                       
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,Intercept]"~ "gamma[americas]",                                                   
                               #Parameter ==  "r_collection_regionname[asia,median_anseriformes_wild_prop]"~ "gamma[asia-anseriformes]",                                                     
                               #Parameter ==  "r_collection_regionname[africa,median_anseriformes_wild_prop]"~ "gamma[africa-anseriformes]",                                                   
                               #Parameter ==  "r_collection_regionname[europe,median_anseriformes_wild_prop]"~ "gamma[europe-anseriformes]",                                                   
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,median_anseriformes_wild_prop]"~ "gamma[americas-anseriformes]",                               
                               #Parameter ==  "r_collection_regionname[asia,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[asia-step-anseriformes]",                                           
                               #Parameter ==  "r_collection_regionname[africa,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-anseriformes]",                                         
                               #Parameter ==  "r_collection_regionname[europe,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-anseriformes]",                                         
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,Imedian_anseriformes_wild_propEQEQ0TRUE]"~ "gamma[americas-step-anseriformes]",                     
                               #Parameter ==  "r_collection_regionname[asia,median_charadriiformes_wild_prop]"~ "gamma[asia-charadriiformes]",                                                  
                               #Parameter ==  "r_collection_regionname[africa,median_charadriiformes_wild_prop]"~ "gamma[africa-charadriiformes]",                                                
                               #Parameter ==  "r_collection_regionname[europe,median_charadriiformes_wild_prop]"~ "gamma[europe-charadriiformes]",                                                
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,median_charadriiformes_wild_prop]"~ "gamma[americas-charadriiformes]",                            
                               #Parameter ==  "r_collection_regionname[asia,Imedian_charadriiformes_wild_propEQEQ0TRUE]" ~ "gamma[asia-step-charadriiformes]",                                       
                               #Parameter ==  "r_collection_regionname[africa,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[africa-step-charadriiformes]",                                      
                               #Parameter ==  "r_collection_regionname[europe,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~ "gamma[europe-step-charadriiformes]",                                      
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,Imedian_charadriiformes_wild_propEQEQ0TRUE]"~  "gamma[americas-step-charadriiformes]",                  
                               #Parameter ==  "r_collection_regionname[asia,persist.time_log1p]" ~ "gamma[asia-persist]",                                                               
                               #Parameter ==  "r_collection_regionname[africa,persist.time_log1p]"~ "gamma[africa-persist]",                                                              
                               #Parameter ==  "r_collection_regionname[europe,persist.time_log1p]"~ "gamma[europe-persiste]",                                                              
                               #Parameter ==  "r_collection_regionname[central.&.northern.america,persist.time_log1p]" ~ "gamma[americas-persist]",                                         
                               #Parameter ==  "r_collection_year[2017,Intercept]" ~ "delta[2017]",                                                                              
                               #Parameter ==  "r_collection_year[2018,Intercept]"~ "delta[2018]",                                                                               
                               #Parameter ==  "r_collection_year[2019,Intercept]"~ "delta[2019]",                                                                               
                               #Parameter ==  "r_collection_year[2020,Intercept]"~ "delta[2020]",                                                                               
                               #Parameter ==  "r_collection_year[2021,Intercept]"~ "delta[2021]",                                                                               
                               #Parameter ==  "r_collection_year[2022,Intercept]"~ "delta[2022]",                                                                               
                               #Parameter ==  "r_collection_year[2023,Intercept]"~ "delta[2023]",      
                               
                               
                               Parameter == 'b_shape_collection_regionnameasia'~ "rho['asia']",
                               Parameter == 'b_shape_collection_regionnameafrica'~ "rho['africa']",
                               Parameter == 'b_shape_collection_regionnameeurope'~ "rho['europe']",
                               Parameter == 'b_shape_collection_regionnamecentral&northernamerica'~ "rho['americas']",
                               Parameter == 'lprior' ~ 'prior',
                               Parameter == 'lp__' ~ 'log~probability')) %>%
  drop_na(Parameter) %>%
  ggplot(aes(y = AC, x = Lag, colour = as.factor(Chain))) + 
  geom_path() + 
  facet_wrap(~Parameter, labeller = label_parsed, ncol = 4) + 
  theme_minimal()  + 
  scale_colour_brewer('Chain') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/diffusion_autocorrelation.jpeg',
       dpi = 360,
       height = 20,
       width = 16,
       units = 'cm')


# Check Ratio of Effective Population Size to Total Sample Size 
# values <0.1 should raise concerns about autocorrelation
neff_ratio(diffusion_model) %>% 
  as_tibble(rownames = 'param') %>%
  mutate(param = fct_reorder(param, desc(value))) %>%
  ggplot() + 
  geom_segment(aes(yend = value,
                   xend=param, 
                   y=0,
                   x = param,
                   colour = value > 0.1)) +
  geom_point(aes(y = value,
                 x = param,
                 colour = value > 0.1)) + 
  geom_hline(aes(yintercept = 0.1), linetype = 'dashed') + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dashed') + 
  scale_colour_manual(values = c( '#0047AB',  'red')) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0),
                     expression(N["eff"]/N)) + 
  scale_x_discrete(expand= c(0.1,0), 'Fitted Parameter') + 
  theme_classic() + 
  coord_flip()  + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none') 



###################################### Posterior Checks ############################################

#Posterior predictive check
pp_check(diffusion_model, 
         ndraws = 100, 
         type = 'dens_overlay_grouped', 
         group = 'collection_regionname' ) %>%
  .$data %>%
  mutate(value = log1p(value)) %>%
  #filter(!is_y) %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group= rep_id, 
                   alpha = is_y,  
                   linewidth = is_y, 
                   colour = is_y), 
               key_glyph = draw_key_path) + 
  scale_alpha_manual(values = c('FALSE' = 0.00001, 
                                'TRUE' = 1)) + 
  scale_linewidth_manual('Data', values = c('FALSE' = 0.1, 
                                            'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual('Data',values = c('FALSE' = '#cbc9e2', 
                                        'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  guides(alpha= 'none', 
         colour=guide_legend()) + facet_wrap(~group, labeller  = as_labeller(str_to_title)) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km^{2}, ' ', year^{-1}, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 8, by = 2)))),
                     labels = expression(0,  1%*%10^{2},  1%*%10^{4},1%*%10^{6},1%*%10^{8}),
                     limits = c(-0.01, log1p(10^8.5)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Density',expand = c(0,0)) +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/diffusion_ppc.png',
       dpi = 360,
       device = 'png' ,
       height = 17,
       width = 17, 
       units = 'cm')

#Step 3B. Summarize posterior of variables
t <- get_variables(diffusion_model)

beta_draws <- diffusion_model %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  #filter(grepl('^b_', .variable)) %>%
  
  mutate(label = case_when(.variable == 'b_Intercept'~ "alpha['0']",
                           .variable == 'b_shape_Intercept'~ "alpha['0']^{shape}",
                           .variable == 'b_n_seq'~ "beta['3']",
                           .variable == "b_count_cross_species_log1p"~ "beta['1']",
                           .variable == 'b_V1'~ "beta['4']",
                           .variable == 'b_V2'~ "beta['5']",
                           .variable == 'b_V3'~ "beta['6']",
                           .variable == 'b_persistence'~ "beta['2']",
                           
                           .variable == 'b_n_seq:count_cross_species_log1p'~ "beta['7']",
                           
                           .variable == "sd_collection_regionname__Intercept" ~ "sigma[continent]",                                                                           
                           .variable == "sd_collection_regionname__V2" ~ "sigma[lr2]",  
                           .variable == "sd_collection_regionname__V3" ~ "sigma[lr3]",  
                           .variable == "sd_collection_year__Intercept" ~ "sigma[year]",  
                           .variable == "sd_collection_regionname__shape_Intercept" ~ "sigma[continent]^{shape}"))
  
ggplot() + 
  geom_histogram(data = beta_draws %>% drop_na(label), 
                 aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 bins = 70,
                 fill = '#1b9e77') + 
  
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['0']"),
                args = list(mean = 5, sd = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['0']^{shape}"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data =tibble(label = "beta['1']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['2']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +   
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['3']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['4']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label ="beta['5']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['6']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
   
   stat_function(fun = dnorm,
                 data = tibble(label = "beta['7']"),
                 args = list(mean = 0, sd = 2),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
   
   stat_function(fun = dexp,
                 data = tibble(label = "sigma[continent]"),
                 args = list(rate = 0.5),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
   
   stat_function(fun = dexp,
                 data = tibble(label = "sigma[lr2]"),
                 args = list(rate = 0.5),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
   
   stat_function(fun = dexp,
                 data = tibble(label = "sigma[lr3]"),
                 args = list(rate = 0.5),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
   
   stat_function(fun = dexp,
                 data =tibble(label = "sigma[year]"),
                 args = list(rate = 0.5),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
   
   stat_function(fun = dstudent_t,
                 data = tibble(label = "sigma[continent]^{shape}"),
                 args = list(df = 3, mu= 0, sigma =0.25),
                 fill = '#d95f02',
                 geom = 'area',
                 alpha = 0.5) +
 
   
   scale_y_continuous('Probability Density') + 
   scale_x_continuous('Parameter Value') + 
  facet_wrap(~label, scales = 'free',  ncol = 3, labeller = label_parsed) +
  theme_minimal() +
   theme(axis.text = element_text(size = 8),
         axis.title = element_text(size = 10),
         legend.text = element_text(size = 8))
 
 ggsave('~/Downloads/flu_plots/diffusion_identifiability.jpeg',
        dpi = 360,
        height = 29,
        width = 20,
        units = 'cm')
 


###################################### Residuals Checks ############################################
# Check Residuals using DHARMA 
# sample from the Posterior Predictive Distribution
preds <- posterior_predict(diffusion_model, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = t(posterior_predict(diffusion_model)),
  observedResponse = diffusion_data$weighted_diff_coeff,
  fittedPredictedResponse = apply(t(posterior_epred(diffusion_model)), 1, mean),
  integerResponse = FALSE)


# QQ plot
qq_data <- data.frame(
  sample = sort(res$scaledResiduals),
  theoretical = sort(ppoints(length(res$scaledResiduals)))
)

ggplot(qq_data, aes(sample = sample)) +
  stat_qq(distribution = stats::qunif) +
  stat_qq_line(distribution = stats::qunif) +
  scale_y_continuous('Observed')+
  scale_x_continuous('Expected')+
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/diffusion_qq.jpeg',
       dpi = 360,
       height = 12,
       width = 12,
       units = 'cm')



############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
###################################################################################################