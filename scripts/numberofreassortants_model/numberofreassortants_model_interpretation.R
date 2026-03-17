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
countmodel_mv_fit <- readRDS('./saved_models/count_model.rds')
count_data <- read_csv('./saved_models/count_data.csv')


# MCMC chains
mcmc_countmodel <- ggs(countmodel_mv_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion


# Posterior predictive distribution
posteriorpredictive_reassortantclass <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'reassortantclass')
posteriorpredictive_nreassortants <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'nreassortants')

countmodel_posteriorpredictive <- bind_rows(posteriorpredictive_nreassortants$data %>%
  mutate(resp = 'nreassortants'),
  posteriorpredictive_reassortantclass$data %>%
    mutate(resp = 'reassortantclass'))
  

                                             

############################################## MAIN ################################################
# Expected values with tidybayes

#avg_prediction
diffusionmodel1_fit_gamma_19 %>% 
  epred_draws(newdata = diffusion_data) %>%  
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)




plt_c <- 
  
countmodel_mv_fit %>% 
  epred_draws(newdata = count_data %>%
                select(collection_season) %>%
                distinct(),
              resp = "nreassortants",
              re_formula = NA)

countmodel_mv_fit %>% 
  emmeans(~collection_season,
          re_formula = NA,
          dpar = 'nreassortants',
          epred = TRUE)

countmodel_mv_fit %>% 
  add_predicted_draws(newdata = count_data %>%
                        select(collection_regionname) %>%
                        distinct(),
                      re_formula = NA)


countmodel_mv_fit %>% 
  emmeans(~ collection_regionname,
          at = list(collection_regionname = count_data %>%
                      select(collection_regionname)),
                      re_formula = NA)


countmodel_mv_fit %>% 
  add_predicted_draws()

joint_preds <- expand_grid(collection_regionname = count_data %>%
              pull(collection_regionname) %>% 
              unique(),
            
            collection_season = count_data %>%
              pull(collection_season) %>% 
              unique(),
            
            n_cases = 20,
            
            previous_reassortant_class = count_data %>%
              pull(previous_reassortant_class) %>% 
              unique()
        
            ) %>%
  add_predicted_draws(countmodel_mv_fit,re_formula = NA)%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  summarise(n = n(),.by = c(collection_regionname, nreassortants, reassortantclass))%>%
  mutate(freq=n/sum(n))



  

  

# Conditional prediction draws  + Conditional marginal means stratified by continent

# Required outputs: Pre/post epizootic, marginal effect of region ,
# 


# Posterior prediction of the number and type of reassortant at 'average' parameters (global and stratified by region)


# marginal effect of of continent on number and type of reassortant (ie how many more reassortants in a vs b)


# reassortant 'transition' probabilities (ie based on lagged reassortant type), global and stratified by region
predict_season <- countmodel_mv_fit %>%
  emmeans( ~ collection_season,
          var = "nreassortants",
          at = list(collection_season = unique(count_data$collection_season)),
          epred = TRUE, 
          dpar = 'mu_nreassortants',
          re_formula = NA, 
          regrid = "response",
          #allow_new_levels = TRUE
          ) %>% 
  gather_emmeans_draws()

ggplot(predict_season, aes(x = .value, y= collection_season)) +
  stat_slabinterval()
  


countmodel_mv_fit %>%
  emmeans( ~ collection_regionname,
           #var = "nreassortants",
           at = list(collection_regionname= unique(count_data$collection_regionname)),
           epred = TRUE, 
           #dpar = 'mu_nreassortants',
           re_formula = NA, 
           regrid = "response",
           #allow_new_levels = TRUE
  ) %>% 
  gather_emmeans_draws()



test <- countmodel_mv_fit %>%
  emmeans( ~ previous_reassortant_class,
           at = list(),
           #epred = TRUE, 
           dpar = 'reassortantclass',
          re_formula = NA,
          dpar ='reassortantclass',
          regrid = "response",
           #allow_new_levels = TRUE
  )  %>% 
  gather_emmeans_draws()

avg_predictions(countmodel_mv_fit)





ggplot(aes(x  = 1-.value,
           y = season,
           slab_colour = season,
           slab_fill = season)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)



countmodel_mv_fit %>%
  conditional_effects(., effects = '')
##############################################

test <- countmodel_mv_fit %>%
  add_epred_draws(newdata = expand_grid(previous_reassortant_class = unique(count_data$previous_reassortant_class),
                                        collection_regionname = unique(count_data$collection_regionname)),
                  re_formula = NA,
                  resp = 'reassortantclass')

tidy_dag <- test %>%
  filter(previous_reassortant_class != 'none') %>%
  filter(.category != 'none') %>%
  median_hdci() %>%
  rename(name = previous_reassortant_class, 
         to = .category) %>% 
  as_tidy_dagitty()  #%>%
  #mutate(circular = case_when(name == to ~ yend - TRUE, .default = FALSE)) 
  #rename(from = name) %>% 
  #igraph::graph_from_data_frame(., directed = TRUE) %>%
  #plot()
  
test %>%
  filter(previous_reassortant_class != 'none') %>%
  filter(.category != 'none') %>%
  ungroup(previous_reassortant_class) %>%
  median_hdi() 

##########################################




emms <- emmeans(

  countmodel_mv_fit, ~ reassortantclass | previous_reassortant_class,  resp = "reassortantclass"                    
)

# Summarize probabilities
summary(emms, type = "response")
plt_5e <- predict_hu_season %>%
  gather_emmeans_draws()

countmodel_mv_fit %>%
  predicted_draws(newdata = tibble(previous_reassortant_class  = unique(count_data$previous_reassortant_class)),
                  re_formula = NA,
                  dpar ='munone_reassortantclass')

  add_epred_draws(.,
                  )


  
##### PLOTS #####
plt_a <- count_data %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  filter(reassortant_class != 'none') %>%
  ggplot() + 
  geom_histogram(aes(x = collection_monthyear, 
                     fill = collection_regionname, 
                     colour = collection_regionname,
                     alpha = reassortant_class))  +
  scale_x_date('Reassortant Emergence',expand = c(0,0))  +
  scale_y_continuous('Count',expand = c(0,0)) +
    scale_fill_manual('Continent' ,values = region_colours,
                      labels = str_to_title) +
    scale_colour_manual('Continent' , values = region_colours,
                        labels = str_to_title) +
    scale_alpha_manual('Reassortant Class' , values = c('dominant' = 1, 'major' = 0.7, 'minor' = 0.5),
                       labels = str_to_title) +
  global_theme + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.2,0.7))


plt_b <- count_data %>%
  filter(reassortant_class != 'none') %>%
  summarise(n = n(), .by = c(collection_regionname, reassortant_class)) %>%
  group_by(collection_regionname) %>%
  mutate(freq = n/sum(n)) %>%
  ggplot() + 
  geom_bar(stat = 'identity',
           aes(x = collection_regionname,
               y = freq,
               fill = collection_regionname, 
               colour = collection_regionname,
               alpha = reassortant_class))  + 
  #scale_fill_brewer('Reassortant Class', palette = 'Set1') + 
  #scale_colour_brewer('Reassortant Class', palette = 'Set1') + 
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  
  scale_alpha_manual(values = c('dominant' = 1, 'major' = 0.7, 'minor' = 0.5)) +
  scale_x_discrete('Continent',expand = c(0,0), labels = str_to_title)  +
  scale_y_continuous('Proportion',expand = c(0,0)) +
  global_theme


# joint predictions
plt_c <- joint_preds %>% 
  #mutate(freq=n/sum(n)) %>%
  mutate(reassortantclass = case_when(reassortantclass == 3 ~ 'minor', 
                                      reassortantclass == 1 ~ 'dominant',
                                      reassortantclass == 2 ~ 'major')) %>%
  filter(nreassortants != 0) %>% 
  filter(freq > 0.001) %>% 
  drop_na(reassortantclass) %>%
  ggplot() + 
  geom_point(aes(y = nreassortants, 
                 x = as.character(reassortantclass), 
                 size = freq,
                 fill = collection_regionname,
                 colour = collection_regionname,
                 alpha = reassortantclass)) +
  facet_grid(cols= vars(collection_regionname), labeller = labeller(collection_regionname = str_to_title)) + 
  scale_y_continuous('N' , breaks= c(1,3,5,7), limits = c(1,7)) +
  scale_x_discrete('Reassortant Class', labels = str_to_title) +
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  
  scale_alpha_manual(values = c('dominant' = 1, 'major' = 0.7, 'minor' = 0.5)) +
  global_theme


# transition probabilities
plt_d <- tidy_dag %>%
  mutate(name = case_when(name == 'dominant'~ 'D',
                          name == 'major' ~ "Mj",
                          name == 'minor' ~ 'Mn')) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(alpha = name, colour = collection_regionname)) +
  geom_dag_text(size = 3) +
  geom_dag_edges_fan(aes(edge_width = .epred*2,
                         label = paste0(round(.epred, digits = 3), 
                                        ' (', round(.lower, digits = 3),
                                        '-', round(.upper, digits = 3), ')')), 
                     spread = 5,
                     label_size = 3,
                     label_dodge	= 0.5) + 
  #geom_dag_edges_arc(aes(edge_width = .epred))+ 
  facet_grid(cols = vars(collection_regionname),
             labeller = labeller(collection_regionname= str_to_title)) + 
  scale_alpha_manual(values = c('D' = 1, 'Mj' = 0.7, 'Mn' = 0.5)) +
  scale_colour_manual(values = region_colours) +
  theme_dag(legend.position = 'none', base_size = 8)



plots <- align_plots(plt_a, plt_c, plt_d, align = 'v', axis = 'l')


plot_grid(plot_grid(plots[[1]], plt_b, ncol = 2, labels = 'AUTO'), plots[[2]], plots[[3]], labels = c('', 'C', 'D'), nrow = 3)


ggsave('~/Downloads/figure6.jpeg',
       height = 30,
       width = 35,
       dpi= 360,
       units= 'cm')
# number of reassortants ~ season

##### Plot Posterior Predictive Check #####
pp_a <- countmodel_posteriorpredictive  %>%
  filter(resp == 'nreassortants') %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group= rep_id, 
                   alpha = is_y,  
                   linewidth = is_y, 
                   colour = is_y), 
               key_glyph = draw_key_path) + 
  scale_alpha_manual(values = c('FALSE' = 0.00001, 
                                'TRUE' = 1)) + 
  scale_linewidth_manual(values = c('FALSE' = 0.1, 
                                    'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual(values = c('FALSE' = '#cbc9e2', 
                                 'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  guides(alpha= 'none', 
         colour=guide_legend()) + 
  scale_x_continuous('Number of Reassortants',
                     expand = c(0,0),
  ) +
  scale_y_continuous('Density',expand = c(0,0)) + 
  global_theme+ 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))


pp_b <- countmodel_posteriorpredictive  %>%
  filter(resp == 'reassortantclass') %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group= rep_id, 
                   alpha = is_y,  
                   linewidth = is_y, 
                   colour = is_y), 
               key_glyph = draw_key_path) + 
  scale_alpha_manual(values = c('FALSE' = 0.00001, 
                                'TRUE' = 1)) + 
  scale_linewidth_manual(values = c('FALSE' = 0.1, 
                                    'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual(values = c('FALSE' = '#cbc9e2', 
                                 'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  guides(alpha= 'none', 
         colour=guide_legend()) + 
  scale_x_continuous('Reassortant Class',
                     expand = c(0,0),
  ) +
  scale_y_continuous('Density',expand = c(0,0)) + 
  global_theme+ 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))

cowplot::plot_grid(pp_a,  pp_b, nrow=1,align='h',axis='tb',labels='AUTO')


#### Plot MCMC chains for Key parameters #####
plt_mcmc <- mcmc_diffusion %>% 
  filter(Iteration > 400) %>%
  mutate(Parameter = factor(Parameter,
                            levels = c("b_Intercept", "b_hu_Intercept", "b_median_anseriformes_wild",
                                       "b_median_charadriiformes_wild", "b_hu_seasonmigrating_spring", "b_hu_seasonmigrating_autumn",
                                       "b_hu_seasonoverwintering", "sd_collection_regionname__Intercept",
                                       "sd_segment__Intercept", "sd_collection_regionname__hu_Intercept", 
                                       "sd_segment__hu_Intercept", "sigma"),
                            labels = c(expression(paste(beta[0])),
                                       expression(paste(beta['hu'])),
                                       expression(paste(beta['anseriformes'])),
                                       expression(paste(beta['charadriiformes'])),
                                       expression(paste(beta['hu_seasonmigrating_spring'])),
                                       expression(paste(beta['hu_seasonmigrating_autumn'])),
                                       expression(paste(beta['hu_seasonoverwinterin'])),
                                       expression(paste(sigma['collection_regionname'])),
                                       expression(paste(sigma[mu])),
                                       expression(paste(sigma['hu_collection_regionname'])),
                                       expression(paste(sigma['hu_segment'])),
                                       expression(paste(sigma[epsilon]))
                            )
  ))  %>% drop_na(Parameter) %>%
  ggplot(aes(x = Iteration,
             y = value, 
             col = as.factor(Chain)))+
  geom_line(alpha = 0.8)+
  facet_wrap(~ Parameter,
             ncol = 2,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'bottom')


#### Plot Prior & posterior distributions of key parameters ####
posterior_beta_draws <- as.data.frame(countmodel_mv_fit) %>%
  as_tibble() %>%
  pivot_longer(., cols = everything(),
               names_to = 'parameter',
               values_to = 'estimate') %>% 
  filter(!grepl('lp__|lprior|^r_|^cor|^sd', parameter)) %>%
  #filter(grepl('b_', parameter)) %>%
  mutate(draw = 'posterior') 


plt_params <- ggplot() + 
  geom_histogram(data = posterior_beta_draws , 
                 aes(x = estimate,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') + 
  
  # Intercepts
  stat_function(fun = dlogis,
                data = tibble(parameter = "Intercept_zi_nreassortants" ),
                args = list(location = 0, scale = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_nreassortants"),
                args = list(mu = 3, sigma = -2.3, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_mumajor_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_muminor_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_munone_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  # Beta - number of reassortants
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_n_cases" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classmajor" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classminor" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classnone" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonmigrating_autumn" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonmigrating_spring" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonoverwintering" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classnone" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  # Beta - ZI nreassortants
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_n_cases" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  # Beta - reassortant class
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dgamma,
                data = tibble(parameter = "shape_nreassortants"),
                args = list(shape = 0.01, rate = 0.01),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  
  
  xlim(c(-15,20)) + 
  facet_wrap(~parameter, scales = 'free_y',  ncol = 4) +
  theme_minimal(base_size = 8)

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################