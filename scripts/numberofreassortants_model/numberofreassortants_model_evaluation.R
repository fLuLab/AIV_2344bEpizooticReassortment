

###################################### MCMC Diagnostics #############################################
t <- get_variables(numbers_model_2)
# Trace plot
numbers_model_2 %>%
  gather_draws(., !!!syms(t)) %>%
  filter(grepl('^continent|^beta|^theta|^lp|^sigma', .variable)) %>%
  
  mutate(label = case_when(.variable == 'continent_specific_theta[1]'~ "theta['africa']",
                           .variable == 'continent_specific_theta[2]'~ "theta['asia']",
                           .variable == 'continent_specific_theta[3]'~ "theta['americas']",
                           .variable == "continent_specific_theta[4]"~ "theta['europe']",
                           
                           .variable == 'continent_specific_abundance[1]'~ "lambda['africa']",
                           .variable == 'continent_specific_abundance[2]'~ "lambda['asia']",
                           .variable == 'continent_specific_abundance[3]'~ "lambda['americas']",
                           .variable == "continent_specific_abundance[4]"~ "lambda['europe']",
                           
                           .variable == 'continent_specific_detection[1]'~ "p['africa']",
                           .variable == 'continent_specific_detection[2]'~ "p['asia']",
                           .variable == 'continent_specific_detection[3]'~ "p['americas']",
                           .variable == "continent_specific_detection[4]"~ "p['europe']",
                           
                           .variable == "beta_cases"~ "beta['cases']",
                           .variable == "beta_sequences"~ "beta['sequences']",
                           .variable == "sigma_year_abundance"~ "sigma['year-abundance']",
                           .variable == "sigma_year_detection"~ "sigma['year-detection']",
                           .variable == 'lp__' ~ 'log~probability')) %>%
  drop_na(label) %>%
  ggplot(aes(x = .iteration,
             y = .value, 
             col = as.factor(.chain)))+
  geom_line(alpha = 0.8) + 


facet_wrap(~label, labeller = label_parsed, ncol = 2, scale  = 'free_y',
           strip.position = 'left') + 
  theme_minimal()  + 
  scale_colour_brewer('Chain', palette = 'GnBu', 'Chains') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/number_traces.jpeg',
       dpi = 360,
       height = 20,
       width = 17,
       units = 'cm')


# Autocorrelation
stan_acf <- posterior::as_draws_array(numbers_model_2, nchains = 4) %>% 
  mcmc_acf()

stan_acf$data %>% 
  filter(!grepl('y_rep', Parameter)) %>% 
  
  mutate(label = case_when(Parameter == 'continent_specific_theta[1]'~ "theta['africa']",
                           Parameter == 'continent_specific_theta[2]'~ "theta['asia']",
                           Parameter == 'continent_specific_theta[3]'~ "theta['americas']",
                           Parameter == "continent_specific_theta[4]"~ "theta['europe']",
                           
                           Parameter == 'continent_specific_abundance[1]'~ "lambda['africa']",
                           Parameter == 'continent_specific_abundance[2]'~ "lambda['asia']",
                           Parameter == 'continent_specific_abundance[3]'~ "lambda['americas']",
                           Parameter == "continent_specific_abundance[4]"~ "lambda['europe']",
                           
                           Parameter == 'continent_specific_detection[1]'~ "p['africa']",
                           Parameter == 'continent_specific_detection[2]'~ "p['asia']",
                           Parameter == 'continent_specific_detection[3]'~ "p['americas']",
                           Parameter == "continent_specific_detection[4]"~ "p['europe']",
                           
                           Parameter == "beta_cases"~ "beta['cases']",
                           Parameter == "beta_sequences"~ "beta['sequences']",
                           Parameter == "sigma_year_abundance"~ "sigma['year-abundance']",
                           Parameter == "sigma_year_detection"~ "sigma['year-detection']",
                           
                           Parameter == 'lp__' ~ 'log~probability',
                           
                           
                           #Parameter == 'year[1]' ~ 'alpha[2019]',
                          # Parameter == 'year[2]' ~ 'alpha[2020]',
                          # Parameter == 'year[3]' ~ 'alpha[2021]',
                          # Parameter == 'year[4]' ~ 'alpha[2022]',
                          # Parameter == 'year[5]' ~ 'alpha[2023]',
                          # Parameter == 'year[6]' ~ 'alpha[2024]'
                           )) %>% 
  
  drop_na(label) %>%
  ggplot(aes(y = AC, 
             x = Lag,
             colour = as.factor(Chain))) +
  geom_path() + 
  facet_wrap(~label, labeller = label_parsed, ncol = 4) + 
  theme_minimal()  + 
  scale_colour_brewer('Chain') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/number_autocorrelation.jpeg',
       dpi = 360,
       height = 20,
       width = 16,
       units = 'cm')

 # Might be worth checking if any of the params are at the 'worry about' threshold neff/n
neff_ratio(numbers_model) %>% 
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
  scale_y_continuous(expand = c(0,0),
                     expression(N["eff"]/N)) + 
  scale_x_discrete('Fitted Parameter') + 
  theme_classic() + 
  coord_flip()  + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none') 



###################################### Posterior Checks ############################################
y_rep_matrix <- numbers_model_2$draws('y_rep') %>%
  posterior::as_draws_matrix()

# Global
ppc_bars(y =  data_processed_2 %>% pull(n_reassortants),
         yrep = y_rep_matrix[sample(0:16000, 100),]) %>% 
  .$data %>%
  ggplot() + 
  geom_bar(aes(x = x,
               y = y_obs),
           stat = 'identity',
           fill = '#cbc9e2',
           colour = '#cbc9e2',
           alpha = 0.7) + 
  geom_pointinterval(aes(x = x,
                         y=m,
                         ymin = l,
                         ymax = h), orientation = 'x',
                     colour = '#54278f') +
  scale_x_continuous('Reassortants per Month')+
  scale_y_continuous('Count', expand = c(0,0)) +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


# Grouped by Continent
ppc_bars_grouped(y =  data_processed_2 %>% pull(n_reassortants), 
                 group = data_processed_2 %>% pull(collection_regionname) %>% as.factor(),
                 yrep = y_rep_matrix[sample(0:16000, 250),]) %>%
  .$data %>%
  #mutate(group = case_when(group == 1 ~ 'Africa',
                           #group == 2 ~ 'Asia',
                           #group == 3 ~ 'Northern and Central America',
                           #group == 4 ~ 'Europe') ) %>%
  
  ggplot() + 
  geom_bar(aes(x = x,
               y = y_obs),
           stat = 'identity',
           fill = '#cbc9e2',
           colour = '#cbc9e2',
           alpha = 0.7) + 
  geom_pointinterval(aes(x = x,
                         y=m,
                         ymin = l,
                         ymax = h), orientation = 'x',
                     colour = '#54278f') +
  
  facet_wrap(~group, labeller  = as_labeller(str_to_title)) + 
  scale_x_discrete('Reassortants per Month')+
  scale_y_continuous('Count') +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/numbers_ppc.jpeg',
       dpi = 360,
       device = 'jpeg' ,
       height = 12,
       width = 12, 
       units = 'cm')


# Identifiability (Prior and Posterior Plots)
t <- get_variables(numbers_model_2) 

numbers_params <- numbers_model_2 %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^continent|beta|sigma|year', .variable)) %>%
  
  mutate(label = case_when(.variable == 'continent_specific_theta[1]'~ "theta['africa']",
                           .variable == 'continent_specific_theta[2]'~ "theta['asia']",
                           .variable == 'continent_specific_theta[3]'~ "theta['americas']",
                           .variable == "continent_specific_theta[4]"~ "theta['europe']",
                           
                           .variable == 'continent_specific_abundance[1]'~ "lambda['africa']",
                           .variable == 'continent_specific_abundance[2]'~ "lambda['asia']",
                           .variable == 'continent_specific_abundance[3]'~ "lambda['americas']",
                           .variable == "continent_specific_abundance[4]"~ "lambda['europe']",
                           
                           .variable == 'continent_specific_detection[1]'~ "p['africa']",
                           .variable == 'continent_specific_detection[2]'~ "p['asia']",
                           .variable == 'continent_specific_detection[3]'~ "p['americas']",
                           .variable == "continent_specific_detection[4]"~ "p['europe']",
                           
                           .variable == "beta_cases"~ "beta['cases']",
                           .variable == "beta_sequences"~ "beta['sequences']",
                           .variable == "sigma_year_abundance"~ "sigma['year-abundance']",
                           .variable == "sigma_year_detection"~ "sigma['year-detection']",
                        
                           .variable == 'lp__' ~ 'log~probability'))

numbers_params %>%
  drop_na(label) %>%
  ggplot() + 
  geom_histogram(aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 bins = 70, 
                 fill = '#1b9e77')+
  
    # thetas
  stat_function(fun = dbeta,
                data = tibble(label = "theta['africa']"),
                args = list(shape1 = 2, shape2 =5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  data = tibble(label = "theta['asia']"),
                  args = list(shape1 = 2, shape2 =5),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  data = tibble(label = "theta['americas']"),
                  args = list(shape1 = 2, shape2 =5),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  data = tibble(label = "theta['europe']"),
                  args = list(shape1 = 2, shape2 =5),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    # lambdas
    stat_function(fun = dnorm,
                  args = list(mean = 3, sd = 1.5),
                  data = tibble(label = "lambda['africa']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dnorm,
                  args = list(mean = 3, sd = 1.5),
                  data = tibble(label = "lambda['asia']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dnorm,
                  args = list(mean = 3, sd = 1.5),
                  data = tibble(label = "lambda['americas']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dnorm,
                  args = list(mean = 3, sd = 1.5),
                  data = tibble(label = "lambda['europe']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    # p
    stat_function(fun = dbeta,
                  args = list(shape1 = 1.5, shape2 = 2),
                  data = tibble(label = "p['africa']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  args = list(shape1 = 1.5, shape2 = 2),
                  data = tibble(label = "p['asia']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  args = list(shape1 = 1.5, shape2 = 2),
                  data = tibble(label = "p['americas']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dbeta,
                  args = list(shape1 = 1.5, shape2 = 2),
                  data = tibble(label = "p['europe']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = 1),
                  data = tibble(label = "beta['cases']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = 1),
                  data = tibble(label = "beta['sequences']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
    
    stat_function(fun = dexp,
                  args = list(rate = 0.5),
                  data = tibble(label = "sigma['year-abundance']"),
                  fill = '#d95f02',
                  geom = 'area',
                  alpha = 0.5) +
  
  stat_function(fun = dexp,
                args = list(rate = 0.5),
                data = tibble(label = "sigma['year-detection']"),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  facet_wrap(~label, scales = 'free',  ncol = 3, labeller = label_parsed) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/numbers_identifiability.jpeg',
       dpi = 360,
       height = 29,
       width = 20,
       units = 'cm')

  
  
  
numbers_model %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^continent_specific_abundance', .variable)) %>%
  ggplot() + 
  
  geom_density(aes(x = .value#,
                   #y = after_stat(density)
  ),
  inherit.aes = F, 
  fill = '#1b9e77') +
  
  stat_function(fun = dnorm,
                args = list(mean = 3, sd = 1.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  facet_wrap(~`.variable`)



###################################### Residuals Checks ############################################

simulated_residuals <- createDHARMa(
  simulatedResponse = t(y_rep_matrix[sample(0:16000, 250),]),
  observedResponse = data_processed_2 %>% pull(n_reassortants)
)


# QQ plot
qq_data <- data.frame(
  sample = sort(simulated_residuals$scaledResiduals),
  theoretical = sort(ppoints(length(simulated_residuals$scaledResiduals)))
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

ggsave('~/Downloads/flu_plots/numbers_qq.jpeg',
       dpi = 360,
       height = 12,
       width = 12,
       units = 'cm')