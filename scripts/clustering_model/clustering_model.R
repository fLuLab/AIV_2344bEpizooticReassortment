################################################################################
## Script Name:        Reassortant Clustering Modelk
## Purpose:            To classify reassortants according to observational and 
##                     phylodynamic variables
## Author:             James Baxter
## Date Created:       2024-09-14
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
library(recipes)
library(rsample)
library(mclust)
library(purrr)
library(fpc)

################################### DATA #######################################
# Read and inspect data
combined_data <- read_csv('./beast_outputs_summary/2024-09-20_combined_data.csv')

summary_data <- read_csv('./beast_outputs_summary/summary_reassortant_metadata_20240904.csv') %>%
  dplyr:::select(-c(cluster_label,
            clade)) 

meta <- read_csv('./data/2024-09-09_meta.csv')

key <- meta %>%
  dplyr::select(cluster_profile,
                cluster_label) %>%
  distinct() %>%
  drop_na()

reassortant_offspring <- read_csv('./beast_outputs_summary/reassortant_offspring.csv') %>% 
  dplyr::select(-cluster_class) %>%
  rename(cluster_label = name,
         offspring = importance) %>%
  left_join(key)


################################### MAIN #######################################
# Main analysis or transformation steps
observational_data <- summary_data %>%
  
  # select variables of interest
  dplyr:::select(c(
    cluster_profile,
    host_richness,
    if_mammal,
    max_distance_km,
    Num_Sequence
  ))

phylodynamic_data <- combined_data %>%
  
  #inlcude only HA and PB2
  filter(segment %in% c('ha', 'pb2')) %>%
  
  # select variables of interest
  dplyr:::select(
    segment,
    cluster_profile,
    weighted_diff_coeff,
    evoRate,
    persist.time,
    count_cross_species,
    count_to_mammal
  ) 

kclust_data <- phylodynamic_data %>%
  
  # Model pre-processing
  group_by(segment, cluster_profile) %>%
  slice_sample(n=1) %>%
  ungroup() %>% 

  # pivot wider so one row/reassortant
  pivot_wider(names_from = segment,
              values_from = where(is.double)) %>%
  
  left_join(reassortant_offspring ,
            by = join_by(cluster_profile)) %>%
  
  left_join(observational_data ,
            by = join_by(cluster_profile)) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  rename_with(~gsub('-', '_', .x)) %>%
  relocate(starts_with('cluster')) %>%
  
  # start tidymodels
  recipe(~ .) %>% 
  
  #drop na
  step_naomit(-starts_with('cluster')) %>%
  
  # normalise numeric vectors
  step_zv(everything()) %>%
  step_normalize(all_numeric()) %>% 
  
  # create dummy variables for categorical predictors
  #step_dummy(starts_with('collection_regionname'), one_hot = TRUE) %>%
  #step_dummy(starts_with('host_simplifiedhost'), one_hot = TRUE) %>%
  
  # run tidymodels recipe
  prep() %>%
  bake(NULL) 


# Run K means clustering for K \in 1-20
kclusts <- tibble(k = 1:20) %>%
  mutate(
    kclust = purrr::map(k, ~kmeans(kclust_data %>% 
                                     dplyr:::select(-starts_with('cluster')), .x)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
    augmented = purrr::map(kclust, augment, kclust_data))

# tidy results
clusters <- kclusts %>%
  unnest(cols = c(tidied))

assignments <- kclusts %>% 
  unnest(cols = c(augmented))

clusterings <-  kclusts %>%
  unnest(cols = c(glanced))


# Calculate cluster stats
dat <- assignments %>%
  filter(k == 3) %>%
  dplyr:::select(-c( kclust,
             tidied,
             glanced,
             starts_with('cluster'),
             k)) 

d <- dist(dat)
clusters <- as.integer(assignments %>% filter(k == 3) %>% pull(.cluster))
stats <- cluster.stats(d = d, clustering = clusters)
#stats$median.distance

# Permutation analysis - determining key variables invvled in cluster asignation
# data points used for clustering and cluster_profiles
labelled <- assignments %>%
  filter(k == 3) %>%
  dplyr:::select(-c( kclust,
             tidied,
             glanced,
             cluster_label,
             k,
             .cluster)) 

# generate permutations for all columns (except 1 - the cluster profile label)
kclust_permutations <- list()

for (i in 2:ncol(labelled)){
  kclust_permutations[[i]] <- labelled %>%
    permutations(permute = all_of(i), times = 10)
}

names(kclust_permutations) <- colnames(labelled)

# bind together permutaed data in a single dataframe, and nest
kclust_permutations %<>% 
  bind_rows(., .id = 'permuted_var') %>%
  unite(id, permuted_var, id)  %>%
  mutate(data = purrr::map(splits, ~ analysis(.x))) %>%
  dplyr:::select(-splits) %>%
  unnest(data) 

# execute kmeans on nested data (ie one run of kmeans per permuted dataset)
permuted <- kclust_permutations %>%
  nest(., .by = id) %>%
  mutate(
    kclust = purrr::map(data,  ~dplyr:::select(.x , -cluster_profile) %>% kmeans(3)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
    augmented = map2(kclust, data, augment)
  ) 

# extract results
permuted_clusters <- permuted %>%
  unnest(cols = c(tidied)) 

permuted_assignments <- permuted %>% 
  unnest(cols = c(augmented)) 

permuted_clusterings <- permuted %>%
  unnest(cols = c(glanced))

# Format for comparison (by means os adjusted rand index)
lookup_clusters <- assignments %>%
  filter(k == 3) %>%
  dplyr:::select(c(cluster_profile,
           .cluster)) %>%
  rename(original_cluster = .cluster)

# Calculate adjusted rand index
ari_phylo <- permuted_assignments %>%
  dplyr:::select(c(cluster_profile,
           .cluster,
           id)) %>%
  left_join(lookup_clusters, 
            by = join_by(cluster_profile)) %>%
  mutate(across(ends_with('cluster'), .fns = ~as.integer(.x))) %>%
  separate_wider_delim(id, delim = regex('_(?!.*_)'), 
                       names = c('var', 'permutation')) %>%
  group_by(var, permutation) %>%
  mutate(ari = mclust::adjustedRandIndex(original_cluster, .cluster) )%>%
  summarise(ari = mean(ari)) %>%
  ungroup()

ari_phylo_summary <- ari_phylo %>%
  mutate(importance = 1 - ari) %>% # Higher impact means lower ARI score
  summarise(mean_importance = mean(importance),
            lower_ci = mean(importance) - 1.96 * sd(importance) / sqrt(n()),  # Lower 95% CI
            upper_ci = mean(importance) + 1.96 * sd(importance) / sqrt(n()),  # Upper 95% CI
            .by = c(var))


################################### OUTPUT #####################################
# Save output files, plots, or results
plt_2a <- ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 8) +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(0,20, by = 5)) + 
  geom_vline(xintercept = 3, colour= 'darkgreen', linetype = 'dashed')  + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


marked_clusters <- c('H5N1/2022/R7_NAmerica', 
                     'H5N1/2020/R1_Europe', 
                     'H5N8/2019/R7_Africa', 
                     'H5N1/2021/R3_Europe', 
                     'H5N1/2022/R12_Europe',
                     'H5N1/2021/R1_Europe',
                     'H5N1/2023/R29_NAmerica')

plt_2b <-  
  assignments %>%
  filter(k == 3) %>%
  dplyr:::select(-c(kclust, tidied, glanced)) %>%
  recipe(~ .) %>%
  step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
  prep() %>%
  bake(NULL) %>%
  left_join(meta %>% 
              dplyr:::select(cluster_profile, cluster_label) %>%
              distinct() %>% 
              drop_na()) %>%
  ggplot(., aes(x = PC1,
                y = PC2)) +
  geom_point(aes(colour = .cluster,
                 alpha = ifelse(cluster_label %in% marked_clusters,
                                '1',
                                '0')),
             size = 2) + 
  geom_text(aes(label=ifelse(cluster_label %in% marked_clusters,
                             gsub('_.*', '', cluster_label),
                             ""), 
                colour = .cluster), 
            hjust = 1.1,
            nudge_y = -0.1,
            size = 2) +
  
 # scale_colour_brewer('K-Means Clusters',
                     # palette = 'Set1')+
  
  scale_colour_manual(values = c( '#FA9F42', '#A21817','#2B4162')) + 
  
  
  scale_alpha_manual(values = c('1' = 1, '0' = 0.3)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



plt_2c <- ari_phylo_summary %>%
  mutate(var = reorder(var, -mean_importance)) %>%
  ggplot()+
  geom_point(aes(x = var, y = mean_importance)) +
  geom_linerange(aes(x = var, ymin = lower_ci, ymax = upper_ci)) + 
  
  scale_x_discrete('Variable',
                   labels = c("Num_Sequence" = 'Number of sequences',
                              'count_cross_species_pb2' = 'PB2 species jumps',
                              'count_cross_species_ha' = 'HA species jumps',
                              'count_to_mammal_pb2' = 'PB2 mammal jumps',
                              'count_to_mammal_ha' = 'HA mammal jumps',
                              'evoRate_pb2' = 'PB2 evolutionary rate',
                              'evoRate_ha' = 'HA evolutionary rate',
                              'host_richness' = 'Number of host states',
                              'if_mammal' = 'Mammal infection',
                              "max_distance_km" = 'Maximum distance between sequences',
                              "offspring" = 'Number of offspring reassortants',
                              'persist.time_ha' = 'HA Persistence Time' ,
                              'persist.time_pb2' = 'PB2 Persistence Time' ,
                              'weighted_diff_coeff_ha' = 'HA diffusion coefficient',
                              'weighted_diff_coeff_pb2' = 'PB2 diffusion coefficient') %>%
                     str_wrap(., width = 15)
  ) +
  scale_y_continuous(
    'Variable Importance') +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

corr_vars <- ari_phylo_summary %>%
  mutate(var = reorder(var, -mean_importance)) %>%
  pull(var) %>% as.character()


plt_2d <- kclust_data %>% 
  dplyr::select(any_of(corr_vars)) %>%
  cor(method = 'spearman') %>%
  as.data.frame() %>% 
  rownames_to_column(var = "param1") %>%
  gather(key = param2, value = corr, -param1) %>%
  as_tibble() %>%
  filter(param1 != param2) %>%
  distinct(param = paste(pmin(param1, param2), 
                            pmax(param1, param2), sep = "_vs_"), .keep_all = TRUE) %>%
  mutate(param2 = factor(param2, levels = unique(param2))) %>%
  ggplot(aes(x = param1, y = param2, fill = corr)) + 
  geom_tile() + 
  geom_text(aes(label = round(corr, digits=3)), colour = 'black') + 
  scale_fill_distiller(palette = 'RdBu') + 
  theme_minimal() + 
  
  scale_x_discrete(NULL,
                   labels = c(
                              'count_cross_species_pb2' = 'PB2 species jumps',
                              'count_cross_species_ha' = 'HA species jumps',
                              'count_to_mammal_pb2' = 'PB2 mammal jumps',
                              'count_to_mammal_ha' = 'HA mammal jumps',
                              'evoRate_pb2' = 'PB2 evolutionary rate',
                              'evoRate_ha' = 'HA evolutionary rate',
                              'host_richness' = 'Number of host states',
                              'if_mammal' = 'Mammal infection',
                              "max_distance_km" = 'Maximum distance between sequences',
                              "offspring" = 'Number of offspring reassortants',
                              'persist.time_ha' = 'HA Persistence Time' ,
                              'persist.time_pb2' = 'PB2 Persistence Time' ,
                              'weighted_diff_coeff_ha' = 'HA diffusion coefficient',
                              'weighted_diff_coeff_pb2' = 'PB2 diffusion coefficient') %>%
                     str_wrap(., width = 20)
  ) + 
  scale_y_discrete(NULL,
                   labels = c("Num_Sequence" = 'Number of sequences',
                              'count_cross_species_pb2' = 'PB2 species jumps',
                              'count_cross_species_ha' = 'HA species jumps',
                              'count_to_mammal_pb2' = 'PB2 mammal jumps',
                              'count_to_mammal_ha' = 'HA mammal jumps',
                              'evoRate_pb2' = 'PB2 evolutionary rate',
                              'evoRate_ha' = 'HA evolutionary rate',
                              'host_richness' = 'Number of host states',
                              'if_mammal' = 'Mammal infection',
                              "max_distance_km" = 'Maximum distance between sequences',
                              "offspring" = 'Number of offspring reassortants',
                              'persist.time_ha' = 'HA Persistence Time' ,
                              'persist.time_pb2' = 'PB2 Persistence Time' ,
                              'weighted_diff_coeff_ha' = 'HA diffusion coefficient',
                              'weighted_diff_coeff_pb2' = 'PB2 diffusion coefficient') %>%
                     str_wrap(., width = 20)
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = 'inside',
        legend.position.inside = c(0,1),
        legend.justification = c(0,1)) 


cowplot::plot_grid( plt_2a,
                             plt_2b, 
                             ncol = 1,
                             align = 'hv',
                             axis = 'tb', 
                             scale = 0.95,
                             labels = 'AUTO', 
                             label_size = 10)

ggsave('~/Downloads/flu_plots/figure_clustering.jpeg', height = 22 , width = 18, units = 'cm', dpi = 360)

 cowplot::plot_grid( plt_2c,
                             plt_2d  ,
                             ncol = 1,
                             align = 'hv',
                             axis = 'tlrb', 
                     scale = 0.95,
                             labels = 'AUTO', 
                             label_size = 10)
 
 ggsave('~/Downloads/flu_plots/figure_clustering2.jpeg', height = 30 , width = 25, units = 'cm', dpi = 360)
#################################### END #######################################
################################################################################