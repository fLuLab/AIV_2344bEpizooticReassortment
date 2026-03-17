####################################################################################################
####################################################################################################
## Script name: Number of reassortants model data
##
## Purpose of script: to prepare data for the number of reassortant model
##
##
## Date created: 2024-xx-xx
##
##
########################################## SYSTEM OPTIONS ##########################################
#options(scipen = 6, digits = 7) 
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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# User functions

FormatContinent <- function(dataframe){
  dataframe %<>%
    mutate(collection_regionname = case_when(grepl('europe', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'europe',
                                             
                                             grepl('africa',
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'africa',
                                             
                                             grepl('asia', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'asia',
                                             
                                             grepl('(central|northern) america', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'central & northern america',
                                             
                                             .default = NA_character_ ))
  
  return(dataframe)
}

############################################## DATA ################################################
combined_data <- read_csv('./beast_outputs_summary/2024-09-20_combined_data.csv')
summary_data <- read_csv('./beast_outputs_summary/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
            clade)) 

meta <- read_csv('./data/2024-09-09_meta.csv') 
woah_hpai <- read_csv('./data/Quantitative data 2025-04-23.csv')

############################################## MAIN ################################################
# Data pre-processing
all <- meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  dplyr::select(collection_regionname, collection_dateyear) %>%
  tidyr::expand(collection_regionname, collection_dateyear = full_seq(collection_dateyear,1))


# Sequences
sequences_month <- meta %>%  
  drop_na(cluster_profile) %>%
  dplyr::select(starts_with('collection_date'),
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, collection_regionname) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(collection_regionname,collection_datemonth)


# Format WOAH data and estimate the minimum number of cases
ref <- ne_countries() %>% 
  as_tibble() %>%
  dplyr::select(name, continent)
  
woah_minimuminferredcases <- woah_hpai %>%
  dplyr::select(Year,
         Semester,
         `World region`,
         Country,
         `Animal Category`,
         `Serotype/Subtype/Genotype`,
         `New outbreaks`, 
         Susceptible, 
         `Measuring units`,
         Cases,
         Deaths) %>%
  
  # replace dashes with NA
  mutate(across(!Year, ~ case_when(grepl('^-$', .x) ~ NA_character_,
                                   .default = .x))) %>%
  
  # format biannual periods
  mutate(date_start = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-01-01'),
                                grepl('Jul-Dec', Semester) ~ paste0(Year, '-07-01'),
                                .default = NA_character_),
         date_end = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-06-30'),
                              grepl('Jul-Dec', Semester) ~ paste0(Year, '-12-31'),
                              .default = NA_character_)) %>%
  
  # format country for continent assignment
  mutate(Country = case_when(Country ==  "Cote D'Ivoire" ~ "Côte d'Ivoire",
                             Country == "Congo (Dem. Rep. of the)" ~ "Dem. Rep. Congo",
                             Country == "China (People's Rep. of)"~ "China",
                             Country == "Chinese Taipei" ~ "Taiwan",
                             Country == "Korea (Dem People's Rep. of)" ~"North Korea",
                             Country == "Korea (Rep. of)" ~"South Korea",
                             Country == "Dominican (Rep.)" ~ "Dominican Rep.",
                             grepl( "Falkland Islands \\(Malvinas\\)|South Georgia and the South Sandwich Islands", Country) ~ "Falkland Is.",
                             Country ==  "Türkiye (Rep. of)" ~ "Turkey",
                             Country == "Hong Kong" ~ 'China',                                 
                             Country == "Bosnia and Herzegovina" ~  "Bosnia and Herz." ,                 
                             Country == "Czech Republic" ~ "Czechia",                               
                             Country == "Faeroe Islands" ~ "Denmark",                              
                             Country == "Reunion" ~ 'Madagascar',                                      
                             .default = Country)) %>%
  # Group by continent
  left_join(ref, by = join_by(Country == name)) %>%
  mutate(animal_category  = str_to_lower(`Animal Category`)) %>%   
  mutate(across(any_of(c('New outbreaks', 'Susceptible', 'Cases', 'Deaths')), ~as.numeric(.x))) %>%
  replace_na(list('New outbreaks' = 0,
                  'Susceptible' = 0,
                  'Cases' = 0,
                  'Deaths' = 0)) %>%
  
  # Filter only H5NX
  filter(grepl('^[Hh]5', `Serotype/Subtype/Genotype`)) %>%
  
  # Infer the minimum number of 'cases' This is mostly due to USA seemingly incapable of reporting
  # their data in a similar manner to the rest of the world.
  mutate(inferred_minimum_cases = case_when(Cases >= Deaths ~ Cases, 
                                            Cases < Deaths ~ Deaths,
                                            .default = Cases)) %>%
  
  rename(collection_regionname = continent) %>%
  mutate(collection_regionname = str_to_lower(collection_regionname)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern|north) america|caribbean', collection_regionname) ~  'central & northern america',
                                           grepl('south america|southern ocean|antarctica', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia|oceania', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname))


# Group WOAH minimum cases by month
woah_minimuminferredcases_monthly <- woah_minimuminferredcases %>%
  
  summarise(sum_IMC = sum(inferred_minimum_cases, na.rm = TRUE),
            sum_deaths = sum(Deaths),
            sum_susceptibles = sum(Susceptible),.by = c(date_start, 
                                                        date_end,
                                                        collection_regionname))  %>%
  mutate(across(starts_with('date'), ~ymd(.x))) %>%
  #mutate(interval = interval(date_start, date_end),
        # collection_datemonth = date_start) %>%
 # dplyr::select(-c(date_start, date_end)) %>%
  
  rename(woah_cases = sum_IMC,
         woah_deaths = sum_deaths,
         woah_susceptibles = sum_susceptibles)  %>%
 #mutate(end_date = int_end(interval), 
         #start_date = int_start(interval)) %>%
  transmute(collection_regionname,
            woah_cases,woah_deaths,
            woah_susceptibles, 
            collection_datemonth = map2(date_start, date_end, seq, by = "1 month")) %>%
  unnest(collection_datemonth)
                                                                                                                                                                                                                           




# Join everything together and calculated scaled estimates (not used)
#all_casedata_monthly <- 
  #%>%
  #interval_left_join(fao_hpai_monthly, by = c('date_start', 'date_end')) %>%
  #filter(collection_regionname.x == collection_regionname.y) %>%
  #dplyr::select(-ends_with('y')) %>%
  #rename_with(~gsub('\\.x', '', .x)) %>%
  #filter(collection_datemonth >= as_date('2019-01-01')) %>%
  #dplyr::select(-starts_with('date')) %>%
  #mutate()







# Number and TMRCA of reassortants as inferred from our phylodynamic analysis
reassortant_counts <- combined_data %>%
  FormatContinent() %>%
  
  # restrict to HA for now
  filter(segment == 'ha') %>%
  mutate(collection_datemonth = date_decimal(TMRCA) %>% 
           #format(., "%Y-%m-%d") %>%
           floor_date(unit = 'month') %>% 
           as_date()) %>%
  
  filter(!is.na(collection_regionname)) %>%
  dplyr::select(collection_regionname,
         collection_datemonth, 
         cluster_profile,
         group2) %>% # include all month-years over collection period to generate zer countr
  summarise(n_reassortants = n_distinct(cluster_profile), 
            .by = c(collection_datemonth,
                    collection_regionname,
                    group2
            ))  %>%
  
  # within each continent...
  arrange(collection_regionname, collection_datemonth) %>%
  group_by(collection_regionname) %>%
  
  
  # Get the class of t-1 reassortant
  pivot_wider(names_from = group2, values_from = n_reassortants) %>%
  
  mutate(joint_reassortant_class = case_when(
    minor > 0 & is.na(major) & is.na(dominant) ~ 'minor',
    is.na(minor) & major > 0 & is.na(dominant) ~ 'major',
    is.na(minor)  & is.na(major) & dominant > 0 ~ 'dominant',
    minor > 0 & major > 0 & is.na(dominant) ~ 'minor_major',
    is.na(minor) & major > 0 & dominant > 0 ~ 'major_dominant',
    minor > 0 & is.na(major) & dominant > 0 ~ 'minor_dominant',
    minor > 0 & major > 0 & dominant > 0 ~ 'minor_major_dominant',
    .default = NA)) %>%
  
  fill(joint_reassortant_class, .direction = 'down') %>%
  mutate(previous_reassortant_class = lag(joint_reassortant_class)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(n_reassortants = sum(minor, dominant,major, na.rm = T)) %>%
  as_tibble()


# Data for model
count_data <- woah_minimuminferredcases_monthly %>% 
  
  left_join(reassortant_counts) %>%
  left_join(sequences_month) %>%
  
  mutate(collection_year = year(collection_datemonth),
         collection_month = month(collection_datemonth),
         collection_quarter = quarter(collection_datemonth),
         collection_dec = decimal_date(collection_datemonth)) %>%
  
  #separate_wider_delim(collection_monthyear, '-', names = c('collection_year', 'collection_month'),   cols_remove = FALSE) %>%
  #mutate(across(c('collection_month', 'collection_year'), .fns = ~ as.double(.x))) %>%
  #mutate(time = ym(collection_monthyear) %>% decimal_date()) %>%
  
  # within each continent...
  group_by(collection_regionname) %>%
  
  
  # And infer the time since the last dominant reassortant
  mutate(last_dominant = case_when(!is.na(dominant) ~ collection_dec, .default = NA)) %>%
  fill(last_dominant) %>%
  mutate(time_since_last_dominant = as.numeric(collection_dec - last_dominant)) %>%
  
  # Classify tmrca month according to breeding season
  mutate(collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                                       collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                                       collection_month %in% c(6,7,8)  ~ 'breeding', 
                                       collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
  )) %>%
  
  # dplyr::select variables
  ungroup(collection_regionname) %>%
  dplyr::select(collection_regionname,
         collection_year,
         collection_month,
         collection_season, 
         collection_quarter,
         collection_datemonth,
         woah_cases,
         woah_susceptibles,
         woah_deaths,
         n_sequences,
         n_reassortants,
         minor,
         dominant,
         major,
         joint_reassortant_class,
         previous_reassortant_class,
         collection_dec,
         time_since_last_dominant)



# currently missing lpai diversity

write_csv(count_data, './data/countmodeldata_2025Apr23.csv')
