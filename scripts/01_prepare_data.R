library(tidyverse)
library(checkmate)
library(rstan)
library(sp)
library(spdep)
library(igraph)
library(zoo)
library(kableExtra)

# setup ------------------------------------------------------------------------

## set wd
setwd("D:/MA_HMM")

## load functions
source('./0_functions.R')

## load data
lkr_data <- sf::read_sf('./germany_shp/mit_ewz/vg250-ew_ebenen_1231/VG250_KRS.shp')
kh_data <- read.csv('./K-2020-AI014-1--AI1401--2022-10-18.csv', sep = ';')
divi_data <- read.csv('./zeitreihe-tagesdaten.csv', sep = ";")
ils_data <- read.csv('./ILS_AGS_Zuordnung_complete.csv', sep = ";")
rki_data <- read.csv('./RKI_COVID19.csv', sep = ",")

## Epid Bull 2022;38:7-25 | DOI 10.25646/10582
phases_data <- read.csv('./covid_phases_germany.csv', sep = ";")

# prepare data -----------------------------------------------------------------

#'@details Prepare the data for fitting the models. Results can be obtained for
#'Germany as a whole or any of its 16 federal states. We re-use this code to
#'obtain data for Germany, suffix _germany, and Bavaria, suffix _bavaria

## combine the different data sources
merged_data <- merge_raw_data(
  age_group = NULL, 
  federal_state = '09')

## add lagged values of a covariate for prediction later on
data_w_lagged_values <- lag_covariates(
  input_data = merged_data$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 14,
  cov_name = 'incidence'
)

# set up the neighborhood matrices --------------------------------------------

## Three types of matrices, spatial distance (SD), adjacency (B) and distance
## in bed density (BD)
nb_mat_sd <- make_adjacency_matrix(
  input_sp = merged_data$ils_sp,
  weights_style = 'SD'
)

nb_mat_b <- make_adjacency_matrix(
  input_sp = merged_data$ils_sp,
  weights_style = 'B'
)

nb_mat_bd <- make_adjacency_matrix(
  input_sp = merged_data$ils_sp,
  weights_style = 'BD'
)

## turn into matrix objects
nb_b <- as.matrix(as.data.frame(nb_mat_b))
nb_sd <- as.matrix(as.data.frame(nb_mat_sd))
nb_bd <- as.matrix(as.data.frame(nb_mat_bd))

## cut off at 60 km (respectively 80 km when setting up data for Germany)
nb_sd[nb_sd > 80] <- 0

## round and offset nb_bd
nb_bd <- round(nb_bd, digits = 3) + 1

## everything zero here should be zero in the other matrix as well
nb_bd[nb_sd == 0] <- 0

## scale the results
nb_sd <- scale_matrix(nb_sd)
nb_bd <- scale_matrix(nb_bd)

## inverse of the entries
nb_sd <- nb_sd^-1
nb_bd <- nb_bd^-1

## set Inf to zero
nb_sd[nb_sd == Inf] <- 0
nb_bd[nb_bd == Inf] <- 0

# prepare pandemic phases ------------------------------------------------------

## get start and end date per COVID phase in Germany
phases_data <- phases_data %>%
  mutate(
    start_date = lubridate::parse_date_time(
      paste(y_start, 
            cw_start, 1, sep="/"),'Y/W/w'),
    end_date = lubridate::parse_date_time(
      paste(y_end, 
            cw_end, 1, sep="/"),'Y/W/w')
  )

# fitting ----------------------------------------------------------------------

## Load the particular model code for fitting the models
## The Stan code is in located in the "./models/" folder
model_01 <- rstan::stan_model(file = "./models/model_01_rw2_normal.stan")
model_02 <- rstan::stan_model(file = "./models/model_02_rw2_normal.stan")
model_03 <- rstan::stan_model(file = "./models/model_03_rw2_normal.stan")
model_04 <- rstan::stan_model(file = "./models/model_04_rw2_normal.stan")

list_of_dates <- get_seed_dates(phases_data, share = 0.05)

#saveRDS(data_w_lagged_values, 'data_w_lagged_values.rds')
#saveRDS(list_of_dates, 'list_of_dates.rds')

#data_w_lagged_values <- readRDS('data_w_lagged_values.rds')
#list_of_dates <- readRDS('list_of_dates.rds')


# select data and prepare stan import  -----------------------------------------

for (j in 1:length(list_of_dates)) {
  
  ## get training data
  training_data <- extract_time_slice(
    input_data = data_w_lagged_values, 
    max_t = 30,
    from_date = list_of_dates[j]
  )
  
  ## prepare data for import into stan
  stan_data <- prepare_stan_data(
    input_data = training_data,
    distribution = 'lognormal',
    W_input = nb_b)
  
  ## same with covariates
  cov_list <- make_covariates()
  stan_data_cov <- add_covariates_to_stan()
  stan_data_cov <- add_age_percentages_to_stan_data()
  
  
  # fit models -----------------------------------------------------------------
  
  fit_01 <- rstan::vb(
    object = model_01,
    data = stan_data,
    iter = 200000
  )

  stan_data$W <- nb_b

  fit_02 <- rstan::vb(
    object = model_02,
    data = stan_data,
    iter = 200000
  )

  stan_data$W <- nb_sd

  fit_03 <- rstan::vb(
    object = model_02,
    data = stan_data,
    iter = 200000
  )

  stan_data$W <- nb_bd

  fit_04 <- rstan::vb(
    object = model_02,
    data = stan_data,
    iter = 200000
  )
  
  fit_05 <- rstan::vb(
    object = model_03, 
    data = stan_data_cov, 
    iter = 200000
  )
  
  stan_data_cov$W <- nb_b
  
  fit_06 <- rstan::vb(
    object = model_04, 
    data = stan_data_cov, 
    iter = 200000
  )
  
  stan_data_cov$W <- nb_sd
  
  fit_07 <- rstan::vb(
    object = model_04, 
    data = stan_data_cov, 
    iter = 200000
  )
  
  stan_data_cov$W <- nb_bd
  
  fit_08 <- rstan::vb(
    object = model_04, 
    data = stan_data_cov, 
    iter = 200000
  )
  
  ## save to results
  for (k in 1:8) {
    
    saveRDS(
      get(paste0('fit_0', k)), 
      paste0(
        './stan_fits/fit_0', k, '_iter_', j, 
        '_normal_bayern_rw2.rds')
    )
    
  }
  
}

# extract Stan results ---------------------------------------------------------

#'@details In this section, extract the Stan results from the fit objects for
#'further processing. 

final_bay_data <- prepare_eval_data()

# final_bay_data_MULTI <- prepare_eval_data(
#   dates_list = list_of_dates,
#   file_suffix = '_normal_bayern_rw2_MULTI.rds')
# 
# final_german_data <- prepare_eval_data(
#   dates_list = list_of_dates,
#   prepared_data = data_w_lagged_values_germany,
#   n_models = 4,
#   file_suffix = '_normal_germany_65.rds')
# 
# ## refit the models which did not converge
# final_german_data_refit2 <- prepare_eval_data(
#   dates_list = list_of_dates[2:5],
#   prepared_data = data_w_lagged_values_germany,
#   n_models = 4,
#   file_suffix = '_normal_germany_65_new.rds')


# fit benchmark models ---------------------------------------------------------

CAR_results_nb_AR1 <- train_benchmark_models(
  mat_type = nb_b, model_type = 'CAR AR(1)',  suffix = ' 2', nr_ar = 1)

CAR_results_sd_AR1 <- train_benchmark_models(
  mat_type = nb_sd, model_type = 'CAR AR(1)', suffix = ' 3', nr_ar = 1)

CAR_results_bd_AR1 <- train_benchmark_models(
  mat_type = nb_bd, model_type = 'CAR AR(1)', suffix = ' 4', nr_ar = 1)

# saveRDS(CAR_results_nb_AR1, "CAR_results_nb_AR1.rds")
# saveRDS(CAR_results_sd_AR1, "CAR_results_sd_AR1.rds")
# saveRDS(CAR_results_bd_AR1, "CAR_results_bd_AR1.rds")


# train without de-trending ----------------------------------------------------

## train on the original scale of the data (without de-trending first)
CAR_results_nb_AR1_o <- train_benchmark_models(
  mat_type = nb_b, model_type = 'CAR AR(1)',  suffix = ' 2', 
  nr_ar = 1, make_stationary = FALSE)

CAR_results_sd_AR1_o <- train_benchmark_models(
  mat_type = nb_sd, model_type = 'CAR AR(1)', suffix = ' 3', 
  nr_ar = 1, make_stationary = FALSE)

CAR_results_bd_AR1_o <- train_benchmark_models(
  mat_type = nb_bd, model_type = 'CAR AR(1)', suffix = ' 4', 
  nr_ar = 1, make_stationary = FALSE)


## save predicitons
# saveRDS(CAR_results_nb_AR1_o, "CAR_results_nb_AR1_o.rds")
# saveRDS(CAR_results_sd_AR1_o, "CAR_results_sd_AR1_o.rds")
# saveRDS(CAR_results_bd_AR1_o, "CAR_results_bd_AR1_o.rds")


# desriptive statistics --------------------------------------------------------

## data for germany
merged_data_germany <- merge_raw_data(
  input_ils = ils_data,
  age_group = NULL, 
  federal_state = NULL)

# saveRDS(merged_data_germany, 'merged_data_germany.rds')

## add lagged values of a covariate for prediction later on
data_w_lagged_values_germany <- lag_covariates(
  input_data = merged_data_germany$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 14,
  cov_name = 'incidence'
)

## data for bavaria
merged_data_bavaria <- merge_raw_data(
  input_ils = ils_data,
  age_group = NULL, 
  federal_state = '09')

#saveRDS(merged_data_bavaria, 'merged_data_bavaria.rds')

## add lagged values of a covariate for prediction later on
data_w_lagged_values_bavaria <- lag_covariates(
  input_data = merged_data_bavaria$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 14,
  cov_name = 'incidence'
)

## select relevant phases
phases_tmp <- phases_data %>%
  filter(
    type == 'wave'#,
    #!is.na(cw_end)
  )

## insert artificial end date
max_date_rki <- rki_data$Meldedatum %>% as.Date() %>% max()
phases_tmp$end_date[nrow(phases_tmp)] <- max_date_rki

## get the dates from the phases for filtering the data
phase_dates <- NULL

for (j in seq(nrow(phases_tmp))) {
  
  phase_dates_tmp <- as.Date(
    phases_tmp$start_date[j]):as.Date(phases_tmp$end_date[j])
  
  assign(paste0('phase_', j), phase_dates_tmp)
  
  phase_dates <- c(phase_dates, phase_dates_tmp)
  
}

phase_dates <- as.Date(phase_dates)

## make summary statistics for rki data (germany)
rki_table_ger <- rki_tmp %>%
  #filter(Bundesland == 'Bayern') %>%
  mutate(
    phase = case_when(
      date %in% phase_1 ~ 1,
      date %in% phase_2 ~ 2,
      date %in% phase_3 ~ 3,
      date %in% phase_4 ~ 4,
      date %in% phase_5 ~ 5,
      date %in% phase_6 ~ 6
    )
  ) %>%
  filter(!is.na(phase)) %>%
  group_by(date, AgeGroups, phase) %>%
  summarise(
    case_count = sum(AnzahlFall)
  ) %>%
  group_by(phase, AgeGroups) %>%
  summarise(
    mean_value = mean(case_count, na.rm = TRUE),
    median_value = median(case_count, na.rm = TRUE),
    sd_value = sd(case_count, na.rm = TRUE),
    min_value = min(case_count, na.rm = TRUE),
    max__value = max(case_count, na.rm = TRUE)
  )

# rki_table_ger %>%
#   kbl() %>%
#   kable_styling()

## make summary statistics for rki data (bavaria)
rki_table_bay <- rki_tmp %>%
  filter(Bundesland == 'Bayern') %>%
  mutate(
    phase = case_when(
      date %in% phase_1 ~ 1,
      date %in% phase_2 ~ 2,
      date %in% phase_3 ~ 3,
      date %in% phase_4 ~ 4,
      date %in% phase_5 ~ 5,
      date %in% phase_6 ~ 6
    )
  ) %>%
  filter(!is.na(phase)) %>%
  group_by(date, AgeGroups, phase) %>%
  summarise(
    case_count = sum(AnzahlFall)
  ) %>%
  group_by(phase, AgeGroups) %>%
  summarise(
    mean_value = mean(case_count, na.rm = TRUE),
    median_value = median(case_count, na.rm = TRUE),
    sd_value = sd(case_count, na.rm = TRUE),
    min_value = min(case_count, na.rm = TRUE),
    max__value = max(case_count, na.rm = TRUE)
  )

# rki_table_bay %>%
#   kbl() %>%
#   kable_styling()

# prepare data for plotting ----------------------------------------------------

## plot the spatial units of interest as wel as the beds density
merged_data_germany$ils_sp %>% mapview(., zcol = 'beds_density')

## prepare rki data for plotting
age_group <- rki_data$Altersgruppe %>% unique()

rki_tmp <- rki_data %>%
  mutate(AGS_raw = as.character(IdLandkreis), 
         AGS = case_when(
           nchar(AGS_raw) == 4 ~ paste0("0", AGS_raw),
           nchar(AGS_raw) >= 5 ~ AGS_raw, 
         ),
         date = as.Date(Meldedatum)
  ) %>% 
  filter(
    Altersgruppe %in% all_of(age_group),
    NeuerFall %in% c(0, 1)
  ) %>%
  mutate(AgeGroups = as.factor(Altersgruppe)) 

age_dummies <- model.matrix( ~ 0 + rki_tmp[, 'AgeGroups']) %>%
  as.data.frame()

age_dummies <- rki_tmp$AnzahlFall * age_dummies

age_labels <- sapply(str_split(colnames(age_dummies), ']'), '[', 2)
age_labels <- gsub('-', '_', age_labels)
age_labels <- gsub('\\+', '', age_labels)
age_labels <- paste0('Cat_', age_labels)

colnames(age_dummies) <- age_labels
rki_tmp <- cbind(rki_tmp, age_dummies)

# stacked shares of covariates (plot) ------------------------------------------

p_ger <- make_stacked_share_plots(
  input_rki = rki_tmp,
  input_phases = phases_tmp,
  input_bay = FALSE,
  no_legend = TRUE)

p_bay <- make_stacked_share_plots(
  input_rki = rki_tmp,
  input_phases = phases_tmp,
  input_bay = TRUE)


stacked_shares_plot <- cowplot::plot_grid(p_ger, p_bay,
                                          ncol = 1,
                                          labels = c('A', 'B'))

# summary stats icu ------------------------------------------------------------

icu_table_ger <- data_w_lagged_values_germany$long %>% 
  mutate(
    phase = case_when(
      date %in% phase_1 ~ 1,
      date %in% phase_2 ~ 2,
      date %in% phase_3 ~ 3,
      date %in% phase_4 ~ 4,
      date %in% phase_5 ~ 5,
      date %in% phase_6 ~ 6
    )
  ) %>%
  filter(!is.na(phase)) %>%
  group_by(date, phase) %>%
  summarise(
    icu = sum(icu),
    incidence = sum(incidence)
  ) %>%
  group_by(phase) %>%
  summarise(
    mean_icu = mean(icu, na.rm = TRUE),
    median_icu = median(icu, na.rm = TRUE),
    sd_icu = sd(icu, na.rm = TRUE),
    min_icu = min(icu, na.rm = TRUE),
    max_icu = max(icu, na.rm = TRUE),
    
    mean_incidence = mean(incidence, na.rm = TRUE),
    median_incidence = median(incidence, na.rm = TRUE),
    sd_incidence = sd(incidence, na.rm = TRUE),
    min_incidence = min(incidence, na.rm = TRUE),
    max_incidence = max(incidence, na.rm = TRUE)
  )


# icu_table_ger %>%
#   kbl() %>%
#   kable_styling()


icu_table_bay <- data_w_lagged_values_bavaria$long %>% 
  mutate(
    phase = case_when(
      date %in% phase_1 ~ 1,
      date %in% phase_2 ~ 2,
      date %in% phase_3 ~ 3,
      date %in% phase_4 ~ 4,
      date %in% phase_5 ~ 5,
      date %in% phase_6 ~ 6
    )
  ) %>%
  filter(!is.na(phase)) %>%
  group_by(date, phase) %>%
  summarise(
    icu = sum(icu),
    incidence = sum(incidence)
  ) %>%
  group_by(phase) %>%
  summarise(
    mean_icu = mean(icu, na.rm = TRUE),
    median_icu = median(icu, na.rm = TRUE),
    sd_icu = sd(icu, na.rm = TRUE),
    min_icu = min(icu, na.rm = TRUE),
    max_icu = max(icu, na.rm = TRUE),
    
    mean_incidence = mean(incidence, na.rm = TRUE),
    median_incidence = median(incidence, na.rm = TRUE),
    sd_incidence = sd(incidence, na.rm = TRUE),
    min_incidence = min(incidence, na.rm = TRUE),
    max_incidence = max(incidence, na.rm = TRUE)
  )


# icu_table_bay %>%
#   kbl() %>%
#   kable_styling()


# parameters plot --------------------------------------------------------------

## helper function for string editing 
## https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## define prefix and suffix
file_prefix = './stan_fits/'
file_suffix = '_normal_bayern_rw2.rds'
#file_suffix = '_normal_bayern_rw2_MULTI.rds'

## select results from model 8.1
mod_8.1 <- readRDS(
  paste0(file_prefix, 'fit_0', 8, '_iter_', 1, "_normal_lognormal_bayern_rw2.rds")
)

## summarize model
s <- summary(mod_8.1, probs = c(0.025, 0.975))

## get parameter names and index
b <- s$summary %>%
  as.data.frame() %>%
  filter(
    #!grepl('mu_err', rownames(.), fixed = TRUE),
    grepl('15]', substrRight(rownames(.), 3), fixed = TRUE)
  ) %>%
  mutate(
    
    par_index = rownames(.),
    
    par_name = par_index %>%
      strsplit( "\\[" ) %>%
      sapply( "[", 1),
    
    time_index = par_index %>% 
      strsplit( "\\[" ) %>%
      sapply( "[", 2 ) %>%
      strsplit(",") %>%
      sapply("[", 1),
    
    # static_par = case_when(
    #   grepl(']', time_index) ~ 1,
    #   !grepl(']', time_index) ~ 0
    # ),
    # 
    index_count = par_index %>%
      str_count(., ','),
    
    time_index_3 = par_index %>%
      strsplit( "\\[" ) %>%
      sapply( "[", 2 ) %>%
      strsplit(",") %>%
      sapply("[", 2),
    
    model_index = par_index %>%
      strsplit( "\\[" ) %>%
      sapply( "[", 2 ) %>%
      strsplit(",") %>%
      sapply("[", 1),
    
    time_index = case_when(
      index_count == 2 ~ time_index_3,
      index_count == 1 ~ time_index
    ),
    
    model_suffix = case_when(
      (index_count == 2 & model_index == 1) ~ '_CC',
      par_name == 'eta_m' ~ '_CC',
      (index_count == 2 & model_index == 2) ~ '_AUX',
      par_name %in% c('b', 'b_err', 'eta_log') ~ '_AUX',
      (index_count == 2 & model_index == 3) ~ '_ICU',
      par_name == 'eta_obs' ~ '_ICU'
    ),
    
    par_name = case_when(
      index_count == 1 ~ par_name,
      index_count == 2 ~ paste0(par_name, model_suffix),
    ),
    
    par_err = grepl('_err', par_name),
    
    par_name = gsub('_err', '', par_name),
    
    par_eta = par_name %in% c('eta_m', 'eta_log', 'eta_obs'),
    
    par_name = case_when(
     par_err == TRUE ~ paste0(par_name, '_err'),
     par_err == FALSE ~ par_name
    ),
    
    par_name = case_when(
      par_eta == TRUE ~ paste0('a', par_name),
      par_eta == FALSE ~ par_name
    )
    
  ) %>%
  group_by(par_name) %>%
  mutate(
    static_par = case_when(
      max(as.numeric(time_index)) != 44 ~ 1,
      max(as.numeric(time_index)) == 44 ~ 0
    )
  )

## table with static parameter values
tab_8.1 <- b %>%
  filter(static_par == 1) %>%
  dplyr::select(mean, `2.5%`, `97.5%`, model_suffix) %>%
  kable() %>%
  kable_styling()

## plot ICU
p_ICU_81 <- b %>%
  filter(static_par == 0) %>%
  filter(model_suffix == '_ICU') %>%
  mutate(time_index = as.numeric(time_index)) %>%
  filter(par_name != "nu") %>%
  ggplot(., aes(x = time_index, y = mean)) +
  #geom_point() +
  geom_line(aes(group = par_name)) +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  group = par_name,
                  fill = par_name), alpha = 0.3) +
  geom_vline(xintercept = 30, linetype="dashed") +
  facet_wrap( ~ par_name, scales = "free") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    subtitle = 'Model ID: 8.1, Iteration: 1, Unit: Leitstelle München' ,
    x = ' ', y = " ",
    title = ' '
  ) +
  scale_fill_discrete(
    name = 'Parameters',
    labels = c(
      parse(text = paste("eta[", 'ICU', "]")),
      parse(text = paste("mu[", 'ICU', "]")),
      parse(text = paste("epsilon[", 'mu[', 'ICU',']', "]")),
      parse(text = paste("mu[", 'ICU[', 'RW2',']', "]")),
      parse(text = paste("varphi[", 'ICU', "]")),
      parse(text = paste("nu[", 'ICU', "]")),
      parse(text = paste("epsilon[", 'nu[', 'ICU', ']', "]"))
    )
  ) 

## plot AUX
p_AUX_81 <- b %>%
  filter(static_par == 0) %>%
  filter(model_suffix == '_AUX') %>%
  mutate(time_index = as.numeric(time_index)) %>%
  ggplot(., aes(x = time_index, y = mean)) +
  #geom_point() +
  geom_line(aes(group = par_name)) +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  group = par_name,
                  fill = par_name), alpha = 0.3) +
  geom_vline(xintercept = 30, linetype="dashed") +
  facet_wrap( ~ par_name, scales = "free") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    subtitle = 'Model ID: 8.1, Iteration: 1, Unit: Leitstelle München' ,
    x = ' ', y = "Median",
    title = ' '
  ) +
  scale_fill_discrete(
    name = 'Parameters',
    labels = c(
      parse(text = paste("eta[", 'AUX', "]")),
      parse(text = paste("gamma[", 'AUX', "]")),
      parse(text = paste("epsilon[", 'gamma[', 'AUX', ']', "]")),
      parse(text = paste("mu[", 'ICU', "]")),
      parse(text = paste("epsilon[", 'mu[', 'AUX', ']', "]")),
      parse(text = paste("mu[", 'AUX[', 'RW2',']', "]")),
      parse(text = paste("varphi[", 'AUX', "]")),
      parse(text = paste("nu[", 'AUX', "]")),
      parse(text = paste("epsilon[", 'nu[', 'AUx', ']', "]"))
    )
  ) 

## plot CC
p_CC_81 <- b %>%
  filter(static_par == 0) %>%
  filter(model_suffix == '_CC') %>%
  mutate(time_index = as.numeric(time_index)) %>%
  filter(par_name != "nu") %>%
  ggplot(., aes(x = time_index, y = mean)) +
  #geom_point() +
  geom_line(aes(group = par_name)) +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  group = par_name,
                  fill = par_name), alpha = 0.3) +
  facet_wrap( ~ par_name, scales = "free") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    subtitle = 'Model ID: 8.1, Iteration: 1, Unit: Leitstelle München' ,
    x = 'Time index', y = " ",
    title = ' '
  )  +
  scale_fill_discrete(
    name = 'Parameters',
    labels = c(
      parse(text = paste("eta[", 'CC', "]")),
      parse(text = paste("mu[", 'CC', "]")),
      parse(text = paste("epsilon[", 'mu[', 'CC', ']', "]")),
      parse(text = paste("mu[", 'CC[', 'RW2',']', "]")),
      parse(text = paste("varphi[", 'CC', "]")),
      parse(text = paste("nu[", 'CC', "]")),
      parse(text = paste("epsilon[", 'nu[', 'CC', ']', "]"))
    )
  ) 


## unite plots
p_81 <- cowplot::plot_grid(p_ICU_81, p_AUX_81, p_CC_81, 
                   labels = c(
                     'Part: ICU',
                     'Part: AUX',
                     'Part: CC'),
                   ncol = 1)
  