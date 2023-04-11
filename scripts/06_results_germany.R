
#'@details procedure for computing aggregate forecasts (bottom up)

# set up -----------------------------------------------------------------------

ils_raw = ils_data
file_prefix = './stan_fits/'
county_list = c(paste0('0', 1:9), as.character(10:16))
dates_list = list_of_dates
prepared_data = merged_data_germany
n_models = 1:4
file_suffix = '_normal_germany_65.rds'
input_divi = divi_data

## set up empty results data frame of correct dimensionality
d <- length(dates_list)
z <- length(n_models)
o <- length(county_list) + 1

results_out <- matrix(NA, nrow = 44 * d * z * o , ncol = 8) %>%
  as.data.frame()

# functions --------------------------------------------------------------------

## helper function for finding the next index to fill the matrix iteratively
find_min_index <- function(x) { min(which(rowSums(is.na(x)) == 8)) }

# aggregate samples ------------------------------------------------------------

for (j in seq(length(dates_list))) {
  
  for (k in n_models) {
    
    mod_tmp <- NULL
    pred_tmp <- NULL
    
    model_name <- paste0(file_prefix, 'fit_0', k, '_iter_', j, file_suffix)
    
    ## load fitted models
    mod_tmp <- readRDS(model_name)
    
    pred_tmp <- aggregate_samples(
      input_mod = force(mod_tmp),
      length_ts = 44,
      input_merged_data = prepared_data,
      input_ils = ils_raw,
      input_county = NULL,
      input_date = dates_list[j]
    )
    
    pred_tmp$iter <- j
    pred_tmp$model_nr <- k
    pred_tmp$county <- 'Germany'
    
    a <- find_min_index(force(results_out))
    
    ## update results
    results_out[a:(a + 43), ] <- pred_tmp
    
    if (!is.null(county_list)) {
      
      for (c in seq(length(county_list))) {
        
        pred_tmp_c <- NULL
        
        print(paste0('bundesland ', c, ' iter ', j, ' model ', k))
        
        pred_tmp_c <- aggregate_samples(
          input_mod = force(mod_tmp),
          length_ts = 44,
          input_merged_data = prepared_data,
          input_ils = ils_raw,
          input_county = county_list[c],
          input_date = dates_list[j]
        )
        
        pred_tmp_c$iter <- j
        pred_tmp_c$model_nr <- k
        pred_tmp_c$county <- county_list[c]
        
        a <- find_min_index(force(results_out))
        
        ## update results
        results_out[a:(a + 43), ] <- pred_tmp_c
        
      }
      
    }
    
  } 
  
}

saveRDS(results_out, 'german_combined_samples.rds')

# finalize results -------------------------------------------------------------

results_out <- readRDS('german_combined_samples.rds')

colnames(results_out) <- colnames(pred_tmp)
results_out$date <- as.Date(results_out$date)
results_out$lower[results_out$lower < 0] <- 0
results_out$upper[results_out$upper < 0] <- 0
results_out$pred[results_out$pred < 0] <- 0


# add observed values ----------------------------------------------------------

## prepare divi data
obs_tmp <- input_divi %>%
  mutate(
    county = case_when(
      bundesland < 10 ~ paste0('0', bundesland),
      bundesland >= 10 ~ paste0(bundesland)
    )
  )

obs_tmp <- obs_tmp %>%
  group_by(date, county) %>%
  summarise(obs = sum(faelle_covid_aktuell, na.rm = TRUE)) %>%
  rbind(
    .,
    obs_tmp %>%
      group_by(date) %>%
      summarise(obs = sum(faelle_covid_aktuell, na.rm = TRUE)) %>%
      mutate(county = 'Germany')
  ) %>%
  mutate(date = as.Date(date))

## join with our results
joined_tmp <- input_aggregated %>%
  left_join(
    .,
    obs_tmp,
    by = c('date', 'county')
  ) 

long_tmp <- make_results_w_intervals(
  input_loss_data = joined_tmp
)

long_tmp <- long_tmp %>%
  rename(ILS = county)


# make plots -------------------------------------------------------------------

plot_forecasts(
  long_tmp %>%
    filter(iter == 2, model_nr == 2),
  add_intervals = TRUE,
  ILS = 'Germany')


plot_forecasts(
  long_tmp %>%
    filter(iter == 2, model_nr == 2),
  add_intervals = TRUE,
  ILS = 'Germany')


p3_germany <- plot_all_results(
  input_w_intervals =  long_tmp,
  time_ids = seq(20),
  model_ids = 4,
  input_ils = 'Germany',
  model_suffix = '.1'
)

p4_bundeslaender <- plot_all_results(
  input_w_intervals =  long_tmp %>% filter(ILS != 'Germany'),
  time_ids = 1,
  model_ids = 4,
  input_ils = NULL,
  model_suffix = '.1'
)

# make tables ------------------------------------------------------------------

## overall
train_results <- joined_tmp %>%
filter(
  time_index <= 30,
  county == 'Germany'
  ) %>%
filter(!is.na(obs)) %>%
  group_by(model_nr) %>%
  summarise(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
         MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
         RMSE_mean = RMSE / mean(obs),
         MAE_mean = MAE / mean(obs)) 

test_results <- joined_tmp %>%
  filter(
    time_index > 30,
    county == 'Germany'
    ) %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr) %>%
  summarise(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
            MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
            RMSE_mean = RMSE / mean(obs),
            MAE_mean = MAE / mean(obs)) %>%
  magrittr::set_colnames(paste0(colnames(.), '_test'))


aggregated_results_germany <- train_results %>%
  full_join(
    .,
    test_results %>%
      rename(model_nr = model_nr_test),
    by = c('model_nr')
  ) %>%
  kbl() %>%
  kable_styling()

aggregated_results_germany


## aggregated by federal state
train_results <- joined_tmp %>%
  filter(time_index <= 30,
         county != 'Germany') %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr) %>%
  summarise(
    RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
    MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
    RMSE_mean = RMSE / mean(obs),
    MAE_mean = MAE / mean(obs)) 

test_results <- joined_tmp %>%
  filter(time_index > 30,
         county != 'Germany') %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr) %>%
  summarise(
    RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
    MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
    RMSE_mean = RMSE / mean(obs),
    MAE_mean = MAE / mean(obs)) %>%
  magrittr::set_colnames(paste0(colnames(.), '_test'))


aggregated_results_germany <- train_results %>%
  full_join(
    .,
    test_results %>%
      rename(model_nr = model_nr_test),
    by = c('model_nr')
  ) %>%
  kbl() %>%
  kable_styling()

aggregated_results_germany



## aggregated by federal state
train_results <- joined_tmp %>%
  filter(time_index <= 30,
         county != 'Germany') %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr, county) %>%
  summarise(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
            MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
            RMSE_mean = RMSE / mean(obs),
            MAE_mean = MAE / mean(obs)) 

test_results <- joined_tmp %>%
  filter(time_index > 30,
         county != 'Germany') %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr, county) %>%
  summarise(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
            MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
            RMSE_mean = RMSE / mean(obs),
            MAE_mean = MAE / mean(obs)) %>%
  magrittr::set_colnames(paste0(colnames(.), '_test'))


aggregated_results_germany <- train_results %>%
  full_join(
    .,
    test_results %>%
      rename(model_nr = model_nr_test,
             county = county_test),
    by = c('model_nr', 'county')
  ) %>%
  filter(model_nr == 4) %>%
  kbl() %>%
  kable_styling()

aggregated_results_germany


# plot showing the dependence on the number of units ---------------------------

## enable aggregation by county
ILS_ID_list <-  merged_data_germany$icu_data[, c('ILS', 'ILS_ID')] %>% 
  distinct() %>%
  left_join(.,
            ils_data %>%
              mutate(
                AGS = as.character(AGS),
                AGS = case_when(
                  nchar(AGS) < 5 ~ paste0('0', AGS),
                  nchar(AGS) >= 5 ~ paste0(AGS)
                ),
                County = substr(AGS, 1, 2)
              ),
            by = c('ILS')
  )

a <- ILS_ID_list %>%
  distinct(ILS, County) %>%
  group_by(County) %>%
  count(County)

## read states
federal_states <- read.csv("D:/MA_HMM/federal_states.csv", sep=";")

## join data
plot_data_tmp <- train_results %>%
  full_join(
    .,
    test_results %>%
      rename(model_nr = model_nr_test,
             county = county_test),
    by = c('model_nr', 'county')
  ) %>%
  filter(model_nr == 4) %>%
  full_join(
    .,
    ILS_ID_list %>%
      distinct(ILS, County) %>%
      group_by(County) %>%
      count(County) %>%
      rename(
        county = County,
        n_ils = n
      ),
    by = c('county')
  ) %>%
  ungroup() %>%
  dplyr::select(RMSE_mean_test, MAE_mean_test, n_ils, county) %>%
  pivot_longer(., !c('county', 'n_ils'), 
               names_to = 'loss_metric', values_to = 'value') %>%
  right_join(
    .,
    federal_states %>%
      mutate(
        county = gsub('"', '', county),
        ),
    by = c('county')
  )


plot_data_tmp %>%
  ggplot(., aes(x = n_ils, y = value, color = loss_metric)) +
  geom_point() +
  geom_smooth(method = 'gam') +
  geom_text(
    label = plot_data_tmp$fs_name, 
    #nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = TRUE
  ) +
  theme_minimal() +
  scale_colour_discrete(
    name = 'Normalized loss metrics',
    labels = c('nMAE', 'nRMSE')) +
  labs(
    subtitle = 'Results on held out data',
    x = 'number of spatial units per federal state', y = "value of loss metric (nMAE or nRMSE)",
    title = 'Normalized loss by number of spatial units'
  )
  
