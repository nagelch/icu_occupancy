
#'@details procedure for combining the estimates from the state space models to
#'an ensemble. This procedure calculates a naive median ensemble. Furthermore,
#'it uses the state-space models as base learners, pooled via a random forest
#'as the meta learner.


library(ranger)

# median ensemble --------------------------------------------------------------

## add a median forecast
results_median <- final_bay_data %>%
    filter(type == 'fitted') %>%
    group_by(id, time_index, iter, date, ILS) %>%
    summarise(pred = median(case_count),
              lower = median(lower),
              upper = median(upper)) %>%
    mutate(model_nr = 'median ensemble') %>%
  full_join(
    .,
    final_bay_data %>%
      filter(type == 'observed') %>%
      dplyr::select(id, time_index, case_count, iter) %>%
      distinct() %>%
      rename(obs = case_count),
    by = c('id', 'time_index', 'iter')
  )


## ensemble --------------------------------------------------------------------

#'@details we prepare the ensemble data

ensemble_tmp <- final_bay_data %>%
  filter(type == 'fitted') %>%
  full_join(
    .,
    final_bay_data %>%
      filter(type == 'observed') %>%
      rename(obs = case_count) %>%
      dplyr::select(obs, time_index, date, ILS, model_nr, iter),
    by = c('time_index', 'date', 'ILS', 'model_nr', 'iter')
  ) 

ensemble_final <- ensemble_tmp %>%
  dplyr::select(-date, -lower, -upper, -type) %>%
  mutate(model_nr = paste0('mod_pred_', model_nr)) %>%
  group_by(model_nr) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = model_nr, values_from = case_count) %>%
  dplyr::select(-row) %>%
  mutate(hold_out = case_when(
    time_index <= 30 ~ 0,
    time_index > 30 ~ 1
  )) %>% 
  full_join(
    .,
    ensemble_tmp %>%
      dplyr::select(-date, -case_count, -upper, -type) %>%
      mutate(model_nr = paste0('mod_lower_', model_nr)) %>%
      group_by(model_nr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = model_nr, values_from = lower) %>%
      dplyr::select(-row),
    by = c('id', 'time_index', 'ILS', 'iter', 'obs')
  ) %>%
  full_join(
    .,
    ensemble_tmp %>%
      dplyr::select(-date, -lower, -case_count, -type) %>%
      mutate(model_nr = paste0('mod_upper_', model_nr)) %>%
      group_by(model_nr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = model_nr, values_from = upper) %>%
      dplyr::select(-row),
    by = c('id', 'time_index', 'ILS', 'iter', 'obs')
  ) 

# random forest ----------------------------------------------------------------

results_ensemble <- NULL

for (i in seq(max(ensemble_final$iter))) {
  
  if (i > 1) {
    
    ## train on past data only
    rf <- ensemble_final %>%
      filter(iter < i,
             hold_out == 1) %>%
      dplyr::select(-id, -ILS, -iter, -hold_out) %>%
      ranger(obs ~ ., data = ., num.trees =  100, quantreg = TRUE)
    
    ## predict on next data in date_list
    reuslts_tmp <- ensemble_final %>%
      filter(iter == i)
    
    pred <- reuslts_tmp %>%
      dplyr::select(-id, -ILS, -iter, -obs) %>%
      predict(rf, ., type = "quantiles",  quantiles = c(0.05, 0.5, 0.95))
    
    pred_tmp <- pred$predictions %>%
      as.data.frame() %>%
      magrittr::set_colnames(c('lower', 'pred', 'upper'))
    
    ## add variables
    results_tmp <- cbind(
      reuslts_tmp %>%
        dplyr::select(id, time_index, ILS, iter, obs)
      , pred_tmp) %>%
      mutate(model_nr = 'ssm ensemble')  
    
    ## update results
    results_ensemble <- rbind(results_ensemble, results_tmp)
    
  }
  
}

# save ensembles ---------------------------------------------------------------

## combine ensembles
all_ensembles <- rbind(
  results_median,
  results_ensemble
  )

## add missing date variable back to data
all_ensembles <- all_ensembles %>%
  ungroup() %>%
  dplyr::select(-date) %>%
  full_join(
    .,
    full_join(
      .,
      final_bay_data %>%
        filter(type == 'observed') %>%
        dplyr::select(id, time_index, iter, date) %>%
        distinct(),
      by = c('id', 'time_index', 'iter')
    )
  )