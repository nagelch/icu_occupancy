
#'@details Calcualtes lsos metrics and interval coverage from extracted and
#'pre-processed Stan results

# load data --------------------------------------------------------------------

## load the results (_germany, _bay, etc.)
## final_bay_data <- readRDS('final_bay_data.rds')

final_data <- final_bay_data

# prepare data -----------------------------------------------------------------

## prepare data
loss_data <- final_data %>%
  filter(type == 'observed') %>%
  dplyr::select(-lower, -upper, - type) %>%
  rename(obs = case_count) %>%
  full_join(
    final_data %>%
      filter(type == 'fitted') %>%
      dplyr::select(- type) %>%
      rename(pred = case_count),
    by = c('id', 'time_index', 'date', 'ILS', 'model_nr', 'iter')
  ) 

# ensembles --------------------------------------------------------------------

## Run this code for evaluating the results of the ensemble

# source('./05_prepare_ensembles.R')
# 
# ## add ensemble results to the loss data
# loss_data <- rbind(
#   loss_data %>%
#     mutate(model_nr = paste0('ssm ', model_nr)), 
#   all_ensembles) %>%
#   mutate(model_nr = as.factor(model_nr)
#   )

# prepare out of sample and within sample data sets ----------------------------

## load benchmark models (de-trended and original scale)
CAR_results_nb_AR1 <- readRDS("CAR_results_nb_AR1.rds")
CAR_results_sd_AR1 <- readRDS("CAR_results_sd_AR1.rds")
CAR_results_bd_AR1 <- readRDS("CAR_results_bd_AR1.rds")

CAR_results_nb_AR1_o <- readRDS("CAR_results_nb_AR1_o.rds")
CAR_results_sd_AR1_o <- readRDS("CAR_results_sd_AR1_o.rds")
CAR_results_bd_AR1_o <- readRDS("CAR_results_bd_AR1_o.rds")

## add benchmark model results to the data
loss_data <- rbind(
  loss_data,
  CAR_results_nb_AR1,
  CAR_results_sd_AR1,
  CAR_results_bd_AR1,
  CAR_results_nb_AR1_o,
  CAR_results_sd_AR1_o,
  CAR_results_bd_AR1_o
)


## out of sample
oos_data <- loss_data %>%
  filter(time_index > 30) %>%
  mutate(time_index = time_index - 30,
         model_nr = as.factor(model_nr))

## within sample
wis_data <- loss_data %>%
  filter(time_index <= 30) %>%
  mutate(model_nr = as.factor(model_nr))

## ex-post identify models which did not converge and exclude them
## only applies to fits from the bavarian data.
list_of_outliers <- oos_data %>%
  group_by(model_nr) %>%
  mutate(res = abs(obs - pred)) %>%
  mutate(res_scaled = scale(res)[, 1]) %>%
  ungroup() %>%
  mutate(
    outlier = case_when(
      res_scaled > 10 ~ 1,
      res_scaled <= 10 ~ 0
    )
  ) %>%  filter(
    outlier == 1
  ) %>%
  dplyr::select(model_nr, iter) %>%
  mutate(out_tmp = paste0(model_nr, '_', iter)) %>%
  distinct()

wis_data <- wis_data %>%
  mutate(out_tmp = paste0(as.character(model_nr), '_', iter)) %>%
  filter(!(out_tmp %in% list_of_outliers$out_tmp))

oos_data <- oos_data %>%
  mutate(out_tmp = paste0(as.character(model_nr), '_', iter)) %>%
  filter(!(out_tmp %in% list_of_outliers$out_tmp))

## metrics
oos_tmp <- oos_data %>%
  #filter(iter == max(iter)) %>%
  #filter(time_index != 8) %>%
  filter(!is.na(obs)) %>%
  group_by(model_nr) %>%
  mutate(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
         MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2)) 


## when summing up results for Germany set na.rm = TRUE
oos_h <- oos_data %>%
  #filter(iter == max(iter)) %>%
  group_by(model_nr, time_index) %>%
  mutate(RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
         MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2)) 

oos_i <- oos_data %>%
  #filter(iter == max(iter)) %>%
  group_by(model_nr, iter) %>%
  mutate(RMSE = round(sqrt(mean((pred - obs)^2)), digits = 2),
         MAE = round(median(abs(pred - obs)), digits = 2)) 

wis_tmp <- wis_data %>%
  group_by(model_nr, iter) %>%
  mutate(RMSE_fitted = round(sqrt(mean((pred - obs)^2)), digits = 2),
         MAE_fitted = round(median(abs(pred - obs)), digits = 2)) %>%
  dplyr::select(RMSE_fitted, MAE_fitted, model_nr, iter) %>%
  distinct(.)

wis_h <- wis_data %>%
  group_by(model_nr, time_index) %>%
  mutate(RMSE_fitted = round(sqrt(mean((pred - obs)^2)), digits = 2),
         MAE_fitted = round(median(abs(pred - obs)), digits = 2)) %>%
  dplyr::select(RMSE_fitted, MAE_fitted, model_nr, iter) %>%
  distinct(.)


# tables -----------------------------------------------------------------------


## make a table showing the RMSE and MAE on the training and hold out data
train_results <- wis_tmp %>%
  ungroup() %>%  
  rename(RMSE = RMSE_fitted,
         MAE = MAE_fitted) %>%
  dplyr::select(model_nr, RMSE, MAE) %>%
  distinct(model_nr, .keep_all = TRUE)

test_results <- oos_tmp %>%
  dplyr::select(model_nr, RMSE, MAE) %>%
  distinct(model_nr, .keep_all = TRUE) %>%
  magrittr::set_colnames(paste0(colnames(.), '_test'))


calibration_results_bavaria <- train_results %>%
  full_join(
    .,
    test_results %>%
      rename(model_nr = model_nr_test),
    by = c('model_nr')
  ) %>%
  kbl() %>%
  kable_styling()

## models excluded from the evaluation (see above, only bavarian data)
kept_out_models <- list_of_outliers %>%
   dplyr::select(-out_tmp) %>%
  kbl() %>%
  kable_styling()


# calculate interval coverage --------------------------------------------------

## interval coverage
intervals_tmp <- oos_tmp %>%
  full_join(
    ., 
    wis_tmp, 
    by = c('iter', 'model_nr')
  ) %>%
  mutate(
    h = time_index,
    sigma_h = RMSE_fitted * sqrt(h * (1 + h / max(h))),
    upper_f = pred + 1.96 * sigma_h,
    lower_f = pred - 1.96 * sigma_h,
    lower_f = case_when(
      lower_f < 0 ~ 0,
      lower_f >= 0 ~ lower_f
    ),
    inside_bayes = case_when(
      ((lower <= obs) & (upper >= obs)) ~ 1,
      !((lower <= obs) & (upper >= obs)) ~ 0),
    inside_freq = case_when(
      ((lower_f <= obs) & (upper_f >= obs)) ~ 1,
      !((lower_f <= obs) & (upper_f >= obs)) ~ 0)
  )


pred_intervals_data <- rbind(
  final_data %>%
    filter(!(time_index > 30 & type == 'fitted')),
  intervals_tmp %>%
    mutate(
      case_count = pred,
      lower = lower_f,
      upper = upper_f,
      type = 'fitted',
      model_nr = model_nr,
      time_index = time_index + 30
    ) %>%
    dplyr::select(colnames(final_data))) %>%
  arrange(id, type, model_nr, iter, time_index) %>%
  mutate(type = fct_recode(type, "predicted" = "fitted"))


intervals_tmp <- intervals_tmp %>%
  group_by(model_nr, time_index) %>%
  dplyr::summarize(
    coverage_bayes = mean(inside_bayes, na.rm = TRUE),
    coverage_freq = mean(inside_freq, na.rm = TRUE))

## check NAs
map(oos_tmp, ~sum(is.na(.)))

# plot intervals ---------------------------------------------------------------

intervals_plot_tmp <- intervals_tmp %>%
  filter(!grepl('MULTI', model_nr)) %>%
  filter(!grepl('CAR', model_nr)) %>%
  #mutate(model_nr = as.character(gsub('MULTI_', '', model_nr))) %>%
  mutate(
    has_covariates = case_when(
      model_nr %in% c(as.character(seq(4))) ~ TRUE,
      model_nr %in% c(as.character(5:8)) ~ FALSE
      # grepl('1 | 2 | 3 | 4', model_nr) ~ TRUE,
      # grepl('5 | 6 | 7 | 8', model_nr)  ~ FALSE
    )
  ) %>%
  mutate(model_nr = paste0(model_nr, '.1'))

p1_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_bayes)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), breaks = seq(14)) +
  labs(title = 'Interval coverage', subtitle = 'Credible intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') +
  theme(legend.position = "none")


p2_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_freq)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), breaks = seq(14)) +
  labs(title = ' ', subtitle = 'Frequentist prediction intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') 

grobs <- ggplotGrob(p2_coverage)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

p2_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_freq)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), 
                     breaks = seq(14)) +
  labs(title = ' ', subtitle = 'Frequentist prediction intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') +
  theme(legend.position = "none")

pgrid_coverage <- cowplot::plot_grid(p1_coverage, p2_coverage)
p_coverage_1 <- cowplot::plot_grid(
  pgrid_coverage, legend, ncol = 2, rel_widths = c(1, .1))

saveRDS(p_coverage_1, 'p_coverage_1.rds')
p_coverage_1 <- readRDS('p_coverage_1.rds')
#p_coverage_1 <- readRDS('p_coverage_1_GERMANY.rds')


# once again for MULTI ---------------------------------------------------------

intervals_plot_tmp <- intervals_tmp %>%
  filter(grepl('MULTI', model_nr)) %>%
  filter(!grepl('CAR', model_nr)) %>%
  mutate(model_nr = as.character(gsub('MULTI_', '', model_nr))) %>%
  mutate(
    has_covariates = case_when(
      model_nr %in% c(as.character(seq(4))) ~ TRUE,
      model_nr %in% c(as.character(5:8)) ~ FALSE
      # grepl('1 | 2 | 3 | 4', model_nr) ~ TRUE,
      # grepl('5 | 6 | 7 | 8', model_nr)  ~ FALSE
    )
  ) %>%
  mutate(model_nr = paste0(model_nr, '.2'))


p1_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_bayes)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), breaks = seq(14)) +
  labs(title = 'Interval coverage', subtitle = 'Credible intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') +
  theme(legend.position = "none")


p2_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_freq)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), breaks = seq(14)) +
  labs(title = ' ', subtitle = 'Frequentist prediction intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') 

grobs <- ggplotGrob(p2_coverage)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

p2_coverage <- intervals_plot_tmp %>%
  ggplot(., aes(x = time_index, y = coverage_freq)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), 
                     breaks = seq(14)) +
  labs(title = ' ', subtitle = 'Frequentist prediction intervals',
       x = 'steps ahead', y = 'share',
       color = 'Model ID') +
  theme(legend.position = "none")

pgrid_coverage <- cowplot::plot_grid(p1_coverage, p2_coverage)
p_coverage_2 <- cowplot::plot_grid(
  pgrid_coverage, legend, ncol = 2, rel_widths = c(1, .1))

saveRDS(p_coverage_2, 'p_coverage_2.rds')

p_coverage_2 <- readRDS('p_coverage_2.rds')


# loss plots -------------------------------------------------------------------

##https://doi.org/10.1002/bimj.201900351
##https://robjhyndman.com/hyndsight/narrow-pi/
loss_plot_tmp <- oos_h %>%
  filter(!grepl('CAR', model_nr)) %>%
  mutate(
    has_covariates = case_when(
      model_nr %in% c(as.character(seq(4))) ~ TRUE,
      model_nr %in% c(as.character(5:8)) ~ FALSE
    )
  )

loss_plot_tmp <- loss_plot_tmp %>%
  mutate(model_nr = paste0(model_nr, '.1'))

## loss
p1_loss <- loss_plot_tmp %>%
  ggplot(., aes(x = time_index, y = RMSE)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  theme_bw() +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)),
                     breaks = seq(14)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  labs(title = 'Loss per forecast step',
       x = 'steps ahead', y = 'RMSE',
       color = 'Model ID') +
  theme(legend.position = "none")


p2_loss <- loss_plot_tmp %>%
  ggplot(., aes(x = time_index, y = MAE)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  theme_bw() +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)),
                     breaks = seq(14)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  #scale_y_continuous(limits = c(0, 20)) +
  labs(title = ' ',
       x = 'steps ahead', y = 'MAE',
       color = 'Model ID') 


grobs <- ggplotGrob(p2_loss)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

p2_loss <- loss_plot_tmp %>%
  ggplot(., aes(x = time_index, y = MAE)) +
  geom_line(aes(color = model_nr)) +
  geom_point(aes(shape = has_covariates, color = model_nr)) +
  theme_bw() +
  scale_x_continuous("steps ahead", labels = as.character(seq(14)), breaks = seq(14)) +
  scale_shape_discrete(
    name = 'Covariates',
    labels = c('no', 'yes')) +
  #scale_y_continuous(limits = c(0, 20)) +
  labs(title = ' ',
       x = 'steps ahead', y = 'MAE',
       color = 'Model ID') +
  theme(legend.position = "none")

pgrid_loss <- cowplot::plot_grid(p1_loss, p2_loss)
p_loss <- cowplot::plot_grid(
  pgrid_loss, legend, ncol = 2, rel_widths = c(1, .1))


saveRDS(p_loss, 'p_loss_GERMANY.rds')
p_loss <- readRDS('p_loss_GERMANY.rds')


# individual forecasts plots ---------------------------------------------------

## plot indi
interval_results <-  make_results_w_intervals()

p1_all <- plot_all_results(
  input_w_intervals =  interval_results,
  time_ids = seq(20),
  model_ids = 4,
  model_suffix = '.1'
)

p1_all_covariates <- plot_all_results(
  input_w_intervals =  interval_results,
  time_ids = seq(20),
  model_ids = 8,
  model_suffix = '.1'
)

p2_all <- plot_all_results(
  input_w_intervals =  interval_results,
  time_ids = 3,
  model_ids = seq(8),
  model_suffix = '.1'
)

p3_all <- plot_all_results(
  input_w_intervals =  interval_results,
  time_ids = 4,
  model_ids = 4,
  input_ils = NULL,
  model_suffix = '.1'
)