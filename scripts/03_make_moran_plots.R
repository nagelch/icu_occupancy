
#'@details procedure for generating the autocorrelation plots of
#'ICU occupancy and case counts

# set up -----------------------------------------------------------------------

## install ape package
install.packages('ape')

## get data for germany
merged_data_germany <- merge_raw_data(
  age_group = NULL, 
  federal_state = NULL)


## add lagged values of a covariate for prediction later on
data_w_lagged_values_germany <- lag_covariates(
  input_data = merged_data_germany$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 14,
  cov_name = 'incidence'
)


# set up the neighborhood matrices --------------------------------------------

## Three types of matrices, spatial distance (SD), adjacency (B) and distance
## in bed density (BD)
nb_mat_sd_germany <- make_adjacency_matrix(
  input_sp = merged_data_germany$ils_sp,
  weights_style = 'SD'
)

nb_mat_b_germany <- make_adjacency_matrix(
  input_sp = merged_data_germany$ils_sp,
  weights_style = 'B'
)

nb_mat_bd_germany <- make_adjacency_matrix(
  input_sp = merged_data_germany$ils_sp,
  weights_style = 'BD'
)

## turn into matrix objects
nb_b_germany <- as.matrix(as.data.frame(nb_mat_b_germany))
nb_sd_germany <- as.matrix(as.data.frame(nb_mat_sd_germany))
nb_bd_germany <- as.matrix(as.data.frame(nb_mat_bd_germany))

## cut off at 60 km
nb_sd_germany[nb_sd_germany > 60] <- 0

## round and offset nb_bd
nb_bd_germany <- round(nb_bd_germany, digits = 3) + 1

## everything zero here should be zero in the other matrix as well
nb_bd_germany[nb_sd_germany == 0] <- 0

## scale the results
nb_sd_germany <- scale_matrix(nb_sd_germany)
nb_bd_germany <- scale_matrix(nb_bd_germany)

## inverse of the entries
nb_sd_germany <- nb_sd_germany^-1
nb_bd_germany <- nb_bd_germany^-1

## set Inf to zero
nb_sd_germany[nb_sd_germany == Inf] <- 0
nb_bd_germany[nb_bd_germany == Inf] <- 0

# Plot phases of spatial autocorrelation (ICU) ---------------------------------

## prepare Moran's I data
dat_sd <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long,  
  input_mat = nb_sd_germany,
  win_size = 14,
  var_name = 'icu'
)

dat_bd <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long,  
  input_mat = nb_bd_germany,
  win_size = 14,
  var_name = 'icu'
)

dat_nb <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long,  
  input_mat = nb_b_germany,
  win_size = 14,
  var_name = 'icu'
)

## make rectangles
rect_sd <-  extract_sig_phase(input_data = dat_sd)
rect_bd <-  extract_sig_phase(input_data = dat_bd)
rect_nb <-  extract_sig_phase(input_data = dat_nb)


## plots
p_moran_sd <- ggplot(data = dat_sd, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_sd, 
    inherit.aes = FALSE, 
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Spatical distance',
       x = 'Date', y = "Moran's I")

p_moran_bd <- ggplot(data = dat_bd, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_bd, 
    inherit.aes=FALSE, 
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Distance in structural quality',
       x = 'Date', y = "Moran's I")

p_moran_nb <- ggplot(data = dat_nb, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_nb, 
    inherit.aes=FALSE, 
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Adjacency',
       x = 'Date', y = "Moran's I")

dat_icu <- data_w_lagged_values_germany$long %>%
  group_by(date) %>%
  summarise(sum_icu = sum(icu, na.rm = TRUE)) %>%
  mutate(log_sum_icu = log(sum_icu + 0.01))

p_moran_icu <-  ggplot(dat_icu, aes(date, sum_icu)) +
  theme_minimal() +
  geom_line(lty = 1, color = "darkred", lwd = 1.1) +
  geom_point(color = "darkred") +
  
  # ## spatial distance
  # geom_rect(
  #   data = rect_sd, 
  #   inherit.aes=FALSE, 
  #   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
  #       ymax = max(dat_icu$sum_icu), group = group), 
  #   color = "transparent", fill = "orange", alpha = 0.3) +
  # 
  # ## adjacency
  # geom_rect(
#   data = rect_nb,
#   inherit.aes=FALSE,
#   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
#       ymax = max(dat_icu$sum_icu), group = group),
#   color = "transparent", fill = "orange", alpha = 0.3) +
# 
# ## bed density
# geom_rect(
#   data = rect_bd,
#   inherit.aes=FALSE,
#   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
#       ymax = max(dat_icu$sum_icu), group = group),
#   color = "transparent", fill = "orange", alpha = 0.3) +


## spatial distance
geom_rect(
  data = rect_sd, 
  inherit.aes=FALSE, 
  aes(xmin = start, xmax = end, ymin = -Inf,
      ymax = Inf, group = group), 
  color = "transparent", fill = "orange", alpha = 0.3) +
  
  ## adjacency
  geom_rect(
    data = rect_nb,
    inherit.aes=FALSE,
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group),
    color = "transparent", fill = "orange", alpha = 0.3) +
  
  ## bed density
  geom_rect(
    data = rect_bd,
    inherit.aes=FALSE,
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group),
    color = "transparent", fill = "orange", alpha = 0.3) +
  
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'ICU occupancy',
       x = 'Date', y = "ICU beds",
       title = 'Spatial-autocorrelation over time')

moran_plot <- cowplot::plot_grid(
  p_moran_icu, 
  p_moran_nb,
  p_moran_sd,
  p_moran_bd,
  ncol = 1,
  labels = c('A.1', 'B.1', 'C.1', 'D.1'))


# Plot phases of spatial autocorrelation (Incidence)----------------------------

## prepare Moran's I data
dat_sd <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long %>%
    mutate(date = date - 14),  
  input_mat = nb_sd_germany,
  win_size = 14,
  var_name = 'incidence'
)

dat_bd <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long %>%
    mutate(date = date - 14),  
  input_mat = nb_bd_germany,
  win_size = 14,
  var_name = 'incidence'
)

dat_nb <- prepare_moran_data(
  input_data = data_w_lagged_values_germany$long %>%
    mutate(date = date - 14),  
  input_mat = nb_b_germany,
  win_size = 14,
  var_name = 'incidence'
)

## make rectangles
rect_sd <-  extract_sig_phase(input_data = dat_sd)
rect_bd <-  extract_sig_phase(input_data = dat_bd)
rect_nb <-  extract_sig_phase(input_data = dat_nb)


## plots
p_moran_sd <- ggplot(data = dat_sd, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_sd, 
    inherit.aes = FALSE, 
    aes(xmin = start, xmax = end, ymin = min(dat_sd$observed),
        ymax = max(dat_sd$observed), group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Spatical distance',
       x = 'Date', y = "Moran's I")

p_moran_bd <- ggplot(data = dat_bd, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_bd, 
    inherit.aes=FALSE, 
    aes(xmin = start, xmax = end, ymin = min(dat_bd$observed),
        ymax = max(dat_bd$observed), group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Distance in structural quality',
       x = 'Date', y = "Moran's I")

p_moran_nb <- ggplot(data = dat_nb, aes(date, observed)) +
  theme_minimal() +
  geom_line(lty = 1, color = "steelblue", lwd = 1.1) +
  geom_point(color = "steelblue") +
  geom_rect(
    data = rect_nb, 
    inherit.aes=FALSE, 
    aes(xmin = start, xmax = end, ymin = min(dat_nb$observed),
        ymax = max(dat_nb$observed), group = group), 
    color = "transparent", fill = "orange", alpha = 0.3) +
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Type of matrix: Adjacency',
       x = 'Date', y = "Moran's I")

dat_incidence <- rki_data %>%
  mutate(date = as.Date(Meldedatum)) %>%
  filter(
    NeuerFall %in% c(0, 1),
    date >= min(data_w_lagged_values_germany$long$date),
    date <= max(data_w_lagged_values_germany$long$date)
  ) %>%
  group_by(date) %>%
  summarise(sum_incidence = sum(AnzahlFall, na.rm = TRUE))


# dat_incidence <- data_w_lagged_values_germany$long %>%
#   mutate(date = date - 14) %>%
#   group_by(date) %>%
#   summarise(sum_incidence = sum(incidence, na.rm = TRUE))

p_moran_incidence <-  ggplot(dat_incidence, aes(date, sum_incidence)) +
  theme_minimal() +
  geom_line(lty = 1, color = "darkred", lwd = 1.1) +
  geom_point(color = "darkred") +
  
  # ## spatial distance
  # geom_rect(
  #   data = rect_sd, 
  #   inherit.aes=FALSE, 
  #   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
  #       ymax = max(dat_icu$sum_icu), group = group), 
  #   color = "transparent", fill = "orange", alpha = 0.3) +
  # 
  # ## adjacency
  # geom_rect(
#   data = rect_nb,
#   inherit.aes=FALSE,
#   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
#       ymax = max(dat_icu$sum_icu), group = group),
#   color = "transparent", fill = "orange", alpha = 0.3) +
# 
# ## bed density
# geom_rect(
#   data = rect_bd,
#   inherit.aes=FALSE,
#   aes(xmin = start, xmax = end, ymin = min(dat_icu$sum_icu),
#       ymax = max(dat_icu$sum_icu), group = group),
#   color = "transparent", fill = "orange", alpha = 0.3) +


## spatial distance
geom_rect(
  data = rect_sd, 
  inherit.aes=FALSE, 
  aes(xmin = start, xmax = end, ymin = -Inf,
      ymax = Inf, group = group), 
  color = "transparent", fill = "orange", alpha = 0.3) +
  
  ## adjacency
  geom_rect(
    data = rect_nb,
    inherit.aes=FALSE,
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group),
    color = "transparent", fill = "orange", alpha = 0.3) +
  
  ## bed density
  geom_rect(
    data = rect_bd,
    inherit.aes=FALSE,
    aes(xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, group = group),
    color = "transparent", fill = "orange", alpha = 0.3) +
  
  scale_x_date(date_labels="%b %y", date_breaks  ="4 month") +
  labs(subtitle = 'Case counts',
       x = 'Date', y = "Incident cases",
       title = ' '
  )

moran_plot_incidence <- cowplot::plot_grid(
  p_moran_incidence, 
  p_moran_nb,
  p_moran_sd,
  p_moran_bd,
  ncol = 1,
  labels = c('A.2', 'B.2', 'C.2', 'D.2'))


moran_both <- cowplot::plot_grid(
  moran_plot, moran_plot_incidence, nrow = 1, ncol = 2)

#saveRDS(moran_both, "moran_both.rds")
