
# merge_raw_data ---------------------------------------------------------------

#'@details merges all relevenat raw data sets
#'@param input_divi data.frame, divi data containing icu cases
#'@param input_lkr data.frame, spatial districts including shape
#'@param input_ils data.frame, list of ILS units 
#'@param input_rki data.frame, official data on case counts from RKI
#'@param age_group (list of) character(s), age group from RKI data
#'@param federal_state character, selects state, format as in lkr_data: SN_L
#'@return list of data.frames, joined tabular data as well as spatial data

merge_raw_data <- function(
    input_divi = divi_data,
    input_lkr = lkr_data,
    input_ils = ils_data,
    input_rki = rki_data,
    input_kh = kh_data,
    age_group = c('A60-A79', 'A80+'),
    federal_state = '09') {
  
  ## clean up the input data
  input_lkr <- input_lkr %>%
    filter(EWZ > 0)
  
  input_ils <- input_ils %>%
    mutate(AGS_double = duplicated(AGS)) %>%
    filter(AGS_double == FALSE) %>%
    dplyr::select(-AGS_double)
  
  ## if nothing is selected, take all values
  if (is.null(age_group)) {
    
    age_group <- input_rki$Altersgruppe %>% unique()
    
  }
  
  ## change the format of AGS such that it is equal in all datasets
  divi_tmp <- input_divi %>%
    mutate(AGS_raw = as.character(gemeindeschluessel), 
           AGS = case_when(
             nchar(AGS_raw) == 4 ~ paste0("0", AGS_raw),
             nchar(AGS_raw) >= 5 ~ AGS_raw, 
           ),
           date = as.Date(date)) %>%
    dplyr::select(-AGS_raw)
  
  ils_tmp <- input_ils %>%
    mutate(AGS_raw = as.character(AGS), 
           AGS = case_when(
             nchar(AGS_raw) == 4 ~ paste0("0", AGS_raw),
             nchar(AGS_raw) >= 5 ~ AGS_raw, 
           )) %>%
    group_by(ILS) %>%
    mutate(ILS_ID = cur_group_id()) %>%
    ungroup()
  
  rki_tmp <- input_rki %>%
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
  
  ## hospital density
  kh_tmp <- input_kh %>%
    magrittr::set_colnames(c('AGS', 'lkr_name', 'kh_density')) %>%
    filter(nchar(AGS) == 5)
  
  kh_tmp$kh_density <- kh_tmp$kh_density %>%
    as.character() %>%
    as.numeric()
  
  ## mean impute missing kh density information
  mean_kh_density <-  kh_tmp %>% filter(kh_density < 1e+6)
  mean_kh_density <- mean(mean_kh_density$kh_density)
  kh_tmp$kh_density[kh_tmp$kh_density > 1e+6] <- mean_kh_density
  
  kh_tmp <- kh_tmp %>%
    full_join(.,
              input_lkr %>%
                as.data.frame() %>%
                filter(EWZ != 0) %>%
                dplyr::select(AGS, EWZ, - geometry),
              by = c('AGS')) %>%
    mutate(kh_beds = (kh_density * EWZ) / 1000)
  
  
  ## select federal state
  if (is.null(federal_state)) {
    
    lkr_tmp <- input_lkr
    
  } else {
    
    lkr_tmp <- input_lkr %>%
      filter(SN_L %in% all_of(federal_state))
    
  }
  
  ## join data sets
  icu_data <- lkr_tmp %>%
    dplyr::select(AGS, GEN, BEZ) %>%
    left_join(.,
              ils_tmp,
              by = c('AGS')) %>%
    left_join(.,
              divi_tmp,
              by = c('AGS')) %>%
    left_join(.,
              rki_tmp %>%
                group_by(Landkreis, date) %>%
                mutate(
                  case_count = sum(AnzahlFall),
                  Cat_A00_A04 = sum(Cat_A00_A04),
                  Cat_A05_A14 = sum(Cat_A05_A14),
                  Cat_A15_A34 = sum(Cat_A15_A34),
                  Cat_A35_A59 = sum(Cat_A35_A59),
                  Cat_A60_A79 = sum(Cat_A60_A79),
                  Cat_A80 = sum(Cat_A80),
                  Cat_unbekannt = sum(Cat_unbekannt)
                ) %>%
                ungroup() %>%
                dplyr::select(AGS, date, case_count, age_labels) %>%
                distinct(),
              by = c('AGS', 'date')) %>%
    as.data.frame() %>%
    dplyr::select(date, ILS, ILS_ID, case_count, age_labels,
                  faelle_covid_aktuell,
                  -geometry) %>%
    group_by(date, ILS, ILS_ID) %>%
    summarise(
      incidence = sum(case_count, na.rm = TRUE), 
      icu = sum(faelle_covid_aktuell, na.rm = TRUE),
      Cat_A00_A04 = sum(Cat_A00_A04, na.rm = TRUE),
      Cat_A05_A14 = sum(Cat_A05_A14, na.rm = TRUE),
      Cat_A15_A34 = sum(Cat_A15_A34, na.rm = TRUE),
      Cat_A35_A59 = sum(Cat_A35_A59, na.rm = TRUE),
      Cat_A60_A79 = sum(Cat_A60_A79, na.rm = TRUE),
      Cat_A80 = sum(Cat_A80, na.rm = TRUE),
      Cat_unbekannt = sum(Cat_unbekannt, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(ILS, date)
  
  
  ## filter states
  if (!is.null(federal_state)) {
    
    ils_sp <- input_lkr %>%
      filter(SN_L %in% all_of(federal_state))
    
  } else {
    
    ils_sp <- input_lkr 
  }
  
  ## spatial data
  ils_sp <- ils_sp %>%
    dplyr::select(AGS, GEN, BEZ) %>%
    left_join(.,
              ils_tmp,
              by = c('AGS')) %>%
    left_join(.,
              kh_tmp,
              by = c('AGS')) %>%
    group_by(ILS) %>%
    summarize(
      geometry = st_union(geometry),
      beds_density = (
        sum(kh_beds, na.rm = TRUE) / sum(EWZ, na.rm = TRUE)) * 1000
    )
  
  ## remove NAs
  icu_data <- icu_data %>% filter(!is.na(ILS))
  ils_sp <- ils_sp %>% filter(!is.na(ILS))
  
  ## return joined data sets
  out_data <- list()
  out_data$icu_data <- icu_data
  out_data$ils_sp <- ils_sp
  
  return(out_data)
  
}

# get shares -------------------------------------------------------------------

#'@details turns the incidence per age group into shares of the total
#'@param input_merged list of data, merged_data
#'@returns the input with shares instead of integers

get_shares <- function(input_merged = merged_data) {
  
  share_data <- input_merged$icu_data %>% 
    dplyr::select(starts_with('Cat_')) / 
    input_merged$icu_data$incidence
  
  share_data[is.na(share_data)] <- 0
  share_data <- round(share_data, digits = 2)
  
  input_merged$icu_data[, colnames(share_data)] <- share_data
  
  return(input_merged)
  
}

# make covariates --------------------------------------------------------------

#'@details prepares n matrices with covariates pertaining to agegroups
#'@param input_merged list, merged_data
#'@param training_dat, data.frame, the final training data
#'@param k_input, positive integer, time lag of covariates
#'@return a list of n matrices to be imported into Stan

make_covariates <- function(
    input_merged = merged_data,
    training_dat = training_data,
    k_input = 14) {
  
  cov_list <- list()
  
  age_labels <- merged_data$icu_data %>% 
    dplyr::select(starts_with('Cat_')) %>%
    colnames()
  
  for (i in seq(length(age_labels))) {
    
    lag_tmp <- NULL
    
    ## lag covariates
    lag_tmp <- lag_covariates(
      input_data = input_merged$icu_data %>%
        dplyr::select(date, ILS, ILS_ID, icu, age_labels[i]),
      k = k_input,
      cov_name = age_labels[i]
    )
    
    ## wide
    lag_tmp <- lag_tmp$wide %>% 
      select(-starts_with('icu')) %>%
      filter(date >= training_dat$from_date) %>%
      head(nrow(training_dat$wide)) %>%
      dplyr::select(-date)
    
    cov_list[[i]] <- lag_tmp
    
  }
  
  return(cov_list)
  
}


# make_stacked_share_plots -----------------------------------------------------

#'@details makes stacked plots showing the covariates as over time
#'@param input_rki, data.frame, the prepared rki data
#'@param input_phases, data.frame, the prepared phases data
#'@param input_bay, boolean, if TRUE: select bavaria
#'@param no_legend, boolean, if TRUE, no legend is added
#'@return returns ggplots of stacked covariates over time

make_stacked_share_plots <- function(
    input_rki = rki_tmp,
    input_phases = phases_tmp,
    input_bay = FALSE,
    no_legend = FALSE) {
  
  ## filter bayern
  if (input_bay == TRUE) {
    
    rki_tmp <- input_rki %>%
      filter(Bundesland == 'Bayern')
    
  } else {
    
    rki_tmp <- input_rki
    
  }
  
  ## age group
  levels(rki_tmp$AgeGroups)[levels(rki_tmp$AgeGroups) == "unbekannt"] <- 
    "unknown"
  
  rki_tmp$phase <- NA
  
  ## add info pertaining to phase
  for (j in seq(nrow(input_phases))) {
    
    current_phase <- 
      as.Date(input_phases[j, 'start_date']):
      as.Date(input_phases[j, 'end_date'])
    
    phase_indices <- which(rki_tmp$date %in% as.Date(current_phase))
    
    rki_tmp$phase[phase_indices] <- input_phases$phase[j]
    
  }
  
  rki_dat_tmp <- rki_tmp %>%
    filter(!is.na(phase)) %>%
    group_by(phase, date, AgeGroups) %>%
    summarise(case_counts = sum(AnzahlFall)) %>%
    group_by(date) %>%
    mutate(
      total_cases = sum(case_counts),
      group_share = case_counts / total_cases)
  
  age_groups <- rki_dat_tmp$AgeGroups %>% unique() %>% as.character()
  levels(age_groups)[levels(age_groups) == "unbekannt"] <- "unknown"
  
  rki_dates <- rki_dat_tmp$date %>% unique()
  
  k <- length(age_groups)
  j <- length(rki_dates)
  
  rki_dates <- cbind(rep(rki_dates, k), as.character(rep(age_groups, j))) %>%
    as.data.frame()
  
  rki_dates <- rki_dates %>%
    magrittr::set_colnames(c('date', 'AgeGroups')) %>%
    mutate(date = as.Date(as.numeric(date)))
  
  rki_dat_tmp <- rki_dat_tmp %>%
    full_join(
      .,
      rki_dates,
      by = c('date', 'AgeGroups')
    ) %>% 
    arrange(date) %>%
    mutate(
      phase = zoo::na.locf(phase),
      ymw = format(date, format(date, '%Y-%m-%V')),
    ) %>%
    replace_na(list(case_counts = 0, total_cases = 0, group_share = 0))
  
  
  k <- 1
  
  plot_title <- 'Share per age group over time'
  
  if (input_bay == TRUE) {
    
    plot_title <- paste0(plot_title, ' (Bavaria)')
    
  } else {
    
    plot_title <- paste0(plot_title, ' (Germany)')
    
  }
  
  for (j in input_phases$phase) {
    
    dat_tmp <- rki_dat_tmp %>%
      filter(phase == j)
    
    p_tmp <- dat_tmp %>%
      ggplot(., aes(x = date, y = group_share, fill = AgeGroups)) + 
      geom_area(alpha = 0.6, size = 1, colour = "black") +
      theme_bw() +
      guides(fill = guide_legend(title = "Age groups"))
    
    if (j == 1) {
      
      p_tmp <- p_tmp +  
        labs(
          title = plot_title,
          subtitle =  
            paste0('Phase: ', j, ' from: ', 
                   input_phases$start_date[k],
                   ' to: ', 
                   input_phases$end_date[k]),
          x = 'Date', y = "share")
      
    } else {
      
      p_tmp <- p_tmp +  
        labs(
          title = ' ',
          subtitle =  paste0('Phase: ', j, ' from: ', 
                             input_phases$start_date[k],
                             ' to: ', 
                             input_phases$end_date[k]),
          x = 'Date', y = "share")
      
    }
    
    if (no_legend == FALSE) {
      
      if (j < input_phases$phase[length(input_phases$phase)]) {
        
        p_tmp <- p_tmp + theme(legend.position = "none")
        
      }
      
    } else {
      
      p_tmp <- p_tmp + theme(legend.position = "none")
      
    }
    
    
    assign(paste0('pl_', k), p_tmp)
    
    k <- k + 1
    
  }
  
  cp_grid <- cowplot::plot_grid(plotlist = mget(paste0("pl_", 1:(k - 1))))
  
  return(cp_grid)
  
}

# make_summary_tables_rki ------------------------------------------------------

#'@details computes summary statics for the rki data
#'@param input_rki, data.frame, prepared rki data
#'@param input_bay, boolean, TRUE = filter for bavaria
#'@param input_phases, data.frame, prepared data for the phases
#'@return returns a html table with summary statistics by phase

make_summary_tables_rki <- function(
    input_rki = rki_tmp,
    input_bay = FALSE,
    input_phases = phases_tmp) {
  
  
  ## filter bayern
  if (input_bay == TRUE) {
    
    rki_tmp <- input_rki %>%
      filter(Bundesland == 'Bayern')
    
  } else {
    
    rki_tmp <- input_rki
    
  }
  
  rki_tmp$phase <- NA
  
  ## add info pertaining to phase
  for (j in seq(nrow(input_phases))) {
    
    current_phase <- 
      as.Date(input_phases[j, 'start_date']):
      as.Date(input_phases[j, 'end_date'])
    
    phase_indices <- which(rki_tmp$date %in% as.Date(current_phase))
    
    rki_tmp$phase[phase_indices] <- input_phases$phase[j]
    
  }
  
  ## summarise statistics
  rki_table_tmp <- rki_tmp %>%
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
  
  ## make table
  tab_out <- rki_table_tmp %>%
    kbl() %>%
    kable_styling()
  
  return(tab_out)
  
}

# make_icu_summary_stats -------------------------------------------------------

#'@details makes summary stats from icu and incidence data
#'@param input_data, data.frame, data_w_lagged_values$long (bay)
#'@param input_phases, data.frame, prepared phases data phases_tmp
#'@return returns a html table with summary statistics per phase

make_icu_summary_stats <- function(
    input_data = data_w_lagged_values_bavaria$long,
    input_phases = phases_tmp) {
  
  input_data$phase <- NA
  
  ## add info pertaining to phase
  for (j in seq(nrow(input_phases))) {
    
    current_phase <- 
      as.Date(input_phases[j, 'start_date']):
      as.Date(input_phases[j, 'end_date'])
    
    phase_indices <- which(input_data$date %in% as.Date(current_phase))
    
    input_data$phase[phase_indices] <- input_phases$phase[j]
    
  }
  
  icu_table <- input_data %>% 
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
  
  
  out_table <- rbind(
    icu_table %>%
      dplyr::select(
        -mean_incidence,
        -median_incidence,
        -sd_incidence,
        -min_incidence,
        -max_incidence) %>%
      magrittr::set_colnames(c('Phase', 'Mean', 'Median', 'SD', 'Min', 'Max')),
    icu_table %>%
      dplyr::select(
        -mean_icu,
        -median_icu,
        -sd_icu,
        -min_icu,
        -max_icu) %>%
      magrittr::set_colnames(c('Phase', 'Mean', 'Median', 'SD', 'Min', 'Max'))
  ) %>%
    kbl() %>%
    kable_styling()
  
  return(out_table)
  
}


# lag_covariates ---------------------------------------------------------------

#'@details adds lagged values of a covariate for prediction later on
#'@param input_data data.frame, prepared icu_data
#'@param k positive integer, number of lagged days of covariate
#'@param z positive integer, number of days for rolling sum
#'@param cov_name character, name of variable to be shifted
#'@return list of data.frames, data with lagged covariates in long and wide

lag_covariates <- function(input_data = icu_data, 
                           k = 14,
                           z = NULL,
                           cov_name = 'incidence') {
  
  
  ## list of all the dates
  date_list <- min(input_data$date, na.rm = TRUE):
    max(input_data$date, na.rm = TRUE) %>%
    as.Date(., origin = '1970-01-01') %>%
    as.data.frame() %>%
    magrittr::set_colnames('date')
  
  ## full range of dates
  icu_tmp <- input_data %>%
    dplyr::select(-ILS, -icu) %>%
    pivot_wider(names_from = ILS_ID, values_from = all_of(cov_name)) %>%
    full_join(.,
              date_list,
              by = c('date')) %>%
    filter(!is.na(date)) 
  
  
  ## sum cases to represent the number of infected (e.g. 14 days sum)
  if (!is.null(z)) {
    
    icu_tmp <- icu_tmp %>%
      read.zoo(.) %>%
      rollapplyr(., z, sum, partial = TRUE) %>%
      as.data.frame() %>%
      mutate(date = as.Date(rownames(.)) + z)
    
  }
  
  
  ## merge icu cases and pivot to longer (for ggplot)
  icu_long <- icu_tmp %>%
    pivot_longer(!date, names_to = "ILS_ID", values_to = all_of(cov_name)) %>%
    mutate(ILS_ID = as.integer(ILS_ID)) %>%
    full_join(.,
              input_data %>%
                dplyr::select(-all_of(cov_name), -ILS),
              by = c('date', 'ILS_ID')) %>%
    full_join(., 
              input_data %>%
                dplyr::select(ILS_ID, ILS) %>%
                distinct(ILS_ID, ILS),
              by = c('ILS_ID')) %>%
    filter(date >= min(icu_tmp$date)) %>%
    arrange(ILS_ID, date)
  
  ## lag covariates
  icu_long <- icu_long %>%
    dplyr::select(date, icu, ILS, ILS_ID) %>%
    left_join(
      .,
      icu_long %>%
        dplyr::select(date, all_of(cov_name), ILS, ILS_ID) %>%
        mutate(date = date + k),
      by = c('date', 'ILS', 'ILS_ID')
    ) %>%
    drop_na()
  
  ## pivot to wider for passing (y, X) to stan
  icu_wide <- icu_long %>%
    dplyr::select(-ILS) %>%
    pivot_wider(names_from = ILS_ID, values_from = c(icu, all_of(cov_name)))
  
  ## export results
  out_data <- list()
  out_data$long <- icu_long
  out_data$wide <- icu_wide
  out_data$horizon <- k
  
  return(out_data)
  
}


# make_adjacency_matrix --------------------------------------------------------

#'@details sets up an adjacency matrix for the CAR model
#'@param input_sp data.frame, ILS units with geometry
#'@param weights_style character, can be either B, C, W, S (see nb2mat), 
#'SD (spatial distance), or BD (delta in beds density)
#'@return adjacency matrix with spatial weights

make_adjacency_matrix <- function(input_sp, weights_style = 'B') {
  
  ## adjacency matrix
  nb_mat <- input_sp %>%
    st_as_sf() %>%
    spdep::poly2nb(c('ILS'))
  
  w_mat <-  nb_mat %>%
    spdep::nb2mat(zero.policy = TRUE, style = 'B')
  
  if (weights_style %in% c('B', 'W', 'S', 'C')) {
    
    w_mat <- nb_mat %>%
      spdep::nb2mat(zero.policy = TRUE, style = all_of(weights_style))
    
  } else if (weights_style == 'SD') {
    
    centroid_mat <- input_sp %>% 
      st_centroid() %>% 
      st_geometry() %>%
      sf::st_distance()
    
    w_mat <- round(centroid_mat * 0.001)
    
    # w_mat <- centroid_mat * w_mat
    # w_mat <- round(w_mat / max(w_mat), digits = 2)
    
  } else if (weights_style == 'BD') {
    
    bed_density <- input_sp$beds_density
    bd_mat <- w_mat
    
    for (i in seq(nrow(input_sp))) {
      
      for (j in seq(nrow(input_sp))) {
        
        bd_mat[i, j] <-  bed_density[j] - bed_density[i]  
        
      }
      
    } 
    
    w_mat <- abs(bd_mat)
    
  }
  
  ## add column names
  w_mat <- w_mat %>%
    magrittr::set_colnames(seq(ncol(.)))
  
  return(w_mat)
  
}


# scale_matrix ----------------------------------------------------------------

#'@details scales a matrix (distances)
#'@param input_mt, matrix, result of make_adjacency_matrix()
#'@return returns a scaled version of the input

scale_matrix <- function(input_mat = NULL) {
  
  mat_tmp <- input_mat / find_nb_dist(input_mat, min_max = 'min')
  mat_tmp <- mat_tmp %>% as.data.frame() %>% as.matrix
  
  return(mat_tmp)
  
}

# find_min_dist ----------------------------------------------------------------

#'@details finds the min and max distance between cenroids
#'@param mat_input, matrix, nb_matrix from make_adjacency_matrix
#'@param min_max, either 'min' or 'max'
#'@return returns numeric minimum distance

find_nb_dist <- function(mat_input = nb_mat, min_max = 'max') {
  
  dist <- mat_input %>% 
    as.list() %>% 
    unlist() %>%
    as.data.frame() %>%
    magrittr::set_colnames('x') %>%
    filter(x > 0)
  
  if (min_max == 'min') {
    
    dist <- dist %>% min()
    
  } else if (min_max == 'max') {
    
    dist <- dist %>% max()
    
  }
  
  return(dist)
  
}

# prepare_eval_data ------------------------------------------------------------

#'@description loads fitted models and stacks results for evaluation
#'@param dates_list list, list of dates
#'@param n_models integer, number of models in scheme
#'@param file_prefix character, first part of file name
#'@param file_suffix character, last part of file name
#'@param par_eval character, name of parameter to be extracted
#'@param prepared_data data.frame, the edited data
#'@return returns a data.frame with stacked results of all fitted models

prepare_eval_data <- function(dates_list = list_of_dates[1:9], 
                              n_models = 8,
                              prepared_data = data_w_lagged_values,
                              file_prefix = './stan_fits/',
                              file_suffix = '_bayern.rds',
                              par_eval = 'eta_obs',
                              transf_params = FALSE) {
  
  results_combined <- NULL
  
  if (length(n_models) == 1) {
    
    n_models <- seq(n_models)
    
  }
  
  for (j in seq(length(dates_list))) {
    
    for (k in n_models) {
      
      ## load fitted models
      mod_tmp <- readRDS(
        paste0(file_prefix, 'fit_0', k, '_iter_', j, file_suffix)
      )
      
      ## get training data
      training_data <- extract_time_slice(
        input_data = prepared_data, 
        max_t = 30,
        from_date = dates_list[j]
      )
      
      ## extract results
      eval_tmp <- combine_fit_w_observed(
        fitted_model = mod_tmp, 
        parameter = par_eval, 
        observed_data = training_data$long,
        transform_param = transf_params,
        ts_type = 'all'
      )
      
      ## extract data frame from results
      eval_tmp <- eval_tmp[[1]]$long
      
      ## add ids
      eval_tmp$model_nr <- k
      eval_tmp$iter <- j 
      
      ## rectify normal results
      eval_tmp <- eval_tmp %>%
        mutate(
          lower = case_when(
            lower < 0 ~ 0,
            lower >= 0 ~ lower),
          upper = case_when(
            upper < 0 ~ 0,
            upper >= 0 ~ upper),
          case_count = case_when(
            case_count < 0 ~ 0,
            case_count >= 0 ~ case_count)
        )
      
      ## combine results
      results_combined <- rbind(results_combined, eval_tmp)
      
    }
    
  }
  
  print(paste0('k is: ', k, ' and j is: ', j))
  
  return(results_combined)
  
}


# prepare_eval_data_pooled -----------------------------------------------------

#'@description loads fitted models and stacks results for evaluation
#'@param dates_list list, list of dates
#'@param n_models integer, number of models in scheme
#'@param file_prefix character, first part of file name
#'@param file_suffix character, last part of file name
#'@param par_eval character, name of parameter to be extracted
#'@param prepared_data data.frame, the edited data
#'@return returns a data.frame with stacked results of all fitted models

prepare_eval_data <- function(dates_list = list_of_dates[1:9], 
                              n_models = 8,
                              prepared_data = data_w_lagged_values,
                              file_prefix = './stan_fits/',
                              file_suffix = '_bayern.rds',
                              par_eval = 'eta_obs',
                              transf_params = FALSE) {
  
  results_combined <- NULL
  
  if (length(n_models) == 1) {
    
    n_models <- seq(n_models)
    
  }
  
  for (j in seq(length(dates_list))) {
    
    for (k in n_models) {
      
      ## load fitted models
      mod_tmp <- readRDS(
        paste0(file_prefix, 'fit_0', k, '_iter_', j, file_suffix)
      )
      
      ## get training data
      training_data <- extract_time_slice(
        input_data = prepared_data, 
        max_t = 30,
        from_date = dates_list[j]
      )
      
      ## extract results
      eval_tmp <- combine_fit_w_observed(
        fitted_model = mod_tmp, 
        parameter = par_eval, 
        observed_data = training_data$long,
        transform_param = transf_params,
        ts_type = 'all'
      )
      
      ## extract data frame from results
      eval_tmp <- eval_tmp[[1]]$long
      
      ## add ids
      eval_tmp$model_nr <- k
      eval_tmp$iter <- j 
      
      ## rectify normal results
      eval_tmp <- eval_tmp %>%
        mutate(
          lower = case_when(
            lower < 0 ~ 0,
            lower >= 0 ~ lower),
          upper = case_when(
            upper < 0 ~ 0,
            upper >= 0 ~ upper),
          case_count = case_when(
            case_count < 0 ~ 0,
            case_count >= 0 ~ case_count)
        )
      
      ## combine results
      results_combined <- rbind(results_combined, eval_tmp)
      
    }
    
  }
  
  print(paste0('k is: ', k, ' and j is: ', j))
  
  return(results_combined)
  
}

# aggregate_samples ------------------------------------------------------------

#'@details extracts stan samples and aggregates them for Germany
#'@param input_mod, stan fit object
#'@param length_ts, integer, length of test plus training
#'@param input_merged_data, list, merged_data
#'@param input_ils, data.frame, input_ils data
#'@param input_county, character, prefix of county (09 == bavaria)
#'@param input_date, date, for time indexing the results
#'@return returns aggregated pred lower and upper

aggregate_samples <- function(
    input_mod = mod_tmp,
    length_ts = 44,
    input_merged_data = merged_data_germany,
    input_ils = ils_data,
    input_county = NULL,
    input_date = list_of_dates[1]) {
  
  
  ## enable aggregation by county
  ILS_ID_list <-  input_merged_data$icu_data[, c('ILS', 'ILS_ID')] %>% 
    distinct() %>%
    left_join(.,
              input_ils %>%
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
  
  if (!is.null(input_county)) {
    
    ILS_ID_list <- ILS_ID_list %>%
      filter(County == input_county)
    
  }
  
  ILS_ID_list <- ILS_ID_list$ILS_ID %>% unique()
  
  ## combine samples
  samples_combined <- matrix(0, nrow = 1000, ncol = length_ts)
  results_out <- matrix(0, nrow = length_ts, ncol = 3) %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('lower', 'upper', 'pred'))
  
  
  for (k in seq(length_ts)) {
    
    for (i in ILS_ID_list) {
      
      samples_tmp <- mod_tmp@sim$samples[[1]][paste0('eta_obs','.', k, '.', i)]
      samples_combined[, k] <- samples_combined[, k] + unlist(samples_tmp)
      
    }
    
    quantiles_tmp <- quantile(samples_combined[, k], c(0.025, 0.5, 0.975))
    results_out[k, 'lower'] <- quantiles_tmp[1]
    results_out[k, 'pred'] <- quantiles_tmp[2]
    results_out[k, 'upper'] <- quantiles_tmp[3]
    
  }
  
  ## add date and time index
  results_out <- results_out %>%
    mutate(
      time_index = seq(length_ts),
      date = input_date:(input_date + (length_ts - 1)) %>% 
        as.Date()
    )
  
  
  return(results_out)
  
}


# train_benchmark_models -------------------------------------------------------

#'@details trains a model for comparison
#'@param training_input, list of data.frames, the prepared training data
#'@param nb_bd, matrix, distance matrix
#'@param dates_list, list, list of dates for the training
#'@param suffix, character, add a suffix to the model name
#'@param model_type, character, name of the model
#'@param nr_ar, integer, either 1 or 2 for (AR(1) or AR(2) process)
#'@param make_stationary, boolean, de-trend data befor training
#'@return returns a data.frame with pred, lower and upper, similar to loss_data

train_benchmark_models <- function(
    training_input = data_w_lagged_values,
    mat_type = nb_b,
    dates_list = list_of_dates,
    suffix = ' 1',
    nr_ar = 2,
    model_type = 'CAR AR(2)',
    make_stationary = TRUE
) {
  
  
  for (j in 1:length(list_of_dates)) {
    
    ## get training data
    training_data <- extract_time_slice(
      input_data = data_w_lagged_values, 
      max_t = 30,
      from_date = list_of_dates[j]
    )
    
    ## add covariates
    training_tmp <- training_data$long %>%
      left_join(
        .,
        merged_data$icu_data %>%
          mutate(date = date + 14) %>%
          dplyr::select(date, ILS_ID, ILS, contains('Cat_')),
        by = c('date', 'ILS_ID', 'ILS')
      ) %>%
      mutate(
        X1 = round((Cat_A00_A04 / (incidence + 0.01)), digits = 3),
        X2 = round((Cat_A05_A14 / (incidence + 0.01)), digits = 3),
        X3 = round((Cat_A15_A34 / (incidence + 0.01)), digits = 3),
        X4 = round((Cat_A35_A59 / (incidence + 0.01)), digits = 3),
        X5 = round((Cat_A60_A79 / (incidence + 0.01)), digits = 3),
        X7 = round((Cat_A80 / (incidence + 0.01)), digits = 3),
        X8 = round((Cat_unbekannt / (incidence + 0.01)), digits = 3)
      )
    
    if (make_stationary == TRUE) {
      
      ## first differences
      training_tmp <- training_tmp %>% 
        group_by(ILS_ID) %>% 
        filter(n() > 1) %>%
        mutate(
          icu = icu - lag(icu),
          incidence = incidence - lag(incidence)
        ) %>%
        replace_na(list(icu = 0, incidence = 0))
      
    } 
    
    ## prepare training data
    train <- training_tmp
    train$obs <- training_data$long$icu
    
    ## set to NA
    train$icu[which(train$hold_out == 1)] <- NA
    
    ## prepare hold_out data
    hold_out <- train %>% filter(hold_out == 1)
    
    Y <- train$icu
    X <- train$incidence
    X1 <- train$X1
    X2 <- train$X2
    X3 <- train$X3
    X4 <- train$X4
    X5 <- train$X5
    
    ## Fit model
    model <- ST.CARar(
      formula = Y ~ X + X1 + X2 + X3 + X4 + X5, 
      family = "gaussian", 
      W = nb_bd,
      AR = nr_ar,
      burnin = 2000,
      n.sample = 6000) 
    
    # Median column for the average posterior prediction
    preds <- apply(
      model$samples$Y, 
      2 ,  
      quantile, 
      probs =  c(0.05, 0.50, 0.95), 
      na.rm = TRUE
    )
    
    hold_out$pred <- preds[2, ]
    hold_out$lower <- preds[1, ]
    hold_out$upper <- preds[3, ]
    
    ## add predictions to training data
    train$pred <- model$fitted.values
    train$lower <- NA
    train$upper <- NA
    
    joined_data <- rbind(train %>% filter(hold_out == 0), hold_out)
    
    if (make_stationary == TRUE) {
      
      ## go back to original scale
      joined_data <- joined_data %>% 
        full_join(
          .,
          training_data$long %>%
            group_by(ILS_ID) %>%
            filter(date == min(date)) %>%
            mutate(const_integration = icu) %>%
            dplyr::select(ILS_ID, const_integration),
          by = c('ILS_ID')
        ) 
      
      ## adjust intervals
      hold_out_tmp <- joined_data %>%
        filter(hold_out == 1) %>%
        group_by(ILS_ID) %>%
        mutate(
          pred = cumsum(pred) + const_integration,
          lower = cumsum(lower) + const_integration,
          upper = cumsum(upper) + const_integration
        )
      
      joined_data <- rbind(
        joined_data %>%
          filter(hold_out == 0),
        hold_out_tmp
      ) %>%
        #rectify results
        mutate(
          pred = case_when(
            pred < 0 ~ 0,
            pred >= 0 ~ pred
          ),
          lower = case_when(
            lower < 0 ~ 0,
            lower >= 0 ~ lower
          ),
          upper = case_when(
            upper < 0 ~ 0,
            upper >= 0 ~ upper
          )
        )
      
    }
    
    joined_data <- joined_data %>%
      rename(id = ILS_ID) %>%
      mutate(iter = j,
             model_nr = paste0(model_type, suffix)) %>%
      dplyr::select(
        c("obs", "id", "time_index", "date", "ILS", "model_nr", 
          "iter", "pred", "lower", "upper")
      ) 
    
    ## combine results
    results_combined <- rbind(results_combined, joined_data)
    
  }
  
  return(results_combined)
  
}

# get_ssed_dates ---------------------------------------------------------------

#'@details extracts one date per pandemic phase (phases_data) as a seed
#'@param input_phases, data.frame, phases_data
#'@param input_lagged_data, data.frame, data_w_lagged_values$long
#'@param real between 0 and 1, proportion for the proportionate sampling
#'@return returns a vector of dates

get_seed_dates <- function(input_phases = phases_data,
                           input_lagged_data = data_w_lagged_values$long,
                           share = 0.1,
                           fixed_number = 5) {
  
  date_list <- NULL
  
  min_date <- input_lagged_data$date %>% min()
  
  phases_tmp <- phases_data %>%
    filter(as.Date(start_date) > as.Date(min_date),
           type == 'wave')
  
  
  if (nrow(phases_tmp) > 2) {
    
    for(i in 1:(nrow(phases_tmp)-1)) {
      
      dat_tmp <- phases_tmp[i, ]
      interval_tmp <- as.Date(dat_tmp$start_date):as.Date(dat_tmp$end_date)
      
      
      if (is.null(fixed_number)) {
        
        ## get interval length
        l <- length(interval_tmp)
        l <-  round(share * l)
        
        if (l == 0) {l <- 1}
        
      } else {
        
        l <- fixed_number
        
      }
      
      seed_date <- sample(interval_tmp, l)
      
      date_list <- c(date_list, as.Date(seed_date))
      
    }
    
  }
  
  return(as.Date(date_list))
  
}


# extract_time_slice -----------------------------------------------------------

#'@details extracts a time series of size max_t from the complete series
#'@param input_data, list of data.frames, results from lag_covariates
#'@param max_t, positive integer, length of training data
#'@param from_date, date, min date of time series to be extracted
#'@param predict_unknown, boolean, if TRUE, predict into the unknown future
#'@return returns a list of data.fames, wide and long version of sliced input
#'data with time_index, and a dummmy variable signifying training and hold out

extract_time_slice <- function(
    input_data = capa_data, 
    max_t = 30, 
    from_date = NULL,
    predict_unknown = FALSE) {
  
  ## length of horizon
  h <- input_data$horizon
  
  ## window size
  k <- max_t + h
  
  ## max date
  max_date <- max(input_data$long$date, na.rm = TRUE) - k
  
  if (predict_unknown == FALSE) {
    
    ## min date
    min_date <- min(input_data$long$date, na.rm = TRUE)
    
    if (is.null(from_date)) {
      
      from_date <- sample(min_date:max_date, 1) %>% as.Date()
      
    } else {
      
      from_date <- as.Date(from_date)
      assert(class(from_date) == 'Date')
      assert(from_date >= min_date)
      assert(from_date <= max_date)
      
    }
    
  } else {
    
    from_date <- max_date + 1
    
  }
  
  ## list of dates
  dates_data <- from_date:(from_date + k - 1) %>%
    as.Date() %>%
    as.data.frame() %>%
    mutate(time_index = row_number(.)) %>%
    magrittr::set_colnames(c('date', 'time_index')) %>%
    mutate(hold_out = case_when(
      time_index <= max_t ~ 0,
      time_index > max_t ~ 1
    ))
  
  ## long
  long_tmp <-input_data$long %>%
    filter(date %in% from_date:(from_date + k - 1)) %>%
    full_join(.,
              dates_data,
              by = c('date'))
  
  ## wide
  wide_tmp <-input_data$wide %>%
    filter(date %in% from_date:(from_date + k- 1)) %>%
    full_join(.,
              dates_data,
              by = c('date'))
  
  ## out
  out_data <- list()
  out_data$long <- long_tmp
  out_data$wide <- wide_tmp
  out_data$from_date <- from_date
  out_data$horizon <- h
  out_data$n_training <- max_t
  
  return(out_data)
  
}

# prepare_moran_data -----------------------------------------------------------

#'@details prepares the a time series of Moran's I values
#'@param input_data, data.frame, long data
#'@param input_mat, matrix, spatial adjacency matrix
#'@param var_name, name of variable
#'@param win_size, integer, size of rolling window
#'@returns data.frame, Moran's I values with significance tests for plotting

prepare_moran_data <- function(
    input_data = data_w_lagged_values_germany$long,
    input_mat = nb_mat_b_germany,
    var_name = 'icu',
    win_size = 14){
  
  dates_tmp <- input_data$date %>% 
    unique() %>%
    as.data.frame() %>%
    mutate(date_id = row_number()) %>%
    magrittr::set_colnames(c('date', 'date_id'))
  
  ## roll data
  rolled_data <- input_data %>%
    group_by(ILS_ID) %>%
    mutate(roll_var = rollapply(
      get(var_name), win_size, median, partial = TRUE))
  
  ## find partially missing cases
  partially_na_cases <- rolled_data %>%
    group_by(ILS_ID) %>%
    filter(any(!(dates_tmp$date %in% date))) %>%
    dplyr::select(ILS_ID, ILS) %>%
    distinct()
  
  ## preclude from analyses
  rolled_tmp <- rolled_data %>%
    filter(!(ILS_ID %in% partially_na_cases$ILS_ID)) %>%
    group_by(ILS_ID) %>%
    mutate(row_id = row_number())
  
  k <- which(!(colnames(input_mat) %in% 
                 as.character(partially_na_cases$ILS_ID)))
  
  mat_tmp <- input_mat %>%
    as.data.frame() %>%
    mutate(rowID = row_number()) %>%
    filter(!(rowID %in% partially_na_cases$ILS_ID))
  
  mat_tmp <- mat_tmp[, k] %>% as.matrix()
  
  moran_i_out <- NULL
  
  for (i in seq(rolled_tmp$row_id %>% max())) {
    
    icu_tmp <- rolled_tmp %>%
      filter(row_id == i)
    
    moran_i_tmp <- ape::Moran.I(icu_tmp$roll_var, mat_tmp)
    moran_i_tmp$date <- as.character(icu_tmp$date[1])
    
    moran_i_out <- rbind(moran_i_out, as.data.frame(moran_i_tmp))
    
  }
  
  moran_i_out$sig <- 0
  moran_i_out$sig[which(moran_i_out$p.value < 0.05)] <- 1
  
  moran_i_out$date <- as.Date(moran_i_out$date)
  
  return(moran_i_out)
  
}

# extract_sig_phase ------------------------------------------------------------

#'@details extracts phases of significant autocorrelation for plotting
#'@param input_data data.frame, results of prepare_moran_data
#'@return returns intervals for plotting the phases of spatial autocorrelation

extract_sig_phase <- function(input_data = dat) {
  
  dat <- input_data
  
  ## Get the start and end points for highlighted regions
  inds <- diff(c(0, dat$sig))
  start <- dat$date[inds == 1]
  end <- dat$date[inds == -1]
  
  if (length(start) > length(end)) end <- c(end, tail(dat$date, 1))
  
  ## highlight region data
  rects <- data.frame(start = start, end = end, group = seq_along(start)) 
  
  return(rects)
  
}
# prepare_stan_data ------------------------------------------------------------

#'@details Prepares a list of data to be passed to Stan
#'@param input_data, data.frame, results of extract_time_slice
#'@param y_name, character, name of y
#'@param x_name, character, name of predictors
#'@param h, positive integer, length of horizon
#'@param W, matrix, neighborhood matrix
#'@param log_normal, boolean, if TRUE, replace 0 by 0.1
#'@return returns a list with the data ready to be passed to Stan

prepare_stan_data <- function(
    input_data = NULL,
    y_name = 'icu',
    x_name = 'incidence',
    W_input = nb_mat,
    distribution = c('normal', 'lognormal', 'poisson')) {
  
  ## assertions
  assertChoice(distribution, c('normal', 'lognormal', 'poisson'))
  
  ## prepare input
  data_tmp <- input_data$wide %>%
    dplyr::select(-hold_out, -time_index)
  
  data_tmp <- data_tmp[, str_sort(colnames(data_tmp), numeric = TRUE)]
  
  data_tmp <- data_tmp %>%
    arrange(date) %>%
    dplyr::select(-date)
  
  ## forecast horizon
  h <- input_data$horizon
  
  ## set up matrices
  y_mat <- data_tmp %>%
    dplyr::select(., contains(y_name)) %>%
    as.matrix()
  
  ## replace NAs with 0 in y_mat (does not affect estimation, since we only
  ## train the model on up to Z observations)
  y_mat[is.na(y_mat)] <- 0
  
  x_mat <- data_tmp %>%
    dplyr::select(., contains(x_name)) %>%
    as.matrix()
  
  ## assure positive values
  if (distribution == 'lognormal') { 
    
    y_mat[y_mat == 0] <- 0.1 
    x_mat[x_mat == 0] <- 0.1
    
  }
  
  ## set up stan import list
  s_data <- list(
    K = dim(y_mat)[2],
    N = nrow(y_mat),
    Y = y_mat,
    X = x_mat,
    Z = nrow(y_mat) - h,
    W = W_input
  )
  
  return(s_data)
  
}

# add covariates to stan--------------------------------------------------------

#'@detail adds the age group data as covariates to the stan list
#'@param stan_dat, list, the stan_data
#'@param list_cov, list, the list with the age data
#'@return augmented stan data, including the age covariates

add_covariates_to_stan <- function(
    stan_dat = stan_data,
    list_cov = cov_list) {
  
  stan_dat$X1 <- as.matrix(list_cov[[1]])
  stan_dat$X2 <- as.matrix(list_cov[[2]])
  stan_dat$X3 <- as.matrix(list_cov[[3]])
  stan_dat$X4 <- as.matrix(list_cov[[4]])
  stan_dat$X5 <- as.matrix(list_cov[[5]])
  stan_dat$X6 <- as.matrix(list_cov[[6]])
  stan_dat$X7 <- as.matrix(list_cov[[7]])
  
  return(stan_dat)
  
}

# add age percentages to stan_data ---------------------------------------------

#'@detail modifies the 7(!) covariates pertaining to age groups, turns them
#'into percentages
#'@param input_data, data.frame, the result of 'add_covariates_to_stan()'
#'@return returns the data with percentages instead of integers

add_age_percentages_to_stan_data <- function(
    input_data = stan_data_cov) {
  
  ## avoid dividing by zero
  input_data$X[input_data$X == 0] <- 0.01
  
  input_data$X1 <- (input_data[paste0('X', 1)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  input_data$X2 <- (input_data[paste0('X', 2)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  
  input_data$X3 <- (input_data[paste0('X', 3)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  input_data$X4 <- (input_data[paste0('X', 4)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  
  input_data$X5 <- (input_data[paste0('X', 5)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  
  input_data$X6 <- (input_data[paste0('X', 6)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  input_data$X7 <- (input_data[paste0('X', 7)] %>%
                      as.data.frame()) / input_data$X %>% as.data.frame()
  
  
  input_data$X1 <- round(input_data$X1, digits = 2)
  input_data$X2 <- round(input_data$X2, digits = 2)
  input_data$X3 <- round(input_data$X3, digits = 2)
  input_data$X4 <- round(input_data$X4, digits = 2)
  input_data$X5 <- round(input_data$X5, digits = 2)
  input_data$X6 <- round(input_data$X6, digits = 2)
  input_data$X7 <- round(input_data$X7, digits = 2)
  
  
  return(input_data)
  
}


# extract_parameters -----------------------------------------------------------

#'@description extracts parameters from stan fit objects
#'@param x stan fit;
#'@param param character; name of parameter to be extracted
#'@return augmented stan fit data.frame with additional columns (time and k)

extract_parameters <- function(x = stan_sp_fit, param = "phi") {
  
  
  
  out_list <- list()
  
  dat_tmp <- summary(x, pars = param, probs = c(0.1, 0.9))$summary %>%
    as.data.frame(.) %>%
    mutate(pattern_tmp = str_extract(
      string = rownames(.), pattern = "(?<=\\[).*(?=\\])"),
      yhat = mean)
  
  if (any(str_count(dat_tmp$pattern_tmp, "\\,") == 2)) {
    
    ## if it is a parameter vector with 2 dimensions
    dat_tmp <- dat_tmp %>%
      mutate(
        k_index = sapply(strsplit(.$pattern_tmp,"\\,"), `[`, 3),
        par_id = sapply(strsplit(.$pattern_tmp,"\\,"), `[`, 2),
        time_index = sapply(strsplit(.$pattern_tmp,"\\,"), `[`, 1),
        par_name =  sapply(strsplit(rownames(.),"\\["), `[`, 1),
        par_name = paste0(par_name, par_id)
      ) %>%
      select(-par_id)
    
    out_tmp <- dat_tmp %>%
      group_split(par_name)
    
    for (i in seq(length(out_tmp))) {
      
      out_dat <- out_tmp[[i]] %>% as.data.frame()
      rownames(out_dat) <- paste0(
        out_dat$par_name, '[', out_dat$time_index, ',', out_dat$k_index, ']'
      )
      
      out_dat <- out_dat %>%
        select(-par_name, pattern_tmp)
      
      out_list[[i]] <- out_dat
      
    }
    
  } else {
    
    dat_tmp <- dat_tmp %>%
      mutate(k_index = sub("^[^,]*,", "", pattern_tmp) %>%
               as.integer(),
             time_index = gsub("(.*),.*", "\\1", pattern_tmp) %>%
               as.integer()) %>%
      dplyr::select(-pattern_tmp) %>%
      arrange(k_index, time_index)
    
    if (identical(dat_tmp$time_index, dat_tmp$k_index)) {
      
      dat_tmp$time_index <- NULL
      
    }
    
    out_list[[1]] <- dat_tmp
    
  }
  
  return(out_list)
  
}


#'@details makes problem trend stationary
#'@param input_data, data.frame, data to be differentiated
#'@param prefix, character, string to identify the relevant columns
#'@return data.frame, discrete one step differences of the input data

diff_matrix <- function(
    input_data = data_w_lagged_values$wide,
    prefix = 'icu_') {
  
  data_tmp <- input_data %>%
    dplyr::select(contains(prefix)) %>%
    as.data.frame()
  
  diff_tmp <- lapply(data_tmp, diff, lag = 1) %>%
    as.data.frame()
  
  diff_tmp$date <- input_data$date[2:nrow(data_tmp)]
  
  return(diff_tmp)
  
}

#'@details put the results back to the original scale of the data
#'@param input_data, data.frame, the results from the ensemble containing de-
#'trended (diff) values
#'@param xi_data, data.frame, the original data, serves as constant of 
#'integration
#'@param h, integer, length of horizon
#'@param prefix, character, name of xi
#'@date_xi, date, indicates the initial values for integration
#'@return data.frame, results on the original scale

inv_diff_matrix <- function(
    input_data = test_diff,
    xi_data = data_w_lagged_values$wide,
    date_xi = NULL,
    prefix = 'obs',
    h = H) {
  
  diff_tmp <- input_data %>%
    dplyr::select(-time_index)
  
  if (!is.null(xi_data)) {
    
    c_list <-  xi_data %>%
      dplyr::select(date, contains(prefix)) %>%
      pivot_longer(!date, names_to = 'ILS_ID', values_to = 'obs') %>%
      mutate(ILS_ID = gsub('icu_', '', ILS_ID) %>% as.numeric()) %>%
      filter(date == date_xi) %>%
      pivot_wider(., names_from = ILS_ID, values_from = obs) %>%
      dplyr::select(-date)
    
    ## columns in same order
    diff_tmp <- diff_tmp[, colnames(c_list)]
    
    assert(all(colnames(diff_tmp) == colnames(c_list)))
    
    out_data <- diffinv(as.matrix(diff_tmp), xi = as.matrix(c_list)) %>%
      as.data.frame() %>%
      magrittr::set_colnames(colnames(c_list))
    
  } else {
    
    out_data <- diffinv(as.matrix(diff_tmp)) %>%
      as.data.frame() %>%
      magrittr::set_colnames(colnames(diff_tmp))
    
  }
  
  out_data$time_index <- seq(nrow(out_data)) + h
  
  out_data <- out_data[-1, ] %>%
    mutate(time_index = time_index - 1)
  
  return(out_data)
  
}

# na_locf_covariates -----------------------------------------------------------

#'@details performs na.locf on stan data, since there are occasinally missing
#'values in the RKI covariates
#'@param input_data, list, prepared stan_data_cov
#'@param cov_list, list, name of covariates to be edited
#'@returns input data wihtout missing values


na_locf_covariates <- function(
    input_data = stan_data_cov,
    cov_list = c(paste0('X', seq(7)))) {
  
  for (i in seq(length(cov_list))) {
    
    dat_tmp <- input_data[cov_list[i]]
    
    dat_tmp <- dat_tmp[[1]] %>%
      na.locf(., na.rm = FALSE) %>%
      na.locf(., fromLast = TRUE)
    
    input_data[cov_list[i]] <- list(dat_tmp)
    
  }
  
  return(input_data)
  
}

# prepare_ensemble_data --------------------------------------------------------

#'@details makes data either stationary or brings it back to scale
#'@param input_data, data.frame, ensemble_final or stationary version thereof
#'@param original_data, data.frame, ensemble final data
#'@param focal_var, character, either mod_forecast, lower or upper
#'@return depending on the input provided (original_data), the function either
#'returns first differences or discretely integrated data

prepare_ensemble_data <- function(
    input_data, 
    original_data = NULL,
    focal_var = 'pred'
) {
  
  out_data <- NULL
  
  split_list <- input_data %>%
    group_by(id, iter) %>%
    group_split(.)
  
  if (!is.null(original_data)) {
    
    original_list <- original_data %>%
      group_by(id, iter) %>%
      group_split(.)
    
  }
  
  for (i in seq(length(split_list))) {
    
    split_tmp <- split_list[[i]]
    
    if (!is.null(original_data)) {
      
      original_tmp <- original_list[[i]]
      
      assert(split_tmp$id %>% unique() == original_tmp$id %>% unique())
      assert(split_tmp$iter %>% unique() == original_tmp$iter %>% unique())
      
    }
    
    if (is.null(original_data)) {
      
      ## get first differences
      split_diff_tmp <- split_tmp %>%
        mutate(
          mod_obs = lag(obs),
          mod_forecast = obs
        ) %>%
        diff_matrix(., prefix = 'mod_')
      
      ## add meta info back
      split_tmp <- split_tmp %>%
        dplyr::select(id, time_index, ILS, iter, obs, date, hold_out) %>%
        right_join(
          .,
          split_diff_tmp,
          by = c('date')
        )
      
    } else {
      
      if (focal_var == 'pred') {
        
        ## go back to original scale
        ini_val <-  original_tmp %>%
          filter(time_index == 1) %>%
          dplyr::select(obs) %>%
          as.matrix()
        
        pred_tmp <- split_tmp %>%
          dplyr::select(pred) %>%
          as.matrix() %>%
          diffinv(., xi = ini_val)
        
        split_tmp$pred <- pred_tmp[2:nrow(pred_tmp),]
        
      } else {
        
        pred_tmp <- split_tmp %>%
          dplyr::select(focal_var) %>%
          as.matrix() %>%
          diffinv(.)
        
        split_tmp[, focal_var] <- pred_tmp[2:nrow(pred_tmp),]
        
        
      }
      
    }
    
    out_data <- rbind(out_data, split_tmp)
    
  }
  
  return(out_data)
  
}

# combine_fit_w_observed -------------------------------------------------------

#'@details combines fitted values with observed values in a long data format
#'@param fitted_model, stan_model, fitted stan model
#'@param parameter, character, name of parameter to be extracted
#'@param obsered_data, tibble, long data with actually observed values
#'@param transform_param, boolean, exp of parameters
#'@param ts_type, character, focus on training, forecast or everything 
#'@results data.frame, predicted values alongside observed values for evaluation

combine_fit_w_observed <- function(
    fitted_model = stan_fit, 
    parameter = 'phi', 
    observed_data = training_data$long,
    transform_param = FALSE,
    ts_type = c('all', 'training', 'forecast')) {
  
  out_list <- list()
  
  ## first date of prediction
  hold_out_date <- observed_data %>%
    filter(hold_out == 1) %>%
    filter(date == min(date)) %>%
    dplyr::select(date) %>%
    unique() %>%
    as.data.frame() %>% 
    as.character() %>% 
    as.numeric() %>% 
    as.Date()  # this is madness...
  
  ## extract parameters
  in_tmp <- extract_parameters(x = fitted_model, param = parameter) 
  
  for(i in seq(length(in_tmp))) {
    
    dat_tmp <- in_tmp[[i]]
    
    if (transform_param == TRUE) {
      
      dat_tmp$`10%` <- exp(dat_tmp$`10%`)
      dat_tmp$`90%` <- exp(dat_tmp$`90%`)
      dat_tmp$yhat <- exp(dat_tmp$yhat)
      
    }
    
    
    dat_tmp <- dat_tmp %>% 
      full_join(., 
                observed_data %>%
                  rename(id = `ILS_ID`) %>%
                  dplyr::select(id, ILS) %>%
                  distinct() %>%
                  mutate(k_index = seq(nrow(.))) %>%
                  dplyr::select(k_index, id)
                
      ) %>%
      mutate(k_index = id)
    
    
    ## append observed and fitted
    out_data <- rbind(
      dat_tmp %>%
        mutate(case_count = yhat,
               lower = `10%`,
               upper = `90%`,
               id = k_index,
               type = 'fitted',
               time_index = as.integer(time_index),
               id = as.integer(id)) %>%
        dplyr:: select(case_count, lower, upper, id, time_index, type),   
      
      observed_data %>%
        ungroup() %>%
        mutate(lower = NA,
               upper = NA,
               type = 'observed',
               id = ILS_ID,
               case_count = icu) %>%
        dplyr::select(case_count, lower, upper, id, time_index, type)) %>%
      mutate(type = as.factor(type))
    
    ## add date variables and ILS
    out_data <- out_data %>%
      full_join(.,
                observed_data %>%
                  dplyr::select(time_index, date) %>%
                  distinct(), 
                by = c('time_index')
      ) %>%
      arrange(type, id, date) %>%
      full_join(.,
                observed_data %>%
                  rename(id = `ILS_ID`) %>%
                  dplyr::select(id, ILS) %>%
                  distinct(),
                by = c('id')
      )
    
    if (ts_type == 'training') {
      
      out_data <- out_data %>%
        filter(date < hold_out_date)
      
    } else if (ts_type == 'forecast') {
      
      out_data <- out_data %>%
        filter(date >= hold_out_date)
      
    } else if (ts_type == 'all') {
      
      out_data <- out_data
      
    }
    
    ## export
    conbined_data <- list()
    conbined_data$long <- out_data
    conbined_data$hold_out_date <- hold_out_date
    
    out_list[[i]] <- conbined_data
    
  }
  
  return(out_list)
  
}

# plot_forecasts ---------------------------------------------------------------

#'@details plots the observed and fitted data as time series
#'@param input_data, data.frame, result of combine_fit_w_observed()
#'@param id_var, character, name of grouping variable
#'@param add_intervals, boolean, plot intervals
#'@param ILS, character, name of ILS region to be shown
#'@param hold_out, integer, N of training sample
#'@return returns a ggplot2 object with the time series as facets

plot_forecasts <- function(
    input_data = eval_data, 
    id_var = 'ILS', 
    add_intervals = TRUE,
    ILS = NULL,
    hold_out = 30) {
  
  # if (!is.data.frame(input_data) & is.list(input_data)) {
  #   
  #   ## fetch long data
  #   input_tmp <- fetch_long(x = input_data)
  #   
  # } else {
  #   
  #   input_tmp <- input_data
  #   
  # }
  
  input_tmp <- as.data.frame(input_data)
  
  ## assure factor
  input_tmp[, id_var] <- as.factor(input_tmp[, id_var])
  
  ## filter by ILS
  if (!is.null((ILS))) {
    
    ## assertion
    assertChoice(ILS, input_tmp$ILS %>% unique() %>% as.character())
    
    input_tmp <- input_tmp[which(input_tmp$ILS == ILS), ]
    
  }
  
  ## plot fitted and observed
  p <- input_tmp %>%
    ggplot(., aes(
      x = as.Date(date), 
      y = case_count, 
      color = type,
      ymin = 0)) +
    geom_line() +
    geom_point() +
    geom_vline(
      xintercept = input_tmp$date[input_tmp$time_index == 30] %>% 
        unique() %>%
        as.numeric(), 
      linetype = 4, 
      color = 'black')
  
  if (add_intervals == TRUE) {
    
    p <- p + geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.2,       
      colour ="grey70",
      fill = "red",
      na.rm = TRUE
    )
    
    
    ## add frequentist prediction intervals, if available
    if (c('lower_fp') %in% colnames(input_tmp)) {
      
      p <- p + geom_ribbon(
        aes(ymin = lower_fp, ymax = upper_fp),
        alpha = 0.2,       
        colour ="grey70",
        fill = "orange",
        na.rm = TRUE
      )
      
    }
    
  }
  
  p <- p + labs(
    title = 'Forecast evaluation', 
    color = 'Type' , 
    x = 'Date',
    y = 'ICU')
  
  p <- p + facet_wrap(vars(!!rlang::sym(id_var)))
  
  return(p)
  
}

# plot_all_results -------------------------------------------------------------

#'@details plots specified models at all several iterations
#'@param input_data, data.frame, results_w_intervals
#'@param time_ids, integer, number of iterations
#'@param model_ids, list, integer values pertaining to models' ids
#'@param input_ils, list, name of ILS to be plotted
#'@param model_suffix, character, a suffix for the labels
#'@return cowplot with several models side by side

plot_all_results <-  function(
    input_w_intervals = results_w_intervals,
    time_ids = 20,
    model_ids = 1,
    input_ils = 'Leitstelle München',
    model_suffix = NULL) {

  
  z <- 0
  K <- time_ids %>% sort()
  
  for (m in model_ids) {
    
    for (k in K) {
      
      dat_tmp <- input_w_intervals %>%
        filter(iter == k, model_nr == m)
      
      if (nrow(dat_tmp) > 0) {
        
        plot_tmp <- plot_forecasts(
          input_data = dat_tmp,
          add_intervals = TRUE,
          ILS = input_ils)
        
        ## meta info
        model_tmp <- paste0(
          'Model ID: ', 
          plot_tmp$data$model_nr %>% unique(), model_suffix
        )
        
        range_tmp <- paste0(
          'Range: from ', plot_tmp$data$date %>% min(), 
          ' to ', plot_tmp$data$date %>% max()
        ) 
        
        iter_tmp <- paste0(
          'Iteration: ', plot_tmp$data$iter %>% unique()
        )
        
        stitle_tmp <- paste0(model_tmp, '\n', range_tmp, '\n', iter_tmp)  
        
        plot_tmp <- plot_tmp +  
          labs(title = ' ', subtitle = stitle_tmp)
        
        if (k == K[1] & m == model_ids[1]) {
          
          plot_tmp <- plot_tmp +  
            labs(title = 'Model evaluation', subtitle = stitle_tmp)
          
        }
        
        if (k == K[length(K)] & m == model_ids[length(model_ids)]) {
          
          grobs <- ggplotGrob(plot_tmp)$grobs
          legend <- grobs[[
            which(sapply(grobs, function(x) x$name) == "guide-box")]]
          
        }
        
        plot_tmp <- plot_tmp + theme(legend.position = 'none')
        
        z <- z + 1
        
        assign(paste0('pl_', z), plot_tmp)
        
      }
      
    } 
    
  }
  
  c_grid <- cowplot::plot_grid(plotlist = mget(paste0("pl_", seq(z))))
  
  results_plot <- cowplot::plot_grid(
    c_grid, legend, ncol = 2, rel_widths = c(1, .1))
  
  return(results_plot)
  
}

# make_results_w_intervals -----------------------------------------------------

#'@details turns loss_data into a format suitable for the plotting functions
#'@param input_loss_data, data.frame, loss_data
#'@param remove_list, list, list of none converged models

make_results_w_intervals <- function(
    input_loss_data = loss_data,
    remove_list = NULL) {
  
  if (!is.null(remove_list)) {
    
    ## remove none-converged plots
    loss_tmp <- input_loss_data %>%
      mutate(out_tmp = paste0(as.character(model_nr), '_', iter)) %>%
      filter(!(out_tmp %in% remove_list))
    
  } else {
    
    loss_tmp <- input_loss_data
    
  }
  
  ## out of sample
  oos_tmp <- loss_tmp %>%
    filter(time_index > 30)
  
  ## within sample
  wis_tmp <- loss_tmp %>%
    filter(time_index <= 30)
  
  ## within sample loss
  wis_loss <- wis_tmp %>%
    group_by(model_nr, iter) %>%
    mutate(RMSE_fitted = round(sqrt(mean((pred - obs)^2)), digits = 2),
           MAE_fitted = round(median(abs(pred - obs)), digits = 2)) %>%
    dplyr::select(RMSE_fitted, MAE_fitted, model_nr, iter) %>%
    distinct(.)
  
  ## oos sample loss
  oos_loss <- oos_tmp %>%
    filter(!is.na(obs)) %>%
    group_by(model_nr) %>%
    mutate(
      RMSE = round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2),
      MAE = round(median(abs(pred - obs), na.rm = TRUE), digits = 2),
      time_index = time_index - 30) 
  
  ## add frequentist prediciton intervals
  oos_w_intervals <- oos_loss %>%
    full_join(
      ., 
      wis_loss, 
      by = c('iter', 'model_nr')
    ) %>%
    mutate(
      h = time_index,
      sigma_h = RMSE_fitted * sqrt(h * (1 + h / max(h))),
      upper_fp = pred + 1.96 * sigma_h,
      lower_fp = pred - 1.96 * sigma_h,
      lower_fp = case_when(
        lower_fp < 0 ~ 0,
        lower_fp >= 0 ~ lower_fp
      ),
      upper_fp = case_when(
        upper_fp < 0 ~ 0,
        upper_fp >= 0 ~ upper_fp
      ),
      time_index = time_index + 30
    ) %>%
    dplyr::select(-RMSE_fitted, -MAE_fitted)
  
  dat_tmp <- rbind(
    oos_w_intervals,
    wis_tmp %>%
      mutate(
        h = time_index,
        sigma_h = NA,
        upper_fp = NA,
        lower_fp = NA,
      )
  ) %>%
    arrange(iter, model_nr, time_index)
  
  ## reshape for plotting
  oos_w_intervals <- rbind(
    dat_tmp %>%
      rename(case_count = obs) %>%
      mutate(
        lower = NA,
        upper = NA,
        lower_fp = NA,
        upper_fp = NA,
        pred = NA,
        type = 'observed'
      ) %>%
      dplyr::select(-pred),
    dat_tmp %>%
      rename(case_count = pred) %>%
      mutate(
        type = 'fitted'
      ) %>%
      dplyr::select(-obs)
  )
  
  return(oos_w_intervals)
  
}


# model_validation -------------------------------------------------------------

#'@details re-fits models several times for cross validation of results
#'@param input_training, data.frame, training data for the fitted model
#'@param full_data, data.frame, the complete observed data
#'@param refit_max_iter, positive integer, max number of re-fit iterations
#'@param input_model, list of stan objects, the stan models to be validated
#'@param init_seed, integer, set seed for replication
#'@param date_list, list of dates, the dates for the evaluation procedure
#'@param save_models, boolean, if TRUE - save re-fitted models to wd
#'@return returns a list containing a data.frame with several forecast values
#'for evaluation as well as a list with the re-fitted models

model_validation <- function(
    input_training = training_data,
    full_data = data_w_lagged_values,
    refit_max_iter = 1e+05,
    input_model = model_stan,
    init_seed = get_seed(fit_stan),
    date_list = NULL,
    save_models = FALSE) {
  
  checkmate::assert(!is.null(date_list))
  checkmate::assert(!is.null(input_model))
  
  ## containers
  pred_data <- NULL
  fit_list <- list()
  out_list <- list()
  
  ## counter
  a <- 1
  
  ## model list
  model_list <- list(model_stan)
  
  for (j in seq(length(model_list))) {
    
    ## re-fit the data
    for (i in date_list) {
      
      ## get training data
      training_data_tmp <- extract_time_slice(
        input_data = full_data, 
        max_t = input_training$n_training,
        from_date = i)
      
      ## prepare stan data
      stan_data_tmp <- prepare_stan_data(
        input_data = training_data_tmp,
        log_normal = TRUE
      )
      
      ## re-fit the model
      refit_stan <- rstan::vb(
        object = model_list[[j]], 
        data = stan_data_tmp, 
        iter = refit_max_iter,
        seed = init_seed
      )
      
      ## save fitted models
      if (save_models == TRUE) {
        
        saveRDS(refit_stan, paste0("val_fit_", as.Date(i), "_mod_", j)) 
        
      }
      
      ## get predictions
      eval_data_tmp <- combine_fit_w_observed(
        fitted_model = refit_stan, 
        parameter = 'eta', 
        observed_data = training_data_tmp$long,
        transform_param = TRUE,
        ts_type = 'forecast'
      )
      
      ## meta info
      eval_data_tmp$long$iter <- a
      eval_data_tmp$long$from_date <- i
      eval_data_tmp$long$model <- j
      
      ## store predictions
      pred_data <- rbind(eval_data_tmp$long, pred_data)
      fit_list[[a]] <- refit_stan
      
      ## increment counter
      a <- a + 1
      
    }
    
  }
  
  ## post process dates
  pred_data$from_date <- as.Date(pred_data$from_date)
  
  ## export
  out_list$forecast <- pred_data
  out_list$fitted_models <- fit_list
  out_list$input_models <- model_list
  
  return(out_list)
  
}

# select dates -----------------------------------------------------------------

#'@details selects dates for the validation of the models
#'@param input_phases, data.frame, phases_data with the phases and dates
#'@param full_data, data.frame, the complete observed data
#'@param n_training, positive integer, number of training samples
#'@param horizon, positive integer, forecast horizon 
#'@param n_samples, positive integer, number of samples per stratum
#'@param rolling_window, boolean, rolling evaluation or stratified sample
#'@return returns a list with dates for evaluating the models as well as a
#'data.frame with the edited dates of the phases of the pandemic

select_dates <- function(
    input_phases = phases_data,
    full_data = data_w_lagged_values,
    n_training = K,
    horizon = H,
    n_samples = J,
    rolling_window = FALSE) {
  
  K <- n_training
  H <- horizon
  J <- n_samples
  out_list <- list()
  
  if (is.null(J)) {J <- H}
  
  ## extract dates of phases
  phases_tmp <- input_phases %>%
    mutate(
      end_date = case_when(
        phase == max(phase) ~ max(full_data$long$date),
        TRUE ~ as.Date(end_date)),
      start_date = as.Date(start_date),
      n_days = difftime(end_date, start_date, units = 'days')
    ) %>%
    ## smallest date for sampling in the particular phase
    mutate(
      end_sampling = end_date - (K + H),
      n_days = difftime(end_sampling, start_date, units = 'days')
    ) %>%
    ## phase must be at least as long as H
    ## the min date must be observed
    filter(
      n_days >= H,
      start_date >= min(full_data$long$date, na.rm = TRUE)
    )
  
  ## make date list
  eval_dates <- NULL
  
  if (rolling_window == FALSE) {
    
    ## sample randomly from each phase
    for (i in seq(nrow(phases_tmp))) {
      
      ## temp dates
      dates_tmp <- phases_tmp$start_date[i]:phases_tmp$end_sampling[i] %>%
        as.Date()
      
      ## list of starting dates for evaluation 
      eval_dates <- c(eval_dates, sample(dates_tmp, J))
      
      ## number of strata
      out_list$meta$n_strata <- nrow(phases_tmp)
      
    }
    
  } else if (rolling_window == TRUE) {
    
    ## rolling window cross validation
    if (is.null(date_list)) {
      
      Z <- J - (K + H)
      
      ## extract variables from input
      min_date <- input_training$from_date
      date_list <- as.Date(min_date:(min_date + H + Z))
      
      out_list$meta$n_strata <- 1
      
    }
    
  }
  
  eval_dates <- as.Date(eval_dates)
  
  ## export
  out_list$phases <- phases_tmp
  out_list$eval_dates <- eval_dates
  out_list$meta$samples_per_stratum <- J
  out_list$meta$n_dates <- length(eval_dates)
  
  return(out_list)
  
}


# assess_out_of_sample ---------------------------------------------------------

#'@details out of sample fit by model
#'@param input_data, data.frame, results of model_validation
#'@param by_group, character, choose grouping (see default values)
#'@return returns a data.frame with loss metrics

assess_out_of_sample <- function(
    input_data = cv_data$forecast, 
    by_group = c('ILS and steps ahead', 'steps ahead', 'overall')) {
  
  ## assertions
  assertChoice(by_group, c('ILS and steps ahead', 'steps ahead', 'overall'))
  
  ## prepare data for comparison
  data_tmp <- input_data %>%
    dplyr::select(-lower, -upper, -from_date)
  
  data_tmp <- data_tmp %>%
    filter(type == 'observed') %>%
    rename(cases_obs = `case_count`) %>%
    dplyr::select(-type) %>%
    full_join(.,
              data_tmp %>%
                filter(type == 'fitted') %>%
                rename(cases_fit = `case_count`) %>%
                dplyr::select(-type, -model),
              by = c('id','time_index', 'date', 'ILS', 'iter')
    ) %>%
    group_by(model)
  
  
  if (by_group == 'ILS and steps ahead') {
    
    data_tmp <- data_tmp %>%
      group_by( ILS, time_index) %>%
      calculate_loss_metrics(.) %>%
      ungroup() %>%
      dplyr::select(ILS, id, RMSE, MAE, time_index) %>%
      distinct() %>%
      arrange(ILS)
    
  } else if (by_group == 'steps ahead') {
    
    data_tmp <- data_tmp %>%
      group_by(time_index) %>%
      calculate_loss_metrics(.) %>%
      ungroup() %>%
      dplyr::select(RMSE, MAE, time_index) %>%
      distinct() %>%
      arrange(time_index)  
    
  } else if (by_group == 'overall') {
    
    data_tmp <- data_tmp %>%
      calculate_loss_metrics(.) %>%
      dplyr::select(RMSE, MAE) %>%
      distinct()
    
  }
  
  return(data_tmp)
  
}

# calculate loss metrics -------------------------------------------------------

#'@details calculates loss metrics for assessing out of sample performance
#'@param input_data, data.frame, data_tmp object from assess_out_of_sample()
#'@return returns the input with loss metrics as additonal columns

calculate_loss_metrics <- function(input_data) {
  
  input_data %>%
    mutate(
      RMSE = sqrt(mean((cases_obs - cases_fit)^2, na.rm = TRUE)),
      MAE = mean(abs(cases_obs -cases_fit), na.rm = TRUE)
    )
}


# fetch_long -------------------------------------------------------------------

#'@details fetch the long data from list, if it exists
#'@param x list or data.frame, eval_data or cv_data
#'@return returns the long data 

fetch_long <- function(x) {
  
  if (exists('long', where = x)) {
    
    out <- x$long
    
  } else if (exists('forecast', where = x)) {
    
    out <- x$forecast
    
  } else {
    
    out <- x
    
  }
  
  return(out)
  
}

