
#'@details procedure for making the alignemnt plot between the
#'ln(ICU_t) and ln(case_counts_t) series in bavaria

# plot without shifting #-------------------------------------------------------

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## get ggplot colors
ggc <- gg_color_hue(2)

test_lag <- lag_covariates(
  input_data = merged_data$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 0,
  cov_name = 'incidence'
)

sm_data <- test_lag$long %>%
  group_by(date) %>%
  summarise(
    Incidence = log(sum(incidence)),
    #ICU = sum(icu)#,
    ICU_log = log(sum(icu))
  ) %>%
  pivot_longer(!date, names_to = 'type', values_to = 'value') 

sm_ext <- sm_data %>% group_by(type) %>%
  mutate(smooth =predict(loess(value ~ as.numeric(date), span = .15))) %>% 
  group_by(type) %>%
  mutate(
    smooth_lag = smooth - lag(smooth),
    sign_lag = sign(smooth_lag),
    sum_sign = sign_lag + lag(sign_lag)
  ) %>%
  ungroup() %>%
  filter(
    sum_sign == 0
  ) %>%
  dplyr::select(date, type, value, smooth)

sm_ext <- sm_ext %>%
  filter(
    !(
      (date %in% c(
        as.Date('2022-03-11'), 
        as.Date('2022-03-17')#,
        #as.Date('2022-02-13')
        )) &
        type == 'ICU_log'
    )
  )

# sm_ext$grouping <- rep(seq(1:10), 2) %>% sort()
# 
# stats_p1 <- sm_ext %>%
#   group_by(grouping) %>%
#   mutate(
#     diff_time = as.numeric(date - lag(date))
#   )
# 
# median(stats_p1$diff_time, na.rm = TRUE)

p1 <- test_lag$long %>%
  filter(date >= min(date) + 20) %>%
  group_by(date) %>%
  summarise(
    Incidence = log(sum(incidence)),
    #ICU = sum(icu)#,
    ICU_log = log(sum(icu))
  ) %>%
  pivot_longer(!date, names_to = 'type', values_to = 'value') %>%
  ggplot(., aes(x = date, y = value, color = type)) +
  geom_point(alpha = 0.2) +
  geom_line(alpha = 0.2) +
  geom_smooth(span = 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +  
  geom_vline(
    xintercept = sm_ext$date[which(sm_ext$type == 'ICU_log')], 
    #linetype ="dotted", 
    color = ggc[1], 
    size = 1,
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = sm_ext$date[which(sm_ext$type == 'Incidence')], 
    #linetype ="dotted", 
    color = ggc[2], 
    size = 1,
    alpha = 0.5
  ) +
  geom_point(data=sm_ext, aes(y=smooth),color="black", size = 2) +
  scale_colour_discrete(
    name = 'Time series',
    labels = c('ln(ICU)', 'ln(case counts)')) +
  labs(
    subtitle = 'Observed dates',
    x = 'Date', y = "ln(observations)",
    title = 'Time series on logarithmic scale (Bavaria)'
  )

# plot shifted series ----------------------------------------------------------

test_lag <- lag_covariates(
  input_data = merged_data$icu_data %>%
    dplyr::select(-starts_with('Cat_')),
  k = 20,
  cov_name = 'incidence'
)

sm_data2 <- test_lag$long %>%
  group_by(date) %>%
  summarise(
    Incidence = log(sum(incidence)),
    #ICU = sum(icu)#,
    ICU_log = log(sum(icu))
  ) %>%
  pivot_longer(!date, names_to = 'type', values_to = 'value') 

sm_ext2 <- sm_data2 %>% group_by(type) %>%
  mutate(smooth =predict(loess(value ~ as.numeric(date), span = .15))) %>% 
  group_by(type) %>%
  mutate(
    smooth_lag = smooth - lag(smooth),
    sign_lag = sign(smooth_lag),
    sum_sign = sign_lag + lag(sign_lag)
  ) %>%
  ungroup() %>%
  filter(
    sum_sign == 0
  ) %>%
  dplyr::select(date, type, value, smooth)

sm_ext2 <- sm_ext2 %>% 
  filter(
    !(
      (date %in% c(
        as.Date('2020-12-16'), 
        as.Date('2020-12-25'))) & 
        type == 'Incidence'
    )
    # ,
    # !(
    #   (date %in% c(
    #     as.Date('2022-02-16'),
    #     as.Date('2022-02-17'),
    #     as.Date('2022-03-11'), 
    #     as.Date('2022-03-17')
    #     )) &
    #     type == 'ICU_log'
    # )
  ) 

sm_ext2 <- sm_ext2 %>%
  filter(type != 'ICU_log') %>%
  rbind(
    .,
    sm_ext %>%
        filter(type == 'ICU_log') %>%
      dplyr::select(colnames(sm_ext2))
    ) %>%
  arrange(date)
  
p2 <- test_lag$long %>%
  group_by(date) %>%
  summarise(
    Incidence = log(sum(incidence)),
    #ICU = sum(icu)#,
    ICU_log = log(sum(icu))
  ) %>%
  pivot_longer(!date, names_to = 'type', values_to = 'value') %>%
  ggplot(., aes(x = date, y = value, color = type)) +
  geom_point(alpha = 0.2) +
  geom_line(alpha = 0.2) +
  geom_smooth(span = 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(
    xintercept = sm_ext2$date[which(sm_ext2$type == 'ICU_log')], 
    #linetype ="dotted", 
    color = ggc[1], 
    size = 1,
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = sm_ext2$date[which(sm_ext2$type == 'Incidence')], 
    #linetype ="dotted", 
    color = ggc[2], 
    size = 1,
    alpha = 0.5
  ) +
  geom_point(data=sm_ext2, aes(y=smooth),color="black", size = 2) +
  scale_colour_discrete(
    name = 'Time series',
    labels = c('ln(ICU)', 'ln(case counts)')) +
  labs(
    subtitle = 'ln(case counts) shifted by + 20 days',
    x = 'Date', y = "ln(observations)",
    title = ''
  )

log_series_plot <- cowplot::plot_grid(
  p1, p2, ncol = 1, labels = 'AUTO', align = 'h')
#saveRDS(log_series_plot, 'log_series_plot.rds')
