# icu_occupancy

## ./scripts
./scripts contains the R code for this project

- 0_functions includes functions for data wrangling, fitting, plotting, etc.
- 01_prepare data edits the raw data, fits the models and generates summary statistics and plots
- 02_log_ts_alignemt_plots is a procedure for generating the plot showing the alignment between the covariate and the outcome
- 03_make_moran_plot generates the plot showing the spatial-autocorrelation over time
- 04_evaluation generates loss metrics and interval estimates
- 05_prepare_ensembles pools state-space model results (further research)
- 06_results_germany generates plots and tables on the german data

## ./models
./models contains the Stan code for the state space models

The Stan code relates to the model IDs in the text as follows:

- model_01_rw2_normal == model ID 1.1
- model_01_rw2_normal_MULTI == model ID 1.2
- model_01_rw2_normal == models ID 2.1 to 4.1
- model_01_rw2_normal_MULTI == model ID 2.2 to 4.2
- model_03_rw2_normal_lognormal == model ID 5.1 
- model_03_rw2_normal_lognormal_MULTI == model ID 5.2
- model_04_rw2_normal_lognormal == model ID 6.1 to 8.1 
- model_04_rw2_normal_lognormal_MULTI == model ID 6.2 to 8.2
