library(targets)
library(tarchetypes)
library(future)
library(future.callr)

# Use the future package to run the workflow in parallel with multiple
# processes locally
plan(callr)

# Read functions necessary for the workflow. These functions make up what most
# targets, or steps, in the workflow do
source("R/functions.R")

# Set defaults for upcoming targets
tar_option_set(packages = c("tidyverse", "furrr", "future", "lubridate"),
               garbage_collection = TRUE)


# Define targets in workflow ----------------------------------------------

list(
  
  # Input files to track ----------------------------------------------------
  
  # EULI data, source: https://knb.ecoinformatics.org/view/doi:10.5063/F12V2D1V
  tar_file(euli_path,
           command = "data/inputs/SGL.5.1.csv"),
  
  # Berkeley Earth temperature data ("BEST"), source: http://berkeleyearth.org/data/
  tar_file(name = best_path,
           command = "data/inputs/best_daily"),
  
  # CRU files, source: https://data.ceda.ac.uk//badc/cru/data/cru_ts/cru_ts_4.01/data/
  # Air temperature (Note, not used in final modeling)
  tar_file(name = cru_airtemp_path,
           command = "data/inputs/cru_data/cru_ts4.01.1901.2016.tmp.dat.nc"),
  # Precipitation
  tar_file(name = cru_precip_path,
           command = "data/inputs/cru_data/cru_ts4.01.1901.2016.pre.dat.nc"),
  
  # Lake shapefiles, source: https://www.naturalearthdata.com/downloads/50m-physical-vectors/
  tar_file(name = lakes_path,
           command = "data/ne_50m_lakes.shp"),
  
  
  # Targets / steps of the workflow -----------------------------------------
  
  # Open raw EULI
  tar_target(euli,
             read_csv(file = euli_path)),
  
  # Clean EULI
  tar_target(euli_filtered,
             clean_euli(euli),
             packages = c("tidyverse", "furrr", "future", "lubridate", "ncdf4",
                          "zoo", "data.table", "sqldf")),
  
  # Pull BEST dataset
  tar_target(best_data,
             command = pull_best_data(best_path = best_path,
                                      euli_filtered = euli_filtered),
             packages = c("tidyverse", "furrr", "future", "lubridate", "ncdf4",
                          "sqldf")),
  
  # Join BEST and EULI datasets
  tar_target(best_euli_join,
             command = join_best_euli(euli_filtered = euli_filtered,
                                      compiled_best_temps = best_data),
             packages = c("tidyverse", "furrr", "future", "lubridate", "sqldf")),
  
  # Estimate ice-formation dates
  tar_target(iceform_dates,
             command = determine_iceform(best_euli_join,
                                         euli_filtered = euli_filtered),
             packages = c("tidyverse", "furrr", "future", "lubridate", "zoo")),
  
  # Estimate pre-stratification dates
  tar_target(prestrat_dates,
             command = determine_prestrat(iceform_dates$tavg_rollmean),
             packages = c("tidyverse", "furrr", "future", "lubridate", "zoo")),
  
  # Take a subset of euli_filtered columns
  tar_target(euli_slim,
             command = euli_filtered %>%
               select(lakename, poolstation, season, euli_year, lat, long,
                      start_date_min, end_date_max)),
  
  # Calculate BEST temperature averages during each season/period
  tar_target(best_tavg,
             command = calc_best_tavg(tavg_full = best_euli_join,
                                      iceform_period_rollmean = iceform_dates$iceform_period_rollmean,
                                      prestrat_period_rollmean = prestrat_dates)),
  
  # Extract SLP data
  tar_target(raw_slp,
             command = pull_agg_slp(),
             packages = c("targets", "tidyverse", "raster", "lubridate",
                          "furrr", "RNCEP")),
  
  # Match SLP locations with EULI locations
  tar_target(slp_by_location,
             command = extract_slp(slp_file_location = raw_slp,
                                   location_data = euli_slim),
             packages = c("tidyverse", "furrr", "future", "lubridate", "sf")),
  
  # Filter SLP for winter dates
  tar_target(filtered_slp_winter,
             command = filter_slp_winter(extracted_slp = slp_by_location,
                                         euli_slim = euli_slim,
                                         iceform_period_rollmean = iceform_dates$iceform_period_rollmean),
             packages = c("tidyverse", "furrr", "future", "lubridate", "janitor", "sf")),
  
  # Filter SLP for summer dates
  tar_target(filtered_slp_summer,
             command = filter_slp_summer(extracted_slp = slp_by_location,
                                         euli_slim = euli_slim,
                                         prestrat_period_rollmean = prestrat_dates),
             packages = c("tidyverse", "furrr", "future", "lubridate", "janitor", "sf")),
  
  # Process SLP after extraction
  tar_target(processed_slp,
             command = process_slp(daily_slp_iceform = filtered_slp_winter$df_period,
                                   daily_slp_underice = filtered_slp_winter$df_sample,
                                   daily_slp_prestrat = filtered_slp_summer$df_period,
                                   daily_slp_strat = filtered_slp_summer$df_sample),
             packages = c("tidyverse", "furrr", "future", "lubridate")),
  
  # Extract CRU monthly air temperature data
  tar_target(cru_airtemp,
             extract_cru_airtemp(cru_airtemp_path = cru_airtemp_path,
                                 winter_slp_lake_centrs = filtered_slp_winter$winter_slp_lake_centrs, 
                                 summer_slp_lake_centrs = filtered_slp_summer$summer_slp_lake_centrs),
             packages = c("tidyverse", "furrr", "future", "lubridate", "ncdf4")),
  
  # Extract CRU monthly precipitation data
  tar_target(cru_precip,
             command = extract_cru_precip(cru_precip_path = cru_precip_path,
                                          winter_slp_lake_centrs = filtered_slp_winter$winter_slp_lake_centrs, 
                                          summer_slp_lake_centrs = filtered_slp_summer$summer_slp_lake_centrs),
             packages = c("tidyverse", "furrr", "future", "lubridate", "ncdf4")),
  
  
  # Build the modeling datasets ---------------------------------------------
  
  # Winter
  tar_target(winter_compiled_data,
             command = build_winter(precip_iceform = cru_precip$precip_monthly_iceform,
                                    tmp_iceform = cru_airtemp$air_temp_monthly_iceform,
                                    precip_underice = cru_precip$precip_monthly_underice,
                                    tmp_underice = cru_airtemp$air_temp_monthly_underice,
                                    tavg_iceform = best_tavg$best_iceform_period_tavg,
                                    tavg_underice = best_tavg$best_underice_tavg,
                                    slp_iceform = processed_slp$slp_iceform_period,
                                    slp_underice = processed_slp$slp_underice_period,
                                    euli_filtered = euli_filtered),
             packages = c("tidyverse", "furrr", "future", "lubridate",
                          "corrplot")),
  
  # Summer
  tar_target(summer_compiled_data,
             command = build_summer(precip_prestrat = cru_precip$precip_monthly_prestrat,
                                    tmp_prestrat = cru_airtemp$air_temp_monthly_prestrat,
                                    precip_strat = cru_precip$precip_monthly_strat,
                                    tmp_strat = cru_airtemp$air_temp_monthly_strat,
                                    tavg_prestrat = best_tavg$best_prestrat_period_tavg,
                                    tavg_strat = best_tavg$best_strat_tavg,
                                    slp_prestrat = processed_slp$slp_prestrat_period,
                                    slp_strat = processed_slp$slp_strat_period,
                                    euli_filtered = euli_filtered),
             packages = c("tidyverse", "furrr", "future", "lubridate",
                          "corrplot")),
  
  # Take a look at the time series for Lake Erie as a means of quality control
  tar_target(single_lake_check,
             command = viz_single_lake(
               winter_climate_euli = winter_compiled_data$winter_climate_euli,
               summer_climate_euli = summer_compiled_data$summer_climate_euli)),
  
  # Finalize and export the built datasets
  tar_target(modeling_datasets,
             command = export_model_datasets(
               winter_climate_euli = winter_compiled_data$winter_climate_euli,
               summer_climate_euli = summer_compiled_data$summer_climate_euli)),
  
  # Run models:
  
  # Specify model using a vector of vars
  tar_target(winter_explanatory,
             command = c("TP_log10", "TN_log10", "iceform_mean_tavg",
                         "iceform_mean_slp", "iceform_var_slp", "underice_var_slp",
                         "iceform_mean_precip", "underice_mean_precip",
                         "underice_mean_slp")),
  
  tar_target(summer_explanatory,
             command = c("TP_log10", "TN_log10", "watertemp", "strat_mean_precip",
                         "strat_mean_slp", "strat_var_slp", "prestrat_mean_slp",
                         "prestrat_var_slp" , "prestrat_mean_precip",
                         "prestrat_mean_tavg")),
  
  # Define the resampling method
  tar_target(kfold_control,
             command =
               # Resampling parameters
               trainControl(
                 method = "repeatedcv",
                 # Number of folds:
                 number = 5,
                 repeats = 3,
                 savePredictions = "final"),
             packages = c("tidyverse", "furrr", "future", "lubridate", "caret")),
  
  # Winter random forest
  tar_target(winter_rf,
             command = run_and_plot_winter_model(winter_model_data = modeling_datasets$winter_model_data,
                                                 winter_explanatory = winter_explanatory,
                                                 kfold_control = kfold_control),
             packages = c("tidyverse", "furrr", "future", "lubridate",
                          "corrplot", "randomForest", "caret", "pdp", "rpart",
                          "rpart.plot", "ggpubr", "cowplot", "partykit",
                          "janitor")),
  
  # Winter regression tree
  tar_target(winter_tree,
             command = run_winter_tree(winter_model_data = modeling_datasets$winter_model_data,
                                       winter_explanatory = winter_explanatory),
             packages = c("tidyverse", "future", "lubridate", "rpart")),
  
  # Winter linear model
  tar_target(winter_lm,
             command = run_winter_lm(winter_model_data = modeling_datasets$winter_model_data),
             packages = c("tidyverse", "future", "lmerTest", "mgcv", "car")),
  
  # Summer random forest
  tar_target(summer_rf,
             command = run_and_plot_summer_model(summer_model_data = modeling_datasets$summer_model_data,
                                                 summer_explanatory = summer_explanatory,
                                                 kfold_control = kfold_control),
             packages = c("tidyverse", "furrr", "future", "lubridate",
                          "corrplot", "randomForest", "caret", "pdp", "rpart",
                          "rpart.plot", "ggpubr", "cowplot", "partykit",
                          "janitor")),  
  
  # Summer regression tree
  tar_target(summer_tree,
             command = run_summer_tree(summer_model_data = modeling_datasets$summer_model_data,
                                       summer_explanatory = summer_explanatory),
             packages = c("tidyverse", "future", "lubridate", "rpart")),
  
  # Summer linear model
  tar_target(summer_lm,
             command = run_summer_lm(summer_model_data = modeling_datasets$summer_model_data),
             packages = c("tidyverse", "future", "lmerTest", "mgcv", "car")),
  
  # Create and export a map figure of data locations
  tar_target(map_figure,
             command = create_map(winter_model_data = modeling_datasets$winter_model_data,
                                  summer_model_data = modeling_datasets$summer_model_data,
                                  euli_filtered = euli_filtered,
                                  lakes_path = lakes_path),
             packages = c("tidyverse", "furrr", "future", "lubridate",
                          "rgdal", "zoo", "viridisLite", "ggrepel", "sf",
                          "cowplot")),
  
  # Summary stats for linear models
  tar_target(lm_stats_out,
             command = {
               winterStat <- broom::tidy(winter_lm$model)
               summerStat <- broom::tidy(summer_lm$model)
               
               outStats <- rbind(winterStat, summerStat) %>%
                 data.frame()
               
               outStats[, "model"] <- c(rep("winter", nrow(winterStat)),
                                        rep("summer", nrow(summerStat)))
               
               write_csv(x = outStats, file = "data/outputs/StatsOut.csv")
               
               return(list(lm_stats_table = outStats,
                           lm_stats_out_path = "data/outputs/StatsOut.csv"))
               
             }
  ),
  
  # Summary table for lakes
  tar_target(summary_table,
             command = build_summary_table(euli_filtered = euli_filtered,
                                           winter_model_data = modeling_datasets$winter_model_data,
                                           summer_model_data = modeling_datasets$summer_model_data)),
  
  # Generate an R Markdown report of model results
  tar_render(report,
             path = "model_results.Rmd",
             packages = c("targets", "caret", "randomForest", "kableExtra",
                          "rpart.plot", "effects", "car", "broom", "mgcv",
                          "rpart", "ggpubr")),
  
  
  # Output files to track ---------------------------------------------------
  
  # SLP data
  tar_file(name = slp_path,
           command = raw_slp),
  
  # Histogram of first rolling average 0C dates
  tar_file(name = iceform_histogram,
           command = iceform_dates$iceform_hist_out),
  
  # Correlation plot for winter variables
  tar_file(name = winter_cor_plot,
           command = winter_compiled_data$winter_climate_cor_path),
  
  # Correlation plot for summer variables
  tar_file(name = summer_cor_plot,
           command = summer_compiled_data$summer_climate_cor_path),
  
  # Lake Erie time series
  tar_file(name = single_lake_timeline_plot,
           command = single_lake_check),
  
  # Winter modeling dataset
  tar_file(name = winter_dataset,
           command = modeling_datasets$winter_data_out_path),
  
  # Summer modeling dataset
  tar_file(name = summer_dataset,
           command = modeling_datasets$summer_data_out_path),
  
  # Correlation plot for winter modeling dataset
  tar_file(name = winter_model_cor_plot,
           command = winter_rf$winter_model_cor_path),
  
  # Winter RF plots
  tar_file(name = winter_pdp_plot_grid,
           command = winter_rf$winter_pdp_path),
  
  # Correlation plot for summer modeling dataset  
  tar_file(name = summer_model_cor_plot,
           command = summer_rf$summer_model_cor_path),
  
  # Summer RF plots
  tar_file(name = summer_pdp_plot_grid,
           command = summer_rf$summer_pdp_path),
  
  # Map figure
  tar_file(name = map_figure_file,
           command = map_figure),
  
  # Linear model summary stats
  tar_file(name = lm_stats_file,
           command = lm_stats_out$lm_stats_out_path)
)





