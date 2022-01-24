# This script stores functions used for the stormy analysis workflow


# 1. Larger processes -----------------------------------------------------


# 1.1 Data import and prep ------------------------------------------------

# Manipulate the EULI dataset for greater ease of use
clean_euli <- function(euli){
  
  euli <- euli %>%
    # New numeric month cols
    mutate(start_month_num = match(startmonth, month.abb),
           end_month_num = match(endmonth, month.abb),
           start_ymm = as.yearmon(paste(startyear,
                                        start_month_num,
                                        sep = "-")),
           end_ymm = as.yearmon(paste(endyear,
                                      end_month_num,
                                      sep = "-")),
           start_date = as.Date(paste(startyear,
                                      start_month_num,
                                      startday,
                                      sep = "-")),
           end_date = as.Date(paste(endyear,
                                    end_month_num,
                                    endday,
                                    sep = "-")))
  
  
  euli_filtered <- euli %>%
    # Keep lakes with open season
    filter(stationlat > -66) %>%
    # Take an average of all coord pairs for a given station
    # Will make it so that 1 station name = 1 pair of coords
    group_by(lakename, poolstation) %>%
    mutate(lat = mean(stationlat), long = mean(stationlong)) %>%
    # If there are multiple sampling points within the same season*station BUT
    # they have different start and end dates, keep the min/max bounds of these
    # periods and collapse into one averaged record
    group_by(lakename, poolstation, season, year) %>%
    mutate(start_ymm_min = min(start_ymm),
           end_ymm_max = max(end_ymm),
           start_date_min = min(start_date),
           end_date_max = max(end_date)) %>%
    ungroup() %>%
    group_by(lakename, poolstation, season, start_date_min, end_date_max) %>%
    # Summarize the numeric columns at this level...drops most categorical columns
    summarise(across(.cols = where(is.numeric),
                     .fns = mean, na.rm = T),
              lakecountry = unique(lakecountry)) %>%
    ungroup() %>%
    # Round latitudes for joins later
    mutate(lat = round(lat, 5)) %>%
    # Rename for clarity downstream
    rename(euli_year = year)
  
  # Return out path for tracking
  return(euli_filtered)
  
}


# Extract data from Berkeley Earth ("BEST")
pull_best_data <- function(best_path, euli_filtered){
  
  # Store location of BEST netCDF files
  # Source: http://berkeleyearth.org/data/
  tavg_files <- dir(path = best_path,
                    full.names = TRUE,
                    pattern = ".nc")
  
  # A data frame of places for which to pull BEST data
  lake_centrs <- euli_filtered %>%
    ungroup() %>%
    select(lakename,
           poolstation,
           Y = lat, X = long) %>%
    unique()
  
  # Pull the climatology data first:
  euli_climatology <- future_map_dfr(.x = tavg_files,
                                     .f = ~get_climatology(latlong = lake_centrs,
                                                           file_name = .x)) %>%
    distinct()
  
  # Pull anomalies
  euli_anomalies <- future_map_dfr(.x = tavg_files,
                                   .f = ~get_anomaly(latlong = lake_centrs,
                                                     file_name = .x))
  
  # All lakes and stations should have the same counts...
  euli_anomalies %>%
    count(lakename, poolstation) %>%
    pull(n) %>%
    unique()
  
  # Make the daily temperature dataset by adding anomaly + climatology.
  # If you run into a memory error, try restarting R and only loading packages +
  # the two datasets needed for the join
  compiled_best_temps <- full_join(x = euli_anomalies,
                                   y = euli_climatology,
                                   by = c("lakename", "poolstation",
                                          "doy" = "day_num")) %>%
    # The sqldf methods that follow don't like date formats
    mutate(tavg_date = as.character(tavg_date))
  
  # Remove objects for efficiency
  rm(euli_anomalies, euli_climatology)
  gc()
  
  # Memory-friendly way to get unique combinations
  compiled_best_temps <- sqldf(x = "SELECT DISTINCT *
                                    FROM compiled_best_temps")
  
  compiled_best_temps <- compiled_best_temps %>%
    mutate(tavg = anomaly_c + climatology_c)
  
  # Return paths for tracking and temperature data for use
  return(compiled_best_temps)
  
}


# Combine BEST with EULI data
join_best_euli <- function(euli_filtered, compiled_best_temps){
  
  # Remove date format for use with sqldf join
  euli_filtered <- euli_filtered %>%
    mutate(across(.cols = contains("date"), .fns = ~as.character(.)))
  
  # Usually would run sqldf joins outside of R with .csv files, but doing it
  # within R if possible this time, b/c there's a comma in one or more lakename
  # strings that causes sqldf to think lakename is two columns of data
  tavg_full <- sqldf(x = "SELECT DISTINCT tavg_daily.lakename, tavg_daily.poolstation, euli.euli_year,
                                    euli.season, euli.start_date_min,euli.end_date_max,
                                    DATE(tavg_daily.tavg_date) tavg_date,
                                    DATE(euli.end_date_max, '+1 year') int_end,
                                    DATE(euli.start_date_min, '-1 year') int_start,
                                    tavg_daily.tavg
                                  FROM compiled_best_temps tavg_daily
                                    JOIN euli_filtered euli
                                      ON tavg_daily.lakename = euli.lakename
                                        AND ((tavg_daily.poolstation = euli.poolstation) OR
                                            (tavg_daily.poolstation IS NULL AND euli.poolstation IS NULL))
                                        AND tavg_date >= int_start
                                        AND tavg_date <= int_end")
  
  rm(compiled_best_temps)
  gc()
  
  tavg_full <- tavg_full %>%
    select(-contains("int"))
  
  # Make sure that each sampling period has a year of data on either side:
  timeline_check <- tavg_full %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(len_period = ymd(end_date_max) - ymd(start_date_min),
              num_days = n(),
              added_days = num_days - len_period)
  
  # All within the range of leap year error:
  unique(timeline_check$added_days)
  
  return(tavg_full)
  
}


# Acquire dates around under-ice sampling 
determine_iceform <- function(tavg_full, euli_filtered){
  
  # First will calculate the two-week rolling temperature average for each
  # lake*station combination. Recall that data are +/- 1 year around each
  # start_date_min
  tavg_rollmean <- future_map_dfr(
    # Turn df of year*dates into a list by row
    .x = tavg_full %>%
      select(lakename, poolstation, euli_year,
             season, start_date_min, end_date_max) %>%
      unique() %>%
      purrr::transpose(),
    # Then for each list item...
    .f = ~ tavg_full %>%
      filter(euli_year == .x$euli_year,
             lakename == .x$lakename,
             season == .x$season,
             # Filter by lakename, and make sure to catch NA lakenames so they
             # aren't dropped
             (poolstation == .x$poolstation |
                is.na(poolstation) & is.na(.x$poolstation))) %>%
      select(tavg, tavg_date, euli_year) %>%
      arrange(tavg_date) %>%
      mutate(roll_mean = rollapply(data = tavg, width = 14, FUN = mean,
                                   align = "right", fill = NA),
             lakename = .x$lakename,
             poolstation = .x$poolstation,
             season = .x$season,
             start_date_min = as_date(.x$start_date_min),
             end_date_max = as_date(.x$end_date_max),
             tavg_date = as_date(tavg_date)) %>%
      filter(!is.na(roll_mean))
  )
  
  # Make sure that each sampling period has a year of data on either side:
  timeline_check <- tavg_rollmean %>%
    distinct() %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(len_period = ymd(end_date_max) - ymd(start_date_min),
              num_days = n(),
              added_days = num_days - len_period)
  
  # All within the range of leap year error after adding 14 lost days for the 
  # leadup to the first rolling average date:
  unique(timeline_check$added_days)
  
  rm(timeline_check)
  gc()
  
  # Filter to remove data from outside a six-week window leading up to the
  # start_date_min. Basically, look only at data in the six weeks before the
  # earliest recorded sampling for the season
  iceform_period_rollmean <- tavg_rollmean %>%
    filter(tavg_date >= (as_date(start_date_min) - weeks(6)),
           tavg_date <= as_date(start_date_min),
           season == "iceon", roll_mean <= 0) %>%
    group_by(lakename, poolstation, euli_year, season, start_date_min, end_date_max) %>%
    arrange(tavg_date) %>%
    summarize(first_szero_date = first(tavg_date),
              start_date_min = unique(start_date_min),
              end_date_max = unique(end_date_max),
              # Length of time between the first 0C date and the earliest sample date
              szero_diff_start = as.integer(first_szero_date - start_date_min),
              iceform_lower_bound = first_szero_date - weeks(2),
              iceform_upper_bound = first_szero_date + weeks(2))
  
  # Q: What percentage of lake*station*years have a date with rolltemp <= 0 within
  # six weeks of the start of their sampling season?
  (tavg_rollmean %>%
      filter(season == "iceon", roll_mean <= 0,
             tavg_date >= (as_date(start_date_min) - weeks(6)),
             tavg_date <= as_date(start_date_min)) %>%
      select(lakename, poolstation, euli_year) %>%
      unique() %>%
      count(lakename) %>%
      summarize(sum(n))) / 
    (euli_filtered %>%
       filter(season == "iceon") %>%
       select(lakename, poolstation, euli_year) %>%
       unique() %>%
       count(lakename) %>%
       summarize(sum(n)))
  
  # Q: What percentage of lake*station*years have roll_mean <= 0 on day 1 of their
  # search interval?
  (tavg_rollmean %>%
      filter(season == "iceon",
             roll_mean <= 0,
             tavg_date == (as_date(start_date_min) - weeks(6))) %>%
      select(lakename, poolstation, euli_year) %>%
      unique() %>%
      nrow() / 
      (tavg_rollmean %>%
         filter(season == "iceon") %>%
         select(lakename, poolstation, euli_year) %>%
         unique() %>%
         nrow()))
  
  # Export the iceform periods for later use
  export_iceform_periods <- iceform_period_rollmean %>%
    select(lakename, poolstation, euli_year, season, start_date_min,
           end_date_max, iceform_lower_bound, iceform_upper_bound, first_szero_date) %>%
    unique()
  
  # Make a histogram of dates for the first rolling average of 0C
  iceform_hist_path <- "figures/first_day_zero_c_iceform.png"
  
  iceform_hist <- ggplot(export_iceform_periods) +
    geom_histogram(aes(yday(first_szero_date)), color = "black", fill = "white") +
    ggtitle("First day of 0C rolling average, 1.5 month buffer")
  
  ggsave(filename = iceform_hist_path, plot = iceform_hist, device = "png",
         width = 6, height = 6, units = "in")
  
  return(list(
    tavg_rollmean = tavg_rollmean,
    iceform_period_rollmean = export_iceform_periods,
    iceform_hist_out = iceform_hist_path
  ))
  
  
}


# Acquire dates around pre-stratification
determine_prestrat <- function(tavg_rollmean){
  
  # Filter to remove data from outside a one-month window leading up to the
  # start_date_min. Basically, look only at data in the 29 days before the earliest
  # recorded sampling for the season (and including that date)
  prestrat_period_rollmean <- tavg_rollmean %>%
    filter(tavg_date >= (as_date(start_date_min) - days(28)),
           tavg_date <= as_date(start_date_min),
           season == "iceoff") %>%
    group_by(lakename, poolstation, euli_year, season, start_date_min, end_date_max) %>%
    arrange(tavg_date) %>%
    mutate(prestrat_lower_bound = (as_date(start_date_min) - days(28)),
           prestrat_upper_bound = as_date(start_date_min))
  
  export_prestrat_periods <- prestrat_period_rollmean %>%
    select(lakename, poolstation, euli_year, season, prestrat_lower_bound,
           start_date_min, prestrat_upper_bound, end_date_max) %>%
    unique()
  
  return(export_prestrat_periods)
  
}


# 1.2 Climate variable wrangling ------------------------------------------

# Get BEST temperatures averages
calc_best_tavg <- function(tavg_full, iceform_period_rollmean, prestrat_period_rollmean){
  
  # BEST during iceform
  non_cru_iceform_int_temp <- inner_join(x = tavg_full %>%
                                           mutate(across(.cols = contains("date"),
                                                         .fns = as_date)),
                                         y = iceform_period_rollmean,
                                         by = c("lakename", "poolstation",
                                                "euli_year", "season",
                                                "start_date_min",
                                                "end_date_max")) %>%
    filter(as_date(tavg_date) %within% interval(as_date(iceform_lower_bound),
                                                as_date(iceform_upper_bound))) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(iceform_lower_bound = unique(iceform_lower_bound),
              iceform_upper_bound = unique(iceform_upper_bound),
              iceform_mean_tavg = mean(tavg))
  
  
  
  # BEST during underice sampling
  non_cru_underice_temp <- inner_join(x = tavg_full %>%
                                        mutate(across(.cols = contains("date"),
                                                      .fns = as_date)),
                                      y = iceform_period_rollmean,
                                      by = c("lakename", "poolstation",
                                             "euli_year", "season",
                                             "start_date_min",
                                             "end_date_max")) %>%
    filter(as_date(tavg_date) %within% interval(as_date(start_date_min),
                                                as_date(end_date_max))) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(start_date_min = unique(start_date_min),
              end_date_max = unique(end_date_max),
              underice_mean_tavg = mean(tavg))
  
  
  
  
  # BEST during prestrat
  non_cru_prestrat_int_temp <- inner_join(x = tavg_full %>%
                                            mutate(across(.cols = contains("date"),
                                                          .fns = as_date)),
                                          y = prestrat_period_rollmean,
                                          by = c("lakename", "poolstation",
                                                 "euli_year", "season",
                                                 "start_date_min",
                                                 "end_date_max")) %>%
    filter(as_date(tavg_date) %within% interval(as_date(prestrat_lower_bound),
                                                as_date(prestrat_upper_bound))) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(prestrat_lower_bound = unique(prestrat_lower_bound),
              prestrat_upper_bound = unique(prestrat_upper_bound),
              prestrat_mean_tavg = mean(tavg))
  
  
  
  # BEST during stratification
  non_cru_strat_temp <- inner_join(x = tavg_full %>%
                                     mutate(across(.cols = contains("date"),
                                                   .fns = as_date)),
                                   y = prestrat_period_rollmean,
                                   by = c("lakename", "poolstation",
                                          "euli_year", "season",
                                          "start_date_min",
                                          "end_date_max")) %>%
    filter(as_date(tavg_date) %within% interval(as_date(start_date_min),
                                                as_date(end_date_max))) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(start_date_min = unique(start_date_min),
              end_date_max = unique(end_date_max),
              strat_mean_tavg = mean(tavg))
  
  
  return(list(
    best_iceform_period_tavg = non_cru_iceform_int_temp,
    best_underice_tavg = non_cru_underice_temp,
    best_prestrat_period_tavg = non_cru_prestrat_int_temp,
    best_strat_tavg = non_cru_strat_temp
  ))
  
  
}


# Extract surface pressure data
pull_agg_slp <- function(){
  
  # Storage path for downloaded files
  slp_path <- "data/inputs/slp_raw"
  
  # Pull SLP data for years of interest and save as rds files by year
  future_walk(.x = 1948:2016,
              .f = ~ {
                temporary_slp <- NCEP.gather(
                  # Pressure near surface
                  variable = "pres.sfc",
                  level = "surface",
                  years.minmax = c(.x, .x),
                  months.minmax = c(1, 12),
                  lon.westeast = c(-150, 0),
                  lat.southnorth = c(36, 70)) %>%
                  NCEP.array2df()
                
                # Save output for later
                write_rds(x = temporary_slp,
                          file = paste0(slp_path,
                                        "/pres_sfc_",
                                        .x,
                                        ".rds"))
              })
  
  # Iterate over yearly .rds files and combine into single object
  compiled_slp <- future_map(.x = list.files(path = slp_path,
                                             pattern = "sfc_[0-9]{4}.rds",
                                             full.names = TRUE),
                             .f = ~ read_rds(file = .x) %>%
                               rename(pres = variable1) %>%
                               mutate(year = substr(x = datetime, start = 1, stop = 4),
                                      month = substr(x = datetime, start = 6, stop = 7),
                                      day = substr(x = datetime, start = 9, stop = 10)) %>%
                               group_by(latitude, longitude, year, month, day) %>%
                               summarize(mean_slp = mean(pres)) %>%
                               ungroup()
  )
  
  # Export the combined object
  future_walk(.x = compiled_slp,
              .f = ~ write_rds(x = .x,
                               file = paste0(slp_path,
                                             "/pres_sfc_agg_",
                                             unique(.x$year),
                                             ".rds")))
  
  return(slp_out_path = slp_path)
  
}

# Extracts SLP values for all locations, but does not subset by time period
extract_slp <- function(slp_file_location, location_data){
  
  # List of yearly SLP files
  slp_agg_list <- list.files(path = slp_file_location,
                             pattern = "_agg_",
                             full.names = TRUE)
  
  # Extract based on location only, not time
  slp_extracted <- future_map_dfr(.x = slp_agg_list,
                                  .f = ~
                                    
                                    {
                                      # Record the year for the SLP data file
                                      temporary_year <- str_extract(string = .x,
                                                                    pattern = "[0-9]{4}\\.rds") %>%
                                        str_remove(pattern = ".rds") %>%
                                        as.numeric()
                                      
                                      # Read in the file
                                      temporary_agg_slp <- read_rds(.x)
                                      
                                      # Create an sf object from it
                                      slp_sf <- st_as_sf(x = temporary_agg_slp,
                                                         coords = c("longitude", "latitude"))
                                      
                                      # Get unique station locations
                                      unique_station_locs <- location_data %>%
                                        st_as_sf(coords = c("long", "lat")) %>%
                                        select(lakename, poolstation) %>%
                                        distinct()
                                      
                                      # Get unique SLP locations
                                      unique_slp_locs <- slp_sf %>%
                                        select(geometry) %>%
                                        unique.data.frame()
                                      
                                      # A data frame with id/name of poolstation
                                      # or GHCN stations matched to location
                                      # ("geometry.1") of the nearest SLP data
                                      # Source: https://stackoverflow.com/questions/59766153/left-join-based-on-closest-lat-lon-in-r
                                      station_slp_match <- unique_station_locs %>% 
                                        cbind(unique_slp_locs[st_nearest_feature(unique_station_locs, unique_slp_locs), ]) %>%
                                        # The coords of the corresponding SLP data
                                        rename(slp_geometry = geometry.1) %>%
                                        mutate(lakename_poolstation = paste(lakename, poolstation, sep = "_"))
                                      
                                      # Switch the geometry that is "active" so
                                      # it's the one for SLP
                                      st_geometry(station_slp_match) <- "slp_geometry"
                                      
                                      # Now pull the SLP data that is closest to
                                      # each station location
                                      slp_by_station <- map_df(.x = unique(station_slp_match$lakename_poolstation),
                                                               .f = ~ st_filter(x = slp_sf,
                                                                                y = station_slp_match %>%
                                                                                  filter(lakename_poolstation == .x),
                                                                                .predicate = st_equals) %>%
                                                                 mutate(lakename_poolstation = .x)) %>%
                                        separate(col = "lakename_poolstation",
                                                 into = c("lakename", "poolstation"),
                                                 sep = "_",
                                                 remove = TRUE,
                                                 extra = "merge")
                                      
                                    }
                                  
  )
  
  return(slp_extracted)
  
}


# Filter SLP for dates relevant to winter
filter_slp_winter <- function(extracted_slp, euli_slim,
                              iceform_period_rollmean){
  
  extracted_slp <- extracted_slp %>%
    mutate(lat = st_coordinates(.)[, "Y"],
           long = st_coordinates(.)[, "X"]) %>%
    st_set_geometry(NULL) %>%
    # Other parts of the workflow treat poolstations with NA names as actual NAs
    # Needs to match in formatting in order to join properly
    mutate(poolstation = if_else(condition = poolstation == "NA",
                                 true = NA_character_,
                                 false = poolstation))
  
  winter_slp_lake_centrs <- left_join(x = iceform_period_rollmean,
                                      y = euli_slim,
                                      by = c("lakename", "poolstation",
                                             "euli_year", "season",
                                             "start_date_min", "end_date_max")) %>%
    ungroup() %>%
    select(lakename,
           poolstation,
           season,
           lat,
           long,
           year = euli_year,
           contains("bound"),
           contains("_date_")) %>%
    unique()
  
  # Ice-formation period
  walk(.x = unique(winter_slp_lake_centrs$lakename),
       .f = ~ {
         
         slp_df <- inner_join(x = winter_slp_lake_centrs %>%
                                filter(lakename == .x),
                              y = extracted_slp %>%
                                filter(lakename == .x),
                              by = c("lakename", "poolstation"),
                              suffix = c("_season", "_slp")) %>%
           mutate(slp_date = ymd(paste(year_slp, month, day, sep = "-"))) %>%
           select(-lat_slp, -long_slp) %>%
           group_by(lakename, poolstation, season, lat = lat_season, long = long_season,
                    year_season) %>%
           filter(slp_date %within% interval(as_date(iceform_lower_bound),
                                             as_date(iceform_upper_bound ))) %>%
           ungroup()
         
         write_rds(x = slp_df,
                   file = paste0("data/inputs/slp_intermediate/",
                                 make_clean_names(.x),
                                 "_iceform_slp.rds"))
       })
  
  filtered_iceform_df <- map_df(.x = list.files(path = "data/inputs/slp_intermediate/",
                                                full.names = TRUE,
                                                pattern = "iceform"),
                                .f = read_rds)
  
  # Under-ice period
  walk(.x = unique(winter_slp_lake_centrs$lakename),
       .f = ~ {
         
         slp_df <- inner_join(x = winter_slp_lake_centrs %>%
                                filter(lakename == .x),
                              y = extracted_slp %>%
                                filter(lakename == .x),
                              by = c("lakename", "poolstation"),
                              suffix = c("_season", "_slp")) %>%
           mutate(slp_date = ymd(paste(year_slp, month, day, sep = "-")),
                  across(.cols = contains("_date_"), .fns = ~as_date(.)),
                  sample_range = interval(start = start_date_min,
                                          end = end_date_max)) %>%
           select(-lat_slp, -long_slp) %>%
           rename(lat = lat_season, long = long_season) %>%
           filter(slp_date %within% sample_range) %>%
           ungroup()
         
         write_rds(x = slp_df,
                   file = paste0("data/inputs/slp_intermediate/",
                                 make_clean_names(.x),
                                 "_underice_slp.rds"))
       })
  
  
  filtered_sample_df <- map_df(.x = list.files(path = "data/inputs/slp_intermediate/",
                                               full.names = TRUE,
                                               pattern = "underice"),
                               .f = read_rds)
  
  
  return(list(df_period = filtered_iceform_df,
              df_sample = filtered_sample_df,
              winter_slp_lake_centrs = winter_slp_lake_centrs %>%
                rename(Y = lat,
                       X = long)))
  
}

# Filter SLP for dates relevant to summer
filter_slp_summer <- function(extracted_slp, euli_slim, prestrat_period_rollmean){
  
  extracted_slp <- extracted_slp %>%
    mutate(lat = st_coordinates(.)[, "Y"],
           long = st_coordinates(.)[, "X"]) %>%
    st_set_geometry(NULL) %>%
    # Other parts of the workflow treat poolstations with NA names as actual NAs
    # Needs to match in formatting in order to join properly
    mutate(poolstation = if_else(condition = poolstation == "NA",
                                 true = NA_character_,
                                 false = poolstation))
  
  summer_slp_lake_centrs <- left_join(x = prestrat_period_rollmean,
                                      y = euli_slim,
                                      by = c("lakename", "poolstation",
                                             "euli_year", "season",
                                             "start_date_min", "end_date_max")) %>%
    ungroup() %>%
    select(lakename,
           poolstation,
           season,
           lat,
           long,
           year = euli_year,
           contains("bound"),
           contains("_date_")) %>%
    unique()
  
  # Pre-stratification period
  walk(.x = unique(summer_slp_lake_centrs$lakename),
       .f = ~ {
         
         slp_df <- inner_join(x = summer_slp_lake_centrs %>%
                                filter(lakename == .x),
                              y = extracted_slp %>%
                                filter(lakename == .x),
                              by = c("lakename", "poolstation"),
                              suffix = c("_season", "_slp")) %>%
           mutate(slp_date = ymd(paste(year_slp, month, day, sep = "-"))) %>%
           select(-lat_slp, -long_slp) %>%
           group_by(lakename, poolstation, season, lat = lat_season, long = long_season,
                    year_season) %>%
           filter(slp_date %within% interval(as_date(prestrat_lower_bound),
                                             as_date(prestrat_upper_bound))) %>%
           ungroup()
         
         write_rds(x = slp_df,
                   file = paste0("data/inputs/slp_intermediate/",
                                 make_clean_names(.x),
                                 "_prestrat_slp.rds"))
       })
  
  filtered_prestrat_df <- map_df(.x = list.files(path = "data/inputs/slp_intermediate/",
                                                 full.names = TRUE,
                                                 pattern = "prestrat"),
                                 .f = read_rds)
  
  
  # Stratification period
  walk(.x = unique(summer_slp_lake_centrs$lakename),
       .f = ~ {
         
         slp_df <- inner_join(x = summer_slp_lake_centrs %>%
                                filter(lakename == .x),
                              y = extracted_slp %>%
                                filter(lakename == .x),
                              by = c("lakename", "poolstation"),
                              suffix = c("_season", "_slp")) %>%
           mutate(slp_date = ymd(paste(year_slp, month, day, sep = "-")),
                  across(.cols = contains("_date_"), .fns = ~as_date(.)),
                  sample_range = interval(start = start_date_min,
                                          end = end_date_max)) %>%
           select(-lat_slp, -long_slp) %>%
           rename(lat = lat_season, long = long_season) %>%
           filter(slp_date %within% sample_range) %>%
           ungroup()
         
         write_rds(x = slp_df,
                   file = paste0("data/inputs/slp_intermediate/",
                                 make_clean_names(.x),
                                 "_summersample_slp.rds"))
       })
  
  
  filtered_sample_df <- map_df(.x = list.files(path = "data/inputs/slp_intermediate/",
                                               full.names = TRUE,
                                               pattern = "summersample"),
                               .f = read_rds)
  
  
  return(list(df_period = filtered_prestrat_df,
              df_sample = filtered_sample_df,
              summer_slp_lake_centrs = summer_slp_lake_centrs %>%
                rename(Y = lat,
                       X = long)))
  
}

# Clean up and process the outputs of the filtered SLP targets
process_slp <- function(daily_slp_iceform, daily_slp_underice,
                        daily_slp_prestrat, daily_slp_strat){
  
  # Summarize the daily data by time period and export
  slp_iceform_period <- daily_slp_iceform %>%
    rename(euli_year = year_season,
           pres = mean_slp) %>%
    filter(slp_date %within% interval(iceform_lower_bound,
                                      iceform_upper_bound)) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(mean_slp = mean(pres),
              var_slp = var(pres))
  
  slp_underice_period <- daily_slp_underice %>%
    rename(euli_year = year_season,
           pres = mean_slp) %>%
    filter(slp_date %within% interval(start_date_min,
                                      end_date_max)) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(mean_slp = mean(pres),
              var_slp = var(pres))
  
  slp_prestrat_period <- daily_slp_prestrat %>%
    rename(euli_year = year_season,
           pres = mean_slp) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    filter(slp_date %within% interval(prestrat_lower_bound,
                                      prestrat_upper_bound)) %>%
    summarize(mean_slp = mean(pres),
              var_slp = var(pres))
  
  slp_strat_period <- daily_slp_strat %>%
    filter(slp_date %within% interval(start_date_min,
                                      end_date_max)) %>%
    rename(euli_year = year_season,
           pres = mean_slp) %>%
    group_by(lakename, poolstation, euli_year, season) %>%
    summarize(mean_slp = mean(pres),
              var_slp = var(pres))
  
  return(list(
    slp_iceform_period = slp_iceform_period,
    slp_underice_period = slp_underice_period,
    slp_prestrat_period = slp_prestrat_period,
    slp_strat_period = slp_strat_period
  ))
  
}

# Retrieve CRU data from .nc files:

# Air temperature
extract_cru_airtemp <- function(cru_airtemp_path,
                                winter_slp_lake_centrs,
                                summer_slp_lake_centrs){
  
  
  winter_air_temp_monthly <- get_cru(latlong = winter_slp_lake_centrs,
                                     variable_name = "tmp",
                                     file_name = cru_airtemp_path)
  
  summer_air_temp_monthly <- get_cru(latlong = summer_slp_lake_centrs,
                                     variable_name = "tmp",
                                     file_name = cru_airtemp_path)
  
  air_temp_monthly_iceform <- winter_air_temp_monthly[["df_period"]] %>%
    select(lakename, poolstation, euli_year, season, 
           iceform_lower_bound = lower_bound,
           iceform_upper_bound = upper_bound,
           tmp_date = clim_date, tmp) %>%
    unique()
  
  air_temp_monthly_underice <- winter_air_temp_monthly[["df_sample"]] %>%
    select(lakename, poolstation, euli_year, season, 
           contains("_date_"),
           tmp_date = clim_date, tmp) %>%
    unique()
  
  air_temp_monthly_prestrat <- summer_air_temp_monthly[["df_period"]] %>%
    select(lakename, poolstation, euli_year, season, 
           prestrat_lower_bound = lower_bound,
           prestrat_upper_bound = upper_bound,
           tmp_date = clim_date, tmp) %>%
    unique()
  
  air_temp_monthly_strat <- summer_air_temp_monthly[["df_sample"]] %>%
    select(lakename, poolstation, euli_year, season, 
           contains("_date_"),
           tmp_date = clim_date, tmp) %>%
    unique()
  
  return(list(
    air_temp_monthly_iceform = air_temp_monthly_iceform,
    air_temp_monthly_underice = air_temp_monthly_underice,
    air_temp_monthly_prestrat = air_temp_monthly_prestrat,
    air_temp_monthly_strat = air_temp_monthly_strat
  ))
  
}

# Precipitation
extract_cru_precip <- function(cru_precip_path,
                               winter_slp_lake_centrs,
                               summer_slp_lake_centrs){
  
  winter_precip_monthly <- get_cru(latlong = winter_slp_lake_centrs,
                                   variable_name = "pre",
                                   file_name = cru_precip_path)
  
  summer_precip_monthly <- get_cru(latlong = summer_slp_lake_centrs,
                                   variable_name = "pre",
                                   file_name = cru_precip_path)
  
  precip_monthly_iceform <- winter_precip_monthly[["df_period"]] %>%
    select(lakename, poolstation, season, euli_year,
           iceform_lower_bound = lower_bound,
           iceform_upper_bound = upper_bound,
           precip_date = clim_date, precip = pre) %>%
    unique()
  
  precip_monthly_underice <- winter_precip_monthly[["df_sample"]] %>%
    select(lakename, poolstation, season, euli_year,
           contains("_date_"),
           precip_date = clim_date, precip = pre) %>%
    unique()
  
  precip_monthly_prestrat <- summer_precip_monthly[["df_period"]] %>%
    select(lakename, poolstation, season, euli_year,
           prestrat_lower_bound = lower_bound,
           prestrat_upper_bound = upper_bound,
           precip_date = clim_date, precip = pre) %>%
    unique()
  
  precip_monthly_strat <- summer_precip_monthly[["df_sample"]] %>%
    select(lakename, poolstation, season, euli_year,
           contains("_date_"),
           precip_date = clim_date, precip = pre) %>%
    unique()
  
  return(list(
    precip_monthly_iceform = precip_monthly_iceform,
    precip_monthly_underice = precip_monthly_underice,
    precip_monthly_prestrat = precip_monthly_prestrat,
    precip_monthly_strat = precip_monthly_strat
  ))
  
}


# 1.3 Building dataset ----------------------------------------------------

# Winter
build_winter <- function(precip_iceform, tmp_iceform,
                         precip_underice, tmp_underice,
                         tavg_iceform, tavg_underice, slp_iceform,
                         slp_underice, euli_filtered){
  
  slp_iceform <- slp_iceform %>%
    rename(iceform_mean_slp = mean_slp,
           iceform_var_slp = var_slp)
  
  slp_underice <- slp_underice %>%
    rename(underice_mean_slp = mean_slp,
           underice_var_slp = var_slp)
  
  # Join together the datasets from CRU that are aggregated to month level:
  cru_iceform_data <- map2(.x = list(precip_iceform, tmp_iceform),
                           .y = list("precip", "tmp"),
                           .f = ~ .x %>%
                             group_by(lakename, poolstation, euli_year, season) %>%
                             # Take an average of the months that overlap with the iceform interval
                             summarize(!!paste0("iceform_mean_", .y) := mean(.data[[.y]], na.rm = TRUE))) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year", "season"))
  
  # Data for the period during which under-ice sampling took place:
  cru_underice_data <- map2(.x = list(precip_underice, tmp_underice),
                            .y = list("precip", "tmp"),
                            .f = ~ .x %>%
                              group_by(lakename, poolstation, euli_year, season,
                                       start_date_min, end_date_max) %>%
                              # Take an average of the months that overlap with the underice interval
                              summarize(!!paste0("underice_mean_", .y) := mean(.data[[.y]], na.rm = TRUE))) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year",
                  "season", "start_date_min", "end_date_max"))
  
  # Combine all CRU data now
  all_winter_cru <- full_join(x = cru_iceform_data, y = cru_underice_data,
                              by = c("lakename", "poolstation", "euli_year", "season"))
  
  # CRU + BEST + SLP
  all_winter_climate_data <- list(all_winter_cru, tavg_iceform, tavg_underice,
                                  slp_iceform, slp_underice) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year", "season"))
  
  # Do the start / end dates all match? If so, will remove duplicate cols
  all_winter_climate_data %>%
    select(lakename, poolstation, euli_year, contains("start"), contains("end"),
           contains("upper"), contains("lower")) %>%
    mutate(same_start = if_else(
      condition = start_date_min.x == start_date_min.y,
      true = TRUE, false = FALSE),
      same_end = if_else(
        condition = end_date_max.x == end_date_max.y,
        true = TRUE, false = FALSE)) %>%
    ungroup() %>%
    count(same_start, same_end)
  
  all_winter_climate_data <- all_winter_climate_data %>%
    rename(start_date_min = start_date_min.x,
           end_date_max = end_date_max.x) %>%
    select(-c(contains(".y"), contains(".x")))
  
  # Join climate to EULI
  winter_climate_euli <- inner_join(x = euli_filtered %>%
                                      select(lakename, poolstation, lakecountry,
                                             stationlat, stationlong, euli_year, season,
                                             start_date_min, end_date_max, avechla,
                                             avetotphos, avetotdissphos, avetotnitro,
                                             avetotdissnitro, watertemp),
                                    y = all_winter_climate_data,
                                    by = c("lakename", "poolstation", "euli_year",
                                           "season", "start_date_min", "end_date_max")) %>%
    mutate(continent = case_when(
      lakecountry %in% c("USA", "Greenland", "Canada", "USA/Canada", "Canada/USA") ~ "North America",
      lakecountry %in% c("Spain", "Russia", "Sweden", "Germany", "Finland", "Estonia", "Italy") ~ "Eurasia",
      TRUE ~ NA_character_))
  
  
  # Take a look at variable relationships as a brief QC measure
  winter_climate_cor_path <- "figures/winter_build_correlations.png"
  
  png(filename = winter_climate_cor_path, width = 5, height = 5, units = "in", res = 300)
  
  winter_climate_euli %>%
    select_if(.predicate = is.numeric) %>%
    cor(use = "pairwise.complete.obs") %>%
    corrplot(type = "upper")
  
  dev.off()
  
  return(list(
    winter_climate_cor_path = winter_climate_cor_path,
    winter_climate_euli = winter_climate_euli))
  
}

# Summer
build_summer <- function(precip_prestrat, tmp_prestrat,
                         precip_strat, tmp_strat,
                         tavg_prestrat, tavg_strat, slp_prestrat,
                         slp_strat, euli_filtered){
  
  slp_prestrat <- slp_prestrat %>%
    rename(prestrat_mean_slp = mean_slp,
           prestrat_var_slp = var_slp)
  
  slp_strat <- slp_strat %>%
    rename(strat_mean_slp = mean_slp,
           strat_var_slp = var_slp)
  
  # Join together the datasets from CRU that are aggregated to month level
  monthly_prestrat_data <- map2(.x = list(precip_prestrat, tmp_prestrat),
                                .y = list("precip", "tmp"),
                                .f = ~ .x %>%
                                  group_by(lakename, poolstation, euli_year, season) %>%
                                  # Take an average of the months that overlap with the prestrat interval
                                  summarize(!!paste0("prestrat_mean_", .y) := mean(.data[[.y]], na.rm = TRUE))) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year", "season"))
  
  # Data for the period during which under ice sampling took place
  monthly_strat_data <- map2(.x = list(precip_strat, tmp_strat),
                             .y = list("precip", "tmp"),
                             .f = ~ .x %>%
                               group_by(lakename, poolstation, euli_year, season,
                                        start_date_min, end_date_max) %>%
                               # Take an average of the months that overlap with the strat interval
                               summarize(!!paste0("strat_mean_", .y) := mean(.data[[.y]], na.rm = TRUE))) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year",
                  "season", "start_date_min", "end_date_max"))
  
  
  # Combine all CRU data now
  all_summer_cru <- full_join(x = monthly_prestrat_data, y = monthly_strat_data,
                              by = c("lakename", "poolstation", "euli_year", "season"))
  
  # CRU + BEST + SLP
  all_summer_climate_data <- list(all_summer_cru, tavg_prestrat, tavg_strat,
                                  slp_prestrat, slp_strat) %>%
    reduce(.f = inner_join,
           by = c("lakename", "poolstation", "euli_year", "season"))
  
  # Do the start / end dates all match? If so, will remove duplicate cols
  all_summer_climate_data %>%
    select(lakename, poolstation, euli_year, contains("start"), contains("end"),
           contains("upper"), contains("lower")) %>%
    mutate(same_start = if_else(
      condition = start_date_min.x == start_date_min.y,
      true = TRUE, false = FALSE),
      same_end = if_else(
        condition = end_date_max.x == end_date_max.y,
        true = TRUE, false = FALSE)) %>%
    ungroup() %>%
    summarize(mismatched_start = sum(!same_start),
              mismatched_end = sum(!same_end))
  
  all_summer_climate_data <- all_summer_climate_data %>%
    rename(start_date_min = start_date_min.x,
           end_date_max = end_date_max.x) %>%
    select(-c(contains(".y"), contains(".x")))
  
  # Join climate to EULI
  summer_climate_euli <- inner_join(x = euli_filtered %>%
                                      select(lakename, poolstation, lakecountry,
                                             stationlat, stationlong, euli_year,
                                             season, start_date_min, end_date_max,
                                             avechla, avetotphos, avetotdissphos,
                                             avetotnitro, avetotdissnitro, watertemp),
                                    y = all_summer_climate_data,
                                    by = c("lakename", "poolstation", "euli_year",
                                           "season", "start_date_min", "end_date_max")) %>%
    mutate(continent = case_when(
      lakecountry %in% c("USA", "Greenland", "Canada", "USA/Canada", "Canada/USA") ~ "North America",
      lakecountry %in% c("Spain", "Russia", "Sweden", "Germany", "Finland", "Estonia", "Italy") ~ "Eurasia",
      TRUE ~ NA_character_))
  
  # Take a look at variable relationships as a brief QC measure
  summer_climate_cor_path <- "figures/summer_build_correlations.png"
  
  png(filename = summer_climate_cor_path, width = 5, height = 5, units = "in", res = 300)
  
  summer_climate_euli %>%
    select_if(.predicate = is.numeric) %>%
    cor(use = "pairwise.complete.obs") %>%
    corrplot(type = "upper")
  
  dev.off()
  
  return(list(
    summer_climate_cor_path = summer_climate_cor_path,
    summer_climate_euli = summer_climate_euli))
  
}

# Quick review of Lake Erie data for QC
viz_single_lake <- function(winter_climate_euli, summer_climate_euli){
  
  ggplot() +
    geom_segment(data = winter_climate_euli %>%
                   filter(lakename == "Lake Erie"),
                 aes(x = start_date_min, xend = end_date_max,
                     y = underice_mean_tavg, yend = underice_mean_tavg,
                     color = "Under-ice Mean BEST"),
                 size = 3.5) +
    geom_segment(data = summer_climate_euli %>%
                   filter(lakename == "Lake Erie"),
                 aes(x = start_date_min, xend = end_date_max,
                     y = strat_mean_tavg, yend = strat_mean_tavg,
                     color = "Stratification Mean BEST"),
                 size = 3.5) +
    geom_segment(data = winter_climate_euli %>%
                   filter(lakename == "Lake Erie"),
                 aes(x = iceform_lower_bound, xend = iceform_upper_bound,
                     y = iceform_mean_tavg, yend = iceform_mean_tavg,
                     color = "Ice-formation Mean BEST"),
                 size = 3.5) +
    geom_segment(data = summer_climate_euli %>%
                   filter(lakename == "Lake Erie"),
                 aes(x = prestrat_lower_bound, xend = prestrat_upper_bound,
                     y = prestrat_mean_tavg, yend = prestrat_mean_tavg,
                     color = "Pre-stratification Mean BEST"),
                 size = 3.5) +
    xlab("Temperature period") +
    ylab("Mean temp. C (B.E.S.T.)") +
    ggtitle("Lake Erie time series") +
    scale_color_manual(name = "Temperature measurement period",
                       values = c("tomato", "goldenrod1",
                                  "darkorchid1", "turquoise")) +
    theme_bw()
  
  erie_out_path <- "figures/lake_erie_timelines.png"
  
  ggsave(filename = erie_out_path, device = "png",
         width = 8, height = 6, units = "in")
  
  return(erie_out_path)
  
}

# Export the modeling datasets as csv files
export_model_datasets <- function(winter_climate_euli, summer_climate_euli){
  
  # Export files for use in modeling script (the versions from the modeling
  # script are what are shared, too)
  winter_data_out <- "data/outputs/stormy_winter_model_data.csv"
  
  winter_model_data <- winter_climate_euli %>%
    mutate(avechla_log10 = log10(avechla + 0.1),
           TN_log10 = log10(avetotnitro + 0.1),
           TP_log10 = log10(avetotphos + 0.1))
  
  write_csv(x = winter_model_data,
            file = winter_data_out)
  
  summer_data_out <- "data/outputs/stormy_summer_model_data.csv" 
  
  summer_model_data <- summer_climate_euli %>%
    mutate(avechla_log10 = log10(avechla + 0.1),
           TN_log10 = log10(avetotnitro + 0.1),
           TP_log10 = log10(avetotphos + 0.1))
  
  write_csv(x = summer_model_data,
            file = summer_data_out)
  
  return(list(
    winter_data_out_path = winter_data_out,
    summer_data_out_path = summer_data_out,
    winter_model_data = winter_model_data,
    summer_model_data = summer_model_data
  ))
  
}


# 2. Modeling -------------------------------------------------------------


# 2.1 Run and plot RF models ----------------------------------------------

# Run winter season random forest model
run_and_plot_winter_model <- function(winter_model_data, winter_explanatory,
                                      kfold_control){
  
  # Correlation plot: general winter
  winter_model_cor_path <- "figures/winter_model_correlations.png"
  
  png(filename = winter_model_cor_path,
      width = 10, height = 10, units = "in", res = 150)
  
  winter_model_data %>%
    select(where(is.numeric)) %>%
    cor(use = "pairwise.complete.obs") %>%
    # https://stackoverflow.com/questions/40509217/how-to-have-r-corrplot-title-position-correct
    corrplot(type = "upper", title = "winter data", mar = c(0, 0, 1, 0),
             method = "color", addCoef.col = "grey39",
             order = "alphabet", number.cex = .73)
  
  dev.off()
  
  
  # The model that we'll fit:
  winter_model_formula <- reformulate(termlabels = winter_explanatory,
                                      response = "avechla_log10")
  
  # Model training:
  # (Setting seed using a randomly generated integer)
  set.seed(3289)
  winter_rf <- train(form = winter_model_formula, 
                     # Remove rows that are NA for any of the predictors or
                     # the response. Don't filter on watertemp
                     data = winter_model_data %>%
                       filter(across(.cols = any_of(c("avechla_log10",
                                                      winter_explanatory)),
                                     .fns = ~!is.na(.))),
                     # Random forest
                     method = "rf",
                     # importance = TRUE,
                     trControl = kfold_control)
  
  
  plot(winter_rf)
  
  winter_rf$finalModel
  
  winter_rf_imp_df <- varImp(winter_rf) %>%
    pluck("importance") %>% 
    rownames_to_column() %>% 
    rename("variable" = rowname) %>% 
    arrange(Overall) %>%
    mutate(variable = forcats::fct_inorder(variable)) %>%
    rename(imp = Overall)
  
  winter_imp_plot <- ggplot(winter_rf_imp_df) +
    geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
                 size = 1, alpha = 0.7) +
    geom_point(aes(x = variable, y = imp), 
               size = 2, show.legend = F, color = "blue") +
    coord_flip() +
    xlab("Variable") +
    ylab("Importance") +
    theme_bw() +
    theme(text = element_text(size = 22),
          axis.text.y = element_text(margin = ggplot2::margin(r = 7)),
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1))
  
  
  # Single-predictor partial dependence plots for the top 5 predictors of the RF
  
  # Define vars & labels for plots
  single_pdp_vars_winter <- winter_rf_imp_df %>%
    arrange(desc(imp)) %>%
    .[1:5, ] %>%
    pull(variable) %>%
    as.character() %>%
    tibble(vars = .) %>%
    mutate(labels = case_when(
      vars == "iceform_mean_slp" ~ "Mean SLP ice-formation",
      vars == "underice_var_slp" ~ "Var SLP under-ice",
      vars == "TP_log10" ~ "logTP under-ice",
      vars == "TN_log10" ~ "logTN under-ice",      
      vars == "underice_mean_precip" ~ "PPT (mm/month) under-ice",
      vars == "underice_mean_slp" ~ "Mean SLP under-ice",
      vars == "iceform_mean_precip" ~ "PPT (mm/month) ice-formation"
    ))
  
  set.seed(3289)
  winter_pdps <- map2(.x = single_pdp_vars_winter$vars,
                      .y = single_pdp_vars_winter$labels,
                      .f = ~ {
                        
                        # Plot SLP with scientific format
                        if(.x == "underice_mean_slp") {
                          
                          partial(object = winter_rf,
                                  pred.var = .x,
                                  # Includes marks showing the deciles of
                                  # the training dataset
                                  rug = TRUE,
                                  # Create a plot using ggplot
                                  plot = TRUE,
                                  plot.engine = "ggplot2") +
                            scale_x_continuous(
                              label = scales::label_scientific(2)
                            ) +
                            ylab(expression(hat(y))) +
                            xlab(.y) +
                            theme_bw() +
                            theme(text = element_text(size = 18),
                                  axis.text = element_text(size = 16),
                                  axis.text.x = element_text(angle = 25,
                                                             margin = ggplot2::margin(
                                                               t = 10),
                                                             hjust = 1),
                                  panel.border = element_rect(fill = NA,
                                                              colour = "black", 
                                                              size = 1))
                          
                        } else{
                          
                          # Plot other variables the same across the board
                          partial(object = winter_rf,
                                  pred.var = .x,
                                  rug = TRUE,
                                  plot = TRUE,
                                  plot.engine = "ggplot2") +
                            ylab("") +
                            xlab(.y) +
                            theme_bw() +
                            theme(text = element_text(size = 20),
                                  axis.text = element_text(size = 18),
                                  panel.border = element_rect(fill = NA,
                                                              colour = "black", 
                                                              size = 1))
                          
                        }
                      })
  
  # Add names
  names(winter_pdps) <- winter_rf_imp_df %>%
    arrange(desc(imp)) %>%
    .[1:5, ] %>%
    pull(variable) %>%
    as.character()
  
  # Then alphabetize the list
  winter_pdps <- winter_pdps[c(grep(pattern = "TP_", x = names(winter_pdps)),
                               grep(pattern = "TN_", x = names(winter_pdps)),
                               grep(pattern = "TP_|TN_", x = names(winter_pdps),
                                    invert = TRUE))]
  # Multi-predictor partial dependence plot
  
  # Define vars & labels for plot
  mult_pdp_vars_winter <- data.frame(
    row.names = c("x_axis", "y_axis"), 
    var = c("TN_log10", "TP_log10"),
    label = c("logTN under-ice", "logTP under-ice"))
  
  set.seed(3289)
  winter_pn_pdp <- autoplot(partial(object = winter_rf,
                                    pred.var = c(mult_pdp_vars_winter["x_axis", "var"],
                                                 mult_pdp_vars_winter["y_axis", "var"]),
                                    parallel = TRUE,
                                    # Generate a convex hull that outlines
                                    # the part of the plot where training
                                    # data existed
                                    chull = TRUE)) +
    xlab(mult_pdp_vars_winter["x_axis", "label"]) +
    ylab(mult_pdp_vars_winter["y_axis", "label"]) +
    scale_fill_viridis_c(guide = guide_colorbar(
      title = expression(hat(y)),
      title.position = "left",
      frame.colour = "black",
      frame.linewidth = 1.5,
      barwidth = 10)) +
    theme_bw() +
    theme(text = element_text(size = 20),
          legend.position = "bottom",
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1))
  
  # Export figure of only partial density plots
  winter_pdp_grid <- plot_grid(
    winter_pn_pdp,
    plotlist = winter_pdps,
    nrow = 2, ncol = 3)
  
  # Export figure with both importance and pdps. First update the importance
  # plot
  winter_imp_grid <- plot_grid(NULL,
                               winter_imp_plot +
                                 scale_x_discrete(
                                   labels = c("TP_log10" = "logTP under-ice",
                                              "TN_log10" = "logTN under-ice",      
                                              "underice_mean_precip" = "PPT (mm/month) under-ice",
                                              "underice_mean_slp" = "Mean SLP under-ice",
                                              "iceform_mean_precip" = "PPT (mm/month) ice-formation",
                                              "iceform_mean_tavg" = "Temp (\u00B0C) ice-formation,",
                                              "iceform_mean_slp" = "Mean SLP ice-formation",
                                              "underice_var_slp" = "VAR SLP under-ice",
                                              "iceform_var_slp" = "VAR SLP ice-formation")
                                 ),
                               NULL,
                               ncol = 3, rel_widths = c(1, 4, 1))
  
  full_winter_plot <- plot_grid(winter_imp_grid, NULL, winter_pdp_grid,
                                nrow = 3,
                                rel_heights = c(1, 0.1, 2),
                                labels = c("A", "B"),
                                label_x = 0.05,
                                label_size = 23)
  
  full_winter_pdp_path <- "figures/winter_summary_plot.png"
  
  ggsave(filename = full_winter_pdp_path,
         plot = full_winter_plot,
         width = 13.75, height = 13.5, device = "png", units = "in")
  
  
  return(list(
    winter_model_cor_path = winter_model_cor_path,
    winter_pdp_path = full_winter_pdp_path,
    winter_model = winter_rf
  ))  
  
}


# Run summer season random forest model
run_and_plot_summer_model <- function(summer_model_data, summer_explanatory,
                                      kfold_control){
  
  # Correlation plot
  summer_model_cor_path <- "figures/summer_model_correlations.png"
  
  png(filename = summer_model_cor_path,
      width = 10, height = 10, units = "in", res = 150)
  
  summer_model_data %>%
    select(where(is.numeric)) %>%
    cor(use = "pairwise.complete.obs") %>%
    # https://stackoverflow.com/questions/40509217/how-to-have-r-corrplot-title-position-correct
    corrplot(type = "upper", title = "summer data", mar = c(0, 0, 1, 0),
             method = "color", addCoef.col = "grey39",
             order = "alphabet", number.cex = .73)
  
  dev.off()
  
  
  summer_model_formula <- reformulate(termlabels = summer_explanatory,
                                      response = "avechla_log10")
  
  # Model training: this formulation has < 400 samples b/c watertemp has a lot of NAs
  set.seed(3289)
  summer_rf <- train(form = summer_model_formula, 
                     # Remove rows that are NA for any of the predictors or
                     # the response
                     data = summer_model_data %>%
                       filter(across(.cols = any_of(c("avechla_log10",
                                                      summer_explanatory)),
                                     .fns = ~!is.na(.))),
                     # Random forest
                     method = "rf",
                     trControl = kfold_control)
  
  
  plot(summer_rf)
  
  summer_rf$finalModel
  
  summer_rf_imp_df <- varImp(summer_rf) %>%
    pluck("importance") %>% 
    rownames_to_column() %>% 
    rename("variable" = rowname) %>% 
    arrange(Overall) %>%
    mutate(variable = forcats::fct_inorder(variable)) %>%
    rename(imp = Overall)
  
  
  summer_imp_plot <- ggplot(summer_rf_imp_df) +
    geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
                 size = 1, alpha = 0.7) +
    geom_point(aes(x = variable, y = imp), 
               size = 2, show.legend = F, color = "blue") +
    coord_flip() +
    xlab("Variable") +
    ylab("Importance") +
    theme_bw() +
    theme(text = element_text(size = 22),
          axis.text.y = element_text(margin = ggplot2::margin(r = 7)),
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1))
  
  # Single-predictor partial dependence plots for the top 5 predictors of the RF
  
  # Define vars & labels for plots
  single_pdp_vars_summer <- summer_rf_imp_df %>%
    arrange(desc(imp)) %>%
    .[1:5, ] %>%
    pull(variable) %>%
    as.character() %>%
    tibble(vars = .) %>%
    mutate(labels = case_when(
      vars == "TP_log10" ~ "logTP stratification",
      vars == "TN_log10" ~ "logTN stratification",      
      vars == "watertemp" ~ "Water Temp (\u00B0C) stratification",
      vars == "strat_mean_slp" ~ "Mean SLP stratification",
      vars == "strat_var_slp" ~ "Var SLP stratification",
      vars == "prestrat_mean_tavg" ~ "Temp (\u00B0C) pre-stratification" ,
      vars == "prestrat_mean_slp" ~ "Mean SLP pre-stratification"
    ))
  
  set.seed(3289)
  summer_pdps <- map2(.x = single_pdp_vars_summer$vars,
                      .y = single_pdp_vars_summer$labels,
                      .f = ~ {
                        
                        # Plot SLP with scientific format
                        if(grepl(pattern = "slp", ignore.case = TRUE, x = .x)) {
                          
                          partial(object = summer_rf,
                                  pred.var = .x,
                                  rug = TRUE,
                                  plot = TRUE,
                                  plot.engine = "ggplot2") +
                            scale_x_continuous(
                              label = scales::label_scientific(2)
                            ) +
                            ylab("") +
                            xlab(.y) +
                            theme_bw() +
                            theme(text = element_text(size = 18),
                                  axis.text = element_text(size = 16),
                                  axis.text.x = element_text(angle = 25,
                                                             margin = ggplot2::margin(
                                                               t = 10),
                                                             hjust = 1),
                                  panel.border = element_rect(fill = NA,
                                                              colour = "black", 
                                                              size = 1))
                          
                        } else{
                          # Plot other variables the same across the board
                          partial(object = summer_rf,
                                  pred.var = .x,
                                  rug = TRUE,
                                  plot = TRUE,
                                  plot.engine = "ggplot2") +
                            ylab("") +
                            xlab(.y) +
                            theme_bw() +
                            theme(text = element_text(size = 20),
                                  axis.text = element_text(size = 18),
                                  panel.border = element_rect(fill = NA,
                                                              colour = "black", 
                                                              size = 1))
                          
                        }
                      })
  
  # Add names
  names(summer_pdps) <- summer_rf_imp_df %>%
    arrange(desc(imp)) %>%
    .[1:5, ] %>%
    pull(variable) %>%
    as.character()
  
  # Then alphabetize the list
  summer_pdps <- summer_pdps[c(grep(pattern = "TP_", x = names(summer_pdps)),
                               grep(pattern = "TN_", x = names(summer_pdps)),
                               grep(pattern = "TP_|TN_", x = names(summer_pdps),
                                    invert = TRUE))]
  
  # Label third plot differently bc it'll be on the left side of the panel
  summer_pdps[[3]] <- summer_pdps[[3]] +
    ylab(expression(hat(y)))
  
  # Multi-predictor partial dependence plot
  
  # Define vars & labels for plot
  mult_pdp_vars_summer <- data.frame(
    row.names = c("x_axis", "y_axis"), 
    var = c("TN_log10", "TP_log10"),
    label = c("logTN stratification", "logTP stratification"))
  
  set.seed(3289)
  summer_pn_pdp <- autoplot(partial(object = summer_rf,
                                    pred.var = c(mult_pdp_vars_summer["x_axis", "var"],
                                                 mult_pdp_vars_summer["y_axis", "var"]),
                                    parallel = TRUE,
                                    chull = TRUE)) +
    xlab(mult_pdp_vars_summer["x_axis", "label"]) +
    ylab(mult_pdp_vars_summer["y_axis", "label"]) +
    scale_fill_viridis_c(guide = guide_colorbar(
      title = expression(hat(y)),
      title.position = "left",
      frame.colour = "black",
      frame.linewidth = 1.5,
      barwidth = 10)) +
    theme_bw() +
    theme(text = element_text(size = 20),
          legend.position = "bottom",
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1))
  
  # Export figure of only partial density plots
  summer_pdp_grid <- plot_grid(
    summer_pn_pdp,
    plotlist = summer_pdps,
    nrow = 2, ncol = 3)
  
  # Export figure with both importance and pdps. First update the importance
  # plot
  summer_imp_grid <- plot_grid(NULL,
                               summer_imp_plot +
                                 scale_x_discrete(
                                   labels = c("TP_log10" = "logTP stratification",
                                              "TN_log10" = "logTN stratification",      
                                              "strat_mean_slp" = "Mean SLP stratification",
                                              "watertemp" = "Water Temp (\u00B0C) stratification,",
                                              "prestrat_mean_slp" = "Mean SLP pre-stratification",
                                              "prestrat_mean_tavg" = "Temp (\u00B0C) pre-stratification",
                                              "strat_mean_precip" = "PPT (mm/month) stratification",
                                              "strat_var_slp" = "VAR SLP stratification",
                                              "prestrat_var_slp" = "VAR SLP pre-stratification",
                                              "prestrat_mean_precip" = "PPT (mm/month) pre-stratification"
                                   )
                                 ),
                               NULL,
                               ncol = 3, rel_widths = c(1, 4, 1))
  
  full_summer_plot <- plot_grid(summer_imp_grid, NULL, summer_pdp_grid,
                                nrow = 3,
                                rel_heights = c(1, 0.1, 2),
                                labels = c("A", "B"),
                                label_x = 0.05,
                                label_size = 23)
  
  full_summer_pdp_path <- "figures/summer_summary_plot.png"
  
  ggsave(filename = full_summer_pdp_path,
         plot = full_summer_plot,
         width = 13.75, height = 13.5, device = "png", units = "in")
  
  
  return(list(
    summer_model_cor_path = summer_model_cor_path,
    summer_pdp_path = full_summer_pdp_path,
    summer_model = summer_rf
  ))  
  
}


# Build a map figure
create_map <- function(winter_model_data, summer_model_data, euli_filtered,
                       lakes_path){
  
  # 1. Data prep ------------------------------------------------------------
  
  # Join modeling locations to geographic data
  winter_lake_locs <- left_join(
    x = winter_model_data,
    y = euli_filtered,
    by = c("lakename", "poolstation", "euli_year", "season", "lakecountry")
  ) %>%
    select(lakename, lakecountry, poolstation, euli_year, season, lat, long) %>%
    distinct()
  
  summer_lake_locs <- left_join(
    x = summer_model_data,
    y = euli_filtered,
    by = c("lakename", "poolstation", "euli_year", "season", "lakecountry")
  ) %>%
    select(lakename, lakecountry, poolstation, euli_year, season, lat, long) %>%
    distinct()
  
  # So let's take an average of the coordinates at the lake level so we can just
  # have one set of coordinates per lake for global plotting purposes. I'm going
  # to combine seasons because o/w we will possibly end up with different coord
  # pairs for each lake between the two data frames. Can always break out by season
  # later if needed
  
  unique_locs <- bind_rows(winter_lake_locs, summer_lake_locs) %>%
    ungroup() %>%
    group_by(lakename) %>%
    summarize(lakecountry = unique(lakecountry),
              mean_lat = mean(lat),
              mean_long = mean(long))
  
  # Western geography
  west_locs <- unique_locs %>%
    filter(lakecountry %in% c(
      "USA", "Greenland", "Canada", "USA/Canada", "Canada/USA"))
  
  # Eastern geography
  east_locs <- unique_locs %>%
    filter(lakecountry %in% c("Spain", "Russia", "Sweden", "Germany", "Finland",
                              "Estonia", "Italy"))
  
  # Group together nearby lakes so that there's less overplotting of points on
  # the map. Generate number labels to indicate how many lakes a point represents
  west_map_data <- west_locs %>%
    mutate(lake_groups = case_when(
      lakename %in% c("Broderick Reservoir", "Blackstrap Reservoir",
                      "Lake Diefenbaker", "Buffalo Pound Lake", "St. Denis Pond 1",
                      "St. Denis Pond 90", "St. Denis Pond S5338") ~ "n=7_A",
      lakename %in% c("West Long Lake", "Pony Lake") ~ "n=2_A",
      lakename %in% c("Lake 239", "Lake 227") ~ "n=2_B",
      lakename %in% c("Trout Bog", "Trout Lake", "Sparkling Lake",
                      "Allequash Lake", "Big Muskellunge Lake", "Crystal Bog",
                      "Crystal Lake") ~ "n=7_B",
      lakename %in% c("Lake Mendota", "Lake Monona", "Lake Wingra", "Fish Lake") ~ "n=4_A",
      lakename %in% c("Lake Huron", "Lake Erie") ~ "n=2_C",
      lakename %in% c("Upper Cascade, ALSC #020274", "Lower Cascade, ALSC #020270",
                      "Shelburne Pond","Chapel Pond, ALSC #020274") ~ "n=4_B",
      lakename %in% c("Dodge Pond", "Bride Lake", "Rogers Lake",
                      "Pattagansett Lake", "Quonnipaug Lake", "Linsley Pond") ~ "n=6_A",
      lakecountry == "Greenland" ~ "n=17_A",
      TRUE ~ lakename),
      # Labels can be repeated, but groupings can't
      lake_label = case_when(
        grepl(pattern = "_[A-Z]$", x = lake_groups) ~
          gsub(pattern = "_[A-Z]$|n=", replacement = "", x = lake_groups),
        TRUE ~ "")) %>%
    group_by(lake_groups) %>%
    summarize(group_lat = mean(mean_lat),
              group_long = mean(mean_long),
              lake_label = unique(lake_label)) %>%
    ungroup()
  
  east_map_data <- east_locs %>%
    mutate(lake_groups = case_when(
      lakename %in% c("Lake Santo Parmense",  "Lake Scuro Parmense") ~ "n=2_A",
      lakename %in% c("Lake Stechlin", "Lake Grosse Fuchskuhle (SW Basin)",
                      "Lake Grosse Fuchskuhle (NE Basin)", 
                      "Scharmuetzelsee", "Muggelsee") ~ "n=4_A",
      lakename %in% c("Lake Vortsjarv",  "Lake Peipsi") ~ "n=2_B",
      lakename %in% c("Lake Valkea-Kotinen", "Lake Vanajanselka") ~ "n=2_C",
      lakename %in% c("Lake Muddus", "Lake Vastus") ~ "n=2_D",
      lakename %in% c("Lake Kuohkima", "Lake Oiko","Lake Siilas", "Saanajarvi",
                      "Lake Kilpis", "Lake Ropi", "Lake Kivi") ~ "n=7_A",
      lakename %in% c("Cimera Lake", "Penalara Lake") ~ "n=2_E",
      lakename %in% c() ~ "",
      TRUE ~ lakename),
      # Labels can be repeated, but groupings can't
      lake_label = case_when(
        grepl(pattern = "_[A-Z]$", x = lake_groups) ~
          gsub(pattern = "_[A-Z]$|n=", replacement = "", x = lake_groups),
        TRUE ~ "")) %>%
    group_by(lake_groups) %>%
    summarize(group_lat = mean(mean_lat),
              group_long = mean(mean_long),
              lake_label = unique(lake_label)) %>%
    ungroup()
  
  # Read in and prepare spatial data for lake shapes on map
  lakes.ne <- sf::read_sf(lakes_path)
  
  
  # 2. Map setup ------------------------------------------------------------
  
  # Base map of country borders 
  mapWorld <- borders("world", colour = "gray55", fill = "wheat1", alpha = 0.7)
  
  main_map <- ggplot() +
    mapWorld +
    ylab("") +
    xlab("")
  
  # Transform the coords into spatial data, including grabbing the bounding box
  # for each "hemisphere"
  west_sf <- st_as_sf(west_map_data,
                      coords = c("group_long", "group_lat"),
                      crs = 4326)
  bbox_west <- st_bbox(west_sf)
  
  east_sf <- st_as_sf(east_map_data,
                      coords = c("group_long", "group_lat"),
                      crs = 4326)
  bbox_east <- st_bbox(east_sf)
  
  
  # 3. Create three map layouts ---------------------------------------------
  
  # Make maps with bubbles + count labels for lakes
  west_bubble_map <- ggplot() +
    mapWorld +
    geom_sf(data = lakes.ne, fill = "lightblue", alpha = 0.7) + 
    geom_sf(data = west_sf, size = 5.5, color = "black", fill = "salmon",
            pch = 21) +
    coord_sf(xlim = c(bbox_west["xmin"], bbox_west["xmax"]),
             ylim = c(bbox_west["ymin"], bbox_west["ymax"])) +
    geom_text(data = west_map_data,
              aes(x = group_long, y = group_lat, label = lake_label),
              size = 3.5, fontface = "bold") +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none", axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "lightcyan1",
                                          colour = "lightcyan1"),
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1),
          panel.grid = element_blank(),
          aspect.ratio = 0.45)
  
  east_bubble_map <- ggplot() +
    mapWorld +
    geom_sf(data = lakes.ne, fill = "lightblue", alpha = 0.7) + 
    geom_sf(data = east_sf, size = 5.5, color = "black", fill = "salmon",
            pch = 21) +
    coord_sf(xlim = c(bbox_east["xmin"], bbox_east["xmax"]),
             ylim = c(bbox_east["ymin"], bbox_east["ymax"])) +
    geom_text(data = east_map_data,
              aes(x = group_long, y = group_lat, label = lake_label),
              size = 3.5, fontface = "bold") +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none", axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "lightcyan1",
                                          colour = "lightcyan1"),
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1),
          panel.grid = element_blank(),
          aspect.ratio = 0.45)
  
  # Similar to an inset map, want to include a world map for location reference
  ref_map <- ggplot() +
    geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf"),
            fill = "gray90", color = NA) +
    geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf")) +
    geom_sf(data = lakes.ne, fill = "white") + 
    geom_sf(data = st_as_sfc(bbox_west),
            color = "red", fill = NA, size = 1) +
    geom_sf(data = st_as_sfc(bbox_east),
            color = "red", fill = NA, size = 1) +
    ylab("") +
    xlab("") +
    theme_bw() +
    # https://stackoverflow.com/questions/58663607/missing-axis-ticks-and-labels-when-plotting-world-map-with-geom-sf
    coord_sf(expand = FALSE) +
    scale_x_continuous(breaks = seq(from = -180, to = 180, by = 60)) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 50, margin = margin(t = 15)),
          panel.border = element_rect(fill = NA,
                                      colour = "black", 
                                      size = 1),
          panel.grid = element_blank())
  
  
  # 4. Combine maps into one figure file ------------------------------------
  
  point_maps <- plot_grid(west_bubble_map, east_bubble_map, nrow = 2)
  
  full_lake_plot <- plot_grid(point_maps, ref_map,
                              ncol = 2, rel_widths = c(2.25, 1.5))
  
  map_figure_path <- "figures/sample_map.png"
  
  ggsave(filename = map_figure_path, plot = full_lake_plot,
         device = "png", width = 10, height = 6.3, units = "in")
  
  return(map_figure_path)
  
}


# 2.2 Run regression trees ------------------------------------------------

# Winter season regression tree
run_winter_tree <- function(winter_model_data, winter_explanatory){
  
  winter_model_formula <- reformulate(
    termlabels = winter_explanatory,
    response = "avechla_log10")
  
  # (Setting seed using a randomly generated integer)
  set.seed(3289)
  winter_tree <- rpart(formula = winter_model_formula,
                       method = "anova",
                       xval = 100,
                       data = winter_model_data %>%
                         filter(across(.cols = any_of(c("avechla_log10",
                                                        winter_explanatory)),
                                       .fns = ~!is.na(.))))
  
  
  # Prune the Tree
  p_winter_tree <- prune(winter_tree,
                         cp = winter_tree$cptable[which.min(winter_tree$cptable[, "xerror"]), "CP"])
  
  return(p_winter_tree)
  
}

# Summer season regression tree
run_summer_tree <- function(summer_model_data, summer_explanatory){
  
  summer_model_formula <- reformulate(termlabels = summer_explanatory,
                                      response = "avechla_log10")
  
  # (Setting seed using a randomly generated integer)
  set.seed(3289)
  summer_tree <- rpart(formula = summer_model_formula,
                       method = "anova",
                       xval = 100,
                       data = summer_model_data %>%
                         filter(across(.cols = any_of(c("avechla_log10",
                                                        summer_explanatory)),
                                       .fns = ~!is.na(.))))
  
  #Prune the Tree
  p_summer_tree <- prune(summer_tree,
                         cp = summer_tree$cptable[which.min(summer_tree$cptable[, "xerror"]), "CP"])
  
  return(p_summer_tree)
  
}


# 2.3 Run linear models ---------------------------------------------------

# Winter season linear model
run_winter_lm <- function(winter_model_data){
  
  ## Fit linear model
  winterNA <- winter_model_data %>%
    filter(!is.na(TP_log10)) %>%
    filter(!is.na(TN_log10)) %>%
    filter(!is.na(avechla_log10))
  
  ## Fit linear model with interactions
  # (Setting seed using a randomly generated integer)
  set.seed(3289)
  m1 <-  lm(avechla_log10 ~ (TP_log10 + TN_log10) *  underice_mean_precip,
            data = winterNA)
  
  return(list(model = m1,
              winterNA = winterNA))
  
}

# Summer season linear model
run_summer_lm <- function(summer_model_data){
  
  ## Fit linear model
  summerNA <- summer_model_data %>% 
    filter(!is.na(watertemp)) %>%
    filter(!is.na(TP_log10)) %>%
    filter(!is.na(TN_log10))
  
  ## Fit linear model with interactions
  # (Setting seed using a randomly generated integer)
  set.seed(3289)
  m1 <- lm(avechla_log10 ~ (TP_log10 + TN_log10) * (watertemp + scale(strat_mean_slp)),
           data = summerNA)
  
  return(list(model = m1,
              summerNA = summerNA))
  
}

# A function to build a summary table for the manuscript
build_summary_table <- function(euli_filtered, winter_model_data, summer_model_data){
  
  # Average station coordinates within each lake
  avg_lake_coords <- euli_filtered %>%
    select(lakename, stationlat, stationlong, lakeelevation) %>%
    group_by(lakename) %>%
    summarize(avg_lat = mean(stationlat),
              avg_long = mean(stationlong),
              lakeelevation = unique(lakeelevation)) %>%
    arrange(lakename)
  
  # Get lake continent designations
  lake_continents <- bind_rows(
    tar_read(modeling_datasets)$winter_model_data %>%
      select(lakename, continent),
    tar_read(modeling_datasets)$summer_model_data %>%
      select(lakename, continent)) %>%
    distinct()
  
  # Pull out the winter data needed
  lake_winter_vals <- winter_model_data %>%
    group_by(lakename) %>%
    summarize(across(.cols = c("TP_log10", "TN_log10", "iceform_mean_tavg",
                               "iceform_mean_slp", "iceform_var_slp", "underice_var_slp",
                               "iceform_mean_precip", "underice_mean_precip",
                               "underice_mean_slp"),
                     .fns = list(mean = mean, min = min, max = max),
                     na.rm = TRUE))
  
  # Text format summary table with mean and range of each variable
  winter_table <- lake_winter_vals %>%
    transmute(
      lakename,
      TP_log10_winter = make_range(mean = TP_log10_mean,
                                   min = TP_log10_min,
                                   max = TP_log10_max),
      TN_log10_winter = make_range(mean = TN_log10_mean,
                                   min = TN_log10_min,
                                   max = TN_log10_max),
      
      iceform_mean_tavg = make_range(mean = iceform_mean_tavg_mean,
                                     min = iceform_mean_tavg_min,
                                     max = iceform_mean_tavg_max),
      
      iceform_mean_slp = make_range(mean = iceform_mean_slp_mean,
                                    min = iceform_mean_slp_min,
                                    max = iceform_mean_slp_max),
      iceform_var_slp = make_range(mean = iceform_var_slp_mean,
                                   min = iceform_var_slp_min,
                                   max = iceform_var_slp_max),
      underice_mean_slp = make_range(mean = underice_mean_slp_mean,
                                     min = underice_mean_slp_min,
                                     max = underice_mean_slp_max),
      underice_var_slp = make_range(mean = underice_var_slp_mean,
                                    min = underice_var_slp_min,
                                    max = underice_var_slp_max),
      iceform_mean_precip = make_range(mean = iceform_mean_precip_mean,
                                       min = iceform_mean_precip_min,
                                       max = iceform_mean_precip_max),
      underice_mean_precip = make_range(mean = underice_mean_precip_mean,
                                        min = underice_mean_precip_min,
                                        max = underice_mean_precip_max))
  
  
  lake_summer_vals <- summer_model_data %>%
    group_by(lakename) %>%
    summarize(across(.cols = c("TP_log10", "TN_log10", "watertemp", "strat_mean_precip",
                               "strat_mean_slp", "strat_var_slp", "prestrat_mean_slp",
                               "prestrat_var_slp" , "prestrat_mean_precip",
                               "prestrat_mean_tavg"),
                     .fns = list(mean = mean, min = min, max = max),
                     na.rm = TRUE))
  
  summer_table <- lake_summer_vals %>%
    transmute(
      lakename,
      TP_log10_summer = make_range(mean = TP_log10_mean,
                                   min = TP_log10_min,
                                   max = TP_log10_max),
      TN_log10_summer = make_range(mean = TN_log10_mean,
                                   min = TN_log10_min,
                                   max = TN_log10_max),
      
      prestrat_mean_tavg = make_range(mean = prestrat_mean_tavg_mean,
                                      min = prestrat_mean_tavg_min,
                                      max = prestrat_mean_tavg_max),
      
      prestrat_mean_slp = make_range(mean = prestrat_mean_slp_mean,
                                     min = prestrat_mean_slp_min,
                                     max = prestrat_mean_slp_max),
      prestrat_var_slp = make_range(mean = prestrat_var_slp_mean,
                                    min = prestrat_var_slp_min,
                                    max = prestrat_var_slp_max),
      strat_mean_slp = make_range(mean = strat_mean_slp_mean,
                                  min = strat_mean_slp_min,
                                  max = strat_mean_slp_max),
      strat_var_slp = make_range(mean = strat_var_slp_mean,
                                 min = strat_var_slp_min,
                                 max = strat_var_slp_max),
      prestrat_mean_precip = make_range(mean = prestrat_mean_precip_mean,
                                        min = prestrat_mean_precip_min,
                                        max = prestrat_mean_precip_max),
      strat_mean_precip = make_range(mean = strat_mean_precip_mean,
                                     min = strat_mean_precip_min,
                                     max = strat_mean_precip_max))
  
  # Combine previous tables
  full_table <- reduce(.x = list(lake_continents, avg_lake_coords,
                                 winter_table, summer_table),
                       .f = full_join,
                       by = "lakename")
  
  # Export the table
  out_path <- "documents/summary_table.csv"
  
  write_csv(x = full_table, file = out_path)
  
  return(list(summary_table = full_table,
              summary_table_out = out_path))
  
}


# 3. General purpose functions --------------------------------------------

# A function to pull anomaly data from the Berkeley Earth Daily Land
# (Experimental; 1880 - Recent) dataset. Retrieves all dates of data for
# each coordinate set provided
get_anomaly <- function(latlong, file_name){
  
  # Latlong to be formatted like:
  #   lakename        poolstation season     Y     X  year
  #   <chr>           <chr>       <chr>  <dbl> <dbl> <dbl>
  # 1 Allequash Lake  NA          iceoff  46.0 -89.6  1982
  # 2 Allequash Lake  NA          iceoff  46.0 -89.6  1983
  # 3 Allequash Lake  NA          iceoff  46.0 -89.6  1984
  
  # file_name = a character string containing a path to a .nc file
  
  # For a single file, pull its lat/lon values and match them to the closest
  # lat/long in our lake dataset. Then use those matches to pull all anomaly
  # data from the file relevant to those locations.
  
  # Get coords from netCDF
  
  nc <- nc_open(file_name)
  lat.vals <- ncvar_get(nc, varid = "latitude")
  lon.vals <- ncvar_get(nc, varid =  "longitude")
  
  # Match lake coords to their closest locations in the netCDF. Returned values
  # are NOT the actual coordinates themselves in the file, but the INDEX
  # of the value that's closest
  coord_locs <- latlong %>%
    rowwise() %>%
    transmute(lakename,
              poolstation,
              latitude = which.min(abs(lat.vals - Y))[1],
              longitude = which.min(abs(lon.vals - X))[1]) %>%
    distinct() %>%
    purrr::transpose()
  
  # Pull data for all dates for a location for several vars
  df_anomaly <- map_df(.x = coord_locs,
                       .f = ~tibble(lakename = .x$lakename,
                                    poolstation = .x$poolstation,
                                    year = ncvar_get(nc = nc, varid = "year"),
                                    month = ncvar_get(nc = nc, varid = "month"),
                                    day = ncvar_get(nc = nc, varid = "day"),
                                    doy = ncvar_get(nc = nc, varid = "day_of_year"),
                                    anomaly_c = ncvar_get(nc = nc, varid = "temperature",
                                                          start = c(
                                                            # Long start position
                                                            .x$longitude,
                                                            # Lat start position
                                                            .x$latitude,
                                                            # Time start position
                                                            1),
                                                          count = c(
                                                            # Read one long
                                                            1,
                                                            # Read one lat
                                                            1,
                                                            # Read all times
                                                            -1)))) %>%
    distinct() %>%
    transmute(lakename,
              poolstation,
              tavg_date = as_date(paste(year, month, day, sep = "-")),
              doy,
              anomaly_c)
  
  
  nc_close(nc)
  
  return(df_anomaly)
}


# A function to pull anomaly data from the Berkeley Earth Daily Land
# (Experimental; 1880 - Recent) dataset
get_climatology <- function(latlong, file_name){
  
  # Same input arguments as get_anomaly()
  
  # For a single file, pull its lat/lon values and match them to the closest
  # lat/long of our lake dataset. Then use those matches to pull all climatology
  # data from the file relevant to those locations.
  
  # Get coords from netCDF
  nc <- nc_open(file_name)
  lat.vals <- ncvar_get(nc, varid = "latitude")
  lon.vals <- ncvar_get(nc, varid =  "longitude")
  
  # Match lake coords to their closest locations in the netCDF. Returned values
  # are NOT the actual coordinates themselves in the file, but the INDEX
  # of the value that's closest
  coord_locs <- latlong %>%
    rowwise() %>%
    transmute(lakename,
              poolstation,
              latitude = which.min(abs(lat.vals - Y))[1],
              longitude = which.min(abs(lon.vals - X))[1]) %>%
    distinct() %>%
    purrr::transpose()
  
  df_clim <- map_df(.x = coord_locs,
                    .f = ~tibble(lakename = .x$lakename,
                                 poolstation = .x$poolstation,
                                 climatology_c = ncvar_get(
                                   nc = nc,
                                   varid = "climatology",
                                   start = c(
                                     # Long start position
                                     .x$longitude,
                                     # Lat start position
                                     .x$latitude,
                                     # Time start position
                                     1),
                                   count = c(
                                     # Read one long
                                     1,
                                     # Read one lat
                                     1,
                                     # Read all times
                                     -1))) %>%
                      mutate(day_num = row_number())
  )
  
  nc_close(nc)
  
  return(df_clim)
}


# A function to pull climate data from the Climatic Research Unit (CRU) dataset
get_cru <- function(latlong, variable_name, file_name){
  
  # Latlong to be formatted same as get_slp()
  
  # variable_name = a character string containing the name of the desired
  #                 climate variable as written in the CRU netCDF file
  
  # file_name = a character string containing a path to a .nc file
  
  # For a single file, pull its lat/lon values and match them to the closest
  # lat/long of our lake dataset. Then use those matches to pull all climate
  # data from the file relevant to those locations.
  
  # Get coords from netCDF
  nc <- nc_open(file_name)
  lat.vals <- ncvar_get(nc, varid = "lat")
  lon.vals <- ncvar_get(nc, varid =  "lon")
  
  nc_time <- as_date("1900-01-01") + ncvar_get(nc = nc, varid = "time")
  
  min_nc_time <- min(nc_time)
  max_nc_time <- max(nc_time)
  
  # Match lake coords to their closest locations in the netCDF. NOT the actual
  # coordinates themselves in the file...the POSITION of the value that's closest
  coord_locs <- latlong %>%
    # Accept either iceform or prestrat bounds
    rename_with(.fn = ~str_remove(pattern = "iceform_|prestrat_", string = .),
                .cols = contains("bound")) %>%
    rowwise() %>%
    mutate(latitude = which.min(abs(lat.vals - Y))[1],
           longitude = which.min(abs(lon.vals - X))[1]) %>%
    distinct() %>%
    purrr::transpose()
  
  # Get data for the ice-formation or pre-stratification period
  df_period <- map_df(.x = coord_locs,
                      .f = ~tibble(
                        # Convert the list back into dataframe form, then combine with...
                        map_df(.x = .,
                               .f = ~.x),
                        # ...the climate data
                        clim_date = as_date("1900-01-01") + days(ncvar_get(nc = nc,
                                                                           varid = "time")),
                        # Name column based on input to function
                        !!variable_name := ncvar_get(nc = nc, varid = variable_name,
                                                     start = c(.x$longitude,
                                                               .x$latitude,
                                                               # 'time' dimension
                                                               1),
                                                     count = c(1, 1, -1)))) %>%
    mutate(across(.cols = c(contains("bound"), contains("date")), .fns = ~as_date(.)),
           # How much does the climate month overlap with the ice-formation/pre-stratification period?
           # https://stackoverflow.com/questions/58517015/lubridate-find-overlap-time-between-interval-and-a-date
           season_overlap = as.duration(lubridate::intersect(
             # The month of the climate date
             x = interval(start = rollback(ymd(clim_date),
                                           roll_to_first = TRUE),
                          end = rollforward(dates = ymd(clim_date))),
             # The ice-formation/pre-stratification period
             y = interval(start = lower_bound,
                          end = upper_bound))), 
           # Format overlap time in weeks
           overlap_weeks = as.numeric(season_overlap, "weeks")) %>%
    # Keep any months that overlap with the interval by two weeks or more
    filter(overlap_weeks >= 2) %>%
    rename(euli_year = year,
           lat = Y,
           lon = X)
  
  # Get data for the full min/max sampling period (under-ice or stratification)
  df_sample <- map_df(.x = coord_locs,
                      .f = ~tibble(
                        # Convert the list back into dataframe form, then combine with...
                        map_df(.x = .,
                               .f = ~.x),
                        # ...the climate data
                        clim_date = as_date("1900-01-01") + days(ncvar_get(nc = nc, varid = "time")),
                        # Name column based on input to function
                        !!variable_name := ncvar_get(nc = nc, varid = variable_name,
                                                     start = c(.x$longitude,
                                                               .x$latitude,
                                                               # 'time' dimension
                                                               1),
                                                     count = c(1, 1, -1)))) %>%
    
    mutate(across(.cols = c(contains("bound"), contains("date")), .fns = ~as_date(.)),
           # Establish time period that CRU month covers
           climate_month = interval(start = rollback(ymd(clim_date),
                                                     roll_to_first = TRUE),
                                    end = rollforward(dates = ymd(clim_date))),
           # How much does the climate month overlap with the sampling period?
           sampling_overlap = as.duration(lubridate::intersect(
             # The month of the climate date
             x = climate_month,
             # The sampling period
             y = interval(start = start_date_min,
                          end = end_date_max))), 
           # Format overlap time in weeks
           overlap_weeks = as.numeric(sampling_overlap, "weeks")) %>%
    group_by(lakename, poolstation, year, start_date_min, end_date_max) %>%
    mutate(overlap_max = max(overlap_weeks, na.rm = TRUE)) %>%
    # How many records are there that overlap at all?
    # Keep any months that overlap with the interval by two weeks or more
    # OR if this would result in no data, then keep the month with max overlap
    # OR (in case single-day periods are treated differently...) keep the month
    # that single-day sample periods fall within
    # Want to double check how this handles periods that, e.g., are two days long
    # and have one day in each of two months. I think they would just be averaged,
    # which is probably fine
    filter( (overlap_weeks >= 2) |
              ((overlap_max < 2) & (overlap_weeks == overlap_max)) |
              ( (start_date_min == end_date_max) & (start_date_min %within% climate_month) ) ) %>%
    rename(euli_year = year,
           lat = Y,
           lon = X)
  
  nc_close(nc)
  
  return(list("df_period" = df_period, "df_sample" = df_sample))
}

# A function to paste together values for output summary table
make_range <- function(mean, min, max){
  
  range_text <- paste0(round(mean, digits = 2),
                       " [",
                       round(min, digits = 2),
                       ", ",
                       round(max, digits = 2),
                       "]")
  
  return(range_text)
  
}

