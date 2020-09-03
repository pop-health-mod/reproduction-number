# ---- Adjustment of the time series ----
ts_adjustment_deconv_incub <- function(inci_dat, delays, variable, min_chi_squared = 1,
                                       max_iter = 30, max_first_delay = 30, 
                                       smooth_incidence = T, verbose = TRUE,
                                       days_incl = 21, degree = 1) {
  
  length_initial_time_series <- inci_dat %>%
    select(dates)%>%
    mutate(dates = as.Date(dates)) %>%
    tidyr::complete(dates = seq.Date(min(dates), max(dates), by = "days")) %>%
    pull() %>%
    length()
  
  minimal_date <- min(as.Date(inci_dat$dates)) - max_first_delay
  maximal_date <- max(as.Date(inci_dat$dates))
  
  all_dates <- seq(minimal_date, maximal_date, by = "days")
  
  # Smoothing the time series observed using local regression
  time_series <- inci_dat
  if (smooth_incidence == T) {
    smoothed_incidence_data <- time_series %>%
      mutate(dates = as.Date(dates)) %>%
      mutate(value = getLOESSCases(dates = dates, count_data = value,
                                   days_incl = days_incl, degree = degree)) %>%
      tidyr::complete(dates = seq.Date(minimal_date, maximal_date, by = "days"),
                      fill = list(value = 0, variable = unique(variable)))
    
    raw_total_incidence <- sum(time_series$value, na.rm = TRUE)
    smoothed_total_incidence <- sum(smoothed_incidence_data$value, na.rm = T)
    
    if (smoothed_total_incidence > 0) {
      smoothed_incidence_data <- smoothed_incidence_data %>%
        mutate(value = value * raw_total_incidence / smoothed_total_incidence)
    }
    
  } else {
    smoothed_incidence_data <- time_series  %>%
      complete(date = seq.Date(minimal_date, maximal_date, by = "days"),
               fill = list(value = 0))
  }
  
  first_incid <-  with(filter(smoothed_incidence_data, cumsum(value) > 0),
                       value[which.min(dates)])
  last_incid <- with(smoothed_incidence_data, value[which.max(dates)])

  
  # Shape and scale parameters of the incubation periods
  shape_incub <- (delays[, "Incubation"][1]^2) / (delays[, "Incubation"][2]^2)
  scale_incub <- (delays[, "Incubation"][2]^2) /  delays[, "Incubation"][1]
  
  # Shape and scale parameters of the events (case, hosp, death)
  shape_event <- (delays[, variable][1]^2) / (delays[, variable][2]^2)
  scale_event <- (delays[, variable][2]^2) /  delays[, variable][1]
  
  # Sample the the delays matrices to adjust the time series
  
  # Constant delays between infection and event
  constant_delay <- get_constant_waiting_time(shape_incub = shape_incub,
                                              scale_incub = scale_incub,
                                              shape_event = shape_event,
                                              scale_event = scale_event)
  
  # Delay distribution matrix based on constant delays
  delay_distribution_matrix <- get_matrix_constant_waiting_time(constant_delay,
                                                                all_dates)
  
  
  Q_matrix <- apply(delay_distribution_matrix, MARGIN = 2, sum)
  
  # First guess to start the deconvolution algorithm and use mode of the 
  # constant delay distribution. Because indices are offset
  # by one as the delay can be 0.
  first_guess_delay <- which.max(constant_delay) - 1
  truncation_delay <-  0
  
  first_guess <- smoothed_incidence_data %>%
    mutate(dates = dates - first_guess_delay) %>%
    filter(cumsum(value) > 0) %>% # remove leading zeroes
    tidyr::complete(dates = seq.Date(minimal_date, min(dates), by = "days"),
                    fill = list(value = first_incid, variable = unique(variable))) %>%
    tidyr::complete(dates = seq.Date(max(dates), maximal_date, by = "days"),
                    fill = list(value = last_incid)) %>% # pad with last recorded value
    arrange(dates) %>% 
    filter(dates >=  minimal_date)
  
  final_estimate <- iterate_RL(
    first_guess$value,
    smoothed_incidence_data$value,
    delay_distribution_matrix = delay_distribution_matrix,
    Q_matrix = Q_matrix,
    threshold_chi_squared = min_chi_squared,
    max_iterations = max_iter,
    max_delay = max_first_delay,
    verbose = verbose)
  
  ## right-truncate trailing zeroes induced by initial shift by 'first_guess_delay'
  infection_dates <- first_guess %>%
    filter(dates <= maximal_date - truncation_delay) %>%
    pull(dates)
  
  infection_ts <- round(final_estimate[seq_len(length(final_estimate) - truncation_delay)])
  
  ## dataframe containing results
  if (variable %in% c("Cas", "Épidémiologique")) {
    df_infections <- data.frame(dates = infection_dates,
                                local = infection_ts,
                                stringsAsFactors = F)
  } else {
    df_infections <- data.frame(dates = infection_dates,
                                I = infection_ts,
                                stringsAsFactors = F)
  }
  
  return(df_infections)
}

