# Adjustment by deconvolution
# August 2020
# Adapted from the Swisss methodology 
# covid-19-Re/shiny-dailyRe app/otherScripts/2_utils_getInfectionIncidence.R
library(dplyr)
library(tidyr)

# ---- Smoothing the time series ----
getLOESSCases <- function(dates, count_data, days_incl = 21, degree = 1,
                          truncation = 0) {
  
  if (truncation != 0) {
    dates <- dates[1:(length(dates) - truncation)]
    count_data <- count_data[1:(length(count_data) - truncation)]
  }
  
  n_points <- length(unique(dates))
  sel_span <- days_incl / n_points
  
  n_pad <- round(length(count_data) * sel_span * 0.5)
  
  c_data <- data.frame(value = c(rep(0, n_pad), count_data),
                       date_num = c(seq(as.numeric(dates[1]) - n_pad, as.numeric(dates[1]) - 1),
                                    as.numeric(dates)))
  c_data.lo <- loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  smoothed <- predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <- round(
    raw_smoothed_counts * sum(count_data, na.rm = T) / sum(raw_smoothed_counts, na.rm = T))
  
  if (truncation != 0) {
    normalized_smoothed_counts <- append(normalized_smoothed_counts, rep(NA, truncation))
  }
  return(normalized_smoothed_counts)
}

# ---- Empirical CDF to draw the reporting probability ----
make_empirical_cdf <- function(sh, sc, numSample = 1e6) {
  draws <- round(rgamma(numSample, shape = sh[1], scale = sc[1]) +
                   rgamma(numSample, shape = sh[2], scale = sc[2]))
  return(Vectorize(ecdf(draws)))
}

# ---- Sampling the waiting time distribution ----
get_constant_waiting_time <- function(shape_incub, scale_incub, shape_event,
                                      scale_event, length_out = 28,
                                      n_samples = 1E6) {
  
  F_h <- make_empirical_cdf(sh = c(shape_incub, shape_event),
                            sc = c(scale_incub, scale_event),
                            numSample = n_samples)
  
  f <- Vectorize(function(x){
    if(x < 0) {
      return(0)
    } else if(x < 0.5) {
      return(F_h(0.5))
    } else {
      # Reprend le PDF
      return(F_h(round(x + 1E-8) + 0.5) - F_h(round(x + 1E-8) - 0.5))
    }
  })
  
  # Returning the mean delay for the given distribution
  pdf_mean <- ceiling(mean(rgamma(n_samples,
                            shape = shape_incub,
                            scale = scale_incub)))
  
  x <- 0:(length_out - 1)
  
  # Normalizing the pdf and mean of 
  val <- f(x) / sum(f(x))
  
  return(list(val, pdf_mean))
}

# ---- Matrix format of the previously computed vector of waiting time ----
get_matrix_constant_waiting_time <- function(waiting_time_distr,
                                            all_dates) {
  N <- length(all_dates)
  
  if(length(all_dates) >= length(waiting_time_distr)) {
    waiting_time_distr <- c(waiting_time_distr, rep(0, times = N - length(waiting_time_distr)))
  }
  
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ),
                                          waiting_time_distr[1:(N - i + 1)])
  }
  
  return(delay_distribution_matrix)
}

# ---- Main iteration of the deconvolution ----
iterate_RL <- function(
  initial_estimate,
  original_incidence,
  delay_distribution_matrix,
  Q_matrix,
  threshold_chi_squared = 1,
  max_iterations = 20,
  max_delay,
  verbose = FALSE) {
  
  current_estimate <- initial_estimate
  N <- length(current_estimate)
  N0 <- N - max_delay
  chi_squared <- Inf
  count <- 1
  
  while(chi_squared > threshold_chi_squared & count <= max_iterations) {
    
    if (verbose) {
      cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    }
    
    E <- as.vector(delay_distribution_matrix %*% current_estimate)
    B <- tidyr::replace_na(original_incidence/E, 0)
    
    current_estimate <- current_estimate / Q_matrix * as.vector(crossprod(B, delay_distribution_matrix))
    current_estimate <- tidyr::replace_na(current_estimate, 0)
    
    chi_squared <- 1/N0 * sum((E - original_incidence)^2/E, na.rm = T)
    count <- count + 1
  }
  
  return(current_estimate)
}

# ---- Function to account for imported cases ----
adjust_imported <- function(imported_ts){
  # We first make it in the appropriate format
  adj_import <- imported_ts %>%
    select(dates, value) %>%
    mutate(dates = as.Date(dates)) %>%
    rename(imported = value) 
  
  return(adj_import)
}
