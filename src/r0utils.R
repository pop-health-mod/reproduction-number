# Functions for producing R0 estimates
# Spring 2020

require(dplyr)
require(heavy)
require(EpiNow)
require(EpiSoon)
require(surveillance)
require(coga)

# Mostly reproducing method from Sciré et al. 2020
# Empirical distribution for time from infection to confirmation

# ---- Utility function to find maximum dates in a list ----
find_max_date_list <- function(x, n_sim) {   
  vec <- c()
  for (i in 1:n_sim) {
    vec <- c(vec, as.character(max(x[[i]]$dates)))
  }
  return(max(as.Date(vec)))
}

# ---- Process data ----
load_df <- function(path_to_data, rss_file, t_start = as.Date("2020-02-28"),
                    t_end = Sys.Date(), plot = T, days_exclude, case_only = F) {
  rss_dat <- get(load(path_to_data))
  # TODO One ifelse statement to remove unecessary data loads.
  rss_dat <- as.data.frame(rss_dat) %>%
    filter(dates >= t_start & dates < t_end)  %>%
    mutate(importation = nc_travel_impute + nc_travel_impute_after_quarantine_date_fin_voyage) %>%
    mutate(overall = nc_lab + importation + nc_epi) %>%
    select(dates, "Cas" = overall, "Épidémiologique" = nc_epi, "Local" = nc_lab,
           "Importé" = importation, "Hospitalisations" = nh_inspq, "Décès" = dx_INSPQ_autre)
  # Coding missing value to 0
  rss_dat[is.na(rss_dat)] <- 0
  # Reshaping to long format
  rss_long <- reshape2::melt(rss_dat, .id = "dates")
  # Plot the TS
  if (plot == T) {
    rss_long_plot <- rss_dat %>%
      reshape2::melt(id.vars = "dates") %>%
      mutate(imported = ifelse(variable == "Local", 1,
                               ifelse(variable == "Épidémiologique", 2,
                                      ifelse(variable == "Importé", 3,
                                             ifelse(variable == "Hospitalisations", 4, 5))))) %>%
      filter(variable != "Cas") %>%
      mutate(variable = ifelse(variable %in% c("Importé", "Local", "Épidémiologique"), "Cas",
                               as.character(variable))) %>%
      mutate(imported = factor(imported, labels = c("Local",
                                                    "Épidémiologique",
                                                    "Importé",
                                                    "Hospitalisations",
                                                    "Décès"))) %>%
      mutate(variable = factor(variable, levels = c("Cas",
                                                    "Hospitalisations",
                                                    "Décès")))
    
    plot_ts_unadj(rss_long_plot, c("Cas", "Hospitalisations", "Décès"), rss_file,
                  days_exclude = days_exclude, case_only = case_only)
  } else {
    return(rss_long)
  }
}


# ---- Empirical CDF to draw the reporting probability ----
make_empirical_cdf <- function(sh, sc, numSample = 1e6) {
  draws <- round(rgamma(numSample, shape = sh[1], scale = sc[1]) +
                   rgamma(numSample, shape = sh[2], scale = sc[2]))
  return(Vectorize(ecdf(draws)))
}


# ---- Balance input delay and incubation samples ----
balance_dfs <- function(df1, df2) {
  if (nrow(df1) > nrow(df2)) {
    df2 <- data.table::rbindlist(list(
      df2,
      df2[sample(1:nrow(df2), (nrow(df1) - nrow(df2)), replace = TRUE), ]
    ))
  }
  return(df2)
}

# ---- Estimation of R0 ----
calc_r0 <- function(inci_dat, results, window, estimated_offset,
                    rt_prior = c(2.6, 2), 
                    si = c(5.2, 1.7), 
                    si_dist = "G",
                    start_date = as.Date("2020-02-21") - 5,
                    step_rt = FALSE) {
  # Last date of report (for NA)
  
  last_date <- max(as.Date(inci_dat$dates)) - estimated_offset

  # Preliminary test Using Cori
  incid_sim <- results
  incid_sim <- incid_sim[incid_sim$dates >= start_date & incid_sim$dates <= last_date, ]
  t_start <- seq(2, as.numeric(max(as.Date(incid_sim$dates)) - min(as.Date(incid_sim$dates))) - window + 1)
  t_end <- t_start + window 
  
  # if we want weekly time step instead
  if (step_rt) {
    incid_sim <- incid_sim[incid_sim$dates >= start_date & incid_sim$dates <= last_date, ]
    
    # Cases - we align start date
    if (any(names(incid_sim) %in% "local")) {
      if (min(incid_sim$dates) > start_date) {
        incid_to_add_first <- data.frame(dates = seq(start_date, as.Date(min(incid_sim$dates) - 1), by = "days"),
                                         local = 0, imported = 0)
      } else { incid_to_add_first <- NULL }
      # Cases - we align end date  
      if (max(incid_sim$dates) < last_date) {
        incid_to_add_last <- data.frame(dates = seq(as.Date(max(incid_sim$dates) + 1), last_date, by = "days"),
                                        local = 0, imported = 0)
      } else { incid_to_add_last <- NULL }
      incid_sim <- rbind(incid_to_add_first, incid_sim, incid_to_add_last)
    }
    
    # Hospit + Deaths - we align start date
    if (!any(names(incid_sim) %in% "local")) {
      if (min(incid_sim$dates) > start_date) {
        incid_to_add_first <- data.frame(dates = seq(start_date, as.Date(min(incid_sim$dates) - 1), by = "days"),
                                         I = 0)
      } else { incid_to_add_first <- NULL }
      #  Hospit + Deaths  - we align end date  
      if (max(incid_sim$dates) < last_date) {
        incid_to_add_last <- data.frame(dates = seq(as.Date(max(incid_sim$dates) + 1), last_date, by = "days"),
                                        I = 0)
      } else { incid_to_add_last <- NULL }
      incid_sim <- rbind(incid_to_add_first, incid_sim, incid_to_add_last)
    }
    
    t_start <- seq(from = 2, to = as.numeric(max(as.Date(incid_sim$dates)) - min(as.Date(incid_sim$dates))) - 7, by = 7)
    t_end <- c(t_start[-1] - 1, t_start[length(t_start)] + 6) 
  }
  
  res_i <- estimate_R(incid_sim, method = "parametric_si",
                      config = make_config(list(
                        t_start = t_start,
                        t_end = t_end,
                        mean_si = si[1],
                        std_si = si[2],
                        mean_prior = rt_prior[1],
                        std_prior = rt_prior[2])))
  
  # Store mean/median and 95% CI for the given simulation
  last_estim <- max(start_date + t_end)
  output_dates <- seq(start_date + 1, last_date, by = 'days')
  # epiestim plots at the end of the sliding window
  obs <- data.frame(dates = incid_sim$dates[t_end])
  obs$mean <- res_i$R$`Mean(R)`
  obs$uci <- res_i$R$`Quantile.0.975(R)`
  obs$lci <- res_i$R$`Quantile.0.025(R)`
  obs$se <- res_i$R$`Std(R)`
  obs <- obs[obs$dates <= last_date, ]
  obs <- na.omit(obs)
  obs_na <- data.frame(dates = output_dates, mean = NA, uci = NA, lci = NA, se = NA)
  df_cori <- obs
  if ((min(output_dates) < min(obs$dates)) | (max(output_dates) > max(obs$dates))) {
    sel <- !(output_dates %in% obs$dates)
    df_cori <- rbind(obs, obs_na[sel, ])
    df_cori <- df_cori[order(df_cori$dates), ]
  }
  if (step_rt) {
    df_cori <- tidyr::fill(df_cori, mean:se, .direction = "down")
  }
  
  return(df_cori)
}