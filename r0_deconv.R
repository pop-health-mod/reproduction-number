# Main script to calculate the effective reproduction number (Rt)

# ---- Set up ----
library(EpiEstim)
# Requires dplyr 1.0.0 or greater.
library(dplyr)
library(ggplot2)

# Requires an R project to use here::here()
source(here::here("src", "r0utils.R"))
source(here::here("src", "plot.R"))
source(here::here("src", "deconv.R"))
source(here::here("src", "deconv_uncer.R"))

# Declare the path to the repository hosting the data and where to output res.
pathData <- "../data/"
pathFile <- here::here("output", "out")

# ---- Parameter set-up----
# Provincial estimates
rss <- data_frame(filename = c("qc"), name = c("Québec"), RSS_code = c(NA))
rss$cases_7days <- NA

# Simulations and sliding window period (days) for continuous Rt estimation
window <- 5

step_rt <- FALSE

# Number of days to exclude to remove reporting delays
days_exclude <- data.frame(variable = c("Cas", "Hospitalisations", "Décès"),
                           number = c(2, 2, 5))
date_data_cut <- as.Date("2020-09-01")

# ---- Around parameters ----
#  Delays and SD from normal approximation
delays <- data.frame(
  "Incubation" = c(3.2, 5.3),
  "Cas" = c(4.19, 3.54),
  "Épidémiologique" = c(6.03, 3.74),
  "Hospitalisations" = c(4.77, 4.98),
  "Décès" = c(8.33, 11.91)
)

# Plotting the time series
for (i in 1:dim(rss)[1]) {
  load_df(sprintf("%s%s.rda", pathData, rss$filename[i]),
                  rss$name[i], t_end = date_data_cut, plot = T,
          days_exclude = days_exclude, case_only = T)
}


for (i in 1:dim(rss)[1]) {
# ---- Loading data ----
  long <- load_df(sprintf("%s%s.rda", pathData, rss$filename[i]),
                  rss$name[i], t_end = date_data_cut, plot = F,
                  days_exclude = days_exclude)
  long_c <- long[long$dates < (date_data_cut - days_exclude$number[days_exclude$variable == "Cas"]), ]
  rss$cases_7days[i] <- sum(long_c[long_c$variable == "Cas" & long_c$dates > (as.Date(max(long_c$dates)) - 7), ]$value)
}

start <- Sys.time()
i <- 1
rt_last <- NULL
for (i in 1:dim(rss)[1]) {
# ---- Loading data ----
  long <- load_df(sprintf("%s%s.rda", pathData, rss$filename[i]),
                  rss$name[i], t_end = date_data_cut, plot = F,
                  days_exclude = days_exclude)
  
  filename <- sprintf("%s/r0_%s_%s.rds", pathFile, rss$name[i], date_data_cut)
  figname <- sprintf("rss-%s-%s", as.character(rss$RSS_code[i]), rss$name[i])
  print(filename)
  print(figname)
  
  # ---- Adjusting the time series ----
  # Cases, either lab-confirmed of via epidemiological links
  long_c <- long[long$dates < (date_data_cut - days_exclude$number[days_exclude$variable == "Cas"]), ]

  # Lab-confirmed cases
  local <- ts_adjustment_deconv_incub(inci_dat = long_c[long_c$variable == "Local",], delays,
                                      variable = "Cas", days_incl = 21, degree = 1)

  # Epidemiological link
  epi <- ts_adjustment_deconv_incub(inci_dat = long_c[long_c$variable == "Épidémiologique",],
                                    delays, variable = "Épidémiologique",
                                    days_incl = 21, degree = 1)

  # Imported cases
  import <- adjust_imported(imported_ts = long_c[long_c$variable == "Importé",])
  
  # Merge the lab-confirmed and epidemiological
  local_epi <- merge(x = local, y = epi, by = "dates", all = T) %>%
    rowwise() %>%
    mutate(local = sum(local.x, local.y, na.rm = T)) %>%
    select(dates, local) %>%
    ungroup()
  
  # Include imported cases
  cas_ts <- merge(x = local_epi, y = import[import$dates <= max(local_epi$dates),],
                  by = "dates", all = T)
  cas_ts[is.na(cas_ts)] <- 0
  cas_ts <- cas_ts %>% filter(dates < max(dates))
  
  # Hospitalisations
  long_h <- long[long$dates < (date_data_cut - days_exclude$number[days_exclude$variable == "Hospitalisations"]), ]
  hosp_adj <- ts_adjustment_deconv_incub(long_h[long_h$variable == "Hospitalisations",],
                                         delays, "Hospitalisations", days_incl = 21, degree = 1)
  hosp_ts <- hosp_adj
  
  # Morts
  long_d <- long[long$dates < (date_data_cut - days_exclude$number[days_exclude$variable == "Décès"]), ]
  mort_adj <- ts_adjustment_deconv_incub(long_d[long_d$variable == "Décès",], delays, "Décès",
                                         days_incl = 21, degree = 1)
  mort_ts <- mort_adj
  
  # ---- Truncation of the series using quantiles of the delays distributions ----
  # The adjustment for right truncation is imperfect... we should remove the days
  # for which we have less than 1/3 right-truncated
  shape_inc <- (delays[, "Incubation"][2]^2) / (delays[, "Incubation"][1]^2)
  scale_inc <- (delays[, "Incubation"][1]^2) /  delays[, "Incubation"][2]
  shape_cas <- (delays[, "Cas"][2]^2) / (delays[, "Cas"][1]^2)
  scale_cas <- (delays[, "Cas"][1]^2) /  delays[, "Cas"][2]
  shape_hop <- (delays[, "Hospitalisations"][2]^2) / (delays[, "Hospitalisations"][1]^2)
  scale_hop <- (delays[, "Hospitalisations"][1]^2) /  delays[, "Hospitalisations"][2]
  shape_dec <- (delays[, "Décès"][2]^2) / (delays[, "Décès"][1]^2)
  scale_dec <- (delays[, "Décès"][1]^2) /  delays[, "Décès"][2]
  
  offset_cases <-  round(quantile(rgamma(1E6, shape = shape_inc, scale = scale_inc) +
                                    rgamma(1E6, shape = shape_cas, scale = scale_cas),
                                  probs = 1 / 3)) 
  offset_hospit <- round(quantile(rgamma(1E6, shape = shape_inc, scale = scale_inc) +
                                    rgamma(1E6, shape = shape_hop, scale = scale_hop),
                                  probs = 1 / 3)) 
  offset_death <-  round(quantile(rgamma(1E6, shape = shape_inc, scale = scale_inc) +
                                    rgamma(1E6, shape = shape_dec, scale = scale_dec),
                                  probs = 1 / 3))
  
  date_cases  <- max(cas_ts$dates)
  date_hospit <- max(hosp_ts$dates)
  date_death  <- max(mort_ts$dates)

  # ---- Cas - Rt estimation with Cori ----
  print("Rt for cases")  
  r0_cas <- calc_r0(inci_dat = long_c[long_c$variable == "Local", ],
                    results = cas_ts, window, 
                    estimated_offset = offset_cases, step_rt = step_rt)
  cas <- r0_cas %>%
    mutate(variable = "Cas") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  # Latest Rt by case (mean and sd)
  r0_cas_last <- r0_cas %>%
    filter(dates %in% as.Date(c(date_death, date_hospit, date_cases))) %>%
    mutate(variable = "Cas") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  # ---- Hospits - Rt estimation with Cori  ----
  print("Rt for hospitalisation")
  r0_hosp <- calc_r0(inci_dat = long_h[long_h$variable == "Hospitalisations", ],
                     results = hosp_ts, window, 
                     estimated_offset = offset_hospit, step_rt = step_rt)
  hosp <- r0_hosp %>%
    mutate(variable = "Hospitalisations") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  r0_hosp_last <- r0_hosp %>%
    filter(dates %in% as.Date(c(date_death, date_hospit))) %>%
    mutate(variable = "Hospitalisations") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  # ---- Décès - Rt estimation with Cori ----
  print("Rt for death")
  r0_mort <- calc_r0(inci_dat = long_d[long_d$variable == "Décès", ],
                     results = mort_ts, window, 
                     estimated_offset = offset_death, step_rt = step_rt)
  
  r0_mort_last <- r0_mort %>%
    filter(dates %in% as.Date(date_death))  %>%
    mutate(variable = "Décès") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  mort <- r0_mort %>%
    mutate(variable = "Décès") %>%
    mutate(region = rss$name[i]) %>%
    mutate(RSS_code = rss$RSS_code[i])
  
  # All R0 estimations in a single dataframe
  r0_long <- do.call(rbind, list(cas, hosp, mort)) %>%
    mutate(variable = factor(variable, levels = c("Cas",
                                                  "Hospitalisations",
                                                  "Décès")))
  # Rt on the last day of the interval
  r0_last <- do.call(rbind, list(r0_cas_last, r0_hosp_last, r0_mort_last)) %>%
    mutate(variable = factor(variable, levels = c("Cas",
                                                  "Hospitalisations",
                                                  "Décès")))
  
  # Saving raw output
  saveRDS(r0_cas, file = sprintf("%s/r0_raw_cas_%s_%s.rds",
                                 pathFile, rss$name[i], date_data_cut))
  saveRDS(r0_hosp, file = sprintf("%s/r0_raw_hosp_%s_%s.rds",
                                  pathFile, rss$name[i], date_data_cut))
  saveRDS(r0_mort, file = sprintf("%s/r0_raw_mort_%s_%s.rds",
                                  pathFile, rss$name[i], date_data_cut))
  
  # Saving last day Rt
  rt_last <- rbind(rt_last, r0_last)
  saveRDS(r0_last, file = sprintf("%s/r0_last_%s_%s.rds",
                                  pathFile, rss$name[i], date_data_cut)) 
  saveRDS(rt_last, file = sprintf("%s/rt_last_%s.rds", pathFile, date_data_cut))
  
  # Saving datafile for daily estimation
  saveRDS(r0_long, file = filename)
}
end <- Sys.time() - start
print(end)



# ---- Plotting the results ----
# ---- Plotting the results ----
for (i in 1:dim(rss)[1]) {
  filename <- sprintf("%s/r0_%s_%s.rds", pathFile, rss$name[i], date_data_cut)
  figname <- sprintf("deconv-rss-%s-%s", as.character(rss$RSS_code[i]), rss$name[i])
  print(figname)
  
  r0_long <- readRDS(file = filename)
  
  if (rss$name[i] == "Québec") {
    r0_inspq <- as.data.frame(r0_long) %>%
      filter(dates >= as.Date("2020-03-13") & variable == "Cas") %>%
      select(dates, mean, lci, uci)
    
    write.csv(r0_inspq,
              file = sprintf("%s/r0_inspq_%s_%s.csv", pathFile, rss$name[i], date_data_cut),
              row.names = F)
  }
  # plot_r0(r0_long, figname)
  plot_r0_one(df = r0_long, rss_name = figname, cv = T, case_only = T)
}
