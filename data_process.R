# ---- Libraries ----
library(tidyverse)

# ---- Load datasets ----
# download from https://www.inspq.qc.ca/covid-19/donnees
data_case <- read_csv("./data-public/graph_1-2_page_principale.csv")
data_hosp <- read_csv("./data-public/graph_3-2_page_principale.csv")
data_mort <- read_csv("./data-public/graph_2-2_page_principale.csv")

names(data_case) <- c("dates", "nc_epi", "nc_lab", "mean_7wk")
names(data_hosp) <- c("dates", "nh_reg", "nh_icu", "mean_7wk")
names(data_mort) <- c("dates", "dx_autre", "dx_domicile", "dx_rpa", "dx_chsld_ch", "mean_7wk")

# ---- Merge and format ----
# cases
data_qc <- data_case[, c("dates", "nc_lab", "nc_epi")]

# hospitalizations
data_qc <- left_join(data_qc, data_hosp[, c("dates", "nh_reg", "nh_icu")], by = "dates")
data_qc$nh_inspq <- data_qc$nh_reg + data_qc$nh_icu

# deaths
data_qc <- left_join(data_qc, data_mort[, c("dates", "dx_autre", "dx_domicile")], by = "dates")
data_qc$dx_INSPQ_autre <- data_qc$dx_autre + data_qc$dx_domicile

data_qc$dx_INSPQ_autre[is.na(data_qc$dx_INSPQ_autre)] <- 0

# fix dates and reorder
data_qc$dates <- as.Date(data_qc$dates)
data_qc$dates <- as.character(data_qc$dates)

data_qc <- data_qc[, c("dates", "nc_lab", "nc_epi", "dx_INSPQ_autre", "nh_inspq")]

data_qc$nc_travel_impute <- 0
data_qc$nc_travel_impute_after_quarantine_date_fin_voyage <- 0

save(data_qc, file = "./data/qc_noise.rda")
