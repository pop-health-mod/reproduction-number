# Plotting functions time series and R0
# Spring 2020
library(ggplot2)
library(ggridges)

# ---- Plot the unadjusted time series ----
plot_ts_unadj <- function(inci_dat_long, vars, rss_name, days_exclude,
                          case_only = F) {
  if (rss_name == "Montréal") {
    # Setting up limit to each variables
    inci_dat_long$min <- NA
    inci_dat_long$max <- NA
    inci_dat_long[inci_dat_long$variable == "Cas",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Cas",]$max <- 650
    inci_dat_long[inci_dat_long$variable == "Hospitalisations",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Hospitalisations",]$max <- 90
    inci_dat_long[inci_dat_long$variable == "Décès",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Décès",]$max <- 30
  } else {
    inci_dat_long$min <- NA
    inci_dat_long$max <- NA
    inci_dat_long[inci_dat_long$variable == "Cas",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Cas",]$max <- 170
    inci_dat_long[inci_dat_long$variable == "Hospitalisations",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Hospitalisations",]$max <- 45
    inci_dat_long[inci_dat_long$variable == "Décès",]$min <- 0
    inci_dat_long[inci_dat_long$variable == "Décès",]$max <- 15
  }

  # For plotting the data in french  
  old_locale <- Sys.getlocale("LC_TIME")
  if (.Platform$OS.type == "windows") {
    Sys.setlocale("LC_TIME", "French_France") 
  } else {
    Sys.setlocale("LC_TIME", "fr_FR")
  }
  
  days_excl <- days_exclude %>%
    mutate(date_up = max(as.Date(inci_dat_long$dates) + 1)) %>%
    group_by(variable) %>%
    mutate(date_low = date_up - number - 0.5) %>%
    ungroup() %>%
    mutate(y_low = 0) %>%
    mutate(y_up = Inf)
  
  # Height and width of the figure
  hgt <- 4
  wth <- 12
  
  # Plotting the time series alone
  if (case_only == F) {
    # Name to save the file
    savename <- sprintf("ts_%s.png", rss_name)
    ylab <- "Nombre par jour"
    showleg <- T
    
    # Initiate plotting instance
    p_inci <- ggplot(data = inci_dat_long %>% filter(variable %in% vars)) +
      facet_wrap(~variable, scale =  "free")
  } else {
    # Creating a dataframe for data exclusion + other parameters
    savename <- sprintf("ts_inspq_%s.png", rss_name)
    ylab <- "Cas par jour"
    showleg <- F
    
    # Only excluded days for cases
    days_excl <- days_excl %>%
      filter(variable == "Cas")
    
    # Initiate plotting instance
    p_inci <- ggplot(data = inci_dat_long %>% filter(variable == "Cas")) 
  }
  p_inci <-  p_inci + 
    theme_minimal(base_family = "Raleway") +
    geom_bar(aes(x = as.Date(dates), y = value, fill = imported, col = imported),
             stat = "identity", show.legend = showleg, size = 0.2, width = 1) +
    geom_rect(data = days_excl,
              aes(xmin = date_low, xmax = date_up, ymin = y_low, ymax = y_up, group = variable),
              fill = "snow3", alpha = 0.6) +
    geom_blank(aes(y = min)) +
    geom_blank(aes(y = max)) +
    scale_x_date(breaks = "2 weeks", date_labels = "%b %d") +
    scale_fill_manual(name = "Évènement", values = c("firebrick3", "firebrick1", "royalblue3")) +
    scale_color_manual(values = c("firebrick4", "firebrick2", "royalblue4",
                                  "cadetblue4", "cornsilk4"), guide=F) +
    labs(x = "Temps (jours)",
         y = ylab) +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.grid.minor = element_blank(),
          legend.position = "right") 
  
  ggsave(savename, path = here::here("output", "fig"),
         device = "png", height = hgt, width = wth, units = "in",
         dpi = 320)
  Sys.setlocale("LC_TIME", old_locale)
  return(p_inci)
}

# ---- Plotting R0 one panel ----
plot_r0_one <- function(df, rss_name, cv = TRUE, case_only = FALSE) {
  
  old_locale <- Sys.getlocale("LC_TIME")
  if (.Platform$OS.type == "windows") {
    Sys.setlocale("LC_TIME", "French_France") 
  } else {
    Sys.setlocale("LC_TIME", "fr_FR")
  }
  
  # conditioning on CV
  if (cv == TRUE) {
    df$RSS_code <- ifelse(is.na(df$RSS_code), 999, df$RSS_code)
    df <- na.omit(df)
    df_c <- df[df$variable == "Cas", ]
    df_cv <- df_c[df_c$dates == max(df_c$dates), ]
    df_h <- df[df$variable == "Hospitalisations", ]
    df_hv <- df_h[df_h$dates == max(df_h$dates), ]
    df_d <- df[df$variable == "Décès", ]
    df_dv <- df_d[df_d$dates == max(df_d$dates), ]  
    
    if (case_only == FALSE){
      # Name to save the figure
      savename <- sprintf("r0_all_%s_%s.png", rss_name, Sys.Date())
      
      # Compute the coefficient of variation to assess whether we include or not
      cv_c <- df_cv$se / df_cv$mean
      cv_h <- df_hv$se / df_hv$mean
      cv_d <- df_dv$se / df_dv$mean
    } else {
      # Name to save the figure
      savename <- sprintf("r0_inspq_%s_%s.png", rss_name, Sys.Date())
      
      # Set coefficient so we only have cases
      cv_c <- df_cv$se / df_cv$mean
      cv_h <- 0.3 + 1
      cv_d <- 0.3 + 1
    }
    
    p <- ggplot() +
      theme_minimal(base_family = "Raleway") 
    
    if (cv_d < 0.3) {
      p <- p + geom_ribbon(data = df_d, aes(x = dates, ymin = lci,
                                            ymax = uci, fill = variable), alpha = 0.4) +
        geom_line(data = df_d, aes(x = dates, y = mean, col = variable), size = 0.75) +
        annotate("segment", x = max(df_d$dates),
                 xend = max(df_d$dates),
                 y = df_dv$mean,
                 yend =  2.25,
                 col = "cornsilk4", lty = 3) +
        annotate("text", x = max(df_d$dates),
                 y = 2.25 + 0.1,
                 label = sprintf("%.2f (%.2f-%.2f)",
                                 round(df_dv$mean, 2),
                                 round(df_dv$lci, 2),
                                 round(df_dv$uci, 2)),
                 col = "cornsilk4", fontface = 2, size = 3.5, hjust = 1) 
    }
    
    if (cv_h < 0.3) {
      p <- p + geom_ribbon(data = df_h, aes(x = dates, ymin = lci,
                                            ymax = uci, fill = variable), alpha = 0.4) +
        annotate("segment", x = max(df_h$dates),
                 xend = max(df_h$dates),
                 y = df_hv$mean,
                 yend = 2.5,
                 col = "cadetblue4", lty = 3) +
        annotate("text", x = max(df_h$dates),
                 y = 2.5 + 0.1,
                 label = sprintf("%.2f (%.2f-%.2f)",
                                 round(df_hv$mean, 2),
                                 round(df_hv$lci, 2),
                                 round(df_hv$uci, 2)),
                 col = "cadetblue4", fontface = 2, size = 3.5, hjust = 1)
    }
    
    # we always add cases
    p <- p + geom_ribbon(data = df_c, aes(x = dates, ymin = lci,
                                          ymax = uci, fill = variable), alpha = 0.4) +
      geom_hline(yintercept = 1, col = "darkgrey", lty = 2) +
      annotate("segment", x = max(df_c$dates),
               xend = max(df_c$dates),
               y = df_cv$mean,
               yend = 2.75,
               col = "firebrick4", lty = 3) +
      annotate("text", x = max(df_c$dates),
               y = 2.75 + 0.1,
               label = sprintf("%.2f (%.2f-%.2f)",
                               round(df_cv$mean, 2),
                               round(df_cv$lci, 2),
                               round(df_cv$uci, 2)),
               col = "firebrick4", fontface = 2, size = 3.5, hjust = 1) 
    
    if (cv_h < 0.3) {
      p <- p + geom_line(data = df_h, aes(x = dates, y = mean, col = variable), size = 0.75)
    }
    if (cv_d < 0.3) {
      p <- p + geom_line(data = df_d, aes(x = dates, y = mean, col = variable), size = 0.75)
    }
    
    p <- p +  geom_line(data = df_c, aes(x = dates, y = mean, col = variable), size = 0.75) +
      labs(x = "Temps (jours)",
           y = expression('Taux de reproduction effectif (R'[t]*')')) +
      scale_x_date(limits = as.Date(c("2020-03-13", max(df$dates) + 1)),
                   breaks = "2 weeks", date_labels = "%b %d") +
      scale_y_continuous(breaks = c(0, 0.5, 0.75, 1, 2, 3)) +
      coord_cartesian(ylim = c(0, 3.5))
    
    if (cv_d < 0.3 & cv_h < 0.3) {
      p <- p +  scale_fill_manual(name = "Évènement", values = c("firebrick3",
                                                                 "cornsilk3",
                                                                 "cadetblue3")) +
        scale_color_manual(name = "Évènement", values = c("firebrick4",
                                                          "cornsilk4",
                                                          "cadetblue4")) 
    }
    if (cv_d >= 0.3 & cv_h < 0.3) {
      p <- p +  scale_fill_manual(name = "Évènement", values = c("firebrick3",
                                                                 "cadetblue3")) +
        scale_color_manual(name = "Évènement", values = c("firebrick4",
                                                          "cadetblue4")) 
    }
    if (cv_d < 0.3 & cv_h >= 0.3) {
      p <- p +  scale_fill_manual(name = "Évènement", values = c("firebrick3",
                                                                 "cornsilk3")) +
        scale_color_manual(name = "Évènement", values = c("firebrick4",
                                                          "cornsilk4")) 
    }
    if (cv_d >= 0.3 & cv_h >= 0.3) {
      p <- p +  scale_fill_manual(name = "Évènement", values = c("firebrick3")) +
        scale_color_manual(name = "Évènement", values = c("firebrick4")) 
    }
    
    p <- p +  theme(strip.text = element_text(face = "bold"),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10),
                    panel.grid.minor = element_blank(),
                    legend.position = "right") 
  } else {
    df$RSS_code <- ifelse(is.na(df$RSS_code), 999, df$RSS_code)
    df <- na.omit(df)
    
    p <- ggplot() +
      theme_minimal(base_family = "Raleway") +
      geom_ribbon(data = df, aes(x = dates, ymin = lci, 
                                 ymax = uci, fill = variable), alpha = 0.4) +
      geom_line(data = df, aes(x = dates, y = mean, col = variable),
                size = 0.75) +
      geom_hline(yintercept = 1, col = "darkgrey", lty = 2) +
      annotate("segment", x = max(df[df$variable == "Décès",]$dates),
               xend = max(df[df$variable == "Décès",]$dates),
               y = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$mean,
               yend = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 0.25,
               col = "cornsilk4", lty = 3) +
      annotate("segment", x = max(df[df$variable == "Hospitalisations",]$dates),
               xend = max(df[df$variable == "Hospitalisations",]$dates),
               y = df[df$variable == "Hospitalisations" & df$dates == max(df[df$variable == "Hospitalisations",]$dates),]$mean,
               yend = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 0.5,
               col = "cadetblue4", lty = 3) +
      annotate("segment", x = max(df[df$variable == "Cas",]$dates),
               xend = max(df[df$variable == "Cas",]$dates),
               y = df[df$variable == "Cas" & df$dates == max(df[df$variable == "Cas",]$dates),]$mean,
               yend = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 1,
               col = "firebrick4", lty = 3) +
      annotate("text", x = max(df[df$variable == "Décès",]$dates),
               y = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 0.25, 
               label = sprintf("%.2f (%.2f-%.2f)",
                               round(df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$mean, 2),
                               round(df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$lci, 2),
                               round(df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci, 2)),
               col = "cornsilk4", fontface = 2, size = 2.5, hjust = 1) +
      annotate("text", x = max(df[df$variable == "Cas",]$dates),
               y = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 1,
               label = sprintf("%.2f (%.2f-%.2f)",
                               round(df[df$variable == "Cas" & df$dates == max(df[df$variable == "Cas",]$dates),]$mean, 2),
                               round(df[df$variable == "Cas" & df$dates == max(df[df$variable == "Cas",]$dates),]$lci, 2),
                               round(df[df$variable == "Cas" & df$dates == max(df[df$variable == "Cas",]$dates),]$uci, 2)),
               col = "firebrick4", fontface = 2, size = 2.5, hjust = 1) +
      annotate("text", x = max(df[df$variable == "Hospitalisations",]$dates),
               y = df[df$variable == "Décès" & df$dates == max(df[df$variable == "Décès",]$dates),]$uci + 0.5,
               label = sprintf("%.2f (%.2f-%.2f)",
                               round(df[df$variable == "Hospitalisations" & df$dates == max(df[df$variable == "Hospitalisations",]$dates),]$mean, 2),
                               round(df[df$variable == "Hospitalisations" & df$dates == max(df[df$variable == "Hospitalisations",]$dates),]$lci, 2),
                               round(df[df$variable == "Hospitalisations" & df$dates == max(df[df$variable == "Hospitalisations",]$dates),]$uci, 2)),
               col = "cadetblue4", fontface = 2, size = 2.5, hjust = 1) +
      labs(x = "Temps (jours)", y = "R(t)") +
      scale_x_date(limits = as.Date(c("2020-03-13", max(df$dates) + 1)),
                   breaks = "2 weeks", date_labels = "%b %d") +
      scale_y_continuous(breaks = c(0, 0.5, 0.75, 1, 2, 3)) +
      coord_cartesian(ylim = c(0, 3.5)) +
      scale_fill_manual(name = "Évènement", values = c("firebrick3",
                                                       "cadetblue3",
                                                       "cornsilk3")) +
      scale_color_manual(name = "Évènement", values = c("firebrick4",
                                                        "cadetblue4",
                                                        "cornsilk4")) +
      theme(strip.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            panel.grid.minor = element_blank(),
            legend.position = "right") 
  }

  
  ggsave(savename,
         path = here::here("output", "fig"),
         device = "png", height = 3.5, width = 12, units = "in",
         dpi = 320)
  Sys.setlocale("LC_TIME", old_locale)
  return(p)
}
