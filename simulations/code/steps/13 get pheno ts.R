pheno_ts_list <-
  foreach(s = 1:nrow(coord_df)) %dopar% {
    set.seed(42)

    param_site <- param_all %>% filter(site == s)

    pheno_list <- vector(mode = "list", length = nrow(param_site))
    for (y in 1:nrow(param_site)) {
      date_list_year <- ts_all %>%
        filter(site == s) %>%
        filter(year == param_site$year[y]) %>%
        dplyr::select(date) %>%
        mutate(date = as.character(date)) %>%
        unlist() %>%
        as.Date()

      pheno <- pheno_mis <- rep(NA, length(date_list_year))
      for (i in 1:length(date_list_year)) {
        date <- date_list_year[i]
        if (param == "m2") {
          pheno[i] <- double_logistics(date,
            m2 = param_site$param[y]
          )
          pheno_mis[i] <- double_logistics(date,
            m2 = param_site$param_mis[y]
          )
        }
        if (param == "m3") {
          pheno[i] <- double_logistics(date,
            m3 = param_site$param[y]
          )
          pheno_mis[i] <- double_logistics(date,
            m3 = param_site$param_mis[y]
          )
        }
        if (param == "m4") {
          pheno[i] <- double_logistics(date,
            m4 = param_site$param[y]
          )
          pheno_mis[i] <- double_logistics(date,
            m4 = param_site$param_mis[y]
          )
        }
        if (param == "m8") {
          pheno[i] <- double_logistics(date,
            m8 = param_site$param[y]
          )
          pheno_mis[i] <- double_logistics(date,
            m8 = param_site$param_mis[y]
          )
        }
      }
      pheno_year <- data.frame(date = date_list_year, pheno = pheno, pheno_mis = pheno_mis, pheno_sd = pheno_sd)
      pheno_list[[y]] <- pheno_year
    }
    pheno_site <- bind_rows(pheno_list) %>%
      mutate(site = s)
    pheno_site
  }
pheno_ts <- bind_rows(pheno_ts_list)

ts <- env_ts %>%
  left_join(pheno_ts, by = c("date", "site"))
