param_all <-
  foreach(s = 1:nrow(coord_df)) %dopar% {
    if (param == "m2") {
      param_site <- ts_all %>%
        filter(site == s) %>%
        dplyr::select(temp, year) %>%
        group_by(year) %>%
        slice_head(n = 90) %>%
        summarize(temp_summ = mean(temp)) %>%
        ungroup() %>%
        mutate(param = env_to_param(
          env = temp_summ,
          lower = 0.2,
          upper = 1,
          steepness = -3,
          midpoint = 1.5
        )) %>%
        mutate(param_mis = case_when((year > midyear) ~ param + 0.1)) %>%
        mutate(site = s)
    }
    if (param == "m3") {
      param_site <- ts_all %>%
        filter(site == s) %>%
        dplyr::select(temp, year) %>%
        group_by(year) %>%
        slice_tail(n = 90) %>% # last 90 days
        summarize(temp_summ = mean(temp)) %>%
        ungroup() %>%
        mutate(temp_summ = lag(temp_summ)) %>%
        mutate(param = env_to_param(
          env = temp_summ,
          lower = 40,
          upper = 120,
          steepness = -2,
          midpoint = 1
        )) %>%
        mutate(param_mis = case_when((year > midyear) ~ param + 20)) %>%
        mutate(site = s)
    }
    if (param == "m4") {
      param_site <- ts_all %>%
        filter(site == s) %>%
        dplyr::select(temp, year) %>%
        group_by(year) %>%
        slice_head(n = 14) %>%
        summarize(temp_summ = mean(temp)) %>%
        ungroup() %>%
        mutate(param = env_to_param(
          env = temp_summ,
          lower = 10,
          upper = 80,
          steepness = 3.5,
          midpoint = 1.2
        )) %>%
        mutate(param_mis = case_when((year > midyear) ~ param / 2)) %>%
        mutate(site = s)
    }
    if (param == "m8") {
      param_site <- ts_all %>%
        filter(site == s) %>%
        dplyr::select(temp, year) %>%
        group_by(year) %>%
        summarize(temp_summ = mean(temp)) %>%
        ungroup() %>%
        mutate(param = env_to_param(
          env = temp_summ,
          lower = 1,
          upper = 1.5,
          steepness = 2,
          midpoint = 1
        )) %>%
        mutate(param_mis = case_when((year > midyear) ~ param - 0.1)) %>%
        mutate(site = s)
    }
    param_site
  }
param_all <- bind_rows(param_all)
