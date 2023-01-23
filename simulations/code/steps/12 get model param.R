# names used in figures
param_name <- list(
  "m2" = "Summer-winter difference",
  "m3" = "Timing of spring onset",
  "m4" = "Slope of curve in spring",
  "m8" = "Number of life cycles"
)
pheno_name <- list(
  "m2" = "Net primary productivity",
  "m3" = "Bird breeding activity",
  "m4" = "Enhanced vegetation index",
  "m8" = "Adult insect abundance"
)
env_name <- list(
  "m2" = expression(italic(T)["1:90"]),
  "m3" = expression(italic(T)["-90:-1"]),
  "m4" = expression(italic(T)["1:14"]),
  "m8" = expression(italic(T)["1:365"])
)

param_all <- vector(mode = "list")
for (s in 1:nrow(coord_df)) {
  if (param == "m2") {
    param_site <- env_ts %>%
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
    param_site <- env_ts %>%
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
    param_site <- env_ts %>%
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
    param_site <- env_ts %>%
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
  param_all[[s]] <- param_site
}
param_all <- bind_rows(param_all)
