path_scaling <- paste0(path_sub, "scaling/")
dir.create(path_scaling, recursive = T)

# convert data and variance into 3d arrays
x <- array(NA,
  dim = c(nrow(coord_df), length(date_list), length(var_list)),
  dimnames = list(
    as.character(1:nrow(coord_df)),
    as.character(date_list),
    var_list
  )
)
for (j in 1:length(date_list)) { # date
  for (v in 1:length(var_list)) { # covariate
    if (var_list[v] %in% c("pheno", "temp", "prcp")) {
      ts_date <- ts %>% filter(date == date_list[j])
      if (nrow(ts_date) > 0) {
        x[, j, v] <- ts_date[var_list[v]] %>%
          unlist() %>%
          as.numeric()
      } else {
        x[, j, v] <- rep(NA, nrow(coord_df))
      }
    } else if (var_list[v] == "doy") {
      x[, j, v] <- rep(sin(as.numeric(format(date_list[j], "%j")) * 2 * pi), nrow(coord_df))
    }
  }
  print(date_list[j])
}

Sigma <- array(NA,
  dim = c(nrow(coord_df), length(date_list), length(var_list)),
  dimnames = list(
    as.character(1:nrow(coord_df)),
    as.character(date_list),
    var_list
  )
)
for (j in 1:length(date_list)) { # date
  for (v in 1:length(var_list)) { # covariate
    if (var_list[v] %in% c("pheno")) {
      ts_date <- ts %>% filter(date == date_list[j])
      if (nrow(ts_date) > 0) {
        Sigma[, j, v] <- (ts_date[paste0(var_list[v], "_sd")] %>% unlist() %>% as.numeric())^2
      } else {
        Sigma[, j, v] <- rep(NA, nrow(coord_df))
      }
    } else {
      Sigma[, j, v] <- rep(0, nrow(coord_df))
    }
  }
  print(date_list[j])
}

# scale data to be roughly between 0 and 1
df_upper_lower <- vector(mode = "list")
for (j in 1:length(var_list)) {
  if (var_list[j] %in% c("pheno")) { # scale by each site
    df_upper_lower[[j]] <- data.frame(x[, , j, drop = F]) %>%
      mutate(site = row_number()) %>%
      gather(key = "date", value = "value", -site) %>%
      drop_na() %>%
      group_by(site) %>%
      dplyr::summarize(
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975)
      ) %>%
      mutate(range = upper - lower)
  } else { # scale for all sites
    all_upper_lower <- data.frame(x[, , j, drop = F]) %>%
      mutate(site = row_number()) %>%
      gather(key = "date", value = "value", -site) %>%
      drop_na() %>%
      dplyr::summarize(
        lower = quantile(value, 0),
        upper = quantile(value, 1)
      ) %>%
      mutate(range = upper - lower)
    df_upper_lower[[j]] <- data.frame(x[, , j, drop = F]) %>%
      mutate(site = row_number()) %>%
      gather(key = "date", value = "value", -site) %>%
      drop_na() %>%
      distinct(site) %>%
      mutate(
        lower = all_upper_lower$lower,
        upper = all_upper_lower$upper,
        range = all_upper_lower$range
      )
  }

  lower <- matrix(df_upper_lower[[j]]$lower) %*% matrix(1, nrow = 1, ncol = ncol(x[, , j, drop = F]))
  range <- matrix(df_upper_lower[[j]]$range) %*% matrix(1, nrow = 1, ncol = ncol(x[, , j, drop = F]))

  x[, , j] <- (x[, , j] - lower) / range - 0.5
}

for (j in 1:length(var_list)) {
  Sigma[, , j] <- Sigma[, , j, drop = F] / (df_upper_lower[[j]]$range)^2
}

for (j in 1:length(var_list)) {
  write_csv(df_upper_lower[[j]], paste0(path_scaling, j, ".csv"))
  print(j)
}

# linear interpolation
for (j in 1:length(var_list)) {
  for (i in 1:nrow(coord_df)) {
    min_id <- min(which(!is.na(x[i, , j])))
    max_id <- max(which(!is.na(x[i, , j])))
    x[i, min_id:max_id, j] <- zoo::na.approx(object = x[i, min_id:max_id, j], x = as.Date(names(x[i, min_id:max_id, j])), maxgap = 14)
  }
}

for (j in 1:length(var_list)) {
  for (i in 1:nrow(coord_df)) {
    min_id <- min(which(!is.na(Sigma[i, , j])))
    max_id <- max(which(!is.na(Sigma[i, , j])))
    Sigma[i, min_id:max_id, j] <- zoo::na.approx(object = Sigma[i, min_id:max_id, j], x = as.Date(names(Sigma[i, min_id:max_id, j])), maxgap = 14)
  }
}

x_raw <- x
Simga_raw <- Sigma

# whittaker smoothing
for (j in 1:length(var_list)) {
  for (i in 1:nrow(coord_df)) {
    max_id <- 0
    done <- F
    while (!done) {
      min_id <- min(which(!is.na(x[i, (max_id + 1):length(date_list), j]))) + (max_id)
      if (min_id == Inf) {
        done <- T
      } else {
        max_id <- min(which(is.na(x[i, min_id:length(date_list), j]))) - 1 + (min_id - 1)
        if (max_id == Inf) {
          max_id <- length(date_list)
          done <- T
        }
        x[i, min_id:max_id, j] <- ptw::whit1(x[i, min_id:max_id, j], 5)
      }
    }
  }
}
