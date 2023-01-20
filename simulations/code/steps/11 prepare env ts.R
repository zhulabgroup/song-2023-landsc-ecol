date_list <- seq(as.Date("2021-01-01"), as.Date("2040-12-31"), by = 1)
midyear <- floor(mean(c(as.numeric(format(min(date_list), "%Y")), as.numeric(format(max(date_list), "%Y")))))
coord_df <- data.frame(lon = 0, lat = 1:5)

distMat <- matrix(
  c(
    0, 1, 2, 3, 4,
    1, 0, 1, 2, 3,
    2, 1, 0, 1, 2,
    3, 2, 1, 0, 1,
    4, 3, 2, 1, 0
  ),
  nrow = 5
)

env_ts_list <-
  foreach(
    s = 1:nrow(coord_df),
    .packages = c("lubridate", "tidyverse")
  ) %dopar% {
    set.seed(42)
    temp <- rep(NA, length(date_list))
    for (i in 1:length(date_list)) {
      date <- date_list[i]
      temp[i] <- waves(
        t = date,
        t_start = date_list[1],
        intercept = 0.3 + 0.2 * coord_df$lat[s],
        slope = 0.0001,
        amplitude1 = 0.8,
        phase1 = 0,
        period1 = 1,
        amplitude2 = 0.5,
        phase2 = 0,
        period2 = 5,
        sd = temp_sd
      )
    }
    ts_site <- data.frame(date = date_list, temp = temp, temp_sd = temp_sd) %>%
      mutate(year = as.numeric(format(date, "%Y"))) %>%
      mutate(site = s)
    ts_site
  }
env_ts <- bind_rows(env_ts_list)

p_temp <-
  ggplot(env_ts) +
  geom_line(aes(x = date, y = temp, group = site, col = site)) +
  theme_classic() +
  theme(legend.position = "right") +
  ylab("Temperature (°C)") +
  xlab("Time (date)") +
  labs(col = "Site")

if (!file.exists("output/temp.pdf")) {
  cairo_pdf(paste0(path, "output/temp.pdf"), width = 11, height = 4)
  print(p_temp)
  dev.off()
}
