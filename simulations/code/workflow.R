path <- "./simulations/"
path_output <- paste0(path, "output/")
dir.create(path_output, recursive = T)

set.seed(1)
source(paste0(path, "code/steps/01 utils.R"))
source(paste0(path, "code/steps/02 settings.R"))

param_list <- c("m2", "m3", "m4", "m8")

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

# fit model and save output for each parameter
stats_list <- vector(mode = "list", length = length(param_list))
for (p in 1:length(param_list)) {
  param <- param_list[p]
  path_sub <- paste0(path, "archive/", param, "/")
  dir.create(path_sub, recursive = T)

  if (!file.exists(paste0(path_output, param, ".csv"))) {
    cl <- makeCluster(num_part, outfile = "")
    registerDoSNOW(cl)

    # prepare environment time series
    source(paste0(path, "code/steps/11 prepare env ts.R"))

    # get phenology model parameter in each year (match and mismatch)
    source(paste0(path, "code/steps/12 get model param.R"))

    # get phenology time series (match and mismatch)
    source(paste0(path, "code/steps/13 get pheno ts.R"))

    # preprocess data
    source(paste0(path, "code/steps/21 preprocess data.R"))

    # use first half to train model
    source(paste0(path, "code/steps/22 prepare embeddings.R"))
    source(paste0(path, "code/steps/23 train GP model.R"))

    # predict for whole duration
    source(paste0(path, "code/steps/24 predict.R"))

    stopCluster(cl)
  }

  # output table and figure
  source(paste0(path, "code/steps/35 output table and figure.R"))
}

# plot all phenological mismatch time series
ts_df_list <- vector(mode = "list", length = length(param_list))
for (p in 1:length(param_list)) {
  param <- param_list[p]
  ts_df_list[[p]] <- read_csv(paste0(path_output, param, ".csv"),
    col_types = list(
      pheno_mis = col_double(),
      mismatch_actual = col_double(),
      mismatch_model = col_double(),
      mismatch_model_upper = col_double(),
      mismatch_model_lower = col_double()
    )
  ) %>%
    mutate(param = param) %>%
    rowwise() %>%
    mutate(param_v = param_name[[param]])
}
ts_df <- bind_rows(ts_df_list) %>%
  mutate(param_v = factor(param_v, levels = c("Summer-winter difference", "Timing of spring onset", "Slope of curve in spring", "Number of life cycles")))

colors <- c(
  "simulated mismatch" = "purple",
  "estimated mismatch" = "dark red",
  "predictive error" = "dark blue"
)
p_ts <-
  ggplot(ts_df %>% filter(site == 3) %>% filter(date >= as.Date(paste0(midyear + 1, "-01-01")))) +
  geom_line(aes(x = date, y = mismatch_actual, col = "simulated mismatch"), alpha = 0.5) +
  geom_line(aes(x = date, y = mismatch_model, col = "estimated mismatch"), alpha = 0.5) +
  geom_ribbon(aes(x = date, ymin = mismatch_model_lower, ymax = mismatch_model_upper, fill = "estimated mismatch"), alpha = 0.25) +
  theme_classic() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill = F) +
  labs(
    x = "Date",
    y = "Phenological mismatch",
    color = ""
  ) +
  theme(legend.position = "top") +
  facet_wrap(. ~ param_v, ncol = 1, scales = "free_y")
p_ts

# summarize phenological mismatch metrics
mismatch_df <- bind_rows(stats_list) %>%
  mutate(stats = factor(stats, levels = c("corr", "R2", "RMSE", "nRMSE"))) %>%
  rowwise() %>%
  mutate(param_v = param_name[[param]]) %>%
  mutate(param_v = factor(param_v, levels = c("Summer-winter difference", "Timing of spring onset", "Slope of curve in spring", "Number of life cycles") %>% rev()))
mismatch_df %>% filter(stats == "nRMSE")

write_csv(mismatch_df, paste0(path_output, "mismatch metrics.csv"))

# plot phenological mismatch metrics
colors <- c(
  "simulated mismatch" = "purple",
  "estimated mismatch" = "dark red",
  "predictive error" = "dark blue"
)

p_stat <- ggplot(mismatch_df %>% filter(stats == "nRMSE")) +
  geom_point(aes(x = param_v, y = theo_mismatch, col = "simulated mismatch")) +
  geom_point(aes(x = param_v, y = est_mismatch, col = "estimated mismatch")) +
  scale_color_manual(values = colors) +
  labs(
    x = "type of mismatch",
    y = "nRMSE",
    color = ""
  ) +
  ylim(0, 0.2) +
  coord_flip() +
  theme_classic()
p_stat

cairo_pdf(paste0(path_output, "simulated and estimated mismatch.pdf"), height = 8, width = 8)
grid.arrange(
  annotate_figure(p_ts, fig.lab = "A"),
  annotate_figure(p_stat, fig.lab = "B"),
  heights = c(0.75, 0.25), ncol = 1
) %>%
  print()
dev.off()
