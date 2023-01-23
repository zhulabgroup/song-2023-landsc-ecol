# read prediction data
combine_df_ori <- read_csv(paste0(path_output, param, ".csv"), col_types = str_c(c("D", rep("d", 15)), collapse = ""))
combine_df_ori_fit <- combine_df_ori %>% filter(year <= midyear)
combine_df_ori_fore <- combine_df_ori %>% filter(year > midyear)

# retrieve scaling info
path_scaling <- paste0(path_sub, "scaling/")
df_upper_lower <- vector(mode = "list")
for (j in 1:length(var_list)) {
  df_upper_lower[[j]] <- read_csv(paste0(path_scaling, j, ".csv"))
}

# calculate statistics for predictive skills
theo_mismatch <- data.frame(compare_stats(obs_ori = combine_df_ori_fore$pheno, pred_ori = combine_df_ori_fore$pheno_mis, range = df_upper_lower[[1]]$range[1])) %>%
  mutate(cat = "theo_mismatch")
est_mismatch <- data.frame(compare_stats(obs_ori = combine_df_ori_fore$value, pred_ori = combine_df_ori_fore$pheno_mis, range = df_upper_lower[[1]]$range[1])) %>%
  mutate(cat = "est_mismatch")
model_predskill <- data.frame(compare_stats(obs_ori = combine_df_ori_fore$value, pred_ori = combine_df_ori_fore$pheno, range = df_upper_lower[[1]]$range[1])) %>%
  mutate(cat = "model_predskill")

# store in list
stats_list[[p]] <- bind_rows(theo_mismatch, est_mismatch, model_predskill) %>%
  mutate(cat = factor(cat, levels = c("theo_mismatch", "est_mismatch", "model_predskill"))) %>%
  gather(key = "stats", value = "value", -cat) %>%
  spread(key = "cat", value = "value") %>%
  mutate(stats = factor(stats, levels = c("corr", "R2", "RMSE", "nRMSE"))) %>%
  arrange(stats) %>%
  mutate(param = param) %>%
  rowwise() %>%
  mutate(param_v = param_name[[param]])

# plot environment-parameter relationship
p1 <-
  ggplot(param_all) +
  geom_line(aes(x = temp_summ, y = param), col = "blue") +
  theme_classic() +
  labs(
    x = env_name[[param]],
    y = param_name[[param]]
  )

# plot simulates and predicted phenology
pheno_doy <- combine_df_ori %>%
  dplyr::select(date, site, simulated = pheno, predicted = value) %>%
  gather(key = "cat", value = "value", -date, -site) %>%
  mutate(year = as.numeric(format(date, "%Y"))) %>%
  mutate(doy = as.numeric(format(date, "%j"))) %>%
  left_join(param_all %>% dplyr::select(-param_mis), by = c("year", "site")) %>%
  mutate(cat = factor(cat, levels = c("simulated", "predicted"))) %>%
  arrange(cat)

p2 <-
  ggplot(pheno_doy %>%
    filter(site == 3)) +
  geom_line(aes(x = doy, y = value, group = year, col = temp_summ), alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Time (day of year)",
    y = pheno_name[[param]],
    color = env_name[[param]]
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ cat, ncol = 1) +
  scale_color_viridis_c(direction = -1)

# plot parameter change over time
colors <- c(
  "simulated" = "blue",
  "mismatched" = "red"
)
p3 <-
  ggplot(data = param_all %>% filter(site == 3)) +
  geom_line(aes(x = year, y = param, group = site, col = "simulated"), linewidth = 2, alpha = 0.5) +
  geom_line(aes(x = year, y = param_mis, group = site, col = "mismatched"), linewidth = 2, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors) +
  labs(
    x = "Time (year)",
    y = param_name[[param]],
    color = ""
  ) +
  theme(legend.position = "top")

# plot predicted and simulated phenology time series
colors <- c(
  "simulated phenology" = "blue",
  "simulated phenology with mismatch" = "red",
  "predicted phenology" = "black"
)
p4 <-
  ggplot(combine_df_ori %>% filter(site == 3) %>% filter(date >= as.Date(paste0(midyear + 1, "-01-01")))) +
  geom_line(aes(x = date, y = pheno, col = "simulated phenology"), alpha = 0.5) +
  geom_line(aes(x = date, y = pheno_mis, col = "simulated phenology with mismatch"), alpha = 0.5) +
  geom_line(aes(x = date, y = value, col = "predicted phenology")) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = "predicted phenology"), alpha = 0.25) +
  theme_classic() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill = "none") +
  labs(
    x = "Time (date)",
    y = pheno_name[[param]],
    color = ""
  ) +
  theme(legend.position = "top")

# plot predicted and simulated phenological mismatch
colors <- c(
  "simulated mismatch" = "purple",
  "estimated mismatch" = "dark red",
  "predictive error" = "dark blue"
)
p5 <-
  ggplot(combine_df_ori %>% filter(site == 3) %>% filter(date >= as.Date(paste0(midyear + 1, "-01-01")))) +
  geom_line(aes(x = date, y = mismatch_actual, col = "simulated mismatch"), alpha = 0.5) +
  geom_line(aes(x = date, y = mismatch_model, col = "estimated mismatch"), alpha = 0.5) +
  geom_ribbon(aes(x = date, ymin = mismatch_model_lower, ymax = mismatch_model_upper, fill = "estimated mismatch"), alpha = 0.25) +
  theme_classic() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill = "none") +
  labs(
    x = "Time (date)",
    y = "Phenological mismatch",
    color = ""
  ) +
  theme(legend.position = "top")

# assemble figure for manuscript
cairo_pdf(paste0(path_output, param, ".pdf"), height = 8.5, width = 11)
grid.arrange(annotate_figure(p1, fig.lab = "a", fig.lab.pos = "top.left"),
  annotate_figure(p2, fig.lab = "c", fig.lab.pos = "top.left"),
  annotate_figure(p3, fig.lab = "b", fig.lab.pos = "top.left"),
  annotate_figure(p4, fig.lab = "d", fig.lab.pos = "top.left"),
  annotate_figure(p5, fig.lab = "e", fig.lab.pos = "top.left"),
  layout_matrix = rbind(
    c(1, 3, 3, 3),
    c(2, 4, 4, 4),
    c(2, 5, 5, 5)
  )
)
dev.off()
