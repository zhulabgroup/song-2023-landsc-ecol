# load library
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(maps)
library(maptools)
library(raster)
library(parallel)
library(doSNOW)

# read and process data
# source: Park et al. (2018)
# https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF309
# pfl: peak flower
# ffr: first fruit
# pfr: peak fruit
# efl: first flower
data_df <- read_csv("./empirical/herbarium/data/hf309-01-crowdcurio.csv") %>%
  filter(metric == "mean") %>%
  filter(phenophase == "efl" | phenophase == "pfl") %>%
  filter(!(state %in% c("Florida"))) %>% # removed data from Florida following Park et al. (2018)
  dplyr::select(
    species = name,
    year,
    fl = doy,
    lat = county.lat,
    lon = county.lon,
    temp = spring,
    maxp,
    link
  ) %>%
  filter(maxp > 0) %>% # select for record with the greatest crowdsourcer consistency (maxp)
  arrange(desc(maxp)) %>%
  distinct(link, .keep_all = T) %>%
  dplyr::select(-maxp, -link) %>%
  mutate(period = case_when(
    year >= 1950 ~ "late",
    TRUE ~ "early"
  )) %>%
  mutate(species = str_replace(species, "_", " "))

# select species
length(data_df %>% pull(species) %>% unique()) # 30 species in total

sp_list <- data_df %>%
  group_by(species, period) %>%
  summarise(n = n()) %>%
  spread(key = "period", value = "n") %>%
  filter(early >= 30 & late >= 30) %>% # at least 30 records in both early and late periods
  pull(species)
length(sp_list) # 19 species

data_df <- data_df %>% filter(species %in% sp_list)

# settings
sp_vis <- "Aquilegia canadensis" # species to plot
userandom <- F # generate results using random grouping (for SI)

if (userandom) {
  set.seed(1)
  data_df <- data_df %>%
    group_by(species) %>%
    mutate(random = rbernoulli(n())) %>%
    mutate(random = case_when(
      random == 1 ~ "in-sample",
      TRUE ~ "out-of-sample"
    )) %>%
    ungroup()
}

# map records
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))

area <- extent(min(data_df$lon) - 0.5, max(data_df$lon) + 0.5, min(data_df$lat) - 0.5, max(data_df$lat) + 0.5)
usa_crop <- crop(usa, area)

p_map <- ggplot() +
  geom_polygon(
    data = usa_crop, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "lightblue", linewidth = .1
  ) +
  geom_point(data = data_df, aes(x = lon, y = lat, col = species), alpha = 0.3, pch = 1) +
  geom_point(data = data_df %>% filter(species == sp_vis), aes(x = lon, y = lat, col = species), alpha = 1, shape = 19) +
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(col = "none") +
  coord_equal()
p_map

# summarize niche
data_df %>%
  group_by(species) %>%
  summarise(
    temp = median(temp),
    fl = median(fl)
  ) %>%
  summary()

# correlation between temperature and phenology
data_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(fl ~ temp, .))) %>%
  filter(term %in% c("temp")) %>%
  dplyr::select(-statistic) %>%
  filter(p.value < 0.05) %>%
  summary()

# change in functional relationship
if (!userandom) {
  data_df %>%
    group_by(species, period) %>%
    do(broom::tidy(lm(fl ~ temp, .))) %>%
    ungroup() %>%
    dplyr::select(species, period, term, estimate) %>%
    spread(key = "period", value = "estimate") %>%
    mutate(diff = late - early) %>%
    group_by(term) %>%
    summarise(
      median = median(diff),
      lower = quantile(diff, 0.025),
      upper = quantile(diff, 0.975)
    )
} else {
  data_df %>%
    group_by(species, random) %>%
    do(broom::tidy(lm(fl ~ temp, .))) %>%
    ungroup() %>%
    dplyr::select(species, random, term, estimate) %>%
    spread(key = "random", value = "estimate") %>%
    mutate(diff = `out-of-sample` - `in-sample`) %>%
    group_by(term) %>%
    summarise(
      median = median(diff),
      lower = quantile(diff, 0.025),
      upper = quantile(diff, 0.975)
    )
}

if (!userandom) {
  ggplot(data_df) +
    geom_point(aes(x = temp, y = fl, col = period), alpha = 0.3) +
    geom_smooth(aes(x = temp, y = fl, group = period, col = period), method = "lm") +
    theme_classic() +
    facet_wrap(. ~ species) +
    theme(legend.position = "bottom") +
    labs(color = "Time period") +
    xlab("Spring mean temperature (°C)") +
    ylab("Flowering time (day of year)")
} else {
  ggplot(data_df) +
    geom_point(aes(x = temp, y = fl, col = random), alpha = 0.3) +
    geom_smooth(aes(x = temp, y = fl, group = random, col = random), method = "lm") +
    theme_classic() +
    facet_wrap(. ~ species) +
    theme(legend.position = "bottom") +
    labs(color = "Group") +
    xlab("Spring mean temperature (°C)") +
    ylab("Flowering time (day of year)")
}

if (!userandom) {
  p_func <- ggplot(data_df %>% filter(species == sp_vis)) +
    geom_point(aes(x = temp, y = fl, col = period), alpha = 0.3) +
    geom_smooth(aes(x = temp, y = fl, group = period, col = period), method = "lm") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "Time period") +
    xlab("Spring mean temperature (°C)") +
    ylab("Flowering time (day of year)")
} else {
  p_func <- ggplot(data_df %>% filter(species == sp_vis)) +
    geom_point(aes(x = temp, y = fl, col = random), alpha = 0.3) +
    geom_smooth(aes(x = temp, y = fl, group = random, col = random), method = "lm") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "Group") +
    xlab("Spring mean temperature (°C)") +
    ylab("Flowering time (day of year)")
}
p_func

# calculate mismatch
if (!userandom) {
  outfile <- "./empirical/herbarium/data/mismatch.rds"
} else {
  outfile <- "./empirical/herbarium/data/mismatch_SI.rds"
}

if (!file.exists(outfile)) {
  cl <- makeCluster(20, outfile = "")
  registerDoSNOW(cl)
  data_df_mis_list <-
    foreach(
      sp = sp_list,
      .packages = c("tidyverse", "spBayes", "gstat", "raster")
    ) %dopar% {
      set.seed(42)

      # filter for species of interest
      data_sp <- data_df %>%
        filter(species == sp) %>%
        drop_na() 

      # use early period as training data
      if (!userandom) {
        data_train <- data_sp %>%
          filter(period == "early")
      } else {
        data_train <- data_sp %>%
          filter(random == "in-sample")
      }

      # prepare to fit Gaussian spatial model
      z <- resid(lm(fl ~ temp, data_train)) # residuals of non-spatial model
      data_train$z <- z
      lm.vgm <- variogram(z ~ 1, data = data_train %>% `coordinates<-`(data_train[, c("x", "y")])) # get variogram
      # plot(lm.vgm)
      lm.fit <- fit.variogram(lm.vgm, model = vgm("Exp"), fit.method = 6, fit.sills = F) # fit variogram with an exponential model without nugget
      d <- lm.fit[2, ]$range # extract effective spatial range (d) from fitted variogram

      # set priors
      priors <- list(
        "beta.Norm" = list(rep(0, 2), diag(100, 2)), # loose normal priors on beta
        "phi.Unif" = c(-log(0.05) / (d * 100), -log(0.05) / (d / 100)), # loose uniform prior on phi (spatial decay parameter)
        "sigma.sq.IG" = c(2, 2), # loose inverse gamma prior on spatial variance parameter (sigma.sq) (shape and scale for IG)
        "tau.sq.IG" = c(2, 0.1) # tight inverse gamma prior on residual error variance (tausq) (shape and scale for IG)
      )

      # set starting and tuning values
      starting <- list(
        "phi" = -log(0.05) / d,
        "sigma.sq" = 50,
        "tau.sq" = 1
      )
      tuning <- list(
        "phi" = (log(0.05) / (d * 100) - log(0.05) / (d / 100)) / 10,
        "sigma.sq" = 0.01,
        "tau.sq" = 0.01
      )

      # fit spatial model
      splm <- spLM(fl ~ temp,
        data = data_train,
        coords = data_train[, c("x", "y")] %>% as.matrix(),
        cov.model = "exponential",
        priors = priors,
        tuning = tuning,
        starting = starting,
        n.samples = 10000,
        n.report = 1000
      )

      # prediction
      splm.pred <- spPredict(splm,
        pred.covars = cbind(1, data_sp[, "temp"]) %>% as.matrix(), pred.coords = data_sp[, c("x", "y")] %>% as.matrix(),
        start = 5001,
        thin = 10
      )
      quant <- function(x) {
        quantile(x, prob = c(0.025, 0.5, 0.975))
      }
      y.hat <- apply(splm.pred$p.y.predictive.samples, 1, quant) %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`(c("lower", "predict", "upper"))

      data_sp %>%
        bind_cols(y.hat)
    }

  data_mis <- bind_rows(data_df_mis_list) %>%
    mutate(resid = fl - predict)
  write_rds(data_mis, outfile)
  stopCluster(cl)
} else {
  data_mis <- read_rds(outfile)
}

# compare predictive skills in early and late periods
if (!userandom) {
  data_mis %>%
    group_by(species, period) %>%
    summarise(
      R2 = cor(predict, fl)^2,
      rmse = Metrics::rmse(predict, fl)
    ) %>%
    gather(key = "stat", value = "value", -species, -period) %>%
    spread(key = "period", value = "value") %>%
    mutate(diff = late - early) %>%
    gather(key = "report", value = "value", -species, -stat) %>%
    group_by(stat, report) %>%
    summarise(
      median = quantile(value, 0.5),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    )
} else {
  data_mis %>%
    group_by(species, random) %>%
    summarise(
      R2 = cor(predict, fl)^2,
      rmse = Metrics::rmse(predict, fl)
    ) %>%
    gather(key = "stat", value = "value", -species, -random) %>%
    spread(key = "random", value = "value") %>%
    mutate(diff = `out-of-sample` - `in-sample`) %>%
    gather(key = "report", value = "value", -species, -stat) %>%
    group_by(stat, report) %>%
    summarise(
      median = quantile(value, 0.5),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    )
}

if (!userandom) {
  p_pred <- ggplot(data_mis %>% filter(species == sp_vis)) +
    geom_point(aes(x = predict, y = fl), alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    geom_text(
      data = data_mis %>%
        filter(species == sp_vis) %>%
        group_by(period) %>%
        summarise(
          median = quantile(resid, 0.5),
          lower = quantile(resid, 0.025),
          upper = quantile(resid, 0.975)
        ),
      aes(
        x = 150, y = 250,
        label = paste0("Ymis=", round(median, 2), "(", round(lower, 2), ", ", round(upper, 2), ")")
      ),
      col = "blue"
    ) +
    theme_classic() +
    facet_wrap(. ~ period, nrow = 1) +
    guides(col = "none") +
    xlab("Predicted flowering time (day of year)") +
    ylab("Observed flowering time (day of year)")
} else {
  p_pred <- ggplot(data_mis %>% filter(species == sp_vis)) +
    geom_point(aes(x = predict, y = fl), alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    geom_text(
      data = data_mis %>%
        filter(species == sp_vis) %>%
        group_by(random) %>%
        summarise(
          median = quantile(resid, 0.5),
          lower = quantile(resid, 0.025),
          upper = quantile(resid, 0.975)
        ),
      aes(
        x = 150, y = 250,
        label = paste0("Ymis=", round(median, 2), "(", round(lower, 2), ", ", round(upper, 2), ")")
      ),
      col = "blue"
    ) +
    theme_classic() +
    facet_wrap(. ~ random, nrow = 1) +
    guides(col = "none") +
    xlab("Predicted flowering time (day of year)") +
    ylab("Observed flowering time (day of year)")
}
p_pred

# summarize and test phenological mismatch
if (!userandom) {
  t_df <- data_mis %>%
    filter(period == "late") %>%
    group_by(species) %>%
    summarise(
      median = quantile(resid, 0.5),
      lower = quantile(resid, 0.025),
      upper = quantile(resid, 0.975),
      p = t.test(x = resid)$p.value,
      estimate = t.test(x = resid)$estimate,
      med_fl = median(fl)
    )
} else {
  t_df <- data_mis %>%
    filter(random == "out-of-sample") %>%
    group_by(species) %>%
    summarise(
      median = quantile(resid, 0.5),
      lower = quantile(resid, 0.025),
      upper = quantile(resid, 0.975),
      p = t.test(x = resid)$p.value,
      estimate = t.test(x = resid)$estimate,
      med_fl = median(fl)
    )
}

t_df %>%
  filter(p < 0.05) %>%
  group_by(median < 0) %>%
  summarise(
    n = n(),
    lower = min(median),
    upper = max(median)
  )

if (!userandom) {
  p_mis <- ggplot() +
    geom_violin(
      data = data_mis %>% filter(period == "late"),
      aes(x = reorder(species, desc(species)), y = resid, col = species),
      draw_quantiles = 0.5
    ) +
    geom_violin(
      data = data_mis %>% filter(period == "late") %>%
        left_join(t_df, by = "species") %>% filter(p < 0.05),
      aes(x = reorder(species, desc(species)), y = resid, fill = species),
      draw_quantiles = 0.5
    ) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    guides(col = "none", fill = "none") +
    coord_flip() +
    ylab("Deviation of observed flowering time \n from predicted flowering time (day)") +
    xlab("Species") +
    theme(axis.text.y = element_text(face = "italic"))
} else {
  p_mis <- ggplot() +
    geom_violin(
      data = data_mis %>% filter(random == "out-of-sample"),
      aes(x = reorder(species, desc(species)), y = resid, col = species),
      draw_quantiles = 0.5
    ) +
    geom_violin(
      data = data_mis %>% filter(random == "out-of-sample") %>%
        left_join(t_df, by = "species") %>% filter(p < 0.05),
      aes(x = reorder(species, desc(species)), y = resid, fill = species),
      draw_quantiles = 0.5
    ) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    guides(col = "none", fill = "none") +
    coord_flip() +
    ylab("Deviation of observed flowering time \n from predicted flowering time (day)") +
    xlab("Species") +
    theme(axis.text.y = element_text(face = "italic"))
}
p_mis

if (!userandom) {
  data_mis %>%
    filter(period == "late") %>%
    group_by(species) %>%
    summarise(median = median(resid)) %>% # median of mismatch among all species
    pull(median) %>%
    t.test() # test if the overall mismatch differs from 0
} else {
  data_mis %>%
    filter(random == "out-of-sample") %>%
    group_by(species) %>%
    summarise(median = median(resid)) %>% # median of mismatch among all species
    pull(median) %>%
    t.test() # test if the overall mismatch differs from 0
}

# make figure for manuscript
if (!userandom) {
  figname <- "./empirical/herbarium/figure/figure.pdf"
} else {
  figname <- "./empirical/herbarium/figure/figure_SI.pdf"
}
cairo_pdf(figname, width = 10, height = 10)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
  annotate_figure(p_func, fig.lab = "b"),
  annotate_figure(p_pred, fig.lab = "c"),
  annotate_figure(p_mis, fig.lab = "d"),
  layout_matrix = rbind(
    c(1, 3),
    c(2, 4)
  ),
  widths = 2:3,
  heights = 2:3
)
dev.off()
