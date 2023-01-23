# load library
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(maps)
library(maptools)
library(sf)
library(raster)
library(parallel)
library(doSNOW)

# read and process data
# source: Hällfors et al. (2020)
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.wstqjq2ht
bird_df <- read_csv("./empirical/bird/data/73_species.csv") %>%
  dplyr::select(
    nestid = NestID,
    X = XEUREF,
    Y = YEUREF,
    species_short = Species,
    year = Year,
    bh = Dayofyear,
    area = BZ
  ) %>%
  mutate(area = fct_relevel(area, levels = c("HB", "SB", "MB", "NB")))

# join with trait data
sp_df <- read_csv("./empirical/bird/data/Traits_73_species.csv") %>%
  dplyr::select(
    species_short = Abbreviation,
    species = `Scientific name`,
    mig = Mig
  )

bird_df <- bird_df %>%
  left_join(sp_df, by = "species_short")

# project to longitude latitude
coord_df <- bird_df %>%
  distinct(X, Y)

coord_sp <- SpatialPointsDataFrame(
  coords = coord_df[, c("X", "Y")], data = coord_df,
  proj4string = CRS("EPSG:3067")
)

coord_sp_reproj <- spTransform(coord_sp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coord_df_reproj <- as.data.frame(coord_sp_reproj) %>%
  rename(lon = X.1, lat = Y.1)

bird_df <- bird_df %>%
  left_join(coord_df_reproj, by = c("X", "Y"))

# read mean annual temperature data at all nest locations
# previously extracted from TerraClimate dataset
climfile <- "./empirical/bird/data/mat.rds"
mat_df <- read_rds(climfile)

# join bird data with climate data
data_df <- bird_df %>%
  dplyr::select(nestid, X, Y, year, bh, area, species, mig, lon, lat) %>%
  spread(key = "year", value = "bh") %>%
  gather(key = "year", value = "bh", -nestid, -X, -Y, -area, -species, -mig, -lon, -lat) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(period = case_when(
    year >= 1995 ~ "late",
    TRUE ~ "early"
  )) %>%
  left_join(mat_df, by = c("lat", "lon", "year"))

# map records
world <- map("world", fill = TRUE)
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world <- map2SpatialPolygons(world, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
finland <- world["Finland"]
finland_andmore <- world[c("Finland", "Russia", "Norway", "Sweden")]
finland_reproj <- spTransform(finland, CRS("EPSG:3067"))
finland_andmore_reproj <- spTransform(finland_andmore, CRS("EPSG:3067"))

area <- extent(min(coord_df$X) - 200000, max(coord_df$X) + 200000, min(coord_df$Y) - 200000, max(coord_df$Y) + 200000)
finland_andmore_reproj_crop <- crop(finland_andmore_reproj, area)

# prepare hexagon for aggregation
size <- 100000

set.seed(1)
hex_points <- spsample(coord_sp, type = "hexagonal", cellsize = size, offset = c(0, 0))
hex_points$id <- 1:length(hex_points)
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
hex_grid$id <- 1:length(hex_grid)

plot(finland_reproj, col = "grey50", bg = "light blue", axes = TRUE)
plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
plot(hex_grid, border = "orange", add = T)

plot_hex <- st_intersection(coord_sp %>% st_as_sf(), hex_grid %>% st_as_sf()) # get hexagon id for each location
hex_df <- plot_hex %>%
  as_tibble() %>%
  dplyr::select(X, Y, hexagon = id) %>%
  left_join(
    hex_points %>%
      as_tibble() %>%
      rename(hex.x = x, hex.y = y, hexagon = id), # join with hexagon centroid coordinates
    by = "hexagon"
  ) %>%
  mutate(hexagon = as.factor(hexagon)) %>%
  mutate(hexagon = fct_shuffle(hexagon))

ggplot() +
  geom_point(data = hex_df, aes(x = X, y = Y, col = hexagon)) +
  geom_polygon(
    data = hex_grid, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = NA, linewidth = .1
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_equal()

data_agg_df <- data_df %>%
  left_join(hex_df, by = c("X", "Y")) %>%
  group_by(species, mig, area, year, period, hexagon, hex.x, hex.y) %>%
  summarise(
    bh = median(bh, na.rm = T),
    mat = median(mat, na.rm = T),
    n = n()
  ) %>%
  ungroup() %>%
  filter(n >= 50) %>%
  rename(x = hex.x, y = hex.y)

# select species
length(bird_df %>% pull(species) %>% unique()) # 73 species in total

sp_list <- data_agg_df %>%
  group_by(species, period) %>%
  drop_na() %>%
  summarise(n = n()) %>%
  spread(key = "period", value = "n") %>%
  filter(early >= 100 & late >= 100) %>% # at least 100 records in both early and late periods
  pull(species)
length(sp_list) # 38 species

data_agg_df <- data_agg_df %>% filter(species %in% sp_list)

# settings
sp_vis <- "Phoenicurus phoenicurus" # species to plot
userandom <- F # generate results using random grouping (for SI)

if (userandom) {
  set.seed(1)
  data_agg_df <- data_agg_df %>%
    group_by(species) %>%
    mutate(random = rbernoulli(n())) %>%
    mutate(random = case_when(
      random == 1 ~ "in-sample",
      TRUE ~ "out-of-sample"
    )) %>%
    ungroup()
}

# map records
p_map <- ggplot() +
  geom_polygon(
    data = finland_andmore_reproj_crop, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "white", linewidth = .1
  ) +
  geom_polygon(
    data = finland_reproj, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "lightblue", linewidth = .1
  ) +
  geom_polygon(
    data = hex_grid, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = NA, linewidth = .1
  ) +
  geom_jitter(
    data = bird_df %>% filter(species == sp_vis) %>%
      mutate(rank = rank(bh)),
    aes(x = X, y = Y, col = rank), alpha = 0.5, cex = 0.5
  ) +
  theme_void() +
  scale_color_viridis_c(
    direction = -1,
    breaks = bird_df %>%
      filter(species == sp_vis) %>%
      mutate(rank = rank(bh)) %>%
      arrange(bh) %>%
      mutate(cut = cut(bh,
        breaks = c(min(bh), c(170, 180, 190), max(bh)),
        include.lowest = T
      )) %>%
      group_by(cut) %>%
      summarise(rank = max(rank) + 0.5) %>%
      slice(-n()) %>%
      pull(rank),
    label = c(170, 180, 190)
  ) +
  coord_equal() +
  labs(col = "Nestling ringing time \n (day of year)") +
  theme(legend.position = "bottom")
p_map

# summarize niche and change over time
data_agg_df %>%
  group_by(species) %>%
  summarise(
    mat = median(mat, na.rm = T),
    bh = median(bh, na.rm = T)
  ) %>%
  summary()

data_agg_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(mat ~ year, .))) %>%
  filter(term %in% c("year")) %>%
  filter(p.value < 0.05) %>%
  summary()

data_agg_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(bh ~ year, .))) %>%
  filter(term %in% c("year")) %>%
  filter(
    p.value < 0.05,
    estimate < 0
  ) %>%
  summary()

# correlation between temperature and phenology
data_agg_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(bh ~ mat, .))) %>%
  filter(term %in% c("mat")) %>%
  dplyr::select(-statistic) %>%
  filter(p.value < 0.05) %>%
  summary()

# change in functional relationship
if (!userandom) {
  data_agg_df %>%
    group_by(species, period) %>%
    do(broom::tidy(lm(bh ~ mat, .))) %>%
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
  data_agg_df %>%
    group_by(species, random) %>%
    do(broom::tidy(lm(bh ~ mat, .))) %>%
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
  ggplot(data_agg_df) +
    geom_point(aes(x = mat, y = bh, col = period), alpha = 0.1) +
    geom_smooth(aes(x = mat, y = bh, group = period, col = period), method = "lm") +
    theme_classic() +
    facet_wrap(. ~ paste0(mig, ": ", species), scales = "free")
} else {
  ggplot(data_agg_df) +
    geom_point(aes(x = mat, y = bh, col = random), alpha = 0.1) +
    geom_smooth(aes(x = mat, y = bh, group = random, col = random), method = "lm") +
    theme_classic() +
    facet_wrap(. ~ paste0(mig, ": ", species), scales = "free")
}

if (!userandom) {
  p_func <- ggplot(data_agg_df %>% filter(species == sp_vis)) +
    geom_point(aes(x = mat, y = bh, col = period), alpha = 0.1) +
    geom_smooth(aes(x = mat, y = bh, group = period, col = period), method = "lm") +
    theme_classic() +
    labs(color = "Time period") +
    xlab("Mean annual temperature (°C)") +
    ylab("Nestling ringing time (day of year)")
} else {
  p_func <- ggplot(data_agg_df %>% filter(species == sp_vis)) +
    geom_point(aes(x = mat, y = bh, col = random), alpha = 0.1) +
    geom_smooth(aes(x = mat, y = bh, group = random, col = random), method = "lm") +
    theme_classic() +
    labs(color = "Group") +
    xlab("Mean annual temperature (°C)") +
    ylab("Nestling ringing time (day of year)")
}
p_func

# calculate mismatch
if (!userandom) {
  outfile <- "./empirical/bird/data/mismatch.rds"
} else {
  outfile <- "./empirical/bird/data/mismatch_SI.rds"
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
      data_sp <- data_agg_df %>%
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
      z <- resid(lm(bh ~ mat, data_train)) # residuals of non-spatial model
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
        "sigma.sq" = 0.1,
        "tau.sq" = 0.1
      )

      # set knots for predictive process
      knots <- c(3, 3)

      # fit spatial model
      splm <- spLM(bh ~ mat,
        data = data_train,
        coords = data_train[, c("x", "y")] %>% as.matrix(),
        cov.model = "exponential",
        priors = priors,
        tuning = tuning,
        starting = starting,
        n.samples = 1000,
        n.report = 100,
        knots = knots
      )

      # prediction
      splm.pred <- spPredict(splm,
        pred.covars = cbind(1, data_sp[, "mat"]) %>% as.matrix(), pred.coords = data_sp[, c("x", "y")] %>% as.matrix(),
        start = 501
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
    mutate(resid = bh - predict)
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
      R2 = cor(predict, bh)^2,
      rmse = Metrics::rmse(predict, bh)
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
      R2 = cor(predict, bh)^2,
      rmse = Metrics::rmse(predict, bh)
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
    geom_point(aes(x = predict, y = bh), alpha = 0.2) +
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
        x = 182, y = 220,
        label = paste0("Ymis=", round(median, 2), "(", round(lower, 2), ", ", round(upper, 2), ")")
      ),
      col = "blue"
    ) +
    theme_classic() +
    facet_wrap(. ~ period, nrow = 1) +
    guides(col = "none") +
    xlab("Predicted nestling ringing time (day of year)") +
    ylab("Observed nestling ringing time (day of year)")
} else {
  p_pred <- ggplot(data_mis %>% filter(species == sp_vis)) +
    geom_point(aes(x = predict, y = bh), alpha = 0.2) +
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
        x = 182, y = 220,
        label = paste0("Ymis=", round(median, 2), "(", round(lower, 2), ", ", round(upper, 2), ")")
      ),
      col = "blue"
    ) +
    theme_classic() +
    facet_wrap(. ~ random, nrow = 1) +
    guides(col = "none") +
    xlab("Predicted nestling ringing time (day of year)") +
    ylab("Observed nestling ringing time (day of year)")
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
      estimate = t.test(x = resid)$estimate
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
      estimate = t.test(x = resid)$estimate
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
      data = data_mis %>% filter(period == "late") %>% left_join(t_df, by = "species") %>% filter(p < 0.05),
      aes(x = reorder(species, desc(species)), y = resid, fill = species),
      draw_quantiles = 0.5
    ) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    guides(col = "none", fill = "none") +
    coord_flip() +
    ylab("Deviation of observed nestling ringing time \n from predicted nestling ringing time (day)") +
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
    ylab("Deviation of observed nestling ringing time \n from predicted nestling ringing time (day)") +
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
  figname <- "./empirical/bird/figure/figure.pdf"
} else {
  figname <- "./empirical/bird/figure/figure_SI.pdf"
}

cairo_pdf(figname, width = 10, height = 12)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
  annotate_figure(p_func, fig.lab = "b"),
  annotate_figure(p_pred, fig.lab = "c"),
  annotate_figure(p_mis, fig.lab = "d"),
  layout_matrix = rbind(
    c(1, 4),
    c(2, 4),
    c(3, 4)
  ),
  widths = c(1, 1),
  heights = c(3, 2, 2)
)
dev.off()
