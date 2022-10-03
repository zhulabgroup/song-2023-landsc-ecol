# source: https://datadryad.org/stash/dataset/doi:10.5061/dryad.wstqjq2ht
library(raster)
library(gdalUtils)
library(maptools)
bird_df<-read_csv("./empirical/bird/data/73_species.csv")  %>% 
  dplyr::select(
    nestid=NestID,
    X=XEUREF,
    Y=YEUREF,
    species_short=Species,
    year=Year,
    bh=Dayofyear,
    area=BZ) %>% 
  # mutate(date=as.Date(paste0(year, "-01-01"))+bh-1) %>% 
  mutate(area=fct_relevel(area, levels=c("HB", "SB", "MB", "NB"))) %>% 
  # group_by(X, Y, species_short, year, area) %>% 
  # summarise(bh=median(bh, na.rm=T)) %>% 
  # ungroup() %>% 
  mutate(period=case_when(year>=1995~"late",
                          TRUE~"early"))


# join with trait data
sp_df<-read_csv("./empirical/bird/data/Traits_73_species.csv") %>% 
  dplyr::select(
    species_short=Abbreviation,
    species=`Scientific name`,
    mig=Mig
  )

bird_df<-bird_df %>% 
  left_join(sp_df, by="species_short")

# project to lon lat
coord_df<-bird_df %>% 
  distinct(X, Y)# %>% 
  # mutate(coordid=row_number())
coord_sp<-SpatialPointsDataFrame(coords=coord_df[,c("X", "Y")],data=coord_df,
                                 proj4string = CRS("EPSG:3067"))

coord_sp_reproj<-spTransform(coord_sp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coord_df_reproj<-as.data.frame(coord_sp_reproj) %>% 
  rename(lon=X.1, lat=Y.1)

bird_df<-bird_df %>% 
  left_join(coord_df_reproj, by=c("X", "Y"))


###
sp_vis<-"Phoenicurus phoenicurus"

library(maps)
world <- map("world", fill = TRUE)
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world <- map2SpatialPolygons(world, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
finland<-world["Finland"]
finland_andmore<-world[c("Finland","Russia", "Norway", "Sweden")]
finland_reproj<-spTransform(finland,CRS("EPSG:3067") )
finland_andmore_reproj<-spTransform(finland_andmore,CRS("EPSG:3067") )

area <- extent(min(coord_df$X) - 200000, max(coord_df$X) + 200000, min(coord_df$Y) - 200000, max(coord_df$Y) + 200000)
finland_andmore_reproj_crop <- crop(finland_andmore_reproj, area)

### aggregate
size <- 100000
# size <- 1.5 # used for supplementary analysis. When doing so, add "_1.5" to file names.
# size <- 0.5 # used for supplementary analysis. When doing so, add "_0.5" to file names.
hex_points <- spsample(coord_sp, type = "hexagonal", cellsize = size, offset = c(0, 0))
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
plot(finland_reproj, col = "grey50", bg = "light blue", axes = TRUE)
plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
plot(hex_grid, border = "orange", add = T)


plot_hex <- spatialEco::point.in.poly(coord_sp, hex_grid)
hex_df<-plot_hex %>% 
  data.frame() %>% 
  dplyr::select(X, Y,hexagon=poly.ids) %>% 
  left_join(hex_points %>% 
              data.frame() %>% 
              rename(hex.x=x, hex.y=y) %>% 
              mutate(hexagon=row_number()),
            by="hexagon") %>% 
  mutate(hexagon = as.factor(hexagon)) %>% 
  mutate(hexagon=fct_shuffle(hexagon))

ggplot() +
  geom_point(data=hex_df ,aes(x = X, y = Y, col = hexagon)) +
  geom_polygon(data=hex_grid,aes(x = long, y = lat, group=group),
               color = "darkblue", fill = NA, size = .1
  )+
  theme_classic() +
  theme(legend.position = "none")

p_map<-ggplot()+
  geom_polygon(
    data = finland_andmore_reproj_crop, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "white", size = .1
  ) +
  geom_polygon(
    data = finland_reproj, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "lightblue", size = .1
  ) +
  geom_polygon(data=hex_grid,aes(x = long, y = lat, group=group),
               color = "darkblue", fill = NA, size = .1
  )+
  # geom_hex(data=bird_df,
  #          aes(x=X, y=Y,fill=area,alpha=log(..count..)),bins=50)+
  geom_jitter(data=bird_df %>% filter(species==sp_vis) %>% 
                mutate(rank=rank(bh)),
             aes(x=X, y=Y, col=rank), alpha=0.5, cex=0.5)+
  theme_void()+
  scale_color_viridis_c(direction=-1,
                        breaks = bird_df %>% 
                          filter(species==sp_vis) %>% 
                          mutate(rank=rank(bh)) %>% 
                          arrange(bh) %>%
                          mutate(cut = cut(bh,
                                           breaks = c(min(bh), c( 170, 180, 190), max(bh)),
                                           include.lowest = T
                          )) %>%
                          group_by(cut) %>%
                          summarise(rank = max(rank) + 0.5) %>%
                          slice(-n()) %>%
                          pull(rank),
                        label = c( 170, 180, 190))+
  # facet_wrap(.~period)+
  coord_equal()+
  # xlab("Latitude")+
  # ylab("Longitude")+
  labs(col="Nestling ringing time \n (day of year)")+
  theme(legend.position="bottom")

p_map

# ### vipphen SOS
# extent_sp<-extent(coord_df_reproj$lon %>% min()-1,
#                   coord_df_reproj$lon %>% max()+1,
#                   coord_df_reproj$lat %>% min()-1,
#                   coord_df_reproj$lat %>% max()+1) %>% 
#   as('SpatialPolygons') %>% 
#   `projection<-` (CRS("+proj=longlat +datum=WGS84"))
# 
# files<-list.files("/data/ZHULAB/phenology/VIPPHEN", pattern =".hdf", full.names = T) %>% sort()
# 
# cl <- makeCluster(20, outfile = "")
# registerDoSNOW(cl)
# sos_df_list<-
# foreach (file = files,
#          .packages = c("gdalUtils", "raster", "tidyverse")) %dopar% {
#   year<-file %>% 
#     str_split(pattern = "/") %>% 
#     unlist() %>% 
#     tail(1) %>% 
#     str_split(pattern = "\\.") %>% 
#     unlist() %>% 
#     .[2] %>% 
#     str_replace("A", "") %>% 
#     as.numeric()
#   sds <- get_subdatasets(file)
#   sos<-sds[1] %>% 
#     raster() %>% 
#     flip(direction = "y") %>% 
#     `extent<-` (c(-180,180,-90,90)) %>% 
#     `projection<-` (CRS("+proj=longlat +datum=WGS84")) %>% 
#     crop( extent_sp)
#   sos[sos<=0]<-NA
#   sos[sos>365]<-NA
#   mask<-sds[26] %>% 
#     raster() %>% 
#     flip(direction = "y") %>% 
#     `extent<-` (c(-180,180,-90,90)) %>% 
#     `projection<-` (CRS("+proj=longlat +datum=WGS84"))%>% 
#     crop( extent_sp)
#   mask[mask>3]<-NA
#   mask[!is.na(mask)]<-1
#   sos_mask<-sos*mask
#   
#   coord_df_reproj %>% 
#     dplyr::select(lon, lat) %>% 
#     mutate(sos=raster::extract(sos_mask,coord_sp_reproj)) %>% 
#     mutate(year=year)
# }
# sos_df<-bind_rows(sos_df_list)
# sos_df
# write_rds(sos_df, "./empirical/bird/data/sos.rds")
# sos_df_goodcoord<-sos_df %>% 
#   drop_na() %>% 
#   group_by(lon, lat) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   filter(n>=25)
# 
# ggplot()+
#   geom_polygon(
#     data = finland, aes(x = long, y = lat, group = group),
#     color = "darkblue", fill = "lightblue", size = .1
#   ) +
#   geom_point(data=sos_df_goodcoord,aes(x=lon, y=lat, col=n))+
#   
#   
# sos_df<-sos_df %>% 
#   right_join(sos_df_goodcoord %>% dplyr::select(lon, lat), by=c("lon", "lat"))
# 
# ggplot(sos_df %>% filter(lon==sos_df$lon[1], lat==sos_df$lat[1]))+
#   geom_point(aes(x=year, y=sos))

###
extent_sp<-extent(coord_df_reproj$lon %>% min()-1,
                  coord_df_reproj$lon %>% max()+1,
                  coord_df_reproj$lat %>% min()-1,
                  coord_df_reproj$lat %>% max()+1) %>% 
  as('SpatialPolygons') %>% 
  `projection<-` (CRS("+proj=longlat +datum=WGS84"))

bri<-brick("/data/ZHULAB/phenology/Velocity/TerraClimate/MAT_0.05degree.nc")

cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)

mat_df_list<-
  foreach (year = (bird_df$year %>% min()):2015,
           .packages = c("gdalUtils", "raster", "tidyverse")) %dopar% {
             
             mat_ras<-bri[[year-1958+1]] %>% 
               crop( extent_sp)
             
             coord_df_reproj %>% 
               dplyr::select(lon, lat) %>% 
               mutate(mat=raster::extract(mat_ras,coord_sp_reproj)) %>% 
               mutate(year=year)
           }
mat_df<-bind_rows(mat_df_list)
mat_df
write_rds(mat_df, "./empirical/bird/data/mat.rds")

# sos_df<-read_rds("./empirical/bird/data/sos.rds")
mat_df<-read_rds("./empirical/bird/data/mat.rds")
# year_list<-bird_df$year %>% unique() %>% sort()
data_df<-bird_df %>% 
  # group_by(X, Y,species_short,area, period,species, mig, coordid, lon, lat, year) %>% 
  # summarise(bh=mean(bh)) %>% 
  dplyr::select(nestid, X, Y, year, bh, area, species, mig, lon, lat) %>% 
  spread(key="year", value="bh") %>% 
  gather(key="year", value="bh", -nestid, -X, -Y, -area, -species, -mig, -lon, -lat) %>% 
  mutate(year=as.numeric(year)) %>% 
  mutate(period=case_when(year>=1995~"late",
                          TRUE~"early")) %>% 
  # left_join(sos_df, by=c("lat", "lon", "year")) %>% 
  left_join(mat_df, by=c("lat", "lon", "year")) #%>% 
  # group_by(species,area, year, period, mig) %>% 
  # summarise(sos=median(sos, na.rm=T),
  #           mat=median(mat, na.rm=T),
  #           bh=quantile(bh, 0.05, na.rm=T),
  #           X=median(X),
  #           Y=median(Y),
  #           n=n()) %>% 
  # ungroup()

data_agg_df<-data_df %>% 
  left_join(hex_df, by=c("X", "Y")) %>% 
  group_by(species, mig, year, period, hexagon, hex.x, hex.y) %>% 
  summarise(bh=median(bh, na.rm=T),
            mat=median(mat, na.rm=T),
            n=n()) %>% 
  ungroup() %>% 
  filter(n>=50)

sp_list<-data_agg_df %>% 
  group_by(species, period) %>% 
  drop_na() %>% 
  summarise(samplesize=n()) %>% 
  spread(key="period", value="samplesize") %>% 
  filter(early>=100&late>=100) %>% 
  pull(species)
data_agg_df<-data_agg_df %>% filter(species %in% sp_list)

### niche
data_agg_df %>% pull(year) %>% summary()
data_agg_df %>% 
  group_by(species) %>% 
  do(broom::tidy(lm(mat ~ year, .))) %>%
  filter(term %in% c("year")) %>% 
  filter(p.value<0.05) %>% 
  summary()

data_agg_df %>% 
  group_by(species) %>% 
  do(broom::tidy(lm(bh ~ year, .))) %>%
  filter(term %in% c("year")) %>%
  filter(p.value<0.05) %>% 
  summary()

data_agg_df %>% 
  group_by(species) %>% 
  summarise(mat=median(mat, na.rm=T),
            bh=median(bh, na.rm=T)) %>% 
  summary()

# ggplot(data_df%>% filter(species==sp_vis) )+
#   geom_point(aes(x=year, y=bh), alpha=0.3)+
#   geom_smooth(aes(x=year, y=bh), method="loess")+
#   theme_classic()+
#   facet_wrap(.~species)

# ggplot(data_df %>% filter(species==sp_vis) %>% 
#          filter(sos<=180))+
#   geom_point(aes(x=year, y=sos), alpha=0.3)+
#   geom_smooth(aes(x=year, y=sos), method="loess")+
#   theme_classic()+
#   facet_wrap(.~species)

ggplot(data_agg_df)+
  geom_point(aes(x=year, y=bh, col=period), alpha=0.1)+
  geom_smooth(aes(x=year, y=bh), method="lm")+
  theme_classic()+
  facet_wrap(.~paste0(mig,": ",species), scales = "free")

### relationship
data_agg_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(bh ~ mat, .))) %>%
  filter(term %in% c("mat")) %>%
  dplyr::select(-statistic) %>% 
  filter(p.value<0.05) %>% 
  summary()

data_agg_df %>%
  group_by(species, period) %>%
  do(broom::tidy(lm(bh ~ mat, .))) %>%
  ungroup() %>% 
  dplyr::select(species, period, term, estimate) %>% 
  spread(key="period", value="estimate") %>% 
  mutate(diff=late-early) %>% 
  group_by(term) %>%
  summarise(median=median(diff),
            lower=quantile(diff, 0.025),
            upper=quantile(diff, 0.975))

ggplot(data_agg_df )+
  geom_point(aes(x=mat, y=bh, col=period), alpha=0.1)+
  geom_smooth(aes(x=mat, y=bh,group=period, col=period), method="lm")+
  theme_classic()+
  facet_wrap(.~paste0(mig, ": ",species), scales = "free")
# ggplot(data_df)+
#   geom_point(aes(x=mat, y=sos, col=period), alpha=0.3)+
#   geom_smooth(aes(x=mat, y=sos,group=period, col=period), method="lm")+
#   theme_classic()+
#   facet_wrap(.~paste0(mig, ": ",species), scales = "free")
# ggplot(data_df)+
#   geom_point(aes(x=sos, y=bh, col=period), alpha=0.3)+
#   geom_smooth(aes(x=sos, y=bh,group=period, col=period), method="lm")+
#   theme_classic()+
#   facet_wrap(.~paste0(mig, ": ",species), scales = "free")

p_func<-ggplot(data_agg_df%>% filter(species==sp_vis))+
  geom_point(aes(x=mat, y=bh, col=period), alpha=0.1)+
  geom_smooth(aes(x=mat, y=bh,group=period, col=period), method="lm")+
  theme_classic()+
  xlab("Mean annual temperature (°C)")+
  ylab("Nestling ringing time (day of year)")+
  labs(color='Time period')
p_func

data_df_new_list<-
foreach (sp = sp_list,
         .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
           # spatial regression
           data_sp<-data_agg_df %>% 
             filter(species  == sp) %>% 
             drop_na() %>% 
             dplyr::select(x=hex.x, y=hex.y, mat, bh, period, species, mig, year)
           ### Fit nonspatial and spatial model
           
           # Prepare data frame
           data_train <-data_sp %>% 
             filter(period=="early")
           
             # Preparing to run Gaussian spatial model; get model phi (spatial decay parameter) estimate
             z <- resid(lm(bh ~ mat, data_train)) # residuals of non-spatial model
             data_train$z <- z
             lm.vgm <- variogram(z ~ 1, data = data_train %>% `coordinates<-`(data_train[, c("x", "y")]))
             # plot(lm.vgm)
             lm.fit <- fit.variogram(lm.vgm, model = vgm("Exp"), fit.method = 6, fit.sills = F)
             phi <- lm.fit[2, ]$range # extract phi value from fitted variogram
             
             # Set priors: loose priors on beta and residual error variance (tausq), and on spatial variance parameter (sigma.sq), but very tight on Phi (spatial decay parameter).
             priors <- list("beta.Norm" = list(rep(0, 2), diag(100, 2)), "phi.Unif" = c(-log(0.05) / (phi * 100), -log(0.05) / (phi / 100)), "sigma.sq.IG" = c(2, 2), "tau.sq.IG" = c(2, 0.1)) # shape and scale for IG
             # Set starting and tuning values
             starting <- list("phi" = -log(0.05) / phi, "sigma.sq" = 50, "tau.sq" = 1)
             tuning <- list("phi" = (log(0.05) / (phi * 100) - log(0.05) / (phi / 100)) / 10, "sigma.sq" = 0.1, "tau.sq" = 0.1)
             
             # Knots for Gaussian models
             # knots = kmeans(data_train[, c("x", "y")] %>% as.matrix(), 10,iter.max=100)$centers
             knots = c(3,3)
             
             # Run Gaussian spatial model
             t1<-Sys.time()
             splm <- spLM(bh ~ mat, data = data_train, coords = data_train[, c("x", "y")] %>% as.matrix(), cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 1000, n.report = 100
                          , knots=knots
             )
             t2<-Sys.time()
             t2-t1
             
             # splm <- spRecover(splm, get.beta = TRUE, get.w = TRUE, start = 1, n.report = 10)
             # library(ggmcmc)
             # splm_mcmc1 <- ggs(splm$p.beta.recover.samples)
             # ggs_traceplot(splm_mcmc1) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0))
             
             ### Prediction
             t1<-Sys.time()
             splm.pred <- spPredict(splm, pred.covars=cbind(1, data_sp[, "mat"]) %>% as.matrix(), pred.coords=data_sp[, c("x", "y")] %>% as.matrix(),
                                    start=501)
             t2<-Sys.time()
             t2-t1
             
             quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
             y.hat <- apply(splm.pred$p.y.predictive.samples, 1, quant) %>% 
               t() %>% 
               as.data.frame() %>% 
               `colnames<-`(c("lower", "predict", "upper"))
             
             data_sp %>% 
               bind_cols(y.hat)
           
}
data_mis<-bind_rows(data_df_new_list) %>% 
  mutate(resid=bh-predict) 
write_rds(data_mis, "./empirical/bird/data/mismatch_mat.rds")
stopCluster(cl)

data_mis<-read_rds("./empirical/bird/data/mismatch_mat.rds")

data_mis %>% 
  group_by(species, period) %>% 
  summarise (R2=cor(predict, bh)^2,
             rmse=Metrics::rmse(predict, bh)) %>% 
  gather(key="stat", value="value", -species, -period) %>% 
  spread(key="period", value="value") %>% 
  mutate(diff=late-early) %>% 
  gather(key="report", value="value", -species, -stat) %>% 
  group_by(stat, report ) %>% 
  summarise(median=quantile(value, 0.5),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))
# plot mismatch

p_pred<-ggplot(data_mis  %>% filter(species==sp_vis))+
  geom_point(aes(x=predict, y=bh), alpha=0.2)+
  # geom_errorbarh(aes(y=bh, xmin=lower, xmax=upper), alpha=0.5)+
  # geom_smooth(aes(x=predict, y=bh), method="lm")+
  geom_abline(intercept = 0, slope=1, col="red")+
  geom_text(data=data_mis  %>% 
              filter(species==sp_vis) %>% 
              group_by(period) %>% 
              summarise(median=quantile(resid, 0.5),
                        lower=quantile(resid, 0.025),
                        upper=quantile(resid, 0.975)),
            aes(x=182, y=220, label=paste0("Ymis=",round(median,2),"(", round(lower,2),", ", round(upper,2),")")), col="blue")+
  # ggpubr::stat_cor(
  #   aes(
  #     x=predict, y=bh,
  #     label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
  #   ),
  #   p.accuracy = 0.05,
  #   label.x.npc = "left",
  #   label.y.npc = "top",
  #   show.legend = F,
  #   col = "blue"
  # ) +
  theme_classic()+
  facet_wrap(.~period, nrow=1)+
  guides(col="none")+
  xlab("Predicted nestling ringing time (day of year)")+
  ylab("Observed nestling ringing time (day of year)")
p_pred

# test if mismatch is different from 0
data_mis %>% 
  filter(period=="late") %>% 
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975))

t.test(data_mis %>% 
         filter(period=="late") %>% 
         pull(resid),
       alternative = "less")

t_df<-data_mis %>% 
  filter(period=="late") %>% 
  group_by(species) %>%
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975),
            p=t.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate) 

t_df %>% 
  filter(p<0.05) %>% 
  group_by(median<0) %>% 
  summarise(n=n(),
            lower=min(median),
            upper=max(median))

p_mis<-ggplot()+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=as.factor("all species"),y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late"),
              aes(x=reorder(species, desc(species)),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late") %>% left_join(t_df, by="species") %>% filter(p<0.05),
              aes(x=reorder(species, desc(species)),y=resid, fill=species), draw_quantiles = 0.5)+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=mig,y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Deviation of observed nestling ringing time \n from predicted nestling ringing time (day)")+
  xlab ("Species")+
  theme(axis.text.y =element_text(face="italic")) 
p_mis


ggplot(data_mis %>% filter(period=="late") %>%
         left_join(t_df, by="species") %>%
         filter(estimate<0 & p<0.05) )+
  # geom_point(aes(x=y, y=resid, col=species), alpha=0.1)+
  geom_smooth(aes(x=y, y=resid, col=species), method="lm")+
  ggpubr::stat_cor(
    aes(
      x=y, y=resid,
      label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
    ),
    p.accuracy = 0.05,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F,
    col = "blue"
  ) +
  geom_hline(yintercept = 0, lty=2)+
  theme_classic()+
  # facet_wrap(.~species, scales="free")+
  guides(col="none")

library(nlme)
# lme.fit <- lme(resid ~ y, random = ~ 1 + y | species, data = data_mis %>% filter(period=="late"))
# summary(lme.fit)
lme.fit <- lme(abs(resid) ~ y, random = ~ 1+y | species,
               # corr = corSpatial(form = ~x + y, type ="exponential", nugget = T),
               data = data_mis %>% filter(period=="late") %>%
                 # left_join(t_df, by="species") %>%
                 # filter(estimate<0 & p<0.05) %>% 
                 mutate(x=jitter(x, amount=0.000001), y=jitter(y, amount=0.000001)) ,
               control = lmeControl(opt = "optim", optimMethod = "SANN"))
summary(lme.fit)

# data_mis_summ<-data_mis %>% 
#   group_by(species, period, mig, area) %>% 
#   summarise(resid=mean(resid)) %>% 
#   ungroup()
# 
# ggplot(data_mis %>% filter(period=="late"))+
#   # geom_violin(aes(x=interaction(mig,area), y=resid, fill=mig))+
#   # geom_violin(aes(x=area, y=resid, col=area))+
#   geom_violin(aes(x=mig, y=resid, col=mig))+
#   geom_jitter(aes(x=mig, y=resid, col=mig), alpha=0.3)+
#   geom_hline(yintercept = 0)+
#   guides(col="none")+
#   theme_classic()

t_df_area<-data_mis %>% 
  filter(period=="late") %>% 
  group_by( area) %>%
  summarise(p=t.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate) 

data_mis %>%
  filter(period=="late") %>% 
  group_by(area) %>% 
  summarise(resid=median(resid))

p_area<-ggplot()+
  geom_violin(data=data_mis %>% filter(period=="late"),
              aes(x=area, y=resid, col=area), draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late") %>% left_join(t_df_area, by="area") %>% filter(p<0.05),
              aes(x=area, y=resid, fill=area), draw_quantiles = 0.5)+
  # geom_jitter(aes(x=area, y=resid, col=area), alpha=0.3)+
  geom_hline(yintercept = 0)+
  guides(fill="none")+
  theme_classic()+
  coord_flip()+
  xlab("Bioclimatic zone")+
  ylab("Difference between observed and predicted BBS (day)")
p_area

cairo_pdf("./empirical/bird/figure.pdf", width = 10, height = 12)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
             annotate_figure(p_func, fig.lab = "b"),
             annotate_figure(p_pred, fig.lab = "c"),
             annotate_figure(p_mis, fig.lab = "d"),
             # annotate_figure(p_area, fig.lab = "(e)"),
             layout_matrix = rbind(
               c(1, 4),
               c(2, 4),
               c(3, 4)
             ),
             widths=c(1,1),
             heights=c(3,2,2)
)
dev.off()
