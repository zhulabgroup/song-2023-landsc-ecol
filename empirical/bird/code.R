# source: https://datadryad.org/stash/dataset/doi:10.5061/dryad.wstqjq2ht
library(raster)
library(gdalUtils)
bird_df<-read_csv("./empirical/bird/data/73_species.csv")  %>% 
  dplyr::select(
    nestid=NestID,
    X=XEUREF,
    Y=YEUREF,
    species_short=Species,
    year=Year,
    doy=Dayofyear,
    area=BZ) %>% 
  mutate(period=case_when(year>=2000~"late",
                          TRUE~"early")) %>% 
  mutate(date=as.Date(paste0(year, "-01-01"))+doy-1) %>% 
  mutate(area=fct_relevel(area, levels=c("HB", "SB", "MB", "NB")))


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
  distinct(X, Y) %>% 
  mutate(coordid=row_number())
coord_sp<-SpatialPointsDataFrame(coords=coord_df[,c("X", "Y")],data=coord_df,
                                 proj4string = CRS("EPSG:3067"))

coord_sp_reproj<-spTransform(coord_sp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coord_df_reproj<-as.data.frame(coord_sp_reproj) %>% 
  rename(lon=X.1, lat=Y.1)

bird_df<-bird_df %>% 
  left_join(coord_df_reproj, by=c("X", "Y"))

sp_vis<-"Phoenicurus phoenicurus"

library(maps)
world <- map("world", fill = TRUE)
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world <- map2SpatialPolygons(world, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
finland<-world["Finland"]
finland_reproj<-spTransform(finland,CRS("EPSG:3067") )

p_map<-ggplot()+
  geom_polygon(
    data = finland_reproj, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "lightblue", size = .1
  ) +
  geom_hex(data=bird_df,
              aes(x=X, y=Y,fill=area,alpha=log(..count..)),bins=50)+
  theme_void()+
  # facet_wrap(.~period)+
  coord_equal()+
  # xlab("Latitude")+
  # ylab("Longitude")+
  labs(alpha="Log of sample size")
  
p_map

### vipphen SOS
extent_sp<-extent(coord_df_reproj$lon %>% min()-1,
                  coord_df_reproj$lon %>% max()+1,
                  coord_df_reproj$lat %>% min()-1,
                  coord_df_reproj$lat %>% max()+1) %>% 
  as('SpatialPolygons') %>% 
  `projection<-` (CRS("+proj=longlat +datum=WGS84"))

files<-list.files("/data/ZHULAB/phenology/VIPPHEN", pattern =".hdf", full.names = T) %>% sort()

cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)
sos_df_list<-
foreach (file = files,
         .packages = c("gdalUtils", "raster", "tidyverse")) %dopar% {
  year<-file %>% 
    str_split(pattern = "/") %>% 
    unlist() %>% 
    tail(1) %>% 
    str_split(pattern = "\\.") %>% 
    unlist() %>% 
    .[2] %>% 
    str_replace("A", "") %>% 
    as.numeric()
  sds <- get_subdatasets(file)
  sos<-sds[1] %>% 
    raster() %>% 
    flip(direction = "y") %>% 
    `extent<-` (c(-180,180,-90,90)) %>% 
    `projection<-` (CRS("+proj=longlat +datum=WGS84")) %>% 
    crop( extent_sp)
  sos[sos<=0]<-NA
  sos[sos>365]<-NA
  mask<-sds[26] %>% 
    raster() %>% 
    flip(direction = "y") %>% 
    `extent<-` (c(-180,180,-90,90)) %>% 
    `projection<-` (CRS("+proj=longlat +datum=WGS84"))%>% 
    crop( extent_sp)
  mask[mask>3]<-NA
  mask[!is.na(mask)]<-1
  sos_mask<-sos*mask
  
  coord_df_reproj %>% 
    dplyr::select(lon, lat) %>% 
    mutate(sos=raster::extract(sos_mask,coord_sp_reproj)) %>% 
    mutate(year=year)
}
sos_df<-bind_rows(sos_df_list)
sos_df
write_rds(sos_df, "./empirical/bird/data/sos.rds")

sos_df<-read_rds("./empirical/bird/data/sos.rds")
data_df<-bird_df %>% 
  left_join(sos_df, by=c("lat", "lon", "year")) %>% 
  group_by(species,area, year, period, mig) %>% 
  summarise(sos=median(sos, na.rm=T),
            bh=quantile(doy, 0.05, na.rm=T),
            X=median(X),
            Y=median(Y),
            n=n()) %>% 
  ungroup()

sp_list<-data_df %>% 
  group_by(species, period) %>% 
  drop_na() %>% 
  summarise(siteyear=n()) %>% 
  spread(key="period", value="siteyear") %>% 
  filter(early>=30&late>=30) %>% 
  pull(species)
data_df<-data_df %>% filter(species %in% sp_list)
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

ggplot(data_df)+
  geom_point(aes(x=year, y=bh, col=period), alpha=0.3)+
  geom_smooth(aes(x=year, y=bh), method="lm")+
  theme_classic()+
  facet_wrap(.~mig*species, scales = "free")


ggplot(data_df)+
  geom_point(aes(x=sos, y=bh, col=period), alpha=0.3)+
  geom_smooth(aes(x=sos, y=bh,group=period, col=period), method="lm")+
  theme_classic()+
  facet_wrap(.~paste0(mig, ": ",species), scales = "free")

p_func<-ggplot(data_df %>% filter(species==sp_vis))+
  geom_point(aes(x=sos, y=bh, col=period), alpha=0.3)+
  geom_smooth(aes(x=sos, y=bh,group=period, col=period), method="lm")+
  theme_classic()+
  xlab("PGS (day of year)")+
  ylab("BBS (day of year)")
p_func

data_df_new_list<-
foreach (sp = sp_list,
         .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
           # spatial regression
           data_sp<-data_df %>% 
             filter(species==sp) %>% 
             drop_na() %>% 
             dplyr::select(x=X, y=Y, sos, bh, period, area, species, mig, year)
           ### Fit nonspatial and spatial model
           
           # Prepare data frame
           data_train <-data_sp %>% 
             filter(period=="early")
           
             # Preparing to run Gaussian spatial model; get model phi (spatial decay parameter) estimate
             z <- resid(lm(bh ~ sos, data_train)) # residuals of non-spatial model
             data_train$z <- z
             lm.vgm <- variogram(z ~ 1, data = data_train %>% `coordinates<-`(data_train[, c("x", "y")]))
             # plot(lm.vgm)
             lm.fit <- fit.variogram(lm.vgm, model = vgm("Exp"), fit.method = 6, fit.sills = F)
             phi <- lm.fit[2, ]$range # extract phi value from fitted variogram
             
             # Set priors: loose priors on beta and residual error variance (tausq), and on spatial variance parameter (sigma.sq), but very tight on Phi (spatial decay parameter).
             priors <- list("beta.Norm" = list(rep(0, 2), diag(100, 2)), "phi.Unif" = c(-log(0.05) / (phi * 100), -log(0.05) / (phi / 100)), "sigma.sq.IG" = c(2, 2), "tau.sq.IG" = c(2, 0.1)) # shape and scale for IG
             # Set starting and tuning values
             starting <- list("phi" = -log(0.05) / phi, "sigma.sq" = 50, "tau.sq" = 1)
             tuning <- list("phi" = (log(0.05) / (phi * 100) - log(0.05) / (phi / 100)) / 10, "sigma.sq" = 0.01, "tau.sq" = 0.01)
             
             # Knots for Gaussian models
             # knots = kmeans(coords, 20,iter.max=100)$centers
             
             # Run Gaussian spatial model
             splm <- spLM(bh ~ sos, data = data_train, coords = data_train[, c("x", "y")] %>% as.matrix(), cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 10000, n.report = 1000
                          # , knots=knots
             )
             
             
             # ### Parameter inference
             #   # Recover parameter samples, ignoring the first half of the run as burn-in.
             #   splm <- spRecover(splm, get.beta = TRUE, get.w = TRUE, start = 5001, n.report = 1000)
             #   beta.hat <- splm$p.beta.recover.samples
             #   theta.hat <- splm$p.theta.recover.samples
             #   
             #   beta.theta.hat <- cbind(beta.hat, theta.hat)
             #   # Coefficient summaries
             #   median <- apply(beta.theta.hat, 2, median)
             #   lower <- apply(beta.theta.hat, 2, quantile, probs = 0.025)
             #   upper <- apply(beta.theta.hat, 2, quantile, probs = 0.975)
             #   
             #   coef.summary <- data.frame(median, lower, upper)
             
             ### Prediction
             splm.pred <- spPredict(splm, pred.covars=cbind(1, data_sp[, "sos"]) %>% as.matrix(), pred.coords=data_sp[, c("x", "y")] %>% as.matrix(),
                                    start=5001, thin = 10)
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
write_rds(data_mis, "./empirical/bird/data/mismatch.rds")
stopCluster(cl)

data_mis<-read_rds("./empirical/bird/data/mismatch.rds")
# ggplot(data_mis %>% filter(period=="early"))+
#   geom_point(aes(x=predict, y=bh, col=area))+
#   geom_errorbarh(aes(y=bh, xmin=lower, xmax=upper, col=area))+
#   geom_abline(intercept = 0, slope=1, col="red")+
#   theme_classic()+
#   facet_wrap(.~mig*species, scales = "free")
# 
# ggplot(data_mis %>% filter(period=="late"))+
#   geom_point(aes(x=predict, y=bh, col=area))+
#   geom_errorbarh(aes(y=bh, xmin=lower, xmax=upper, col=area))+
#   geom_abline(intercept = 0, slope=1, col="red")+
#   theme_classic()+
#   facet_wrap(.~mig*species, scales = "free")

# plot mismatch
p_pred<-ggplot(data_mis  %>% filter(species==sp_vis))+
  geom_point(aes(x=predict, y=bh))+
  geom_errorbarh(aes(y=bh, xmin=lower, xmax=upper), alpha=0.5)+
  geom_abline(intercept = 0, slope=1, col="red")+
  ggpubr::stat_cor(
    aes(
      x=predict, y=bh,
      label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
    ),
    p.accuracy = 0.05,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F,
    col = "blue"
  ) +
  theme_classic()+
  facet_wrap(.~period, nrow=1)+
  guides(col="none")+
  xlab("Predicted BBS (day of year)")+
  ylab("Observed BBS (day of year)")+
  coord_equal()
p_pred

# test if mismatch is different from 0
hist(data_mis %>% 
       filter(period=="late") %>% 
       pull(resid))
# cor_res<-cor.test(data_mis %>%  filter(period=="late") %>% pull(bh),
#          data_mis %>%  filter(period=="late") %>% pull(predict)
#          )
# cor_res
# cor_res$estimate^2

t.test(data_mis %>% 
         filter(period=="late") %>% 
         pull(resid),
       alternative = "greater")
t_df<-data_mis %>% 
  filter(period=="late") %>% 
  group_by(species) %>%
  summarise(p=t.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate) 

p_mis<-ggplot()+
  geom_violin(data=data_mis %>% filter(period=="late"),
              aes(x=as.factor("all species"),y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late"),
              aes(x=paste0(mig, ": ",species),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late") %>% left_join(t_df, by="species") %>% filter(p<0.05),
              aes(x=paste0(mig, ": ",species),y=resid, fill=species), draw_quantiles = 0.5)+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=mig,y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Difference between observed and predicted BBS (day)")+
  xlab ("Species")
p_mis

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

cairo_pdf("./empirical/bird/figure.pdf", width = 10, height = 10)
grid.arrange(annotate_figure(p_map, fig.lab = "(a)"),
             annotate_figure(p_func, fig.lab = "(b)"),
             annotate_figure(p_pred, fig.lab = "(c)"),
             annotate_figure(p_mis, fig.lab = "(d)"),
             annotate_figure(p_area, fig.lab = "(e)"),
             layout_matrix = rbind(
               c(1, 4),
               c(2, 4),
               c(3, 5)
             ),
             widths=c(2,3),
             heights=c(3,3,2)
)
dev.off()
