library(tidyverse)
library(gridExtra)
library(ggpubr)
### Plants and climate
# source: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF309
# pfl: peak flower
# ffr: first fruit
# pfr: peak fruit
# efl: first flower

data_df<-read_csv("./empirical/herbarium/data/hf309-01-crowdcurio.csv") %>% 
  filter(metric=="mean") %>%
  filter(phenophase=="efl"|phenophase=="pfl") %>% 
  filter(!(state %in% c("Florida") )) %>% 
  dplyr::select(species=name,
                # date, 
                year,
                fl=doy,
                lat=county.lat, 
                lon=county.lon,
                temp=spring,
                maxp,
                link
  ) %>%
  filter(maxp>0) %>% 
  arrange(desc(maxp)) %>% 
  distinct(link, .keep_all = T) %>% 
  dplyr::select(-maxp, -link) %>% 
  mutate(period=case_when(year>=1950~"late",
                          TRUE~"early")) %>% 
  mutate(species=str_replace(species, "_", " "))

# select species
sp_list<-data_df %>% 
  group_by(species, period) %>% 
  summarise(siteyear=n()) %>% 
  spread(key="period", value="siteyear") %>% 
  filter(early>=30&late>=30) %>% 
  pull(species)

data_df<-data_df %>% filter(species %in% sp_list)
# visualize data
library(maps)
library(maptools)
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs = IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))

area <- raster::extent(min(data_df$lon) - 0.5, max(data_df$lon) + 0.5, min(data_df$lat) - 0.5, max(data_df$lat) + 0.5)
usa_crop <- raster::crop(usa, area)

sp_vis<-"Aquilegia canadensis"
p_map <- ggplot() +
  geom_polygon(
    data = usa_crop, aes(x = long, y = lat, group = group),
    color = "darkblue", fill = "lightblue", size = .1
  ) +
  geom_point(data=data_df, aes(x=lon, y=lat, col=species), alpha=0.3, pch=1)+
  geom_point(data=data_df %>% filter(species==sp_vis), aes(x=lon, y=lat, col=species), alpha=1, shape=19)+
  theme_minimal() +
  xlab("Longitude") +
  ylab("Latitude") +
  # guides(col = guide_legend(title = "Species")) +
  guides(col = "none") +
  coord_equal()#+
  # facet_wrap(.~period)
p_map

data_df %>% 
  group_by(species) %>% 
  summarise(temp=median(temp),
            fl=median(fl)) %>% 
  summary()
# change in functional relationship 
data_df %>%
  group_by(species) %>%
  do(broom::tidy(lm(fl ~ temp, .))) %>%
  filter(term %in% c("temp")) %>%
  dplyr::select(-statistic) %>% 
  filter(p.value<0.05) %>% 
  summary()

data_df %>%
  group_by(species, period) %>%
  do(broom::tidy(lm(fl ~ temp, .))) %>%
  # filter(term %in% c("temp")) %>%
  ungroup() %>% 
  dplyr::select(species, period, term, estimate) %>% 
  spread(key="period", value="estimate") %>% 
  mutate(diff=late-early) %>% 
  group_by(term) %>%
  summarise(median=median(diff),
            lower=quantile(diff, 0.025),
            upper=quantile(diff, 0.975))


ggplot(data_df)+
  geom_point(aes(x=temp, y=fl, col=period), alpha=0.3)+
  geom_smooth(aes(x=temp, y=fl, group=period, col=period), method="lm")+
  theme_classic()+
  facet_wrap(.~species)+
  theme(legend.position="bottom")+
  labs(color='Time period') +
  xlab("Spring mean temperature (°C)")+
  ylab("Flowering time (day of year)")

p_func<-ggplot(data_df %>% filter(species==sp_vis))+
  geom_point(aes(x=temp, y=fl, col=period), alpha=0.3)+
  geom_smooth(aes(x=temp, y=fl, group=period, col=period), method="lm")+
  theme_classic()+
  # facet_wrap(.~species)+
  theme(legend.position="bottom")+
  labs(color='Time period') +
  xlab("Spring mean temperature (°C)")+
  ylab("Flowering time (day of year)")
p_func

# calculate mismatch
library(parallel)
library(doSNOW)
cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)
data_df_new_list<-
  foreach (sp = sp_list,
           .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
             
             set.seed(42)
             data_sp<-data_df %>% 
               filter(species==sp) %>% 
               drop_na() %>% 
               dplyr::select(x=lon, y=lat, temp, fl, period,  species, year)

             data_train <-data_sp %>% 
               filter(period=="early")
             
             # Preparing to run Gaussian spatial model; get model phi (spatial decay parameter) estimate
             z <- resid(lm(fl ~ temp, data_train)) # residuals of non-spatial model
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
             splm <- spLM(fl ~ temp, data = data_train, coords = data_train[, c("x", "y")] %>% as.matrix(), cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 10000, n.report = 1000
                          # , knots=knots
             )
             
             ### Prediction
             splm.pred <- spPredict(splm, pred.covars=cbind(1, data_sp[, "temp"]) %>% as.matrix(), pred.coords=data_sp[, c("x", "y")] %>% as.matrix(),
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
  mutate(resid=fl-predict) %>% 
  mutate(temp_bin = cut(temp, breaks=5)) %>% 
  mutate(lat_bin = cut(y, breaks=5))
write_rds(data_mis, "./empirical/herbarium/data/mismatch.rds")
stopCluster(cl)

data_mis<-read_rds("./empirical/herbarium/data/mismatch.rds") 
data_mis %>% 
  group_by(species, period) %>% 
  summarise (R2=cor(predict, fl)^2,
             rmse=Metrics::rmse(predict, fl)) %>% 
  gather(key="stat", value="value", -species, -period) %>% 
  spread(key="period", value="value") %>% 
  mutate(diff=late-early) %>% 
  gather(key="report", value="value", -species, -stat) %>% 
  group_by(stat, report ) %>% 
  summarise(median=quantile(value, 0.5),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))
  
data_mis  %>% 
  filter(species==sp_vis) %>% 
  group_by(period) %>% 
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975))

# plot mismatch
p_pred<-ggplot(data_mis  %>% filter(species==sp_vis))+
  geom_point(aes(x=predict, y=fl), alpha=0.5)+
  # geom_smooth(aes(x=predict, y=fl), lty=2, method="lm")+
  # geom_errorbarh(aes(y=fl, xmin=lower, xmax=upper), alpha=0.5)+
  geom_abline(intercept = 0, slope=1, col="red")+
  geom_text(data=data_mis  %>% 
              filter(species==sp_vis) %>% 
              group_by(period) %>% 
              summarise(median=quantile(resid, 0.5),
                        lower=quantile(resid, 0.025),
                        upper=quantile(resid, 0.975)),
            # aes(x=120, y=250, label=paste0("Ymis=",round(median,2),"(", round(lower,2),", ", round(upper,2),")")), col="blue")+
  aes(x=150, y=250, label=paste0("Ymis=",round(median,2),"(", round(lower,2),", ", round(upper,2),")")), col="blue")+
  # ggpubr::stat_cor(
  #   aes(
  #     x=predict, y=fl,
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
  xlab("Predicted flowering time (day of year)")+
  ylab("Observed flowering time (day of year)")
p_pred

# ggplot(data_mis %>% filter(period=="late"))+
#   geom_point(aes(x=x, y=y, col=resid))+
#   theme_classic()+
#   scale_color_viridis_c()+
#   facet_wrap(.~species)

# ggplot(data_mis %>% filter(period=="late"))+
#   geom_violin(aes(x=species, y=resid, col=species))+
#   geom_jitter(aes(x=species, y=resid, col=species), alpha=0.3)+
#   geom_hline(yintercept = 0)+
#   guides(col="none")+
#   theme_classic()

# test if mismatch is different from 0
data_mis %>% 
  filter(period=="late") %>% 
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975))

t.test(data_mis %>% 
         filter(period=="late") %>% 
         pull(resid),
       alternative = "two.sided")

t_df<-data_mis %>% 
  filter(period=="late") %>% 
  group_by(species) %>%
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975),
            p=t.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate,
            med_fl=median(fl)) 
t_df %>% 
  filter(p<0.05) %>% 
  group_by(median<0) %>% 
  summarise(n=n(),
            lower=min(median),
            upper=max(median))

t_df %>% 
  ggplot()+
  geom_point(aes(x=med_fl, y=estimate))+
  geom_smooth(aes(x=med_fl, y=estimate), method="lm")+
  ggpubr::stat_cor(
    aes(
      x=med_fl, y=estimate,
      label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
    ),
    p.accuracy = 0.05,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F,
    col = "blue"
  ) +
  theme_classic()

# sp_list_order<-t_df %>% 
#   arrange(med_fl) %>% 
#   pull(species)
# data_mis<-data_mis %>% 
#   mutate(species=fct_relevel(species, levels=sp_list_order))

p_mis<-ggplot()+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=as.factor("all species"),y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late"),
              aes(x=reorder(species, desc(species)),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(period=="late") %>% left_join(t_df, by="species") %>% filter(p<0.05),
              aes(x=reorder(species, desc(species)),y=resid, fill=species), draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Deviation of observed flowering time \n from predicted flowering time (day)")+
  xlab ("Species")+
  theme(axis.text.y =element_text(face="italic")) 
p_mis

data_mis %>% filter(period=="late") %>% 
  group_by(species) %>% 
  summarise(median=median(resid)) %>% 
  pull(median) %>% 
  t.test()

library(nlme)
lme.fit <- lme((resid) ~ 1, random = ~ 1 | species,
               data = data_mis %>% filter(period=="late") %>%
                 mutate(x=jitter(x, amount=0.000001), y=jitter(y, amount=0.000001)) )
summary(lme.fit)

cairo_pdf("./empirical/herbarium/figure.pdf", width = 10, height = 10)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
             annotate_figure(p_func, fig.lab = "b"),
             annotate_figure(p_pred, fig.lab = "c"),
             annotate_figure(p_mis, fig.lab = "d"),
             layout_matrix = rbind(
               c(1, 3),
               c(2, 4)
             ),
             widths=2:3,
             heights=2:3
)
dev.off()

# pattern in mismatch
reg_df <- data_mis %>%
  filter(period=="late") %>%
  group_by(species) %>%
  do(broom::tidy(lm(resid ~ temp, .))) %>%
  filter(term %in% c("temp")) %>%
  dplyr::select(-statistic)

ggplot()+
  geom_point(data=data_mis %>% filter(period=="late") %>% left_join(t_df, by="species") %>% filter(p<0.05),
             aes(x=temp, y=resid, col=species))+
  geom_hline(yintercept = 0, lty=2)+
  geom_smooth(data=data_mis %>%
                filter(period=="late")  %>% left_join(t_df, by="species") %>% filter(p<0.05)%>%
                left_join(reg_df, by="species") %>%
                filter(p.value<0.05),
              aes(x=temp, y=resid, col=species,group=species), method="lm")+
  # ggpubr::stat_cor(
  #   aes(
  #     x = y, y = resid, group = species,
  #     label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
  #   ),
  #   p.accuracy = 0.05,
  #   label.x.npc = "left",
  #   label.y.npc = "top",
  #   show.legend = F,
  #   col = "blue"
  # ) +
  theme_classic()+
  guides(col="none")+
  facet_wrap(.~species, scales = "free")

library(nlme)
# lme.fit <- lme(resid ~ y, random = ~ 1 + y | species, data = data_mis %>% filter(period=="late"))
# summary(lme.fit)
lme.fit <- lme(abs(resid) ~ temp, random = ~ temp | species,
               corr = corSpatial(form = ~x + y, type ="exponential", nugget = T),
               data = data_mis %>% filter(period=="late") %>% mutate(x=jitter(x, amount=0.000001), y=jitter(y, amount=0.000001)) ,
                control = lmeControl(opt = "optim", optimMethod = "SANN"))
summary(lme.fit)
# 
# # why such pattern
# ggplot(data_mis %>% filter(period=="late"))+
#   geom_point(aes(x=year, y=temp, col=lat_bin,group=lat_bin))+
#   geom_smooth(aes(x=year, y=temp, col=lat_bin, group=lat_bin), method="lm")+
#   theme_classic()
# 
# data_df %>% 
#   filter(period=="late") %>% 
#   mutate(latzone=(lat%/%0.5)*0.5) %>% 
#   group_by(latzone, species) %>%
#   do(broom::tidy(lm(temp ~ year, .))) %>%
#   filter(term %in% c("year")) %>%
#   dplyr::select(-statistic) %>% 
#   ggplot(aes(x=latzone, y=estimate, group=species, col=species))+
#   geom_point()+
#   geom_smooth( method="lm")+
#   theme_classic()+
#   ggpubr::stat_cor(
#     aes(
#       label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
#     ),
#     p.accuracy = 0.05,
#     label.x.npc = "left",
#     label.y.npc = "top",
#     show.legend = F,
#     col = "blue"
#   )+
#   guides(col="none")+
#   facet_wrap(.~species, scales = "free")+
#   geom_hline(yintercept = 0, lty=2)
# 
# data_df %>% 
#   filter(period=="late") %>% 
#   mutate(latzone=(lat%/%0.5)*0.5) %>% 
#   group_by(latzone, species) %>%
#   do(broom::tidy(lm(pfl ~ year, .))) %>%
#   filter(term %in% c("year")) %>%
#   dplyr::select(-statistic) %>% 
#   ggplot(aes(x=latzone, y=estimate, group=species, col=species))+
#   geom_point()+
#   geom_smooth( method="lm")+
#   theme_classic()+
#   ggpubr::stat_cor(
#     aes(
#       label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
#     ),
#     p.accuracy = 0.05,
#     label.x.npc = "left",
#     label.y.npc = "top",
#     show.legend = F,
#     col = "blue"
#   )+
#   guides(col="none")+
#   facet_wrap(.~species, scales = "free")+
#   geom_hline(yintercept = 0, lty=2)