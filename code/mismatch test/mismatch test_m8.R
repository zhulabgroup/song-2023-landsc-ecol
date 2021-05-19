path<-"./archive/m8/"
source("code/utils.R")
source("code/settings.R")

set.seed(1)
# prepare env time series

date_list<-seq(as.Date("2021-01-01"),as.Date("2041-12-31"), by=1)
midyear<-mean(c(as.numeric(format(min(date_list), "%Y")),as.numeric(format(max(date_list), "%Y"))))
coord_df <- data.frame(lon=0, lat=1:5)

distMat<-matrix(c(0,1,2,3,4,
                  1,0,1,2,3,
                  2,1,0,1,2,
                  3,2,1,0,1,
                  4,3,2,1,0), 
                nrow=5)

ts_all<-
foreach (s = 1:nrow(coord_df)#,
         # .packages = c("lubridate")
         ) %dopar% {
           library(lubridate, lib.loc = "/usr/lib64/R/library")
           library(tidyverse, lib.loc = "/usr/lib64/R/library")
  env<-rep(NA, length(date_list))
  for (i in 1:length(date_list)) {
    date<-date_list[i]
    env[i]<-waves(t=date,
                  t_start=date_list[1],
                  intercept=0.3+0.2*coord_df$lat[s],
                  slope = 0.0001,
                  amplitude1 =0.8,
                  phase1 =0,
                  period1=1,
                  amplitude2 = 0.5,
                  phase2 = 0,
                  period2 = 5,
                  sd = 0.05
    )
    # print(i)
  }
  ts_site<-data.frame(date=date_list, env=env, env_sd=0.05)%>% 
    mutate(year=as.numeric(format(date, "%Y"))) %>% 
    mutate(site=s)
  ts_site
}
ts_all<-bind_rows(ts_all)
ggplot(ts_all)+
  geom_line(aes(x=date, y=env, group=site, col=site))+
  theme_classic()

# get phenology model parameter in each year (match and mismatch)

param_all<-
foreach (s = 1:nrow(coord_df)) %dopar% {
  param_site<-ts_all %>% 
    filter(site==s) %>% 
    dplyr::select(env,year) %>% 
    group_by(year) %>% 
    slice_head(n=180) %>% 
    summarize(env_summ=mean(env)) %>% 
    ungroup() %>% 
    mutate(param=env_to_param (env=env_summ,
                               lower=1.2,
                               upper=2.2,
                               steepness=8,
                               midpoint=1.5)) %>% 
    mutate(param_mis = case_when((year>midyear)~param-0.5)) %>%
    mutate(site=s)
  # mutate(param_mis = if_else(param>=25,param,25))
  param_site
}
param_all<-bind_rows(param_all)
param_all

ggplot(param_all)+
  geom_line(aes(x=year, y=env_summ, group=site, col=site))+
  theme_classic()

p1<-
  ggplot(param_all)+
  geom_line(aes(x=env_summ, y=param), col="blue")+
  # geom_line(aes(x=env_summ, y=param_mis), col="red")+
  theme_classic()+
  ylab("m")

colors <- c("m without mismatch" = "blue",
            "m with mismatch" = "red")
p2<-
  ggplot(data=param_all %>% filter(site==3))+
  geom_line(aes(x=year, y=param, group=site, col="m without mismatch"), lwd=2, alpha=0.5)+
  geom_line(aes(x=year, y=param_mis, group=site, col="m with mismatch"), lwd=2, alpha=0.5)+
  # geom_line(data=param_all ,aes(x=year, y=param, group=site, col=site))+
  theme_classic()+
  scale_color_manual(values = colors)+
  labs(x = "year",
       y = "m",
       color = "") +
  theme(legend.position="top") 

# get phenology time series (match and mismatch)
pheno_all<-
foreach (s = 1:nrow(coord_df)) %dopar% {
  param_site<-param_all %>% filter(site==s)
  pheno_list<-vector(mode="list", length=nrow(param_site))
  for (y in 1:nrow(param_site)) {
    date_list_year<-ts_all %>%
      filter(site==s) %>% 
      filter(year==param_site$year[y]) %>% 
      dplyr::select(date) %>% 
      mutate(date=as.character(date)) %>% 
      unlist() %>% 
      as.Date()
    pheno<-pheno_mis<-rep(NA, length(date_list_year))
    for (i in 1:length(date_list_year)) {
      date<-date_list_year[i]
      pheno[i]<-double_logistics(date,
                                 m1=0, # average greenness in winter
                                 m2=1, # difference between summer and winter
                                 m3=100, # spring onset
                                 m4=10, # slope of curve in spring
                                 m5=260, # fall offset
                                 m6=20, # slope of curve in fall
                                 m7=0, # summer greendown
                                 m8=param_site$param[y], #life cycle
                                 sd = 0.05
      )
      pheno_mis[i]<-double_logistics(date,
                                     m1=0, # average greenness in winter
                                     m2=1, # difference between summer and winter
                                     m3=100, # spring onset
                                     m4=10, # slope of curve in spring
                                     m5=260, # fall offset
                                     m6=20, # slope of curve in fall
                                     m7=0, # summer greendown
                                     m8=param_site$param_mis[y], #life cycle
                                     sd = 0.05
      )
      # print(i)
    }
    pheno_year<-data.frame(date=date_list_year, pheno=pheno,pheno_mis=pheno_mis, pheno_sd=0.05)
    pheno_list[[y]]<-pheno_year
  }
  pheno_site<-bind_rows(pheno_list) %>% 
    mutate(site=s)
  pheno_site
}
pheno_all<-bind_rows(pheno_all)

ggplot(pheno_all )+
  geom_line(aes(x=date, y=pheno),col="blue",alpha=0.5)+
  geom_line(aes(x=date, y=pheno_mis), col="red", alpha=0.5)+
  theme_classic()+
  facet_wrap(~site, ncol=1)

p3<-
  ggplot(pheno_all %>%
           filter(site==3) %>% 
           mutate(year=as.numeric(format(date, "%Y")))%>%
           mutate( doy=as.numeric(format(date, "%j"))) %>% 
           left_join(param_all, by=c("year","site")) #%>% 
         # filter(year%in%year_zoom)
  )+
  geom_line(aes(x=doy, y=pheno,group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("pheno")+
  theme(legend.position="bottom")
# facet_wrap(.~site, ncol = 1)

ts_all<-ts_all %>% 
  left_join(pheno_all, by=c("date","site"))

# use first half to train model

ts_train<-ts_all %>% 
  filter(year<=midyear)
ts<-ts_train
source("code/preprocess data.R")
source("code/prepare embeddings.R")
source("code/train GP model.R")

# predict for second half
ts<-ts_all
source("code/preprocess data.R")
source("code/fit.R")

p4<-
  ggplot(fit_df_ori %>%
           filter(site==3) %>% 
           mutate(year=as.numeric(format(date, "%Y")))%>%
           mutate( doy=as.numeric(format(date, "%j"))) %>% 
           left_join(param_all, by=c("year","site")) #%>% 
         # filter(year%in%year_zoom)
  )+
  geom_line(aes(x=doy, y=value,group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("pheno")+
  theme(legend.position="bottom")
  # facet_wrap(.~site, ncol = 1)

combine_df_ori<-ts_all %>% 
  left_join(fit_df_ori %>% dplyr::select(date,site, value, lower, upper), by=c("date","site")) %>% 
  mutate(mismatch_actual=case_when (year>midyear~pheno-pheno_mis),
         mismatch_model=case_when (year<=midyear~ value-pheno,
                                   year>midyear~value-pheno_mis),
         mismatch_model_upper=case_when (year<=midyear~ upper-pheno,
                                   year>midyear~upper-pheno_mis),
         mismatch_model_lower=case_when (year<=midyear~ lower-pheno,
                                         year>midyear~lower-pheno_mis))

colors <- c("phenology without mismatch" = "blue",
            "phenology with mismatch" = "red",
            "predicted phenology" = "black")
p5<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=pheno, col="phenology without mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=pheno_mis, col="phenology with mismatch"), alpha=0.5)+
  geom_line(aes(x=date, y=value, col="predicted phenology"))+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  labs(x = "date",
       y = "phenology",
       color = "") +
  theme(legend.position="top") #+
  # facet_wrap(.~site, ncol = 1)

colors <- c("actual mismatch" = "red",
            "estimated mismatch" = "black")
p6<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=mismatch_actual, col="actual mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=mismatch_model, col="estimated mismatch"), alpha=0.5)+
  geom_ribbon(aes(x=date, ymin=mismatch_model_lower, ymax=mismatch_model_upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  labs(x = "date",
       y = "mismatch",
       color = "") +
  theme(legend.position="top") #+
  # facet_wrap(.~site, ncol=1)

# year_zoom<-2026:2027
# ggplot(combine_df_ori%>%
#          filter(site==3) %>% 
#          filter(as.numeric(format(date, "%Y"))%in%year_zoom))+
#   geom_line(aes(x=date, y=pheno), col="blue",alpha=0.5)+
#   geom_line( aes(x=date, y=pheno_mis), col="red", alpha=0.5)+
#   geom_line(aes(x=date, y=value), col="black")+
#   geom_ribbon( aes(x=date, ymin=lower, ymax=upper), fill="black",alpha=0.25)+
#   theme_classic()

range_df<-combine_df_ori %>% filter(year <=midyear) %>% 
  dplyr::summarize(lower=quantile(pheno, 0.025),
                   upper=quantile(pheno, 0.975)) %>% 
  mutate(range=as.numeric(upper-lower))

combine_df_ori_valid<-combine_df_ori %>% filter(year >midyear)
theo_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$pheno, pred_ori=combine_df_ori_valid$pheno_mis, range=range_df$range))%>% 
  mutate(cat="theo_mismatch")
est_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno_mis,range=range_df$range))%>% 
  mutate(cat="est_mismatch")
model_predskill<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno,range=range_df$range)) %>% 
  mutate(cat="model_predskill")
stats_df<-bind_rows(theo_mismatch,est_mismatch, model_predskill) %>% 
  gather(key="stats", value="value",-cat ) %>% 
  spread(key="cat", value="value") %>% 
  mutate(stats=factor(stats, levels=c("corr", "R2", "RMSE", "nRMSE"))) %>% 
  arrange(stats)
write_csv(stats_df,"./archive/m8.csv")

vis_df<-combine_df_ori %>% filter(site==3) %>% dplyr::select(date, pheno, pheno_mis, value, upper, lower)
pheno_ts <- xts(vis_df$pheno,vis_df$date) 
pheno_mis_ts <- xts(vis_df$pheno_mis,vis_df$date) 
pred_ts <- xts(vis_df$value,vis_df$date) 
pred_ts_upper <- xts(vis_df$upper,vis_df$date) 
pred_ts_lower <- xts(vis_df$lower,vis_df$date) 
vis_ts <- cbind(pheno_ts, pheno_mis_ts,pred_ts, pred_ts_upper, pred_ts_lower)
dygraph(vis_ts) %>% 
  dySeries("pheno_ts", label = "phenology without mismatch", col="blue") %>%
  dySeries("pheno_mis_ts", label = "phenology with mismatch", col='red') %>%
  dySeries(c("pred_ts_lower","pred_ts","pred_ts_upper"), label = "predicted phenology", col="black") %>%
  dyRangeSelector()


vis_df<-combine_df_ori %>% filter(site==3) %>% dplyr::select(date, mismatch_actual, mismatch_model, mismatch_model_upper, mismatch_model_lower)
act_mis_ts <- xts(vis_df$mismatch_actual,vis_df$date) 
est_mis_ts <- xts(vis_df$mismatch_model,vis_df$date) 
est_mis_ts_upper <- xts(vis_df$mismatch_model_upper,vis_df$date) 
est_mis_ts_lower <- xts(vis_df$mismatch_model_lower,vis_df$date) 
vis_ts <- cbind(act_mis_ts, est_mis_ts,est_mis_ts_upper, est_mis_ts_lower)
dygraph(vis_ts) %>% 
  dySeries("act_mis_ts", label = "actual mismatch", col="red") %>%
  dySeries(c("est_mis_ts_lower","est_mis_ts","est_mis_ts_upper"), label = "estimated mismatch", col="black") %>%
  dyRangeSelector()


cairo_pdf("./figure/m8.pdf")
grid.arrange(p1,p2,p3,p4,p5,p6,
             layout_matrix=rbind(c(1,2,2,2),
                                 c(3,5,5,5),
                                 c(4,6,6,6))
)
dev.off()
