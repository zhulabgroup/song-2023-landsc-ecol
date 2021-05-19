path<-"./archive/m2/"
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
      slice_head(n=30) %>% 
      summarize(env_summ=mean(env)) %>% 
      ungroup() %>% 
      mutate(param=env_to_param (env=env_summ,
                                 lower=0.2,
                                 upper=1,
                                 steepness=-3,
                                 midpoint=0.5)) %>% 
      mutate(param_mis = if_else((year<=midyear),param,param+0.2)) %>%
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

p2<-
  ggplot()+
  geom_line(data=param_all %>% filter(site==3), aes(x=year, y=param_mis, group=site), col="red", lwd=2, alpha=0.5)+
  geom_line(data=param_all %>% filter(site==3), aes(x=year, y=param, group=site), col="blue", lwd=2, alpha=0.5)+
  geom_line(data=param_all ,aes(x=year, y=param, group=site, col=site))+
  theme_classic()+
  ylab("m")

# param_all<-ts_all %>% 
#   dplyr::select(env,year) %>% 
#   group_by(year) %>% 
#   slice_head(n=21) %>% # first 21 days
#   summarize(env_summ=mean(env)) %>% 
#   ungroup() %>% 
#   mutate(param=env_to_param (env=env_summ,
#                              vertex_x=0.2,
#                              vertex_y=1,
#                              curvature=-0.3)) %>% 
#   mutate(param_mis = if_else((year<=midyear),param,param-0.3)) #%>%
# # mutate(param_mis = if_else(param>=25,param,25))

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
                                   m2=param_all$param[y], # difference between summer and winter
                                   m3=100, # spring onset
                                   m4=10, # slope of curve in spring
                                   m5=260, # fall offset
                                   m6=20, # slope of curve in fall
                                   m7=0, # summer greendown
                                   sd = 0.05
        )
        pheno_mis[i]<-double_logistics(date,
                                       m1=0, # average greenness in winter
                                       m2=param_all$param_mis[y], # difference between summer and winter
                                       m3=100, # spring onset
                                       m4=10, # slope of curve in spring
                                       m5=260, # fall offset
                                       m6=20, # slope of curve in fall
                                       m7=0, # summer greendown
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

p3<-
  ggplot(fit_df_ori %>%
           filter(site==3) %>% 
           mutate(year=as.numeric(format(date, "%Y")))%>%
           mutate( doy=as.numeric(format(date, "%j"))) %>% 
           left_join(param_all, by=c("year","site")) #%>% 
         # filter(year%in%year_zoom)
  )+
  geom_line(aes(x=doy, y=value,group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("pheno")#+
# facet_wrap(.~site, ncol = 1)

combine_df_ori<-ts_all %>% 
  left_join(fit_df_ori %>% dplyr::select(date,site, value, lower, upper), by=c("date","site")) %>% 
  mutate(mismatch_actual=pheno-pheno_mis,
         mismatch_model=value-pheno_mis,
         mismatch_model_upper=upper-pheno_mis,
         mismatch_model_lower=lower-pheno_mis)

p4<-
  ggplot(combine_df_ori%>% filter(site==3))+
  geom_line(aes(x=date, y=pheno), col="blue",alpha=0.5)+
  geom_line( aes(x=date, y=pheno_mis), col="red", alpha=0.5)+
  geom_line(aes(x=date, y=value), col="black")+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  ylab("pheno")#+
  # facet_wrap(.~site, ncol = 1)

p5<-
  ggplot(combine_df_ori%>% filter(site==3))+
  geom_line(aes(x=date, y=mismatch_actual), col="red",alpha=0.5)+
  geom_line( aes(x=date, y=mismatch_model), col="black", alpha=0.5)+
  geom_ribbon(aes(x=date, ymin=mismatch_model_lower, ymax=mismatch_model_upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  ylab("mismatch")#+
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
  dplyr::summarize(lower=quantile(pheno, 0.025, na.rm = T),
                   upper=quantile(pheno, 0.975, na.rm = T)) %>% 
  mutate(range=as.numeric(upper-lower))

combine_df_ori_valid<-combine_df_ori %>% filter(year >midyear)
compare_stats( obs_ori=combine_df_ori_valid$pheno, pred_ori=combine_df_ori_valid$pheno_mis, range=range_df$range)
compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno_mis,range=range_df$range)
compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno,range=range_df$range)

cairo_pdf("/figure/m2.pdf")
grid.arrange(p1,p2,p3,p4,p5, 
             layout_matrix=rbind(c(1,3,3),
                                 c(2,2,2),
                                 c(4,4,4),
                                 c(5,5,5))
)
dev.off()

