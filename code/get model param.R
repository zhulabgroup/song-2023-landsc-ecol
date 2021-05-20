param_all<-
  foreach (s = 1:nrow(coord_df)) %dopar% {
    if (param == "m2") {
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
        mutate(param_mis = case_when((year>midyear)~param+0.1)) %>%
        mutate(site=s)
    }
    if (param == "m3") {
      param_site<-ts_all %>% 
        filter(site==s) %>% 
        dplyr::select(env,year) %>% 
        group_by(year) %>% 
        slice_tail(n=90) %>% # last 90 days
        summarize(env_summ=mean(env)) %>% 
        ungroup() %>% 
        mutate(env_summ=lag(env_summ)) %>% 
        mutate(param=env_to_param (env=env_summ,
                                   lower=40,
                                   upper=160,
                                   steepness=-2,
                                   midpoint=1)) %>% 
        mutate(param_mis = case_when((year>midyear)~param+20)) %>%
        mutate(site=s)
    }
    if (param == "m4") {
      param_site<-ts_all %>% 
        filter(site==s) %>% 
        dplyr::select(env,year) %>% 
        group_by(year) %>% 
        slice_head(n=14)  %>% 
        summarize(env_summ=mean(env)) %>% 
        ungroup() %>% 
        mutate(env_summ=lag(env_summ)) %>%
        mutate(param=env_to_param (env=env_summ,
                                   lower=0,
                                   upper=80,
                                   steepness=3.5,
                                   midpoint=1.2)) %>% 
        mutate(param_mis = case_when((year>midyear)~param-15)) %>%
        mutate(site=s)
    }
    if (param == "m5") {
      param_site<-ts_all %>% 
        filter(site==s) %>% 
        dplyr::select(env,year) %>% 
        group_by(year) %>% 
        slice_head(n=120) %>% 
        summarize(env_summ=mean(env)) %>% 
        ungroup() %>% 
        mutate(param=env_to_param (env=env_summ,
                                   lower=180,
                                   upper=280,
                                   steepness=-1,
                                   midpoint=0.8)) %>% 
        mutate(param_mis = case_when((year>midyear)~param+0.2)) %>%
        mutate(site=s)
    }
    if (param =="m8") {
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
    }
    param_site
  }
param_all<-bind_rows(param_all)
# param_all

# ggplot(param_all)+
#   geom_line(aes(x=year, y=env_summ, group=site, col=site))+
#   theme_classic()

