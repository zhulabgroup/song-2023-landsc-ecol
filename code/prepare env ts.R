date_list<-seq(as.Date("2021-01-01"),as.Date("2040-12-31"), by=1)
midyear<-floor(mean(c(as.numeric(format(min(date_list), "%Y")),as.numeric(format(max(date_list), "%Y")))))
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
                    sd = env_sd
      )
      # print(i)
    }
    ts_site<-data.frame(date=date_list, env=env, env_sd=env_sd)%>% 
      mutate(year=as.numeric(format(date, "%Y"))) %>% 
      mutate(site=s)
    ts_site
  }
ts_all<-bind_rows(ts_all)

# p<-
#   ggplot(ts_all)+
#   geom_line(aes(x=date, y=env, group=site, col=site))+
#   theme_classic()+
#   theme(legend.position="right")
# cairo_pdf("./figure/env.pdf", width = 11, height=4)
# print (p)
# dev.off()
