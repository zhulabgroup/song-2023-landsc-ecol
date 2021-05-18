double_logistics<-function (t, m1, m2, m3, m4, m5, m6, m7, sd) {
  d<-as.integer(format(t, "%j"))
  
  v <- m1 + (m2-m7*d)*(1/(1+exp((m3-d)/m4))-1/(1+exp((m5-d)/m6))) + rnorm (1,0,sd)
  
  return (v)
}

date_list<-seq(as.Date("2021-01-01"),as.Date("2021-12-31"), by=1)
coord_df <- data.frame(lon=0, lat=0)
distMat<-matrix(0)
level<-rep(NA, length(date_list))
for (i in 1:length(date_list)) {
  date<-date_list[i]
  level[i]<-double_logistics(date,
                             m1=0, # average greenness in winter
                             m2=1, # difference between summer and winter
                             m3=40, # spring onset
                             m4=10, # slope of curve in spring
                             m5=180, # fall offset
                             m6=20, # slope of curve in fall
                             m7=0, # summer greendown
                             sd = 0.1
  )
  # print(i)
}
ts_all<-data.frame(time=date_list, level=level, level_sd=0.1)
ggplot(ts_all)+
  geom_line(aes(x=time, y=level))+
  theme_classic()


path<-"./archive/"
source("code/utils.R")
source("code/settings.R")
source("code/preprocess data.R")
source("code/prepare embeddings.R")
source("code/train GP model.R")



level<-rep(NA, length(date_list))
for (i in 1:length(date_list)) {
  date<-date_list[i]
  level[i]<-double_logistics(date,
                             m1=0, # average greenness in winter
                             m2=1, # difference between summer and winter
                             m3=40, # spring onset
                             m4=20, # slope of curve in spring
                             m5=180, # fall offset
                             m6=20, # slope of curve in fall
                             m7=0, # summer greendown
                             sd = 0.1
  )
  # print(i)
}
ts_all_new<-data.frame(time=date_list, level=level, level_sd=0.1)
ggplot()+
  geom_line(data=ts_all, aes(x=time, y=level), col="blue")+
  geom_line(data=ts_all_new, aes(x=time, y=level), col="red")+
  theme_classic()

source("code/fit.R")
compare_stats( obs_ori=combined_df_ori$y, pred_ori=combined_df_ori$value,obs=combined_df$y,pred=combined_df$value)

source("code/preprocess data_new.R")
source("code/fit_new.R")
compare_stats( obs_ori=combined_df_ori_new$y, pred_ori=combined_df_ori_new$value,obs=combined_df_new$y,pred=combined_df_new$value)
