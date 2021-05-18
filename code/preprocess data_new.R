x<-array(NA, dim= c(nrow(coord_df),length(date_list),length(var_list)),
         dimnames = list(as.character(1:nrow(coord_df)),
                         as.character(date_list),
                         var_list))
for (j in 1:length(date_list)) { #time
  for (v in 1:length(var_list)) { #covariate
    if (var_list[v] %in% c("level")) {
      level_date<-ts_all_new %>% filter(time==date_list[j])
      if (nrow(level_date)>0) {
        x[,j,v]<-level_date$level
      } else {
        x[,j,v]<-rep(NA, nrow(coord_df))
      }
    }
    if (var_list[v]=="doy") {
      x[,j,v]<-rep(sin(as.numeric(format(date_list[j], "%j"))*2*pi), nrow(coord_df))
    }
  }
  print(date_list[j])
}

# export_path<-paste0(path, "raw")
# dir.create(export_path, recursive = T)
# for (i in 1:nrow(coord_df)) {
#   write_csv(as.data.frame(x[i,,]), paste0(export_path, "/", i, ".csv"))
# }

Sigma<-array(NA, dim= c(nrow(coord_df),length(date_list),length(var_list)),
             dimnames = list(as.character(1:nrow(coord_df)),
                             as.character(date_list),
                             var_list))
for (j in 1:length(date_list)) { #time
  for (v in 1:length(var_list)) { #covariate
    if (var_list[v] %in% c("level")) {
      level_date<-ts_all_new %>% filter(time==date_list[j])
      if (nrow(level_date)>0) {
        Sigma[,j,v]<-(level_date$level_sd)^2
      } else {
        Sigma[,j,v]<-rep(NA, nrow(coord_df))
      }
    }
    else {
      Sigma[,j,v]<-rep(0, nrow(coord_df))
    }
  }
  print(date_list[j])
}


# for (i in 1:nrow(coord_df)) {
#   for (j in 1:length(date_list)) { 
#     if(!is.na(Sigma[i,j,1])) {
#       if (Sigma[i,j,1]>0.0001) {
#         Sigma[i,j,1]<-NA
#         x[i,j,1]<-NA
#       }
#     }
#   }
# }

### to percentage###
date_id<-1:length(seq(min(ts_all_new$time),max(ts_all_new$time), by = 1))  #using percentile in NEON data
df_upper_lower<-vector(mode="list")
for(j in 1:length(var_list)) {
  if (var_list[j]%in% c("level")) {
    df_upper_lower[[j]]<-data.frame(x[,,j,drop=F]) %>% 
      mutate(site=row_number()) %>% 
      gather(key="date", value = "value",-site) %>% 
      drop_na() %>% 
      group_by(site) %>% 
      dplyr::summarize(lower=quantile(value, 0.025),
                       upper=quantile(value, 0.975)) %>% 
      mutate(range=upper-lower)
  } else { #scale for all sites
    all_upper_lower<-data.frame(x[,,j,drop=F]) %>% 
      mutate(site=row_number()) %>% 
      gather(key="date", value = "value",-site) %>% 
      drop_na() %>% 
      dplyr::summarize(lower=quantile(value, 0.025),
                       upper=quantile(value, 0.975)) %>% 
      mutate(range=upper-lower)
    df_upper_lower[[j]]<-data.frame(x[,,j]) %>% 
      mutate(site=row_number()) %>% 
      gather(key="date", value = "value",-site) %>% 
      drop_na() %>% 
      distinct(site) %>% 
      mutate(lower=all_upper_lower$lower,
             upper=all_upper_lower$upper,
             range=all_upper_lower$range)
  }
  
  lower<-matrix(df_upper_lower[[j]]$lower)%*%matrix(1, nrow=1, ncol=ncol(x[,,j,drop=F]) )
  range<-matrix(df_upper_lower[[j]]$range)%*%matrix(1, nrow=1, ncol=ncol(x[,,j,drop=F]) )
  
  x[,,j]<-(x[,,j]-lower)/range-0.5
}



### scaling
for(j in 1:length(var_list)) {
  Sigma[,,j]<-Sigma[,,j,drop=F]/(df_upper_lower[[j]]$range)^2
}
# 
# ##### linear interpolation
# for(j in 1:length(var_list)) {
#   for (i in 1:nrow(coord_df)) {
#     min_id<-min(which(!is.na(x[i,,j])))
#     max_id<-max(which(!is.na(x[i,,j])))
#     x[i,min_id:max_id,j]<-zoo::na.approx(object=x[i,min_id:max_id,j], x=as.Date(names(x[i,min_id:max_id,j])),maxgap=14)
#   }
# }
# 
# # linear interpolation
# for(j in 1:length(var_list)) {
#   for (i in 1:nrow(coord_df)) {
#     min_id<-min(which(!is.na(Sigma[i,,j])))
#     max_id<-max(which(!is.na(Sigma[i,,j])))
#     Sigma[i,min_id:max_id,j]<-zoo::na.approx(object=Sigma[i,min_id:max_id,j], x=as.Date(names(Sigma[i,min_id:max_id,j])),maxgap=14)
#   }
# }
# 
# x_raw<-x
# Simga_raw<-Sigma
# 
# # whittaker smoothing
# for(j in 1:length(var_list)) {
#   for (i in 1:nrow(coord_df)) {
#     max_id<-0
#     done<-F
#     while(!done) {
#       min_id<-min(which(!is.na(x[i,(max_id+1):length(date_list),j])))+(max_id)
#       if (min_id==Inf) {
#         done<-T
#       } else {
#         max_id<-min(which(is.na(x[i,min_id:length(date_list),j])))-1+(min_id-1)
#         if (max_id==Inf) {
#           max_id<-length(date_list)
#           done<-T
#         }
#         x[i,min_id:max_id,j]<-ptw::whit1(x[i,min_id:max_id,j],5) #gcc
#       }
#     }
#   }
# }
# # 
# export_path<-paste0(path,"processed")
dir.create(paste0(path, "scaling/"))
for (j in 1:length(var_list)) {
  write_csv(df_upper_lower[[j]], paste0(path, "scaling/", j, "_new.csv"))
  print(j)
}
# for (i in 1:nrow(coord_df)) {
#   write_csv(as.data.frame(x[i,,]), paste0(export_path, "/", i, ".csv"))
#   print(i)
# }

####



# # whittaker smoothing
# for(j in 1:length(var_list)) {
#   for (i in 1:nrow(coord_df)) {
#     max_id<-0
#     done<-F
#     while(!done) {
#       min_id<-min(which(!is.na(Sigma[i,(max_id+1):length(date_list),j])))+(max_id)
#       if (min_id==Inf) {
#         done<-T
#       } else {
#         max_id<-min(which(is.na(Sigma[i,min_id:length(date_list),j])))-1+(min_id-1)
#         if (max_id==Inf) {
#           max_id<-length(date_list)
#           done<-T
#         }
#         Sigma[i,min_id:max_id,j]<-ptw::whit1(Sigma[i,min_id:max_id,j],5) #gcc
#       }
#     }
#   }
# }

#### visualize time series
df_list<-vector(mode="list", length(var_list))
for (i in 1:length(var_list)) {
  df_list[[i]]<-as.data.frame(x[,,i]) %>% 
    mutate(site=1:nrow(coord_df)) %>% 
    gather(key = "date", value="value", -site) %>% 
    mutate(var=var_list[i]) %>% 
    mutate(date=as.Date(date))
}
df<-bind_rows(df_list)

var_df_list<-vector(mode="list", length(var_list))
for (i in 1:length(var_list)) {
  var_df_list[[i]]<-as.data.frame(Sigma[,,i]) %>% 
    mutate(site=1:nrow(coord_df)) %>% 
    gather(key = "date", value="value", -site) %>% 
    mutate(var=var_list[i]) %>% 
    mutate(date=as.Date(date))
}
var_df<-bind_rows(var_df_list)

df<-left_join(df,var_df, by=c("site", "date", "var")) %>% 
  dplyr::rename(value=value.x, variance=value.y) %>% 
  mutate(lower=value-1.96*sqrt(variance),
         upper=value+1.96*sqrt(variance))

ggplot(df %>% filter(var=="gcc") )+
  geom_line(aes(x=date, y=value))+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper), fill="blue", alpha=0.5)+
  theme_classic()+
  facet_wrap(~site*var)