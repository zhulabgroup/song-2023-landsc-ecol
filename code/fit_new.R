P_fit<-as.matrix(1:nrow(coord_df))


obs_df_new<-as.data.frame(x[P_fit,,1, drop=F]) %>%
  as_tibble() %>% 
  mutate(site=P_fit) %>%
  gather(key="date", value="y", -site) %>% 
  mutate(date=as.Date(date))
obs_df_ori_new<-obs_df_new %>% 
  left_join(df_upper_lower[[1]], by="site")%>% 
  mutate(y=y*range+lower)

######
fit_start<-max(unlist(lags))
fit_end<-length(date_list)

steps=fit_end-fit_start
D_fit <- date_list[(fit_start + 1):(fit_start + steps)] %>% 
  as.character() %>% 
  as.matrix()

xnew <- x
Sigmanew <- Sigma

Y_fit<-Var_fit<-matrix(NA, nrow=nrow(P_fit), ncol = steps )      
for (t in 1:steps) {
  
  res<-PrepareEmbedding(x,start=fit_start+t,end=fit_start+t, focalsites = P_fit, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newX<-res$X
  newP<-res$P
  newD<-res$D
  
  missing_id<-which(rowSums(is.na(newX))!=0)
  if (length(missing_id)>0) {
    newX<-newX[-missing_id,,drop=F]
    newP<-newP[-missing_id,,drop=F]
    newD<-newD[-missing_id,,drop=F]
  }
  res<-PrepareEmbedding(Sigma,start=fit_start+t,end=fit_start+t, focalsites = P_fit, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newSigma<-res$X
  if (length(missing_id)>0) {
    newSigma<-newSigma[-missing_id,,drop=F]
  }
  
  Y_fit_all<-
    foreach (i=1:num_part,
             .export = c( 
               "xnew","Sigmanew",
               "GPSDM","res"
             ),
             .packages = c("tidyverse","RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      particle_pick <- pars[, i, drop = F]
      res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis,basisD = D_basis, basisY = Y_basis, newX = newX, newP = newP,newD = newD, newSigma=newSigma, mode = c("predict"))
      cbind(res$mt, res$Ct)
    }
  for(i in 1:nrow(newP)) {
    means<-variances<-c()
    for (j in 1:num_part) {
      means<-cbind(means,Y_fit_all[[j]][i,1])
      variances<-cbind(variances,Y_fit_all[[j]][i,2])
    }
    means<-matrix(means)
    variances<-matrix(variances)
    weighted_mean<-weighted.mean(means,w=exp(log_p-mean(log_p)))
    weighted_variance<-weighted.mean(variances, w=exp(log_p-mean(log_p)))+Hmisc::wtd.var(means,weights=exp(log_p-mean(log_p)), normwt = T)
    #store prediction
    Y_fit[newP[i,1],t]<-weighted_mean
    Var_fit[newP[i,1],t]<-weighted_variance
    
    xnew[newP[i,1],(fit_start+t),1]<-weighted_mean
    Sigmanew[newP[i,1],(fit_start+t),1]<-weighted_variance
  }
  
  print(t)
}



fit_df_new<-as.data.frame(Y_fit) %>% 
  `names<-`(D_fit) %>% 
  mutate(site=P_fit) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))
var_fit_df_new<-as.data.frame(Var_fit) %>% 
  `names<-`(D_fit) %>% 
  mutate(site=P_fit) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))
dir.create(paste0(path,"analyses"))
write_csv(fit_df,paste0(path,"analyses/fit_new.csv"))
write_csv(var_fit_df,paste0(path,"analyses/var_fit_new.csv"))

fit_df_new<-left_join(fit_df_new,var_fit_df_new, by=c("site", "date")) %>% 
  dplyr::rename(value=value.x, variance=value.y) %>% 
  mutate(lower=value-1.96*sqrt(variance),
         upper=value+1.96*sqrt(variance))

fit_df_ori_new<-fit_df_new %>% 
  left_join(df_upper_lower[[1]], by="site",suffix = c("", ".scale"))%>% 
  mutate(value=(value+0.5)*range+lower.scale,
         upper=(upper+0.5)*range+lower.scale,
         lower=(lower+0.5)*range+lower.scale)

combined_df_new<-obs_df_new %>% 
  left_join(fit_df_new, by=c("site", "date")) %>% 
  drop_na()
combined_df_ori_new<-obs_df_ori_new %>% 
  dplyr::select(-lower, -upper,-range)%>% 
  left_join(fit_df_ori_new, by=c("site", "date")) %>% 
  dplyr::select(-lower.scale, -upper.scale,-range) %>% 
  drop_na()

ggplot()+
  # geom_line(data=ts_all, aes(x=time, y=level), col="blue")+
  # geom_line(data=fit_df_ori, aes(x=date, y=value), col="blue")+
  # geom_ribbon(data=fit_df_ori, aes(x=date, ymin=lower, ymax=upper), fill="blue",alpha=0.25)+
  geom_line(data=ts_all_new, aes(x=time, y=level))+
  geom_line(data=fit_df_ori_new, aes(x=date, y=value), col="red")+
  geom_ribbon(data=fit_df_ori_new, aes(x=date, ymin=lower, ymax=upper), fill="red",alpha=0.25)+
  theme_classic()
