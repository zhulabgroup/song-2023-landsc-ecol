P_fit<-as.matrix(1:nrow(coord_df))

# obs_df<-as.data.frame(x[P_fit,,1, drop=F]) %>%
#   as_tibble() %>% 
#   mutate(site=P_fit) %>%
#   gather(key="date", value="y", -site) %>% 
#   mutate(date=as.Date(date))
# obs_df_ori<-obs_df %>% 
#   left_join(df_upper_lower[[1]], by="site")%>% 
#   mutate(y=y*range+lower)

######
fit_start<-max(unlist(lags))
fit_end<-length(date_list)

steps=fit_end-fit_start
D_fit <- date_list[(fit_start + 1):(fit_start + steps)] %>% 
  as.character() %>% 
  as.matrix()

xnew <- x
Sigmanew <- Sigma

step_groups<-split(seq(1:steps) ,ceiling(seq(1:steps) / 100))

Y_fit<-Var_fit<-matrix(NA, nrow=nrow(P_fit), ncol = steps )      
for (g in 1:length(step_groups)) {
  t<-step_groups[[g]]
  res<-PrepareEmbedding(x,start=fit_start+min(t),end=fit_start+max(t), focalsites = P_fit, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newX<-res$X
  newP<-res$P
  newD<-res$D
  
  missing_id<-which(rowSums(is.na(newX))!=0)
  if (length(missing_id)>0) {
    newX<-newX[-missing_id,,drop=F]
    newP<-newP[-missing_id,,drop=F]
    newD<-newD[-missing_id,,drop=F]
    t<-t[-missing_id]
  }
  
  if (nrow(newX)>0) {
    res<-PrepareEmbedding(Sigma,start=fit_start+min(t),end=fit_start+max(t), focalsites = P_fit, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
    newSigma<-res$X
    if (length(missing_id)>0) {
      newSigma<-newSigma[-missing_id,,drop=F]
    }
    
    Y_fit_all<-
      foreach (i=1:num_part,
               .packages = c("tidyverse","RhpcBLASctl")
      ) %dopar% {
        blas_set_num_threads(1)
        omp_set_num_threads(1)
        particle_pick <- pars[, i, drop = F]
        res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis,basisD = D_basis, basisY = Y_basis, newX = newX, newP = newP,newD = newD, newSigma=newSigma, mode = c("predict"))
        cbind(res$mt, res$Ct)
      }
    
    logweights<-log_p_valid-max(log_p_valid)
    for(s in 1:nrow(P_fit)) {
      means<-variances<-c()
      for (i in 1:num_part) {
        means<-cbind(means,Y_fit_all[[i]][newP==s,1])
        variances<-cbind(variances,Y_fit_all[[i]][newP==s,2])
      }
      res_df<-bind_rows(as.data.frame(means) %>% mutate(stat="means") %>% mutate(t=row_number()) ,
      as.data.frame(variances) %>% mutate(stat="variances") %>% mutate(t=row_number()))%>% 
        gather(key="particle", value="value",-t,-stat) %>% 
        spread(key="stat", value="value") %>% 
        group_by(t) %>% 
        summarise(weighted_mean=weighted.mean(means,w=exp(logweights)),
                  weighted_variance=weighted.mean(variances, w=exp(logweights))+Hmisc::wtd.var(means,weights=exp(logweights), normwt = F, method = "ML"))
      
      #store prediction
      Y_fit[P_fit[s,1],t]<-res_df$weighted_mean
      Var_fit[P_fit[s,1],t]<-res_df$weighted_variance
      
      xnew[P_fit[s,1],(fit_start+t),1]<-res_df$weighted_mean
      Sigmanew[P_fit[s,1],(fit_start+t),1]<-res_df$weighted_variance
    }
  }
  
  print(paste0(g, " out of ", length(step_groups)))
}



fit_df<-as.data.frame(Y_fit) %>% 
  `names<-`(D_fit) %>% 
  mutate(site=P_fit) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))
var_fit_df<-as.data.frame(Var_fit) %>% 
  `names<-`(D_fit) %>% 
  mutate(site=P_fit) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))

path_analyses<-paste0(path_sub,"analyses")
dir.create(path_analyses,recursive = T)
write_csv(fit_df,paste0(path_analyses,"/fit.csv"))
write_csv(var_fit_df,paste0(path_analyses,"/var_fit.csv"))

fit_df<-left_join(fit_df,var_fit_df, by=c("site", "date")) %>% 
  dplyr::rename(value=value.x, variance=value.y) %>% 
  mutate(lower=value-1.96*sqrt(variance),
         upper=value+1.96*sqrt(variance))

fit_df_ori<-fit_df %>% 
  left_join(df_upper_lower[[1]], by="site",suffix = c("", ".scale"))%>% 
  mutate(value=(value+0.5)*range+lower.scale,
         upper=(upper+0.5)*range+lower.scale,
         lower=(lower+0.5)*range+lower.scale)

climatology_df<-ts %>% 
  mutate(doy=format( date, "%j") %>% as.numeric(),
         year=format (date, "%Y") %>% as.numeric()) %>% 
  filter(year <=midyear) %>% 
  group_by(doy) %>%
  summarise (pheno_clim=mean(pheno, na.rm=T)) %>% 
  filter(doy!=365)
  
  
  combine_df_ori<-ts %>% 
    mutate(doy=format( date, "%j") %>% as.numeric()) %>% 
    left_join(climatology_df, by="doy") %>% 
    left_join(fit_df_ori %>% dplyr::select(date,site, value, lower, upper), by=c("date","site")) %>% 
    mutate(pred_error=value-pheno) %>% 
    dplyr::select(-doy)



