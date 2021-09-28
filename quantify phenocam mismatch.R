# phenos<-get_phenos()
rois_df<-get_rois() 

source("code/utils.R")
library(daymetr)
source("code/settings.R")

type_list<-c("AG", "GR")
mismatch_list<-vector(mode = "list", length=length(type_list))

for (ty in 2:length(type_list)) {
  type<-type_list[ty]
  
    rois_df_type<-rois_df %>% 
      filter(roitype==type) %>% 
      filter(first_date<=as.Date("2015-12-31")) %>% 
      filter(last_date>=as.Date("2020-12-31")) %>%
      # filter(site_years>=6) %>%
      # filter(active==T) %>% 
    filter(lat< 83.3 & #49.3457868 & 
             lat> 6.6 & #24.7433195 &
             lon< -49 & #-66.9513812 &
             lon> -178.2 ) %>% #-124.7844079) %>% 
    filter(sequence_number%%1000==0)
      
  phenocam_ts<- vector(mode = "list", length = nrow(rois_df_type))
  for (i in 1:nrow(rois_df_type)) {
    try({
      ts <- get_pheno_ts(site = rois_df_type$site[i], vegType =rois_df_type$roitype[i], roiID = rois_df_type$sequence_number[i], type = '1day') %>% 
        dplyr::select(date,midday_gcc) %>% 
        mutate(lon=rois_df_type$lon[i],lat=rois_df_type$lat[i],site = rois_df_type$site[i], roi_name = rois_df_type$roi_name[i]) 
      phenocam_ts[[i]]<-ts
    })
  }
  phenocam_ts<-bind_rows(phenocam_ts) %>% 
    mutate(date=as.Date(date))
  
  coord_df<-data.frame(lon=0, lat=0)
  distMat<-matrix(0)
  
  cl <- makeCluster(num_part, outfile = "")
  registerDoSNOW(cl)
  
  roi_list<-phenocam_ts%>% 
    distinct(roi_name) %>% 
    unlist()
  
  nRMSE_fit<-nRMSE_pre<-rep(NA, length(roi_list))
  for (r in 1:length(roi_list)) {
    roi<-roi_list[[r]]
    basisnumber <-200 
    
    path<-paste0("./archive/phenocam/",type,"/",roi,"/")
    dir.create(path, recursive = T)
    pheno<-phenocam_ts %>% 
      filter(roi_name==roi) %>% 
      dplyr::select(date, pheno=midday_gcc) %>%  
      group_by(date) %>% 
      summarize(pheno=mean(pheno))
    
    daymet <- download_daymet(site = roi,
                              lat = phenocam_ts [phenocam_ts$roi_name==roi,][1,"lat"],
                              lon = phenocam_ts [phenocam_ts$roi_name==roi,][1,"lon"],
                              start = max(phenocam_ts %>% filter(roi_name==roi) %>% dplyr::select(date) %>%arrange(date) %>%  slice(1) %>% format("%Y"),1980),
                              end = min(phenocam_ts %>% filter(roi_name==roi) %>% dplyr::select(date) %>% arrange(desc(date) )%>% slice(1) %>% format("%Y"),2020),
                              internal = TRUE,
                              simplify = TRUE) 
    
    daymet_df<-daymet %>% 
      mutate(date=as.Date(yday, origin = as.Date(paste0(year,"-01-01")))) %>% 
      spread(key="measurement", value="value") %>% 
      dplyr::select(date, tmax=`tmax..deg.c.`, tmin=`tmin..deg.c.`, prcp=`prcp..mm.day.`) 
    
    daymet_df<-daymet_df %>% 
      mutate(env=(tmax+ tmin)/2) %>% 
      dplyr::select(date, env)
    
    ts_all<-pheno %>% 
      left_join(daymet_df, by="date") %>% 
      mutate(env_sd=0.01,
             year=as.integer(format(date, "%Y")),
             site=1,
             pheno_sd=0.01,
             pheno_mis=0)
    ts<-ts_all
    
    date_list<-ts_all$date
    source("./code/preprocess data.R")
    
    # use first half to train model
    midyear=2018
    
    source("./code/prepare embeddings.R")
    source("./code/train GP model.R")
    
    # predict for whole duration
    source("./code/fit.R")
    
    # output table and figure
    p<-
      ggplot(combine_df_ori)+
      geom_line(aes(x=date, y=pheno, col="observed phenology"),alpha=0.5)+
      geom_line(aes(x=date, y=value, col="predicted phenology"))+
      geom_ribbon(aes(x=date, ymin=pheno, ymax=pheno, fill="observed phenology"),alpha=0.25)+
      geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill="predicted phenology"),alpha=0.25)+
      theme_classic()+
      geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
      guides(fill="none")+
      labs(x = "date",
           y = "GCC",
           color = "") +
      theme(legend.position="top") 
    
    write_csv(combine_df_ori, paste0(path,"data.csv"))
    nRMSE_fit[r]<-combine_df_ori %>% 
      filter(year<=midyear) %>%
      drop_na() %>% 
      mutate(sq=(pred_error/(df_upper_lower[[1]]$range))^2)  %>% 
      dplyr::select(sq) %>% 
      unlist() %>% 
      mean() %>% 
      sqrt() 
    
    nRMSE_pre[r]<-combine_df_ori %>% 
      filter(year>midyear) %>%
      drop_na() %>% 
      mutate(sq=(pred_error/(df_upper_lower[[1]]$range))^2)  %>% 
      dplyr::select(sq) %>% 
      unlist() %>% 
      mean() %>% 
      sqrt() 
    
    cairo_pdf(paste0(path,"plot.pdf"))
    print(p)
    dev.off()
  }
  
  mismatch_list[[ty]]<-data.frame(nRMSE_fit, nRMSE_pre,roi=roi_list,type=type) %>% 
    mutate(nRMSE_change=nRMSE_pre-nRMSE_fit)
  
  closeAllConnections()
}

mismatch_df<-bind_rows(mismatch_list)
write_csv(mismatch_df, "./archive/phenocam/mismatch.csv")

ggplot(mismatch_df )+
  geom_segment(aes(x=as.factor(0), xend=as.factor(1), y=nRMSE_fit, yend=nRMSE_pre),
             arrow = arrow(length = unit(0.5, "cm")))+
  geom_label_repel(aes(x=as.factor(0), y=nRMSE_fit, label=roi_list), cex=3, col="red")+
  theme_classic()+
  facet_wrap(.~type)+
  scale_y_reverse()

ggplot(mismatch_df)+
  geom_boxplot(aes(x=type, y=nRMSE_change))+
  geom_point(aes(x=type, y=nRMSE_change), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=type, y=nRMSE_change, label=roi_list), cex=3, col="red")+
  theme_classic()+
  scale_y_reverse()

t.test(mismatch_df %>% filter(type=="AG") %>% pull(nRMSE_change),
       mismatch_df %>% filter(type=="GR") %>% pull(nRMSE_change))

# file <- c("./gpw_v4_population_density_rev11_30_min.nc")
# Human <- brick(file)[[1]]
 
