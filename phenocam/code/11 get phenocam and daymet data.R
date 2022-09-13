
pheno_df<-read_csv(getURL(rois_df_select$one_day_summary[r]), skip=22) %>% 
  # ts <- get_pheno_ts(site = rois_df_select$site[r], vegType =rois_df_select$roitype[r], roiID = rois_df_select$sequence_number[r], type = '1day') %>% 
  dplyr::select(date,pheno=midday_gcc) %>% 
  group_by(date) %>% 
  summarise(pheno=mean(pheno)) %>% 
  # mutate(lon=rois_df_select$lon[r],lat=rois_df_select$lat[r],site = rois_df_select$site[r], roi_name = rois_df_select$roi_name[r]) %>% 
  mutate(date=as.Date(date)) 

daymet_df <- download_daymet(site = roi,
                          lat = rois_df_select$lat[r],
                          lon = rois_df_select$lon[r],
                          start = max(pheno_df %>%pull (date) %>% min() %>%  format("%Y") %>% as.integer()-1,1980),
                          end = min(pheno_df %>%pull (date) %>% max() %>%  format("%Y") %>% as.integer(),2020),
                          internal = TRUE,
                          simplify = TRUE) %>% 
  mutate(date=as.Date(yday, origin = as.Date(paste0(year,"-01-01")))) %>% 
  spread(key="measurement", value="value") %>% 
  dplyr::select(date, tmax=`tmax..deg.c.`, tmin=`tmin..deg.c.`, prcp=`prcp..mm.day.`) %>% 
  mutate(temp=(tmax+ tmin)/2) %>% 
  dplyr::select(date, temp, prcp)

ts<-pheno_df %>% 
  full_join(daymet_df, by="date") %>% 
  mutate(year=as.integer(format(date, "%Y")),
         site=1,
         pheno_sd=0.01) %>% 
  arrange(date)

coord_df<-data.frame(lon=0, lat=0)
distMat<-matrix(0)
date_list<-ts$date