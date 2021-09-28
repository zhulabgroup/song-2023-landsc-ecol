
phenos<-get_phenos()

head(phenos$site)


site_df<-data.frame(site=phenos$site,lat=phenos$lat,lon=phenos$lon, igbp=phenos$landcover_igbp,duration=as.Date(phenos$date_last)-as.Date(phenos$date_first))

site_df_crop<-site_df %>%
  filter(igbp==14) %>% # http://www.eomf.ou.edu/static/IGBP.pdf
  filter(duration>=365*8)

rois_df<-get_rois()

phenocam_ts<- vector(mode = "list", length = nrow(site_df_crop))
for (i in 1:nrow(site_df_crop)) {
  ts <- get_pheno_ts(site = site_df_crop$site[i], vegType =rois_df [rois_df$site==site_df_crop$site[i],][1,"roitype"], roiID = 1000, type = '1day') %>% 
    dplyr::select(date,midday_gcc) %>% 
    mutate(lon=site_df_crop$lon[i],lat=site_df_crop$lat[i],site = site_df_crop$site[i], roi = 1000) 
  phenocam_ts[[i]]<-ts
}
phenocam_ts_crop<-bind_rows(phenocam_ts) %>% 
  mutate(date=as.Date(date))

site_df_dbf<-site_df %>%
  filter(igbp==4) %>% 
  filter(duration>=365*8)

phenocam_ts<- vector(mode = "list", length = nrow(site_df_dbf))
for (i in 1:nrow(site_df_dbf)) {
  ts <- get_pheno_ts(site = site_df_dbf$site[i], vegType =rois_df [rois_df$site==site_df_dbf$site[i],][1,"roitype"], roiID = 1000, type = '1day') %>% 
    dplyr::select(date,midday_gcc) %>% 
    mutate(lon=site_df_dbf$lon[i],lat=site_df_dbf$lat[i],site = site_df_dbf$site[i], roi = 1000) 
  phenocam_ts[[i]]<-ts
}
phenocam_ts_dbf<-bind_rows(phenocam_ts) %>% 
  mutate(date=as.Date(date))


site_df_gs<-site_df %>%
  filter(igbp==10) %>% 
  filter(duration>=365*8)

phenocam_ts<- vector(mode = "list", length = nrow(site_df_gs))
for (i in 1:nrow(site_df_gs)) {
  try({
    ts <- get_pheno_ts(site = site_df_gs$site[i], vegType =rois_df [rois_df$site==site_df_gs$site[i],][1,"roitype"], roiID = 1000, type = '1day') %>% 
      dplyr::select(date,midday_gcc) %>% 
      mutate(lon=site_df_gs$lon[i],lat=site_df_gs$lat[i],site = site_df_gs$site[i], roi = 1000) 
    phenocam_ts[[i]]<-ts
  })
}
phenocam_ts_gs<-bind_rows(phenocam_ts) %>% 
  mutate(date=as.Date(date))
