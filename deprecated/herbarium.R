data<-read_csv("./crowd/hf309-01-crowdcurio.csv") %>% 
  filter(metric=="mean") %>% 
  filter(phenophase=="pfl") %>% 
  dplyr::select(species=name,
                date, 
                year,
                doy,
                lat=county.lat, 
                lon=county.lon,
                temp=spring
) %>%
  mutate(period=case_when(year>=1950~"late",
                          # (year>=1920 & year<1980) ~"middle",
                          TRUE~"early")) %>% 
  mutate(site=row_number())

# visualize data
ggplot (data)+
  geom_point(aes(x=lat, y=lon, col=species))+
  facet_wrap(.~species)+
  theme_classic()

sp_list<-data %>% pull(species) %>% unique()

for (sp in sp_list) {
  
}
data_sp<-data %>% 
  filter(species==sp)



ggplot(data)+
  geom_point(aes(x=temp, y=doy, col=period), alpha=0.3)+
  geom_smooth(aes(x=temp, y=doy, group=period, col=period), method="lm")+
  theme_classic()+
  facet_wrap(.~species)



data<-read_csv("./bird/73_species.csv") %>% 
  mutate(period=case_when(Year>=1995~"late",
                          # (year>=1920 & year<1980) ~"middle",
                          TRUE~"early")) %>% 
  rename(X=XEUREF,
         Y=YEUREF)

data %>% pull(Year) %>% hist()

sp_list<-data %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  pull(Species) 

ggplot(data %>% filter(Species==sp_list[4]))+
  geom_density(aes(x=Dayofyear, group=Year, col=Year))+
  theme_classic()+
  facet_wrap(.~BZ, scales = "free_y")+
  scale_color_viridis_c()

path_npn <- "./data/NPN/"
dir.create(path_npn, recursive = T)
rnpn::npn_phenophases() %>%
  filter(pheno_class_id == 1) %>%
  arrange(phenophase_name)

species_list <- rnpn::npn_species() %>%
  filter(kingdom == "Plantae") %>%
  rowwise() %>%
  filter(sum(species_type$Species_Type == "Calibration") > 0)
# https://www.usanpn.org/files/articles/developing_a_plant_profile.pdf

# run only once to download data
# for (i in 1:nrow(species_list)) {
#   npn_df <- rnpn::npn_download_site_phenometrics(
#     request_source = "YS",
#     years = as.character(1980:2021),
#     species_ids = species_list$species_id[i],
#     pheno_class_ids = 1,
#     climate_data = T
#   )
#   write_csv(npn_df, paste0(path_npn, species_list$species_id[i], ".csv"))
#   print(i)
# }

path_npn<-"/raid/users/ysong67/GitHub/2021-agu-adv/data/NPN/"
data_list <- vector(mode = "list")
for (i in 1:nrow(species_list)) {
  npn_df <- read_csv(paste0(path_npn, species_ids = species_list$species_id[i], ".csv"))
  
  if (nrow(npn_df) > 0) {
    data <- npn_df %>%
      dplyr::select(species_id,
                    genus,
                    species,
                    common_name,
                    lat = latitude,
                    lon = longitude,
                    year = mean_first_yes_year,
                    doy = mean_first_yes_doy,
                    tmax_spring,
                    tmin_spring,
                    prcp_spring
      ) %>%
      filter(
        doy != -9999,
        tmax_spring != -9999,
        tmin_spring != -9999,
        prcp_spring != -9999
      ) %>% 
      filter(doy >= quantile(doy, 0.005),
             doy <= quantile(doy, 0.975)
      ) %>% 
      mutate(period=case_when(year>=2015~"late",
                              TRUE~"early")) %>% 
      mutate(date=as.Date(paste0(year, "-01-01"))) %>% 
      mutate(site=row_number())
    
    ggplot(data %>% 
             mutate(id=row_number()) %>% 
             dplyr::select(id, doy, tmax_spring, tmin_spring, prcp_spring, period) %>% 
             gather(key="covariate", value="value",-id, -doy, -period))+
      geom_point(aes(x=value, y=doy, col=period))+
      geom_smooth(aes(x=value, y=doy, group=period, col=period), method="lm")+
      theme_classic()+
      facet_wrap(.~covariate, nrow=1, scales = "free_x")
    
    data %>% 
      ggplot()+
      geom_point(aes(x=lon, y=lat, col=doy))+
      theme_classic()+
      scale_color_viridis_c()+
      facet_wrap(.~period)
    
    summary(model <- lm(doy ~ tmax_spring + tmin_spring + prcp_spring, data=data %>% filter(period=="early")))
    
    data_predict<-bind_rows(data %>% 
                              filter(period=="early") %>% 
                              mutate(predict=predict(model, newdata=data %>% filter(period=="early"))),
                            data %>% 
                              filter(period=="late") %>% 
                              mutate(predict=predict(model, newdata=data %>% filter(period=="late")))
    ) %>% 
      mutate(residual=predict-doy) 
    
    ggplot(data_predict)+
      geom_point(aes(x=doy, y=predict),alpha=0.5)+
      # geom_smooth(aes(x=doy, y=predict),alpha=0.5, method="lm")+
      geom_abline(intercept = 0, slope=1, col="red")+
      ggpubr::stat_cor(
        aes(
          x = doy, y = predict, group = period,
          # col = interaction( metric, site),
          label = paste(..r.label.., ..rr.label.., sep = "*`,`~")
        ),
        # p.accuracy = 0.05,
        label.x.npc = "left",
        label.y.npc = "top",
        show.legend = F,
        col = "blue"
      ) +
      theme_classic()+
      facet_wrap(.~period)
    
    compare_stats( obs_ori=data_predict %>% filter(period=="early") %>% pull(doy), pred_ori=data_predict %>% filter(period=="early") %>% pull(predict))
    compare_stats( obs_ori=data_predict %>% filter(period=="late") %>% pull(doy), pred_ori=data_predict %>% filter(period=="late") %>% pull(predict))
    
    data_predict %>% 
      ggplot()+
      geom_point(aes(x=lon, y=lat, col=doy))+
      theme_classic()+
      scale_color_viridis_c()+
      facet_wrap(.~period)
    
    data_predict %>% 
      ggplot()+
      geom_histogram(aes(x=(residual)))+
      theme_classic()+
      facet_wrap(.~period)
    
    ggplot(data_predict %>% filter(period=="late") %>% filter(abs(residual)<=30))+
      geom_point(aes(x=lon, y=lat, col=residual), cex=2)+
      theme_classic()+
      colorspace::scale_color_continuous_diverging()#+
    # scale_color_viridis_c()
    
    ggplot(data_predict %>% filter(period=="late") %>% filter(abs(residual)<=30))+
      geom_point(aes(x=doy, y=residual, col=residual), cex=2)+
      theme_classic()+
      colorspace::scale_color_continuous_diverging()
    
  }
}


data

coord_df<-data %>% dplyr::select(lon, lat)
date_list<-data %>% pull(date) %>% unique() %>% sort()

ts<-data %>% 
  full_join(expand(.,date, site), by=c("date", "site")) %>% 
  dplyr::select(date, pheno=doy, tmin=tmin_spring, tmax=tmax_spring, prcp=prcp_spring, site=site) %>% 
  mutate(pheno_sd=0) %>% 
  arrange(site, date)

coord_df<-data %>% dplyr::select(lon, lat)
date_list<-data %>% pull(date) %>% unique()


source(paste0(path, "code/21 preprocess data.R"))

# use first half to train model
midyear=2015

source(paste0(path, "code/22 prepare embeddings.R"))
source(paste0(path, "code/23 train GP model.R"))

# predict for whole duration
source(paste0(path, "code/24 fit.R"))

source(paste0(path, "code/25 output table and figure.R"))
compare_stats( obs_ori=combine_df_ori %>% filter(period=="early") %>% pull(pheno), pred_ori=combine_df_ori %>% filter(period=="early") %>% pull(value))
compare_stats( obs_ori=combine_df_ori %>% filter(period=="late") %>% pull(pheno), pred_ori=combine_df_ori %>% filter(period=="late") %>% pull(value))
