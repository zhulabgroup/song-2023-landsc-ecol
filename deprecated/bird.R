library(raster)
library(daymetr)

whitfun <- function(x, lambda) {
  max_id <- 0
  done <- F
  while (!done) {
    min_id <- min(which(!is.na(x[(max_id + 1):length(x)]))) + (max_id) # first number that is not NA
    if (min_id == Inf) { # all numbers are NA
      done <- T # consider this ts done
    } else {
      max_id <- min(which(is.na(x[min_id:length(x)]))) - 1 + (min_id - 1) # last number in the first consecutive non-NA segment
      if (max_id == Inf) {
        max_id <- length(x) # last non-NA segment is at the end of the whole ts
        done <- T # consider this ts done
      }
      x[min_id:max_id] <- ptw::whit1(x[min_id:max_id], lambda) # whitman smoothing for this non-NA segment
    }
  }
  return(x)
}

# read in data
data_df<-read_csv("./bird/73_species.csv")  %>% 
  dplyr::select(
    nestid=NestID,
    X=XEUREF,
    Y=YEUREF,
    species_short=Species,
    year=Year,
    doy=Dayofyear,
    area=BZ) %>% 
  mutate(period=case_when(year>=1995~"late",
                          # (year>=1920 & year<1980) ~"middle",
                          TRUE~"early")) %>% 
  mutate(date=as.Date(paste0(year, "-01-01"))+doy-1)


# join with trait data
sp_df<-read_csv("./bird/Traits_73_species.csv") %>% 
  dplyr::select(
    species_short=Abbreviation,
    species=`Scientific name`,
    mig=Mig
  )

data_df<-data_df %>% 
  left_join(sp_df, by="species_short")

sp_list<-data_df %>% 
  group_by(species) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  pull(species) 

# project to lon lat
coord_df<-data_df %>% 
  distinct(X, Y) %>% 
  mutate(coordid=row_number())
coord_sp<-SpatialPointsDataFrame(coords=coord_df[,c("X", "Y")],data=coord_df,
                              proj4string = CRS("EPSG:3067"))

coord_sp_reproj<-spTransform(coord_sp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coord_df_reproj<-as.data.frame(coord_sp_reproj) %>% 
  rename(lon=X.1, lat=Y.1)

data_df<-data_df %>% 
  left_join(coord_df_reproj, by=c("X", "Y"))

ggplot(data_df%>% filter(species==sp_list[4]) %>% filter(doy>=season_start, doy<=season_end))+
  geom_point(aes(x=lon, y=lat, col=doy))+
  theme_classic()+
  facet_wrap(.~species)+
  coord_equal()+
  scale_color_viridis_c()

summary(coord_df_reproj)
# filter for species and years

# get CHELSA data
## download

# cl <- makeCluster(12, outfile = "")
# registerDoSNOW(cl)
# 
# chelsa_list<-vector(mode="list")
# for (year in 1979:2016) {
#   files<-list.files("/data/ZHULAB/phenology/CHELSA/daily/", pattern = year %>% as.character(), full.names = T)
#   chelsa_month_list<-foreach (file = files, 
#                        .packages = c("raster", "stringr", "tidyverse")) %dopar% {
#     bri<-brick(file)
#     # plot(bri[[1:3]])
#     raster::extract(bri,coord_sp_reproj) %>% 
#       data.frame() %>% 
#       mutate(coordid=row_number()) %>% 
#       gather(key="date", value="temp", -coordid) %>% 
#       mutate(temp=temp-272.15) %>% 
#       mutate(date=str_replace(date, "X", "")) %>% 
#       mutate(date=str_replace_all(date, "\\.", "-")) %>% 
#       mutate(date=as.Date(date))
#       
#                        }
#   chelsa_list[[year]]<-bind_rows(chelsa_month_list)
#   print(year)
# }
# stopCluster(cl)
# chelsa_df<-bind_rows(chelsa_list) %>% 
#   left_join(coord_df_reproj, by="coordid") %>% 
#   left_join(data_df %>% 
#               distinct(lon, lat, area), 
#             by=c("lat", "lon"))
# write_rds(chelsa_df, "./bird/chelsa.rds")

chelsa_df<-read_rds( "./bird/chelsa.rds")

chelsa_df_area<-chelsa_df  %>% 
  group_by(area, date) %>% 
  summarise(temp=mean(temp, na.rm=T)) %>% 
  mutate(temp=whitfun(temp, lambda = 30)) %>% 
  drop_na(area) 

# ggplot(chelsa_df_area )+
#   geom_line(aes(x=date, y=temp, group=area, col=area))+
#   theme_classic()

ggplot(chelsa_df_area %>% 
         mutate(doy=format(date, "%j") %>% as.numeric(),
                year=format(date, "%Y") %>% as.numeric()))+
  geom_line(aes(x=doy, y=temp, group=year, col=year), alpha=0.3)+
  theme_classic()+
  scale_color_viridis_c()+
  facet_wrap(.~area)


ggplot(chelsa_df_area %>% 
         mutate(doy=format(date, "%j") %>% as.numeric(),
                year=format(date, "%Y") %>% as.numeric()) %>% 
         # filter(area=="SB") %>% 
         filter(doy==110))+
  geom_line(aes(x=year, y=temp))+
  geom_smooth(aes(x=year, y=temp), method="loess")+
  theme_classic()+
  facet_wrap(.~area, scales = "free")

# get bird phenology data
for (sp in sp_list) {
  
}
## summarizing by area
sp=sp_list[2]
pheno_df<-data_df %>%
  filter(species == sp) %>% 
  # filter(area=="HB") %>%
  group_by(area, date) %>% 
  summarise(count=n()) %>% 
  ungroup() %>% 
  right_join(chelsa_df_area %>% distinct( area, date), by=c("area", "date")) %>% 
  arrange(date, area) %>% 
  mutate(count=replace_na(count, 0)) %>% 
  mutate(doy=format(date, "%j") %>% as.numeric(),
         year=format(date, "%Y") %>% as.numeric()) %>% 
  right_join(group_by(.,area, year) %>% 
               summarise(sum=sum(count)) %>% 
               filter(sum>=100) %>% # remove some site and years with low counts
               distinct(area, year),
             by=c("area", "year")) %>% 
  right_join(chelsa_df_area %>% distinct( area, date), by=c("area", "date")) %>% 
  arrange(date, area) %>% 
  group_by(area, year) %>% 
  mutate(count=whitfun(count, lambda = 30)) %>% 
  mutate(freq=count/sum(count)) %>% 
  mutate(freq=freq/max(freq)) %>%
  ungroup() %>% 
  mutate(doy=format(date, "%j") %>% as.numeric())
  
ggplot(pheno_df )+
  geom_line(aes(x=doy,y=freq, group=year, col=year), alpha=0.5)+
  geom_vline(xintercept = season_start)+
  geom_vline(xintercept = season_end)+
  geom_vline(xintercept = 135, lty=2)+
  theme_classic()+
  facet_wrap(.~area, scales = "free_y")+
  scale_color_viridis_c()

ggplot(pheno_df %>% 
         filter(doy==110))+
  geom_line(aes(x=year, y=freq))+
  geom_smooth(aes(x=year, y=freq), method="loess")+
  theme_classic()+
  facet_wrap(.~area, scales = "free_y")


ggplot(pheno_df)+
  geom_line(aes(x=date, y=freq))+
  facet_wrap(.~area)+
  theme_classic()

season_start<-pheno_df %>% filter(freq>=0.01) %>% pull(doy) %>% min()-5
# season_end<-pheno_df %>% filter(freq==1) %>% pull(doy) %>% max()+5
season_end<-pheno_df %>% filter(freq>=0.01) %>% pull(doy) %>% max()+5

ts<-pheno_df %>% 
  right_join(chelsa_df_area, by=c("area", "date")) %>% 
  dplyr::select(date,doy, pheno=freq, temp,site=area) %>% 
  mutate(pheno_sd=0) %>% 
  arrange(site, date) %>% 
  mutate(pheno=case_when((doy>=season_start & doy<=season_end)~pheno)) %>% 
  dplyr::select(-doy)

coord_df<-coord_df_reproj %>% 
  left_join(data_df %>% 
              distinct(lon, lat, area),
            by=c("lat", "lon")) %>% 
  rename(site=area) %>% 
  group_by(site) %>% 
  summarise(lon=mean(lon, na.rm=T), lat=mean(lat, na.rm=T))
date_list<-ts %>% pull(date) %>% unique()

source(paste0(path, "code/21 preprocess data.R"))

# use first half to train model
midyear=1995

path_sub<-paste0(path, "archive/",sp,"/")
source(paste0(path, "code/22 prepare embeddings.R"))

cl <- makeCluster(num_part, outfile = "")
registerDoSNOW(cl)
source(paste0(path, "code/23 train GP model.R"))
stopCluster(cl)
# predict for whole duration
source(paste0(path, "code/24 fit.R"))


# normalize predictions to freq?
combine_df_ori<-combine_df_ori %>% 
  group_by(site, year) %>% 
  mutate(value=value*(1/sum(value, na.rm = T))) %>% 
  mutate(value=value/max(value)) %>%
  ungroup() %>% 
  mutate(doy=format(date, "%j") %>% as.numeric()) %>% 
  mutate(value=case_when((doy>=season_start & doy<=season_end)~value)) 

# # plot
# ggplot(combine_df_ori %>% filter(year>=1980))+
#   geom_line(aes(x=date, y=pheno, col="observed phenology"),alpha=0.5)+
#   geom_line(aes(x=date, y=value, col="predicted phenology"))+
#   # geom_point(aes(x=date, y=pheno_clim, col="climatology"))+
#   # geom_ribbon(aes(x=date, ymin=pheno, ymax=pheno, fill="observed phenology"),alpha=0.25)+
#   # geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill="predicted phenology"),alpha=0.25)+
#   theme_classic()+
#   geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") )-1, alpha=0.5)+
#   guides(fill="none")+
#   labs(x = "date",
#        y = "pheno",
#        color = "") +
#   theme(legend.position="top") +
#   facet_wrap(.~site)

ggplot(combine_df_ori %>% filter(year>=1980))+
  # geom_line(aes(x=doy,y=pheno, group=year, col=year), alpha=0.5)+
  geom_line(aes(x=doy,y=value, group=year, col=year), alpha=0.5)+
  theme_classic()+
  facet_wrap(.~site, scales = "free_y")+
  scale_color_viridis_c()+ 
  xlim(c(season_start, season_end))


ggplot(combine_df_ori %>% 
         filter(year>=1980) %>%  
         filter(doy==150))+
  # geom_line(aes(x=year, y=pheno), col="dark green")+
  # geom_smooth(aes(x=year, y=pheno), method="lm", col="dark green")+
  geom_line(aes(x=year, y=value), col="blue")+
  geom_smooth(aes(x=year, y=value), method="lm", col="blue")+
  theme_classic()+
  facet_wrap(.~site, scales = "free_y")

ggplot(combine_df_ori %>% 
         filter(year>=1980) %>%  
         filter(doy==135))+
  geom_point(aes(x=value, y=pheno))+
  geom_smooth(aes(x=value, y=pheno), method="lm")+
  theme_classic()+
  facet_wrap(.~site, scales = "free")

### calculate mismatch