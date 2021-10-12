library(tidyverse)
library(corrplot)
library(ggrepel)

data<-read_csv("./bethany_bottomtrawl_fall.csv") %>% 
  dplyr::select(-SEASON, -name, -LAT, -LON) %>% 
  gather(key="species", value="CPUE", -YEAR, -STRATUM, -NTOWS, -midlat, -midlon, -group) %>% 
  mutate(abundance=CPUE*NTOWS) %>% 
  mutate(site=paste0(midlon,"_", midlat))
  # group_by (YEAR, group, species) %>% 
  # summarize (abundance=sum(abundance),
  #            lat=mean(midlat),
  #            lon=mean(midlon)) %>% 
  # ungroup()

coord_df<-data %>% 
  distinct(midlon, midlat) %>% 
  arrange(midlat, midlon) %>%
  mutate(id=row_number())

data<-data %>% 
  left_join(coord_df, by=c("midlat", "midlon"))

cairo_pdf("./marine/map.pdf")
p<-ggplot(data)+
  geom_point(aes(x=midlon, y=midlat, col=group))+
  geom_text(aes(x=midlon+0.02, y=midlat+0.02, label=id))+
  theme_minimal()+
  xlab("longitude")+
  ylab("latitude")+
  guides(col=guide_legend(title=""))+
  coord_equal()
print(p)
dev.off()
# test<-data %>% filter(species==unique(data$species)[5]) %>% 
#   pull(abundance)  
# mean(test==0) # check with Casey about 0

sp_list<-unique(data$species)
corr_decay_list<-mean_corr_list<-rep(NA, length(sp_list))
for (i in 1:length(sp_list)) {
  sp<-sp_list[i]
  data_sp<-data %>%
    filter(species==sp_list[i]) %>% 
    arrange(midlat, midlon) %>% 
    dplyr::select(-STRATUM, -NTOWS, -group,-midlat, -midlon, -species, -CPUE) %>% 
    group_by(site,id) %>% 
    filter(sum(abundance)>1000) %>% 
    ungroup() 
  
  cairo_pdf(paste0("./marine/ts/",sp, ".pdf"))
  p<-ggplot(data_sp)+
    geom_line(aes(x=YEAR, y=abundance, col=site, group=site))+
    guides(col="none")+
    theme_classic()+
    xlab("year")
  print(p)
  dev.off()
  
  data_mat<- data_sp %>% 
    dplyr::select(-site) %>% 
    spread(key = "id", value="abundance") %>% 
    dplyr::select(-YEAR)
  
  res<-cor(data_mat,use="complete.obs")
  # variogram
  
  cairo_pdf(paste0("./marine/corr_mat/",sp, ".pdf"))
  corrplot(res, method="circle")
  dev.off()
  
  mean_corr_list[i]<-mean(res-diag(nrow=nrow(res), ncol=ncol(res)), na.rm = T)
  
  corr_df<-as.data.frame(res) %>% 
    rownames_to_column(var="start") %>% 
    gather(key="end", value = "correlation", -start) %>%
    mutate(start=as.integer(start),
           end=as.integer(end)) %>% 
    left_join(coord_df, by=c("start"="id")) %>% 
    rename(start_lat=midlat,
           start_lon=midlon) %>% 
    left_join(coord_df, by=c("end"="id")) %>% 
    rename(end_lat=midlat,
           end_lon=midlon) %>% 
    mutate(distance=sqrt((start_lat-end_lat)^2+(start_lon-end_lon)^2))
  
  cairo_pdf(paste0("./marine/corr_map/",sp, ".pdf"))
  p<-ggplot ()+
    geom_point(data=data,aes(x=midlon, y=midlat))+
    geom_segment(data=corr_df,aes(x=start_lon, y=start_lat, xend=end_lon, yend=end_lat, col=correlation))+
    scale_color_viridis_c(begin=0, end=1)+
    theme_classic()+
    coord_equal()
  print(p)
  dev.off()
  
  cairo_pdf(paste0("./marine/decay/",sp, ".pdf"))
  p<-ggplot(corr_df)+
    geom_point(aes(x=distance, y=correlation))+
    geom_smooth(aes(x=distance, y=correlation),method=loess)+
    geom_smooth(aes(x=distance, y=correlation),method=lm, col="red")+
    xlim(0,12)+
    theme_classic()
  print(p)
  dev.off()
  corr_decay_list[i]<-coef(lm(data=corr_df, correlation~distance))["distance"]
  
  print(i)
}

management=c(1, 1, 1, 1, 0,
          0, 0, 1, 1, 1,
          1, 0, 1, 0, 1,
          1)
# https://media.fisheries.noaa.gov/2021-04/Mid-Atlantic-Managed-Species.pdf
# https://media.fisheries.noaa.gov/2021-04/New-England-Managed-Species.pdf
# https://www.fisheries.noaa.gov/new-england-mid-atlantic/population-assessments/fishery-stock-assessments-new-england-and-mid-atlantic

sync_df<-data.frame(species=sp_list, corr=mean_corr_list,decay=corr_decay_list,
                    management=management) %>% 
  mutate(management=case_when(management==0~"unmanaged",
                           management==1~"managed")) %>% 
  arrange(species)

cairo_pdf("./marine/mean_corr.pdf")
p<-ggplot(sync_df)+
  geom_boxplot(aes(x=management, y=corr))+
  geom_point(aes(x=management, y=corr), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=corr, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("mean correlation")
print(p)
dev.off()

mean(sync_df %>% filter(management=="managed") %>% pull(corr))
sd(sync_df %>% filter(management=="managed") %>% pull(corr))

mean(sync_df %>% filter(management=="unmanaged") %>% pull(corr))
sd(sync_df %>% filter(management=="unmanaged") %>% pull(corr))

t.test(sync_df %>% filter(management=="managed") %>% pull(corr),
       sync_df %>% filter(management=="unmanaged") %>% pull(corr))


p<-ggplot(sync_df)+
  geom_boxplot(aes(x=management, y=decay))+
  geom_point(aes(x=management, y=decay), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=decay, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("correlation decay")

t.test(sync_df %>% filter(management=="managed") %>% pull(decay),
       sync_df %>% filter(management=="unmanaged") %>% pull(decay))


abundance_df<-data %>% 
  group_by(species, YEAR) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  do(broom::tidy(lm(abundance ~ YEAR, .))) %>%
  filter(term == "YEAR") %>%
  dplyr::select(species, roc = estimate, p = p.value) %>% 
  ungroup() %>% 
  mutate(management=management) %>% 
  mutate(management=case_when(management==0~"unmanaged",
                              management==1~"managed")) %>% 
  arrange(species)

cairo_pdf("./marine/abundance_roc.pdf")
p<-ggplot(abundance_df)+
  geom_boxplot(aes(x=management, y=roc))+
  geom_point(aes(x=management, y=roc), cex=2, col="red", pch=1)+
  geom_label_repel(aes(x=management, y=roc, label=species), cex=3, col="red")+
  theme_classic()+
  xlab("")+
  ylab("abundance rate of change")
print(p)
dev.off()

t.test(abundance_df %>% filter(management=="managed") %>% pull(roc),
       abundance_df %>% filter(management=="unmanaged") %>% pull(roc))

# does depth matter