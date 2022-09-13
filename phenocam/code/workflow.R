path<-"./phenocam/"
source(paste0(path, "code/01 utils.R"))
source(paste0(path, "code/02 settings.R"))

# Sys.setenv(CURL_CA_BUNDLE = file.path(Sys.getenv("R_HOME"), "lib/microsoft-r-cacert.pem"))
rois_df<-get_rois() 

type_list<-rois_df %>% 
  group_by(roitype) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  head(10) %>% 
  pull(roitype)

rois_df_select<-rois_df %>% 
  filter(roitype%in%type_list) %>% 
  filter(first_date<=as.Date("2015-12-31")) %>% 
  filter(last_date>=as.Date("2020-12-31")) %>%
  # filter(site_years>=6) %>%
  # filter(active==T) %>% 
  filter(lat< 83.3 & #49.3457868 & 
           lat> 6.6 & #24.7433195 &
           lon< -49 & #-66.9513812 &
           lon> -178.2 ) %>% #-124.7844079) %>% 
  filter(sequence_number%%1000==0)

unique(rois_df_select$roitype)

cl <- makeCluster(num_part, outfile = "")
registerDoSNOW(cl)

stats_list_b<-stats_list
stats_list<-vector(mode="list", length=nrow(rois_df_select))
for (r in 45:nrow(rois_df_select)) {
  type<-rois_df_select$roitype[r]
  roi<-rois_df_select$roi_name[r]
  
  path_sub<-paste0(path, "archive/",type,"/",roi,"/")
  dir.create(path_sub, recursive = T)
  source(paste0(path, "code/11 get phenocam and daymet data.R"))
  
  source(paste0(path, "code/21 preprocess data.R"))
  
  # use first half to train model
  midyear=2018
  
  source(paste0(path, "code/22 prepare embeddings.R"))
  source(paste0(path, "code/23 train GP model.R"))
  
  # predict for whole duration
  source(paste0(path, "code/24 fit.R"))
  # compare with climatology
  
  # output table and figure
  source(paste0(path, "code/25 output table and figure.R"))
}

closeAllConnections()


mismatch_df<-bind_rows(stats_list)
write_csv(mismatch_df, paste0(path, "output/mismatch.csv"))

mismatch_df<-read_csv(paste0(path, "output/mismatch.csv"))
mismatch_df<-mismatch_df %>% 
  left_join(rois_df_select, by=c("roi"="roi_name")) %>% 
  filter(!type %in% c("EB", "WL"))


mismatch_df %>% 
  filter(stats=="nRMSE") %>% 
  mutate(clim_greater=fore_clim-fore) %>%
  arrange(desc(clim_greater)) %>% 
  head()

better_sites<-mismatch_df %>% 
  filter(stats=="nRMSE") %>% 
  mutate(clim_greater=fore_clim-fore) %>% 
  filter(clim_greater>0) %>% 
  pull(roi) %>% 
  unique()

cairo_pdf( paste0(path, "output/map.pdf"))
p<-ggplot(mismatch_df)+
  geom_point(aes(x=lon, y=lat, col=type))+
  theme_minimal()+
  xlab("longitude")+
  ylab("latitude")+
  guides(col=guide_legend(title=""))+
  coord_equal()
print(p)
dev.off()

for (focal_stats in c("corr", "R2", "RMSE", "nRMSE")) {
  cairo_pdf( paste0(path, "output/fit_fore_",focal_stats,".pdf"))
  p<-ggplot(mismatch_df %>% filter(stats==focal_stats))+
    geom_segment(aes(x=as.factor("fitted" ), xend=as.factor("forecasted"), y=fit, yend=fore),
                 arrow = arrow(length = unit(0.5, "cm")))+
    # geom_segment(aes(x=as.factor("fitted" ), xend=as.factor("forecasted"), y=fit_clim, yend=fore_clim),
    #              col="red",
    #              arrow = arrow(length = unit(0.5, "cm")))+
    # geom_label_repel(aes(x=as.factor("fitted"), y=fit, label=roi), cex=3, col="red")+
    theme_classic()+
    facet_wrap(.~type)+
    # scale_y_reverse()+
    xlab("")+
    ylab(focal_stats)
  print(p)
  dev.off()
}


cairo_pdf( paste0(path, "output/stats_change.pdf"))
p<-ggplot(mismatch_df  %>% 
            # filter(type %in% c("GR", "AG")) 
          filter(stats=="nRMSE") %>% 
            filter(roi %in% better_sites) #%>% 
            # mutate(change=fore_clim-fit_clim)
          )+
  geom_boxplot(aes(x=type, y=change))+
  geom_point(aes(x=type, y=change), cex=2, col="red", pch=1)+
  # geom_label_repel(aes(x=type, y=change, label=roi), cex=3, col="red")+
  theme_classic()+
  # scale_y_reverse()+
  xlab("")+
  ylab("stats change")+
  facet_wrap(.~stats, ncol=1, scales = "free_y")
print(p)
dev.off()


anova.model <- lm(change ~ type, data=mismatch_df %>% filter(stats=="nRMSE") %>% filter(roi %in% better_sites))
anova(anova.model)
a1 <- aov(change ~ type, data=mismatch_df %>% filter(stats=="nRMSE") %>% filter(roi %in% better_sites))
posthoc <- TukeyHSD(x=a1, 'type', conf.level=0.90)
plot(posthoc, las=1)

t.test(mismatch_df %>% filter(stats=="nRMSE") %>% filter(type=="AG") %>% pull(change),
       mismatch_df %>% filter(stats=="nRMSE") %>% filter(type=="GR") %>% pull(change))
# t.test(mismatch_df %>% filter(stats=="RMSE") %>% filter(type=="AG") %>% pull(change),
#        mismatch_df %>% filter(stats=="RMSE") %>% filter(type=="GR") %>% pull(change))
# t.test(mismatch_df %>% filter(stats=="corr") %>% filter(type=="AG") %>% pull(change),
#        mismatch_df %>% filter(stats=="corr") %>% filter(type=="GR") %>% pull(change))
# t.test(mismatch_df %>% filter(stats=="R2") %>% filter(type=="AG") %>% pull(change),
#        mismatch_df %>% filter(stats=="R2") %>% filter(type=="GR") %>% pull(change))
 
