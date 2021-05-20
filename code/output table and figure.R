range_df<-combine_df_ori %>% filter(year <=midyear) %>% 
  dplyr::summarize(lower=quantile(pheno, 0.025),
                   upper=quantile(pheno, 0.975)) %>% 
  mutate(range=as.numeric(upper-lower))

combine_df_ori_valid<-combine_df_ori %>% filter(year >midyear)
theo_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$pheno, pred_ori=combine_df_ori_valid$pheno_mis, range=range_df$range))%>% 
  mutate(cat="theo_mismatch")
est_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno_mis,range=range_df$range))%>% 
  mutate(cat="est_mismatch")
model_predskill<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno,range=range_df$range)) %>% 
  mutate(cat="model_predskill")
stats_df<-bind_rows(theo_mismatch,est_mismatch, model_predskill) %>% 
  gather(key="stats", value="value",-cat ) %>% 
  spread(key="cat", value="value") %>% 
  mutate(stats=factor(stats, levels=c("corr", "R2", "RMSE", "nRMSE"))) %>% 
  arrange(stats)
write_csv(stats_df,paste0("./figure/",param,".csv"))

p1<-
  ggplot(param_all)+
  geom_line(aes(x=env_summ, y=param), col="blue")+
  # geom_line(aes(x=env_summ, y=param_mis), col="red")+
  theme_classic()+
  ylab("m")

p2<-
  ggplot(pheno_all %>%
           filter(site==3) %>% 
           mutate(year=as.numeric(format(date, "%Y")))%>%
           mutate( doy=as.numeric(format(date, "%j"))) %>% 
           left_join(param_all, by=c("year","site")) 
  )+
  geom_line(aes(x=doy, y=pheno,group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("pheno")+
  theme(legend.position="bottom")

p3<-
  ggplot(fit_df_ori %>%
           filter(site==3) %>% 
           mutate(year=as.numeric(format(date, "%Y")))%>%
           mutate( doy=as.numeric(format(date, "%j"))) %>% 
           left_join(param_all, by=c("year","site")) 
  )+
  geom_line(aes(x=doy, y=value,group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("pheno")+
  theme(legend.position="bottom")

colors <- c("m without mismatch" = "blue",
            "m with mismatch" = "red")
p4<-
  ggplot(data=param_all %>% filter(site==3))+
  geom_line(aes(x=year, y=param, group=site, col="m without mismatch"), lwd=2, alpha=0.5)+
  geom_line(aes(x=year, y=param_mis, group=site, col="m with mismatch"), lwd=2, alpha=0.5)+
  # geom_line(data=param_all ,aes(x=year, y=param, group=site, col=site))+
  theme_classic()+
  scale_color_manual(values = colors)+
  labs(x = "year",
       y = "m",
       color = "") +
  theme(legend.position="top") 

colors <- c("phenology without mismatch" = "blue",
            "phenology with mismatch" = "red",
            "predicted phenology" = "black")
p5<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=pheno, col="phenology without mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=pheno_mis, col="phenology with mismatch"), alpha=0.5)+
  geom_line(aes(x=date, y=value, col="predicted phenology"))+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  labs(x = "date",
       y = "phenology",
       color = "") +
  theme(legend.position="top") 

colors <- c("actual mismatch" = "red",
            "estimated mismatch" = "black")
p6<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=mismatch_actual, col="actual mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=mismatch_model, col="estimated mismatch"), alpha=0.5)+
  geom_ribbon(aes(x=date, ymin=mismatch_model_lower, ymax=mismatch_model_upper), fill="black",alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  labs(x = "date",
       y = "mismatch",
       color = "") +
  theme(legend.position="top") 

cairo_pdf(paste0("./figure/",param,".pdf"))
grid.arrange(annotate_figure(p1, fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p2, fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p3, fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p4, fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p5, fig.lab = "e", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p6, fig.lab = "f", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             layout_matrix=rbind(c(1,4,4,4),
                                 c(2,5,5,5),
                                 c(3,6,6,6))
)
dev.off()

vis_df<-combine_df_ori %>% filter(site==3) %>% dplyr::select(date, pheno, pheno_mis, value, upper, lower)
pheno_ts <- xts(vis_df$pheno,vis_df$date) 
pheno_mis_ts <- xts(vis_df$pheno_mis,vis_df$date) 
pred_ts <- xts(vis_df$value,vis_df$date) 
pred_ts_upper <- xts(vis_df$upper,vis_df$date) 
pred_ts_lower <- xts(vis_df$lower,vis_df$date) 
vis_ts <- cbind(pheno_ts, pheno_mis_ts,pred_ts, pred_ts_upper, pred_ts_lower)
dygraph(vis_ts) %>% 
  dySeries("pheno_ts", label = "phenology without mismatch", col="blue") %>%
  dySeries("pheno_mis_ts", label = "phenology with mismatch", col='red') %>%
  dySeries(c("pred_ts_lower","pred_ts","pred_ts_upper"), label = "predicted phenology", col="black") %>%
  dyRangeSelector()


vis_df<-combine_df_ori %>% filter(site==3) %>% dplyr::select(date, mismatch_actual, mismatch_model, mismatch_model_upper, mismatch_model_lower)
act_mis_ts <- xts(vis_df$mismatch_actual,vis_df$date) 
est_mis_ts <- xts(vis_df$mismatch_model,vis_df$date) 
est_mis_ts_upper <- xts(vis_df$mismatch_model_upper,vis_df$date) 
est_mis_ts_lower <- xts(vis_df$mismatch_model_lower,vis_df$date) 
vis_ts <- cbind(act_mis_ts, est_mis_ts,est_mis_ts_upper, est_mis_ts_lower)
dygraph(vis_ts) %>% 
  dySeries("act_mis_ts", label = "actual mismatch", col="red") %>%
  dySeries(c("est_mis_ts_lower","est_mis_ts","est_mis_ts_upper"), label = "estimated mismatch", col="black") %>%
  dyRangeSelector()
