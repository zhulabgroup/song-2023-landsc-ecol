range_df<-combine_df_ori %>% filter(year <=midyear) %>% 
  dplyr::summarize(lower=quantile(pheno, 0, na.rm=T),
                   upper=quantile(pheno, 1, na.rm=T)) %>% 
  mutate(range=as.numeric(upper-lower))

combine_df_ori_valid<-combine_df_ori %>% filter(year >midyear)
theo_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$pheno, pred_ori=combine_df_ori_valid$pheno_mis, range=range_df$range))%>% 
  mutate(cat="theo_mismatch")
est_mismatch<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno_mis,range=range_df$range))%>% 
  mutate(cat="est_mismatch")
model_predskill<-data.frame(compare_stats( obs_ori=combine_df_ori_valid$value, pred_ori=combine_df_ori_valid$pheno,range=range_df$range)) %>% 
  mutate(cat="model_predskill")

stats_df<-bind_rows(theo_mismatch,est_mismatch, model_predskill) %>% 
  mutate(cat=factor(cat, levels=c("theo_mismatch", "est_mismatch", "model_predskill"))) %>% 
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

pheno_doy<-combine_df_ori %>% 
  dplyr::select(date, site, simulated=pheno, predicted=value) %>% 
  gather(key="cat", value="value", -date, -site) %>% 
  mutate(year=as.numeric(format(date, "%Y")))%>%
  mutate(doy=as.numeric(format(date, "%j"))) %>% 
  left_join(param_all %>% dplyr::select(-param_mis), by=c("year","site")) %>% 
  mutate(cat=factor(cat, levels=c("simulated", "predicted"))) %>% 
  arrange(cat)
  
p2<-
  ggplot(pheno_doy %>%
           filter(site==3) 
  )+
  geom_line(aes(x=doy, y=value, group=year, col=env_summ), alpha=0.5)+
  theme_classic()+
  ylab("phenology")+
  theme(legend.position="bottom")+
  facet_wrap(.~cat, ncol=1)+
  scale_color_viridis_c(direction=-1)

colors <- c("m" = "blue",
            "m with mismatch" = "red")
p3<-
  ggplot(data=param_all %>% filter(site==3))+
  geom_line(aes(x=year, y=param, group=site, col="m"), lwd=2, alpha=0.5)+
  geom_line(aes(x=year, y=param_mis, group=site, col="m with mismatch"), lwd=2, alpha=0.5)+
  # geom_line(data=param_all ,aes(x=year, y=param, group=site, col=site))+
  theme_classic()+
  scale_color_manual(values = colors)+
  labs(x = "year",
       y = "m",
       color = "") +
  theme(legend.position="top") 

colors <- c("simulated phenology" = "blue",
            "simulated phenology with mismatch" = "red",
            "predicted phenology" = "black")
p4<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=pheno, col="simulated phenology"),alpha=0.5)+
  geom_line( aes(x=date, y=pheno_mis, col="simulated phenology with mismatch"), alpha=0.5)+
  geom_line(aes(x=date, y=value, col="predicted phenology"))+
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill="predicted phenology"),alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  guides(fill=F)+
  labs(x = "date",
       y = "phenology",
       color = "") +
  theme(legend.position="top") 

colors <- c("simulated mismatch" = "purple",
            "estimated mismatch" = "dark red",
            "predictive error" = "dark blue")
p5<-
  ggplot(combine_df_ori %>% filter(site==3))+
  geom_line(aes(x=date, y=mismatch_actual, col="simulated mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=mismatch_model, col="estimated mismatch"), alpha=0.5)+
  geom_line( aes(x=date, y=pred_error, col="predictive error"), alpha=0.5)+
  geom_ribbon(aes(x=date, ymin=mismatch_model_lower, ymax=mismatch_model_upper, fill="estimated mismatch"),alpha=0.25)+
  theme_classic()+
  geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  guides(fill=F)+
  labs(x = "date",
       y = "deviation",
       color = "") +
  theme(legend.position="top") 

cairo_pdf(paste0("./figure/",param,".pdf"), height=8.5, width=11)
grid.arrange(annotate_figure(p1, fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p2, fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p3, fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p4, fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p5, fig.lab = "e", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             layout_matrix=rbind(c(1,3,3,3),
                                 c(2,4,4,4),
                                 c(2,5,5,5))
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
  dySeries("pheno_ts", label = "simulated phenology", col="blue") %>%
  dySeries("pheno_mis_ts", label = "simulated phenology with mismatch", col='red') %>%
  dySeries(c("pred_ts_lower","pred_ts","pred_ts_upper"), label = "predicted phenology", col="black") %>%
  dyRangeSelector()


vis_df<-combine_df_ori %>% filter(site==3) %>% dplyr::select(date, mismatch_actual, mismatch_model, mismatch_model_upper, mismatch_model_lower,pred_error)
act_mis_ts <- xts(vis_df$mismatch_actual,vis_df$date) 
est_mis_ts <- xts(vis_df$mismatch_model,vis_df$date) 
est_mis_ts_upper <- xts(vis_df$mismatch_model_upper,vis_df$date) 
est_mis_ts_lower <- xts(vis_df$mismatch_model_lower,vis_df$date) 
pred_err_ts <- xts(vis_df$pred_error,vis_df$date) 
vis_ts <- cbind(act_mis_ts, est_mis_ts,est_mis_ts_upper, est_mis_ts_lower,pred_err_ts)
dygraph(vis_ts) %>% 
  dySeries("act_mis_ts", label = "simulated mismatch", col="purple") %>%
  dySeries(c("est_mis_ts_lower","est_mis_ts","est_mis_ts_upper"), label = "estimated mismatch", col="darkred") %>%
  dySeries("pred_err_ts", label = "predictive error", col="darkblue") %>%
  dyRangeSelector()

