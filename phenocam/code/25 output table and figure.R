# range_df<-combine_df_ori %>% filter(year <=midyear) %>% 
#   dplyr::summarise(lower=quantile(pheno, 0, na.rm=T),
#                    upper=quantile(pheno, 1, na.rm=T)) %>% 
#   mutate(range=as.numeric(upper-lower))

combine_df_ori_fit<-combine_df_ori %>% filter(year <=midyear)
combine_df_ori_fore<-combine_df_ori %>% filter(year >midyear)

  stats_list[[r]]<-bind_rows(data.frame(compare_stats( obs_ori=combine_df_ori_fit$value, pred_ori=combine_df_ori_fit$pheno,range=df_upper_lower[[1]]$range)) %>% 
    mutate(cat="fit"),
    data.frame(compare_stats( obs_ori=combine_df_ori_fore$value, pred_ori=combine_df_ori_fore$pheno,range=df_upper_lower[[1]]$range)) %>% 
    mutate(cat="fore"), 
    data.frame(compare_stats( obs_ori=combine_df_ori_fit$pheno_clim, pred_ori=combine_df_ori_fit$pheno,range=df_upper_lower[[1]]$range)) %>% 
    mutate(cat="fit_clim"),
    data.frame(compare_stats( obs_ori=combine_df_ori_fore$pheno_clim, pred_ori=combine_df_ori_fore$pheno,range=df_upper_lower[[1]]$range)) %>% 
      mutate(cat="fore_clim"))%>% 
    mutate(cat=factor(cat, levels=c("fit","fit_clim", "fore", "fore_clim"))) %>% 
    gather(key="stats", value="value",-cat ) %>% 
    spread(key="cat", value="value") %>% 
    mutate(stats=factor(stats, levels=c("corr", "R2", "RMSE", "nRMSE"))) %>% 
    arrange(stats) %>% 
    mutate(change=fore-fit) %>% 
  mutate(type=type,
         roi=roi) 
  

  p<-
    ggplot(combine_df_ori)+
    geom_line(aes(x=date, y=pheno, col="observed phenology"),alpha=0.5)+
    geom_line(aes(x=date, y=value, col="predicted phenology"))+
    geom_line(aes(x=date, y=pheno_clim, col="climatology"))+
    geom_ribbon(aes(x=date, ymin=pheno, ymax=pheno, fill="observed phenology"),alpha=0.25)+
    geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill="predicted phenology"),alpha=0.25)+
    theme_classic()+
    geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
    guides(fill="none")+
    labs(x = "date",
         y = "GCC",
         color = "") +
    theme(legend.position="top") 
  
  # ggplot(combine_df_ori %>% drop_na(value) %>% mutate(doy=as.numeric(format(date, "%j"))))+
  #   geom_line(aes(x=doy, y=pheno, lty="observed phenology", group=year, col=year))+
  #   geom_line(aes(x=doy, y=value, lty="predicted phenology", group=year, col=year))+
  #   geom_ribbon(aes(x=doy, ymin=pheno, ymax=pheno, fill=year, group=year),alpha=0.25)+
  #   geom_ribbon(aes(x=doy, ymin=lower, ymax=upper, fill=year, group=year),alpha=0.25)+
  #   theme_classic()+
  #   guides(fill="none")+
  #   labs(x = "date",
  #        y = "GCC",
  #        color = "") +
  #   theme(legend.position="top")+
  #   scale_color_viridis_c()+
  #   facet_wrap(.~year, ncol=1)
  
  dir.create(paste0(path, "output/eachroi"), recursive = T)
  cairo_pdf(paste0(path,"output/eachroi/", type, "_", roi,".pdf"))
  print(p)
  dev.off()

