path<-"./simulations/"
set.seed(1)
source(paste0(path, "code/steps/01 utils.R"))
source(paste0(path, "code/steps/02 settings.R"))

cl <- makeCluster(num_part, outfile = "")
registerDoSNOW(cl)

param_list<-c("m8", "m2","m3","m4")

stats_list<-vector(mode="list", length=length(param_list))
for (p in 1:length(param_list)) {
  param<-param_list[p]
  path_sub<-paste0(path, "archive/",param,"/")
  dir.create(path, recursive = T)
  
  # prepare env time series
  source(paste0(path, "code/steps/11 prepare env ts.R"))
  
  # get phenology model parameter in each year (match and mismatch)
  source(paste0(path, "code/steps/12 get model param.R"))
  
  # get phenology time series (match and mismatch)
  source(paste0(path, "code/steps/13 get pheno ts.R"))
  
  # preprocess data
  ts<-ts_all
  source(paste0(path, "code/steps/21 preprocess data.R"))
  
  # use first half to train model
  source(paste0(path, "code/steps/22 prepare embeddings.R"))
  source(paste0(path, "code/steps/23 train GP model.R")) # make output look better
  
  # predict for whole duration
  source(paste0(path, "code/steps/24 fit.R"))
  
  # output table and figure
  source(paste0(path, "code/steps/25 output table and figure.R"))
}

stopCluster(cl)

ts_df_list<-vector(mode="list", length=length(param_list))
for (p in 1:length(param_list)) {
  param<-param_list[p]
  ts_df_list[[p]]<-read_csv(paste0(path, "output/",param,".csv"),
                            col_types=list(pheno_mis=col_double(),
                                           mismatch_actual=col_double(),
                                           mismatch_model=col_double(),
                                           mismatch_model_upper=col_double(),
                                           mismatch_model_lower=col_double())) %>% 
    mutate(param=param) %>% 
    rowwise() %>% 
    mutate(param_v=param_name[[param]])
}
ts_df<-bind_rows(ts_df_list)%>% 
  mutate(param_v=factor(param_v, levels=c("Summer-winter difference", "Timing of spring onset", "Slope of curve in spring", "Number of life cycles") )) 

colors <- c("simulated mismatch" = "purple",
            "estimated mismatch" = "dark red",
            "predictive error" = "dark blue")
p_ts<-
  ggplot(ts_df %>% filter(site==3)%>% filter(date>=as.Date(paste0(midyear+1,"-01-01") )))+
  geom_line(aes(x=date, y=mismatch_actual, col="simulated mismatch"),alpha=0.5)+
  geom_line( aes(x=date, y=mismatch_model, col="estimated mismatch"), alpha=0.5)+
  # geom_line( aes(x=date, y=pred_error, col="predictive error"), alpha=0.5)+
  geom_ribbon(aes(x=date, ymin=mismatch_model_lower, ymax=mismatch_model_upper, fill="estimated mismatch"),alpha=0.25)+
  theme_classic()+
  # geom_vline(xintercept =as.Date(paste0(midyear+1,"-01-01") ), alpha=0.5)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  guides(fill=F)+
  labs(x = "Date",
       y = "Phenological mismatch",
       color = "") +
  theme(legend.position="top") +
  facet_wrap(.~param_v, ncol=1, scales = "free_y")
p_ts

mismatch_df<-bind_rows(stats_list) %>% 
  mutate(stats=factor(stats, levels=c("corr", "R2", "RMSE", "nRMSE"))) %>% 
  rowwise() %>% 
  mutate(param_v=param_name[[param]]) %>% 
  mutate(param_v=factor(param_v, levels=c("Summer-winter difference", "Timing of spring onset", "Slope of curve in spring", "Number of life cycles") %>% rev())) 

write_csv(mismatch_df, paste0(path, "output/mismatch.csv"))

mismatch_df<-read_csv(paste0(path, "output/mismatch.csv"))
mismatch_df %>% filter(stats=="nRMSE") %>% View()
colors <- c("simulated mismatch" = "purple",
            "estimated mismatch" = "dark red",
            "predictive error" = "dark blue")

p_stat<-ggplot(mismatch_df %>% filter(stats=="nRMSE"))+
  geom_point(aes(x=param_v, y=theo_mismatch, col="simulated mismatch"))+
  geom_point(aes(x=param_v, y=est_mismatch, col="estimated mismatch"))+
  # geom_point(aes(x=param_v, y=model_predskill, col="predictive error"))+
  scale_color_manual(values = colors)+
  # facet_wrap(.~stats, scales="free_y")+
  labs(x = "type of mismatch",
       y = "nRMSE",
       color = "") +
  ylim(0, 0.2)+
  coord_flip()+
  theme_classic()
p_stat

cairo_pdf(paste0(path, "output/simulated and estimated mismatch.pdf"), height=8, width=8)
grid.arrange(
  annotate_figure(p_ts, fig.lab = "A"),
  annotate_figure(p_stat, fig.lab = "B"),
  heights=c(0.75, 0.25),ncol=1) %>% 
  print()
dev.off()
