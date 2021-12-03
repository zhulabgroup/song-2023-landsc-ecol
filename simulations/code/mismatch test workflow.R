path<-"./simulations/"
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
  source(paste0(path, "code/steps/23 train GP model.R"))
  
  # predict for whole duration
  source(paste0(path, "code/steps/24 fit.R"))
  
  # output table and figure
  source(paste0(path, "code/steps/25 output table and figure.R"))
}

closeAllConnections()

mismatch_df<-bind_rows(stats_list) %>% 
  mutate(stats=factor(stats, levels=c("corr", "R2", "RMSE", "nRMSE"))) %>% 
  mutate(param_v=case_when(param=="m2"~ "difference between summer and winter",
                           param=="m3"~ "timing of spring onset",
                           param=="m4"~ "slope of curve in spring",
                           param=="m8"~ "number of life cycles")) %>% 
  mutate(param_v=factor(param_v, levels=c("difference between summer and winter", "timing of spring onset", "slope of curve in spring", "number of life cycles") %>% rev())) 

write_csv(mismatch_df, paste0(path, "output/mismatch.csv"))

colors <- c("simulated mismatch" = "purple",
            "estimated mismatch" = "dark red",
            "predictive error" = "dark blue")

p<-ggplot(mismatch_df %>% filter(stats=="nRMSE"))+
  geom_point(aes(x=param_v, y=theo_mismatch, col="simulated mismatch"))+
  geom_point(aes(x=param_v, y=est_mismatch, col="estimated mismatch"))+
  geom_point(aes(x=param_v, y=model_predskill, col="predictive error"))+
  scale_color_manual(values = colors)+
  # facet_wrap(.~stats, scales="free_y")+
  labs(x = "type of mismatch",
       y = "nRMSE",
       color = "") +
  coord_flip()+
  theme_classic()

cairo_pdf(paste0(path, "output/nRMSE.pdf"), height=4, width=8)
print(p)
dev.off()
