path<-"./simulations/"
source(paste0(path, "code/steps/01 utils.R"))
source(paste0(path, "code/steps/02 settings.R"))

cl <- makeCluster(num_part, outfile = "")
registerDoSNOW(cl)

param_list<-c("m8", "m2","m3","m4")

for (param in param_list) {
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
