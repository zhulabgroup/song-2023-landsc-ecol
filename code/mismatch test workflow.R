source("code/utils.R")
source("code/settings.R")

cl <- makeCluster(num_part, outfile = "")
registerDoSNOW(cl)
set.seed(1)

param_list<-c("m2","m3", "m4","m8")

for (param in param_list) {
  path<-paste0("./archive/",param,"/")
  dir.create(path, recursive = T)
  
  # prepare env time series
  source("./code/prepare env ts.R")
  
  # get phenology model parameter in each year (match and mismatch)
  source("./code/get model param.R")
  
  # get phenology time series (match and mismatch)
  source("./code/get pheno ts.R")
  
  # preprocess data
  ts<-ts_all
  source("./code/preprocess data.R")
  
  # use first half to train model
  source("./code/prepare embeddings.R")
  source("./code/train GP model.R")
  
  # predict for whole duration
  source("./code/fit.R")
  
  # output table and figure
  source("./code/output table and figure.R")
}

closeAllConnections()
