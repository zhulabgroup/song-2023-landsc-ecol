library(parallel)
library(doSNOW)
library(tidyverse)

terraclim_path<-"/nfs/turbo/seas-zhukai/climate/TerraClimate/individual_years/"
scratch_path<-"/scratch/zhukai_root/zhukai0/songyl/"
var_list<-c("tmax", "tmin")

cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)
for (var in var_list) {
  foreach (year = seq(1958, Sys.Date() %>% lubridate::year()) ) %dopar% {
    url<-paste0("http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_",var,"_",year,".nc")
    system(paste0("wget ", url, " -P ", terraclim_path, " -nv"))
    
  }
}

foreach (year=1958:2021,
         .packages = c("raster", "ncdf4")) %dopar%  {
           rasterOptions(tmpdir=paste0(scratch_path,"tmpraster")) #important, otherwise might run out of memory
           tmax_file <- paste0(terraclim_path,"TerraClimate_tmax_",year,".nc")
           monthly_tmax<-stack(tmax_file,varname="tmax")
           tmin_file <- paste0(terraclim_path,"TerraClimate_tmin_",year,".nc")
           monthly_tmin<-stack(tmin_file,varname="tmin")
           monthly_mean<-(monthly_tmax+monthly_tmin)/2
           mat <- mean(monthly_mean)
           
           dir.create(paste0(terraclim_path, "metric/"))
           raster::writeRaster(mat, paste0(terraclim_path,"metric/MAT_1_24degree_",year,".tif"), overwrite=TRUE, format="GTiff")
         }
unlink(paste0(scratch_path,"tmpraster"))

stopCluster(cl)