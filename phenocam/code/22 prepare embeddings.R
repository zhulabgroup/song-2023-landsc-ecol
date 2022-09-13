end_id<-which (date_list==as.Date(paste0(midyear+1,"-01-01"))-1)
dir.create(paste0(path_sub,"train"),recursive = T)
res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=end_id, focalsites =1:nrow(coord_df), lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
X_train<-res$X
Y_train<-res$Y
D_train<-res$D
P_train<-res$P

write_csv(as.data.frame(X_train),paste0(path_sub,"train/X.csv"))
write_csv(as.data.frame(Y_train),paste0(path_sub,"train/Y.csv"))
write_csv(as.data.frame(D_train),paste0(path_sub,"train/D.csv"))
write_csv(as.data.frame(P_train),paste0(path_sub,"train/P.csv"))
