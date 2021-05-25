library(tidyverse)
library(LaplacesDemon)
library(doSNOW)
library(lubridate)
library(gridExtra)
library(xts)
library(ggpubr)

waves<-function (t, t_start,
                 intercept, slope, 
                 amplitude1, phase1, period1,
                 amplitude2, phase2, period2,
                 sd=0.05) {
  # d<-as.integer(format(t, "%j"))
  t_diff<-as.numeric(t-t_start)
  if (leap_year(t)) {d_all<-366} else {d_all<-365}
  v <- intercept + slope* t_diff+
    amplitude1 * sin(2*pi/(period1*d_all) *(t_diff+phase1)) + 
    amplitude2 * sin(2*pi/(period2*d_all) *(t_diff+phase2)) + 
    rnorm (1,0,sd)
  
  return (v)
}

env_to_param<-function (env, lower, upper, steepness, midpoint) {
  param<-(upper-lower)/(1+exp(-steepness*(env-midpoint)))+lower
  return(param)
}
  
double_logistics<-function (t,
                            m1=0,# average greenness in winter
                            m2=1, # difference between summer and winter
                            m3=100, # spring onset
                            m4=10, # slope of curve in spring
                            m5=260, # fall offset
                            m6=20, # slope of curve in fall
                            m7=0, # summer greendown
                            m8=1, #life cycle
                            sd=0.05
                            ) {
  if (!leap_year(t)) {
    d<-(as.integer(format(t, "%j")) %% (365/round(m8)))*round(m8)
  } else {
    d<-(as.integer(format(t, "%j")) %% (366/round(m8)))*round(m8)
  }
  
  
  v <- m1 + (m2-m7*d)*(1/(1+exp((m3-d)/m4))-1/(1+exp((m5-d)/m6))) + rnorm (1,0,sd)
  
  return (v)
}

compare_stats<-function ( obs_ori, pred_ori,obs=NULL,pred=NULL, range=NULL) {
  corr<-cor(obs_ori, pred_ori, use="pairwise.complete.obs")
  R2<-summary(lm(pred_ori~obs_ori))$r.squared
  RMSE<-sqrt(mean((obs_ori-pred_ori)^2, na.rm=T))
  
  if(!is.null(obs)) {
    nRMSE<-sqrt(mean((obs-pred)^2, na.rm=T))
  } 
  else if (!is.null (range)) {
    nRMSE <- RMSE / range
  }
  else {
    nRMSE<-NA
  }
  
  out<-list(corr=corr, R2=R2, nRMSE=nRMSE, RMSE=RMSE)
  return(out)
}

# https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

tdistMat<-matrix(NA, nrow=365, ncol=365) 
for (i in 1:365) {
  for (j in 1:365) {
    tdistMat[i,j]<-min(abs(i-j), 365-abs(i-j))
  }
}
colnames(tdistMat)<-1:365
tdistMat_df<-as.data.frame(tdistMat) %>% 
  mutate(doy1=1:365) %>% 
  gather(key="doy2", value="value", -doy1) %>% 
  mutate(doy1=as.integer(doy1),
         doy2=as.integer(doy2))

# ggplot(tdistMat_df,aes(x=doy1, y=doy2))+
  # geom_raster(aes(fill=value))+
  # scale_fill_viridis_c()

date2doy<-function (x){
  doy_all<-matrix(NA, nrow=nrow(x), ncol=1)
  for (i in 1:nrow (x)) {
    doy<-as.integer(format(as.Date(x[i,]),"%j"))
    if (doy==366) {doy<-365}
    doy_all[i,]<-doy
  }
  return (doy_all)
}

PrepareEmbedding<-function(x,start, end, focalsites=NULL, lags, neighbors, vars,distMat) {
  numPatch<-dim(x)[[1]]
  Amat<-distMat
  
  # vars<-dim(x)[[3]]
  if(is.null(focalsites)) {
    focalsites<-seq(1,numPatch)
  }
  
  Xi_list<-Yi_list<-Pi_list<-Di_list<-vector(mode = "list", length=length(focalsites))
  for (i in focalsites) {
    Ind<-order(Amat[i,])
    Ind<-c(i,setdiff(Ind,i)) #make sure focal site is always the first
    
    Xi<-matrix(NA, nrow=(end-start+1),ncol=0)
    for (v in 1:length(vars)) {
      var<-var_list[vars[v]]
      for (n in 1:length(neighbors[[v]])) {
        neighbor<-neighbors[[v]][n]
        site<-Ind[neighbor]
        if (!is.null(lags[[v]])) {
          for (period in 1:length(lags[[v]])) {
            period_lags<-lags[[v]][[period]]
            values<-matrix(NA, nrow=end-start+1, ncol=length(period_lags))
            for (l in 1:length(period_lags)) {
              period_lag<-period_lags[l]
              values[, l]<-x[site,(start-period_lag):(end-period_lag),var]
            }
            Xi_add<-as.matrix(rowMeans(values, na.rm = T))
            colnames(Xi_add)<-paste0(var,"_",neighbor,"_",period)
            Xi<-cbind(Xi, Xi_add)
          }
        }
      }
    }
    Xi_list[[i]]<-Xi
    
    Yi<-matrix(x[i,start:end,1,drop=F]) #Assuming the first var is always same as the var to predict
    Yi_list[[i]]<-Yi
    
    Pi<-matrix(i, nrow=(end-start+1), ncol=1)
    Pi_list[[i]]<-Pi
    
    Di<-matrix(dimnames(x[i,start:end,1,drop=F])[[2]])
    Di_list[[i]]<-Di
  }
  
  X<-do.call(rbind,Xi_list)
  Y<-do.call(rbind,Yi_list)
  P<-do.call(rbind,Pi_list)
  D<-do.call(rbind,Di_list)
  
  output<-list(X=X,Y=Y,P=P, D=D)
  return(output)
}

fmingrad_Rprop<-function (fun, xinit, sample_n, maxcount) {
  #fun is a handle to a function that has grad as optional second output
  #this uses the sign of the gradient to determine descent directions and an adaptive step size - supposedly smoother convergence than conjugate gradients for GP optimization
  p<-length(xinit)
  
  # optimization parameters for Rprop
  Delta0 <- 0.1*matrix(1, nrow=p, ncol=1)
  Deltamin <- 1e-6*matrix(1, nrow=p, ncol=1)
  Deltamax <- 50*matrix(1, nrow=p, ncol=1)
  eta_minus <- 0.5
  eta_minus <- eta_minus-1
  eta_plus <- 1.2
  eta_plus <- eta_plus-1
  # maxcount <- maxcount
  
  # initialize 
  x <- xinit
  seed=1
  init_res<-fun(xinit, sample_n, seed)
  f<-init_res$neglpost
  g<-init_res$neglgrad
  s <- sqrt(t(g)%*%g)
  
  # loop starts here
  count <- 1
  del <- Delta0
  df <- 10
  
  while ((s>.0001)&(count<=maxcount)&(df>.0000001)) {
    # while (count<maxcount) {
    
    # step 1-move
    xnew<-x-sign(g)*del
    new_res<-fun(xnew, sample_n, seed)
    fnew<-new_res$neglpost
    gnew<-new_res$neglgrad
    s <- sqrt(t(gnew)%*%gnew)
    df <- abs(fnew/f-1)
    
    # step 2 - update step size
    gc <- g*gnew
    del_max<-matrix(apply(cbind(Deltamin,del*(1+eta_plus*(gc>0)+eta_minus*(gc<0))), MARGIN=c(1), max))
    del <- matrix(apply(cbind(Deltamax,del_max), MARGIN=c(1), min))
    
    x<-xnew
    g<-gnew
    f<-fnew
    print (paste0(count,", ",seed, ", ",x[nrow(x),]))
    count<-count+1
    if (count%%10==1) {seed<-seed+1}
  }
  
  res<-fun(x, sample_n,seed)
  output<-list(xopt=x, fopt=res$neglpost, gradopt=res$neglgrad)
  return (output)
}

PrepareInformedPriors<-function (distMat, focalsite=NULL, lags,neighbors, vars) {
  Amat<-distMat
  
  # Calculate distance from focal site to neighboring sites
  if(is.null(focalsite)) {
    focalsite<-1
  }
  Ind<-order(Amat[focalsite,])
  
  dist<-array(NA, dim = c(length(lags), length(neighbors), length(vars)))
  for (v in 1:length(vars)) {
    var<-vars[v]
    for (n in 1:length(neighbors)) {
      neighbor<-neighbors[n]
      site<-Ind[neighbor]
      dist[ ,n,v]<-Amat[focalsite,site]
    }
  }
  dim(dist)<-c(1,length(lags)*length(neighbors)*length(vars))
  
  hnormvars<-matrix(NA,nrow=length(dist) , ncol=1)
  for (j in 1:length(dist)) {
    hnormvars[j,]<-pi/2*(1/sqrt(4*v*(((j-1)%%length(lags))+1))*exp(-dist[j]^2/(2*v*(((j-1)%%length(lags))+1))) )^2
    # hnormvars[j,]<-((j-1)%%lags)+1 # for testing
  }
  
  return(hnormvars)
}


GPSDM<-function (pars, distMat,basisX, basisP,basisD, basisY=NULL, newX=NULL,newSigma=NULL,newP=NULL,newD=NULL,newY=NULL,priors=NULL, priors_tr=NULL,Cinv.fixed=NULL,m.fixed=NULL,xnew=NULL, vnew=NULL, P_forecast=NULL,focalsites=NULL ,start=NULL, steps=NULL, lags=NULL, neighbors=NULL,vars=NULL,C_g=NULL,mode=c("basics", "optimize", "predict", "forecast","update")) {
  X<-basisX
  Y<-basisY
  P<-basisP
  D<-basisD
  J<-date2doy(D)
  ndim<-ncol(X)
  
  # transform parameters from real line to constrained space
  phi<-(phimax-phimin)/(1+exp(-pars[1:ndim,,drop=F]))+phimin
  ve<-(vemax-vemin)/(1+exp(-pars[ndim+1]))+vemin
  tau<-(taumax-taumin)/(1+exp(-pars[ndim+2]))+taumin
  # rho<-(rhomax-rhomin)/(1+exp(-pars[ndim+3]))+rhomin
  gamma1<-(gamma1max-gamma1min)/(1+exp(-pars[ndim+3]))+gamma1min
  gamma2<-(gamma2max-gamma2min)/(1+exp(-pars[ndim+4]))+gamma2min
  
  if (!is.null(Cinv.fixed)&!is.null(m.fixed) & mode!="optimize" ) { #must recalculate when optimizing
    Cinv<-Cinv.fixed
    m<-m.fixed
  }
  else {
    # construct base covariance matrix
    # kXX
    lC0<-0
    DD1<-vector(mode="list")
    for (i in 1:ndim) {
      DD1[[i]]<-abs(X[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(X)) -matrix(1, nrow=nrow(X),ncol=1)%*%t(X[,i, drop=FALSE]))^2
      lC0<-lC0-0.5*phi[i]*DD1[[i]] #I added 0.5
    }
    kXX<-tau*exp(lC0)
    # ps1<-((P%*%matrix(1, nrow=1,ncol=nrow(P))-matrix(1, nrow=nrow(P),ncol=1)%*%t(P))==0)
    Dist1<-matrix(NA, nrow=nrow(X), ncol=nrow(X))
    for (i in 1:nrow(X)) {
      for (j in 1:nrow(X)) {
        Dist1[i,j]<-distMat[P[i],P[j]]
      }
    }
    kXX<-kXX*exp(-Dist1^2*gamma1)
    tDist1<-matrix(NA, nrow=nrow(X), ncol=nrow(X))
    for (i in 1:nrow(X)) {
      for (j in 1:nrow(X)) {
        tDist1[i,j]<-tdistMat[J[i],J[j]]#+abs(rnorm(1,0,0.01))
      }
    }
    kXX<-kXX*exp(-tDist1^2*gamma2)
    # kXX<-tau*exp(lC0)*(ps1+rho*(1-ps1))*exp(-tDist1^2*gamma2)
    
    Id<-diag(1, nrow=nrow(X), ncol=nrow(X))
    C<-kXX+ve*Id
    C<-(C+t(C))/2
    # C<-as.matrix(Matrix::nearPD(C)$mat)
    # chol algorithm from R&W
    L<-chol(C)
    Linv<-solve(L)#%*%Id
    Cinv<-tcrossprod(Linv)
    m<-Cinv%*%Y
  }
  
  if (mode=="basics") {
    output<-list(Cinv=Cinv,m=m)
    return(output)
  }
  
  if (!is.null(newX)&!is.null(newP)) {
    Xt<-newX
    Pt<-newP
    Dt<-newD
    Jt<-date2doy(Dt)
    
    # kXtX
    lC0<-0
    DD2<-vector(mode="list")
    for (i in 1:ndim) {
      DD2[[i]]<-abs(Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(X))-matrix(1, nrow=nrow(Xt),ncol=1)%*%t(X[,i, drop=FALSE]))^2
      lC0<-lC0-0.5*phi[i]*DD2[[i]]
    }
    kXtX<-tau*exp(lC0)
    # ps2<-((Pt%*%matrix(1, nrow=1,ncol=nrow(P))-matrix(1,nrow=nrow(Pt),ncol=1)%*%t(P))==0)
    Dist2<-matrix(NA, nrow=nrow(Xt), ncol=nrow(X))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(X)) {
        Dist2[i,j]<-distMat[Pt[i],P[j]]
      }
    }
    kXtX<-kXtX*exp(-Dist2^2*gamma1)
    tDist2<-matrix(NA, nrow=nrow(Xt), ncol=nrow(X))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(X)) {
        tDist2[i,j]<-tdistMat[Jt[i],J[j]]#+abs(rnorm(1,0,0.01))
      }
    }
    kXtX<-kXtX*exp(-tDist2^2*gamma2)
    # kXtX<-tau*exp(lC0)*(ps2+rho*(1-ps2))*exp(-tDist2^2*gamma2)
    kXXt<-t(kXtX)
    
    # kXtXt
    lC0<-0
    DD3<-vector(mode="list")
    for (i in 1:ndim) {
      DD3[[i]]<-abs(Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(Xt))-matrix(1, nrow=nrow(Xt),ncol=1)%*%t(Xt[,i, drop=FALSE]))^2
      lC0<-lC0-0.5*phi[i]*DD3[[i]]
    }
    kXtXt<-tau*exp(lC0)
    # ps3<-((Pt%*%matrix(1, nrow=1,ncol=nrow(Pt))-matrix(1,nrow=nrow(Pt),ncol=1)%*%t(Pt))==0)
    Dist3<-matrix(NA, nrow=nrow(Xt), ncol=nrow(Xt))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(Xt)) {
        Dist3[i,j]<-distMat[Pt[i],Pt[j]]
      }
    }
    kXtXt<-kXtXt*exp(-Dist3^2*gamma1)
    tDist3<-matrix(NA, nrow=nrow(Xt), ncol=nrow(Xt))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(Xt)) {
        tDist3[i,j]<-tdistMat[Jt[i],Jt[j]]#+abs(rnorm(1,0,0.01))
      }
    }
    kXtXt<-kXtXt*exp(-tDist3^2*gamma2)
    # kXtXt<-tau*exp(lC0)*(ps3+rho*(1-ps3))*exp(-tDist3^2*gamma2)
    
    Id<-diag(1, nrow=nrow(Xt), ncol=nrow(Xt))
    mt<-kXtX%*%m
    Ct<-kXtXt-kXtX%*%Cinv%*%kXXt+ve*Id
    Ct<-(Ct+t(Ct))/2
    # Ct<-as.matrix(Matrix::nearPD(Ct)$mat)
    
    # A<-kXtXt+ve*Id
    # LA<-chol(A)
    # LAinv<-solve(LA)
    # Ainv<-tcrossprod(LAinv)
    # 
    # U<-kXtX
    # V<-kXXt
    # S <- C-V%*%Ainv%*%U
    # LS<-chol(S)
    # LSinv<-solve(LS)
    # Sinv<-tcrossprod(LSinv)
    # 
    # Ctinv<-Ainv+Ainv%*%U%*%Sinv%*%V%*%Ainv
    # logdetCt<-log(det(S))+log(det(Cinv))+log(det(A))
    Lt<-chol(Ct)
    Ltinv<-solve(Lt)#%*%Id
    Ctinv<-tcrossprod(Ltinv)
  }
  
  if (mode=="predict") {
    # # invalid_dim<-which(phi<=1e-50)
    # Winv<-diag(as.numeric(phi), nrow=ndim, ncol=ndim)
    # Sigma_x<-diag(as.numeric(newSigma), nrow=ndim, ncol=ndim)
    # 
    # Deltainv<-diag(diag(Winv)-1/(1/(diag(Winv))+diag(Sigma_x)),nrow=ndim, ncol=ndim)
    # # Deltainv<-Winv-ginv(ginv(Winv)+Sigma_x)
    # # Deltainv[invalid_dim,invalid_dim]<-0
    # lC0<-0
    # DDD2<-vector(mode="list")
    # for (i in 1:ndim) {
    #   DDD2[[i]]<-abs(Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(X))-matrix(1, nrow=nrow(Xt),ncol=1)%*%t(X[,i, drop=FALSE]))^2
    #   lC0<-lC0+0.5*Deltainv[i,i]*DDD2[[i]]
    # }
    # coef<-det(Winv%*%Sigma_x+diag(1, nrow=ndim, ncol=ndim))^(-1/2)
    # Corr<-coef*exp(lC0)
    # 
    # Lambdainv<-diag(2*diag(Winv)-1/(1/diag(Winv)/2+diag(Sigma_x)),nrow=ndim, ncol=ndim)
    # # Lambdainv<-2*Winv-ginv(ginv(Winv)/2+Sigma_x)
    # # Lambdainv[invalid_dim,invalid_dim]<-0
    # X_mu<-matrix(1, nrow=nrow(X), ncol=1)%*%colMeans(X)
    # lC0<-0
    # DDD2<-vector(mode="list")
    # for (i in 1:ndim) {
    #   DDD2[[i]]<-abs(Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(X_mu))-matrix(1, nrow=nrow(Xt),ncol=1)%*%t(X_mu[,i, drop=FALSE]))^2
    #   lC0<-lC0+0.5*Lambdainv[i,i]*DDD2[[i]]
    # }
    # coef<-det(2*Winv%*%Sigma_x+diag(1, nrow=ndim, ncol=ndim))^(-1/2)
    # Corr2<-coef*exp(lC0)
    # 
    # mt_ep<-(kXtX*Corr)%*%m #3.39
    # Ct_ep<-kXtXt-(kXtX*Corr2)%*%(Cinv-m%*%t(m))%*%kXXt-(mt_ep)%*%t(mt_ep)+ve*diag(1, nrow=nrow(Xt), ncol=nrow(Xt))
    # output<-list(mt=matrix(mt_ep),Ct=matrix(diag(Ct_ep)))
    output<-list(mt=matrix(mt),Ct=matrix(diag(Ct)))
    return(output)
  }
  
  if (mode=="optimize") {
    # if (is.null(newY)) {
    #   #####likelihood#####
    #   like<--.5*t(Y)%*%m-sum(log(diag(L)))
    #   
    #   #####gradient#####
    #   dl<-matrix(0, nrow= ndim+3, ncol=1)
    #   vQ<-matrix(tcrossprod(m)-Cinv, nrow = 1)
    #   dC<-vector(mode="list")
    #   for (i in 1:ndim) {
    #     dC[[i]]<--0.5*DD1[[i]]*kXX # I added this 0.5
    #     dl[i]<-.5*vQ%*%matrix(dC[[i]], ncol = 1)
    #   }
    #   dC[[ndim+1]]<-diag(1, nrow=nrow(X), ncol=nrow(X))
    #   dl[ndim+1]<-.5*vQ%*%matrix(dC[[ndim+1]], ncol=1)
    #   dC[[ndim+2]]<-kXX/tau
    #   dl[ndim+2]<-.5*vQ%*%matrix(dC[[ndim+2]], ncol=1)
    #   dC[[ndim+3]]<-kXX*(-Dist1^2)
    #   dl[ndim+3]<-.5*vQ%*%matrix(dC[[ndim+3]], ncol=1)
    # }
    
    if (!is.null(newY)) {
      Yt<-newY
      #####likelihood#####
      like<--.5*t(Yt-mt)%*%Ctinv%*%(Yt-mt)-sum(log(diag(Lt)))#0.5*logdetCt
      
      #####gradient#####
      dl<-matrix(0, nrow= ndim+4, ncol=1)
      vQ<-matrix(tcrossprod(Ctinv%*%(Yt-mt))-Ctinv, nrow = 1)
      dC<-vector(mode="list")
      for (i in 1:ndim) {
        dC[[i]]<--0.5*DD3[[i]]*kXtXt+0.5*(DD2[[i]]*kXtX)%*%Cinv%*%kXXt+0.5*kXtX%*%Cinv%*%(t(DD2[[i]])*kXXt)-0.5*kXtX%*%Cinv%*%(DD1[[i]]*kXX)%*%Cinv%*%kXXt
        dl[i,1]<-.5*vQ%*%matrix(dC[[i]], ncol = 1)
      }
      dC[[ndim+1]]<-diag(1, nrow=nrow(Xt), ncol=nrow(Xt))+kXtX%*%Cinv%*%Cinv%*%kXXt
      dl[ndim+1,1]<-.5*vQ%*%matrix(dC[[ndim+1]], ncol=1)
      dC[[ndim+2]]<-kXtXt/tau-(kXtX/tau)%*%Cinv%*%kXXt-kXtX%*%Cinv%*%(kXXt/tau)+kXtX%*%Cinv%*%(kXX/tau)%*%Cinv%*%kXXt
      dl[ndim+2,1]<-.5*vQ%*%matrix(dC[[ndim+2]], ncol=1)
      # dC[[ndim+3]]<-kXtXt*(1-ps3)/rho-(kXtX*(1-ps2)/rho)%*%Cinv%*%kXXt-kXtX%*%Cinv%*%(kXXt*(1-t(ps2))/rho)+kXtX%*%Cinv%*%(kXX*(1-ps1)/rho)%*%Cinv%*%kXXt
      # dl[ndim+3,1]<-.5*vQ%*%matrix(dC[[ndim+3]], ncol=1)
      dC[[ndim+3]]<-kXtXt*(-Dist3^2)-(kXtX*(-Dist2^2))%*%Cinv%*%kXXt-kXtX%*%Cinv%*%(kXXt*(-t(Dist2)^2))+kXtX%*%Cinv%*%(kXX*(-Dist1^2))%*%Cinv%*%kXXt
      dl[ndim+3, 1]<-.5*vQ%*%matrix(dC[[ndim+3]], ncol=1)
      dC[[ndim+4]]<-kXtXt*(-tDist3^2)-(kXtX*(-tDist2^2))%*%Cinv%*%kXXt-kXtX%*%Cinv%*%(kXXt*(-t(tDist2)^2))+kXtX%*%Cinv%*%(kXX*(-tDist1^2))%*%Cinv%*%kXXt
      dl[ndim+4, 1]<-.5*vQ%*%matrix(dC[[ndim+4]], ncol=1)
    }
    
    # derivative for hyperparameters wrt transformed hyperparameters -- for gradient calculation
    dpars<-matrix(c((phi-phimin)*(1-(phi-phimin)/(phimax-phimin)),
                    (ve-vemin)*(1-(ve-vemin)/(vemax-vemin)),
                    (tau-taumin)*(1-(tau-taumin)/(taumax-taumin)),
                    # (rho-rhomin)*(1-(rho-rhomin)/(rhomax-rhomin)),
                    (gamma1-gamma1min)*(1-(gamma1-gamma1min)/(gamma1max-gamma1min)),
                    (gamma2-gamma2min)*(1-(gamma2-gamma2min)/(gamma2max-gamma2min))
    ))
    
    if (is.null(priors)&is.null(priors_tr)) {
      lpost<-like
      neglpost<--lpost
      
      # J is gradient in parameter space - need gradient in transformed parameters
      GradLpost<-(dl)*dpars
      neglgrad<--GradLpost
    }
    else if (!is.null(priors)&is.null(priors_tr)) {
      # phi
      lp_phi<--.5*sum(((phi-priors$E_phi)^2)/priors$V_phi)
      dlp_phi<--(phi-priors$E_phi)^1/priors$V_phi
      if (!is.null(priors$E_ve)) {
        lp_tau<--.5*sum(((tau-priors$E_tau)^2)/priors$V_tau)
        dlp_tau<--(tau-priors$E_tau)^1/priors$V_tau
        lp_ve<--.5*sum(((ve-priors$E_ve)^2)/priors$V_ve)
        dlp_ve<--(ve-priors$E_ve)^1/priors$V_ve
        # lp_rho<--.5*sum(((rho-priors$E_rho)^2)/priors$V_rho)
        # dlp_rho<--(rho-priors$E_rho)^1/priors$V_rho
        lp_gamma1<--.5*sum(((gamma1-priors$E_gamma1)^2)/priors$V_gamma1)
        dlp_gamma1<--(gamma1-priors$E_gamma1)^1/priors$V_gamma1
        lp_gamma2<--.5*sum(((gamma2-priors$E_gamma2)^2)/priors$V_gamma2)
        dlp_gamma2<--(gamma2-priors$E_gamma2)^1/priors$V_gamma2
      }
      
      lp<-(lp_phi+lp_ve+lp_tau+
             # lp_rho+
             lp_gamma1+
             lp_gamma2
      )
      dlp<-matrix(c(dlp_phi,dlp_ve,dlp_tau,
                    # dlp_rho,
                    dlp_gamma1,
                    dlp_gamma2
      ))
      
      lpost<-like+lp
      neglpost<--lpost
      
      # J is gradient in parameter space - need gradient in transformed parameters
      GradLpost<-(dl+dlp)*dpars
      neglgrad<--GradLpost
    }
    # else if (is.null(priors)&!is.null(priors_tr)) {
    #   
    #   E<-priors_tr$E
    #   V<-priors_tr$V
    #   
    #   lp_theta<- -1/2*t(pars-E)%*%solve(V)%*%(pars-E)
    #   dlp_theta<- -solve(V)%*%(pars-E)
    #   
    #   # derivative for transformed hyperparameters wrt hyperparameters -- for jacobian in likelihood
    #   lj_phi<-0
    #   for (i in 1:ndim) {
    #     # if (phi[i,]<10^-10) {
    #     #   lj_phi<-lj_phi+log(1/(10^-10))
    #     # } else {
    #     #   lj_phi<-lj_phi+log(1/phi[i,])
    #     # }
    #     if ((phimax-phi[i,])<10^-10) {
    #       lj_phi<-lj_phi+log((phimax-phimin)/(10^-10)/(phi[i,]-phimin))
    #     } else if ((phi[i,]-phimin)<10^-10) {
    #       lj_phi<-lj_phi+log((phimax-phimin)/(phimax-phi[i,])/(10^-10))
    #     } else {
    #       lj_phi<-lj_phi+log((phimax-phimin)/(phimax-phi[i,])/(phi[i,]-phimin))
    #     }
    #   }
    #   # lj_phi<-log(sum(1/phi))
    #   if ((vemax-ve)<10^-10) {
    #     lj_ve<-log((vemax-vemin)/(10^-10)/(ve-vemin))
    #   } else if ((ve-vemin)<10^-10) {
    #     lj_ve<-log((vemax-vemin)/(vemax-ve)/(10^-10))
    #   } else {
    #     lj_ve<-log((vemax-vemin)/(vemax-ve)/(ve-vemin))
    #   }
    #   # lj_ve<-log((vemax-vemin)/(vemax-ve)/(ve-vemin))
    #   if ((taumax-tau)<10^-10) {
    #     lj_tau<-log((taumax-taumin)/(10^-10)/(tau-taumin))
    #   } else if ((tau-taumin)<10^-10) {
    #     lj_tau<-log((taumax-taumin)/(taumax-tau)/(10^-10))
    #   } else {
    #     lj_tau<-log((taumax-taumin)/(taumax-tau)/(tau-taumin))
    #   }
    #   # lj_tau<-log((taumax-taumin)/(taumax-tau)/(tau-vemin))
    #   if ((gammamax-gamma)<10^-10) {
    #     lj_gamma<-log((gammamax-gammamin)/(10^-10)/(gamma-gammamin))
    #   } else if ((gamma-gammamin)<10^-10) {
    #     lj_gamma<-log((gammamax-gammamin)/(gammamax-gamma)/(10^-10))
    #   } else {
    #     lj_gamma<-log((gammamax-gammamin)/(gammamax-gamma)/(gamma-gammamin))
    #   }
    #   # lj_gamma<-log((gammamax-gammamin)/(gammamax-gamma)/(gamma-gammamin))
    #   
    #   # derivative for jacobian wrt transformed hyperparameters -- for jacobian in gradients
    #   phi_tr<-pars[1:ndim,,drop=F]
    #   ve_tr<-pars[ndim+1,]
    #   tau_tr<-pars[ndim+2,]
    #   gamma_tr<-pars[ndim+3,]
    #   
    #   # dlj_phi<-matrix(-1, nrow=ndim, ncol=1)
    #   dlj_phi<-(exp(phi_tr)-1)/(exp(phi_tr)+1)
    #   dlj_ve<-(exp(ve_tr)-1)/(exp(ve_tr)+1)
    #   dlj_tau<-(exp(tau_tr)-1)/(exp(tau_tr)+1)
    #   dlj_gamma<-(exp(gamma_tr)-1)/(exp(gamma_tr)+1)
    #   
    #   lp<-lp_theta+(lj_phi+lj_ve+lj_tau+lj_gamma)
    #   dlp<-dlp_theta+matrix(c(dlj_phi,dlj_ve,dlj_tau,dlj_gamma), ncol=1)
    #   
    #   lpost<-like+lp
    #   neglpost<--lpost
    #   
    #   # J is gradient in parameter space - need gradient in transformed parameters
    #   GradLpost<-dl*dpars+dlp
    #   neglgrad<--GradLpost
    # }
    
    # output<-list(neglpost=neglpost, neglgrad=neglgrad, LOO=lnL_LOO, df=df) #LOO and df not used in Rprop?
    output<-list(neglpost=neglpost, neglgrad=neglgrad)
    return(output)
  }
  
  # if (mode =="forecast") {
  #   if(is.null(steps)) {
  #     steps<-sum(is.na(xnew[1,,1]))-1
  #   }
  # 
  #   Y_pred<-Var_pred<-matrix(NA, nrow=length(focalsites), ncol = steps )
  #   for (t in 1:steps) {
  #     print(t)
  #     # prepare for next iteration
  #     # adapted this for spatial-delay embedding
  #     res<-PrepareEmbedding(xnew,start=start+t,end=start+t, focalsites = focalsites, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  #     mu_Xt<-res$X
  #     Pt<-P_forecast
  # 
  #     res<-PrepareEmbedding(vnew,start=start+t,end=start+t, focalsites = focalsites, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  #     Sigma_Xt<-res$X
  # 
  #     # kXtX
  #     lC0<-0
  #     DD2<-vector(mode="list")
  #     for (i in 1:ndim) {
  #       DD2[[i]]<-abs(mu_Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(X))-matrix(1, nrow=nrow(mu_Xt),ncol=1)%*%t(X[,i, drop=FALSE]))^2
  #       # varinv<-1/(Sigma_Xt[site,i]+1/phi[i]) #not the same as the equation proposed but can't find ways to inverse the full matrix
  #       varinv<-phi[i] # for testing
  #       lC0<-lC0-0.5*varinv*DD2[[i]]
  #     }
  #     ps2<-((Pt%*%matrix(1, nrow=1,ncol=nrow(P))-matrix(1,nrow=nrow(Pt),ncol=1)%*%t(P))==0)
  #     kXtX<-tau*exp(lC0)*(ps2+rho*(1-ps2))
  #     # Dist2<-matrix(NA, nrow=nrow(mu_Xt), ncol=nrow(X))
  #     # for (i in 1:nrow(mu_Xt)) {
  #     #   for (j in 1:nrow(X)) {
  #     #     Dist2[i,j]<-distMat[Pt[i],P[j]]
  #     #   }
  #     # }
  #     # # Id<-diag(1, nrow=ndim, ncol=ndim)
  #     # # coef<-(diag(c(phi), nrow=ndim, ncol=ndim))%*%diag(c(Sigma_Xt[site,,drop=FALSE]), nrow=ndim, ncol=ndim)+Id #(3.40)
  #     # # kXtX<-(det(coef))^(-1/2)*tau*exp(lC0)*(ps2+rho*(1-ps2))
  #     # kXtX<-tau*exp(lC0)*exp(-Dist2^2*gamma)
  #     kXXt<-t(kXtX)
  # 
  #     # kXtXt
  #     lC0<-0
  #     DD3<-vector(mode="list")
  #     for (i in 1:ndim) {
  #       DD3[[i]]<-abs(mu_Xt[,i, drop=FALSE]%*%matrix(1, nrow=1,ncol=nrow(mu_Xt))-matrix(1, nrow=nrow(mu_Xt),ncol=1)%*%t(mu_Xt[,i, drop=FALSE]))^2
  #       lC0<-lC0-0.5*phi[i]*DD3[[i]]
  #     }
  #     ps3<-((Pt%*%matrix(1, nrow=1,ncol=nrow(Pt))-matrix(1,nrow=nrow(Pt),ncol=1)%*%t(Pt))==0)
  #     kXtXt<-tau*exp(lC0)*(ps3+rho*(1-ps3))
  #     # Dist3<-matrix(NA, nrow=nrow(mu_Xt), ncol=nrow(mu_Xt))
  #     # for (i in 1:nrow(mu_Xt)) {
  #     #   for (j in 1:nrow(mu_Xt)) {
  #     #     Dist3[i,j]<-distMat[Pt[i],Pt[j]]
  #     #   }
  #     # }
  #     # kXtXt<-tau*exp(lC0)*exp(-Dist3^2*gamma)
  # 
  #     mt<-kXtX%*%m
  #     Id<-matrix(1, nrow=nrow(mu_Xt), ncol=nrow(mu_Xt))
  #     Ct<-kXtXt-kXtX%*%Cinv%*%t(kXtX)+ve*Id #original
  # 
  #     #store prediction
  #     Y_pred[,t]<-as.numeric(mt)
  #     Var_pred[,t]<-as.numeric(diag(Ct))
  # 
  #     xnew[,(start+t),1]<-as.numeric(mt)
  #     vnew[,(start+t),1]<-as.numeric(diag(Ct))
  #   }
  # 
  #   # Y_pred<-matrix(t(Y_pred), ncol=1)
  #   # Var_pred<-matrix(t(Var_pred), ncol=1)
  #   output<-list(Y_pred=Y_pred, Var_pred=Var_pred)
  #   return(output)
  # }
  # 
  # if (mode=="update") {
  #   Yt<-newY
  #   if(is.null(C_g)) {
  #     C_g<-diag(1,nrow=nrow(X), ncol=nrow(X))
  #   }
  #   J_t<-kXtX%*%Cinv #(8)
  #   B<-kXtXt-J_t%*%kXXt #(7)
  #   mu_t_p<-mt #(6)
  #   C_t_p<-B+J_t%*%C_g%*%t(J_t) #below (9)
  #   
  #   Id<-diag(1, nrow=nrow(Xt), ncol=nrow(Xt))
  #   G_t<-C_g%*%t(J_t)%*%solve(C_t_p+ve*Id) #(12)
  #   mu_t_g<-Y+G_t%*%(Yt-mu_t_p)
  #   C_t_g<-C_g-G_t%*%J_t%*%C_g
  #   
  #   output<-list(basisY=mu_t_g,basisY_var=C_t_g)
  #   return (output)
  # }
}
