library(tidyverse)
library(LaplacesDemon)
library(doSNOW)
library(lubridate)
library(foreach)
library(gridExtra)
library(xts)
library(ggpubr)
library(ggrepel)
library(dygraphs)
library(rlist)
library(ptw)
library(RhpcBLASctl)
library(ptw)

waves <- function(t, t_start,
                  intercept, slope,
                  amplitude1, phase1, period1,
                  amplitude2, phase2, period2,
                  sd = 0.05) {
  t_diff <- as.numeric(t - t_start)
  if (leap_year(t)) {
    d_all <- 366
  } else {
    d_all <- 365
  }
  v <- intercept + slope * t_diff +
    amplitude1 * sin(2 * pi / (period1 * d_all) * (t_diff + phase1)) +
    amplitude2 * sin(2 * pi / (period2 * d_all) * (t_diff + phase2)) +
    rnorm(1, 0, sd)

  return(v)
}

env_to_param <- function(env, lower, upper, steepness, midpoint) {
  param <- (upper - lower) / (1 + exp(-steepness * (env - midpoint))) + lower
  return(param)
}

double_logistics <- function(t,
                             m1 = 0, # average greenness in winter
                             m2 = 1, # difference between summer and winter
                             m3 = 100, # spring onset
                             m4 = 10, # slope of curve in spring
                             m5 = 260, # fall offset
                             m6 = 20, # slope of curve in fall
                             m7 = 0, # summer greendown
                             m8 = 1, # life cycle
                             sd = 0.02) {
  if (!leap_year(t)) {
    d <- (as.integer(format(t, "%j")) %% (365 / (m8))) * (m8)
  } else {
    d <- (as.integer(format(t, "%j")) %% (366 / (m8))) * (m8)
  }

  v <- m1 + (m2 - m7 * d) * (1 / (1 + exp((m3 - d) / m4)) - 1 / (1 + exp((m5 - d) / m6))) + rnorm(1, 0, sd)

  return(v)
}

compare_stats <- function(obs_ori, pred_ori, obs = NULL, pred = NULL, range = NULL) {
  corr <- cor(obs_ori, pred_ori, use = "pairwise.complete.obs")
  R2 <- summary(lm(pred_ori ~ obs_ori))$r.squared
  RMSE <- sqrt(mean((obs_ori - pred_ori)^2, na.rm = T))

  if (!is.null(obs)) {
    nRMSE <- sqrt(mean((obs - pred)^2, na.rm = T))
  } else if (!is.null(range)) {
    nRMSE <- RMSE / range
  } else {
    nRMSE <- NA
  }

  out <- list(corr = corr, R2 = R2, nRMSE = nRMSE, RMSE = RMSE)
  return(out)
}

# reference: https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace) {
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if (tolower(triangle.to.replace) == "lower") {
    tri <- lower.tri(m)
  } else if (tolower(triangle.to.replace) == "upper") {
    tri <- upper.tri(m)
  } else {
    stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  }
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints) {
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.

  GeoDistanceInMetres <- function(g1, g2) {
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2) {
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1 = g1$lat, lon.1 = g1$lon, lat.2 = g2$lat, lon.2 = g2$lon, units = "m")))
    }
    return(mapply(DistM, g1, g2))
  }

  n.geopoints <- nrow(df.geopoints)

  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints

  # Create a list of lists
  list.geopoints <- by(df.geopoints[, c("index", "lat", "lon")], 1:n.geopoints, function(x) {
    return(list(x))
  })

  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name

  return(mat.distances)
}

tdistMat <- matrix(NA, nrow = 365, ncol = 365)
for (i in 1:365) {
  for (j in 1:365) {
    tdistMat[i, j] <- min(abs(i - j), 365 - abs(i - j))
  }
}
colnames(tdistMat) <- 1:365
tdistMat_df <- as.data.frame(tdistMat) %>%
  mutate(doy1 = 1:365) %>%
  gather(key = "doy2", value = "value", -doy1) %>%
  mutate(
    doy1 = as.integer(doy1),
    doy2 = as.integer(doy2)
  )

date2doy <- function(x) {
  doy_all <- matrix(NA, nrow = nrow(x), ncol = 1)
  for (i in 1:nrow(x)) {
    doy <- as.integer(format(as.Date(x[i, ]), "%j"))
    if (doy == 366) {
      doy <- 365
    }
    doy_all[i, ] <- doy
  }
  return(doy_all)
}

PrepareEmbedding <- function(x, start, end, focalsites = NULL, lags, neighbors, vars, distMat) {
  numPatch <- dim(x)[[1]]
  Amat <- distMat

  if (is.null(focalsites)) {
    focalsites <- seq(1, numPatch)
  }

  Xi_list <- Yi_list <- Pi_list <- Di_list <- vector(mode = "list", length = length(focalsites))
  for (i in focalsites) {
    Ind <- order(Amat[i, ])
    Ind <- c(i, setdiff(Ind, i)) # make sure focal site is always the first

    Xi <- matrix(NA, nrow = (end - start + 1), ncol = 0)
    for (v in 1:length(vars)) {
      var <- var_list[vars[v]]
      for (n in 1:length(neighbors[[v]])) {
        neighbor <- neighbors[[v]][n]
        site <- Ind[neighbor]
        if (!is.null(lags[[v]])) {
          for (period in 1:length(lags[[v]])) {
            period_lags <- lags[[v]][[period]]
            values <- matrix(NA, nrow = end - start + 1, ncol = length(period_lags))
            for (l in 1:length(period_lags)) {
              period_lag <- period_lags[l]
              values[, l] <- x[site, (start - period_lag):(end - period_lag), var]
            }
            Xi_add <- as.matrix(rowMeans(values, na.rm = T))
            colnames(Xi_add) <- paste0(var, "_", neighbor, "_", period)
            Xi <- cbind(Xi, Xi_add)
          }
        }
      }
    }
    Xi_list[[i]] <- Xi

    Yi <- matrix(x[i, start:end, 1, drop = F]) # Assuming the first var is always same as the var to predict
    Yi_list[[i]] <- Yi

    Pi <- matrix(i, nrow = (end - start + 1), ncol = 1)
    Pi_list[[i]] <- Pi

    Di <- matrix(dimnames(x[i, start:end, 1, drop = F])[[2]])
    Di_list[[i]] <- Di
  }

  X <- do.call(rbind, Xi_list)
  Y <- do.call(rbind, Yi_list)
  P <- do.call(rbind, Pi_list)
  D <- do.call(rbind, Di_list)

  output <- list(X = X, Y = Y, P = P, D = D)
  return(output)
}

fmingrad_Rprop <- function(fun, xinit, sample_n, maxcount) {
  # fun is a handle to a function that has grad as optional second output
  # this uses the sign of the gradient to determine descent directions and an adaptive step size - supposedly smoother convergence than conjugate gradients for GP optimization
  p <- length(xinit)

  # optimization parameters for Rprop
  Delta0 <- 0.1 * matrix(1, nrow = p, ncol = 1)
  Deltamin <- 1e-6 * matrix(1, nrow = p, ncol = 1)
  Deltamax <- 50 * matrix(1, nrow = p, ncol = 1)
  eta_minus <- 0.5
  eta_minus <- eta_minus - 1
  eta_plus <- 1.2
  eta_plus <- eta_plus - 1

  # initialize
  x <- xinit
  seed <- 1
  init_res <- fun(xinit, sample_n, seed)
  f <- init_res$neglpost
  g <- init_res$neglgrad
  s <- sqrt(t(g) %*% g)

  # loop starts here
  count <- 1
  del <- Delta0
  df <- 10

  while ((s > .0001) & (count <= maxcount) & (df > .0000001)) {
    # step 1-move
    xnew <- x - sign(g) * del
    new_res <- fun(xnew, sample_n, seed)
    fnew <- new_res$neglpost
    gnew <- new_res$neglgrad
    s <- sqrt(t(gnew) %*% gnew)
    df <- abs(fnew / f - 1)

    # step 2 - update step size
    gc <- g * gnew
    del_max <- matrix(apply(cbind(Deltamin, del * (1 + eta_plus * (gc > 0) + eta_minus * (gc < 0))), MARGIN = c(1), max))
    del <- matrix(apply(cbind(Deltamax, del_max), MARGIN = c(1), min))

    x <- xnew
    g <- gnew
    f <- fnew
    print(paste0("seed: ", seed, ", iteration: ", count))
    count <- count + 1
    if (count %% 10 == 1) {
      seed <- seed + 1
    }
  }

  res <- fun(x, sample_n, seed)
  output <- list(xopt = x, fopt = res$neglpost, gradopt = res$neglgrad)
  return(output)
}

PrepareInformedPriors <- function(distMat, focalsite = NULL, lags, neighbors, vars) {
  Amat <- distMat

  # Calculate distance from focal site to neighboring sites
  if (is.null(focalsite)) {
    focalsite <- 1
  }
  Ind <- order(Amat[focalsite, ])

  dist <- array(NA, dim = c(length(lags), length(neighbors), length(vars)))
  for (v in 1:length(vars)) {
    var <- vars[v]
    for (n in 1:length(neighbors)) {
      neighbor <- neighbors[n]
      site <- Ind[neighbor]
      dist[, n, v] <- Amat[focalsite, site]
    }
  }
  dim(dist) <- c(1, length(lags) * length(neighbors) * length(vars))

  hnormvars <- matrix(NA, nrow = length(dist), ncol = 1)
  for (j in 1:length(dist)) {
    hnormvars[j, ] <- pi / 2 * (1 / sqrt(4 * v * (((j - 1) %% length(lags)) + 1)) * exp(-dist[j]^2 / (2 * v * (((j - 1) %% length(lags)) + 1))))^2
  }

  return(hnormvars)
}

GPSDM <- function(pars, distMat, basisX, basisP, basisD, basisY = NULL, newX = NULL, newSigma = NULL, newP = NULL, newD = NULL, newY = NULL, priors = NULL, priors_tr = NULL, Cinv.fixed = NULL, m.fixed = NULL, xnew = NULL, vnew = NULL, P_forecast = NULL, focalsites = NULL, start = NULL, steps = NULL, lags = NULL, neighbors = NULL, vars = NULL, C_g = NULL, mode = c("basics", "optimize", "predict", "forecast", "update")) {
  X <- basisX
  Y <- basisY
  P <- basisP
  D <- basisD
  J <- date2doy(D)
  ndim <- ncol(X)

  # transform parameters from real line to constrained space
  phi <- (phimax - phimin) / (1 + exp(-pars[1:ndim, , drop = F])) + phimin
  ve <- (vemax - vemin) / (1 + exp(-pars[ndim + 1])) + vemin
  tau <- (taumax - taumin) / (1 + exp(-pars[ndim + 2])) + taumin
  gamma1 <- (gamma1max - gamma1min) / (1 + exp(-pars[ndim + 3])) + gamma1min
  gamma2 <- (gamma2max - gamma2min) / (1 + exp(-pars[ndim + 4])) + gamma2min

  if (!is.null(Cinv.fixed) & !is.null(m.fixed) & mode != "optimize") { # must recalculate when optimizing
    Cinv <- Cinv.fixed
    m <- m.fixed
  } else {
    # construct base covariance matrix
    # kXX
    lC0 <- 0
    DD1 <- vector(mode = "list")
    for (i in 1:ndim) {
      DD1[[i]] <- abs(X[, i, drop = FALSE] %*% matrix(1, nrow = 1, ncol = nrow(X)) - matrix(1, nrow = nrow(X), ncol = 1) %*% t(X[, i, drop = FALSE]))^2
      lC0 <- lC0 - 0.5 * phi[i] * DD1[[i]]
    }
    kXX <- tau * exp(lC0)
    Dist1 <- matrix(NA, nrow = nrow(X), ncol = nrow(X))
    for (i in 1:nrow(X)) {
      for (j in 1:nrow(X)) {
        Dist1[i, j] <- distMat[P[i], P[j]]
      }
    }
    kXX <- kXX * exp(-Dist1^2 * gamma1)
    tDist1 <- matrix(NA, nrow = nrow(X), ncol = nrow(X))
    for (i in 1:nrow(X)) {
      for (j in 1:nrow(X)) {
        tDist1[i, j] <- tdistMat[J[i], J[j]]
      }
    }
    kXX <- kXX * exp(-tDist1^2 * gamma2)

    Id <- diag(1, nrow = nrow(X), ncol = nrow(X))
    C <- kXX + ve * Id
    C <- (C + t(C)) / 2
    L <- chol(C)
    Linv <- solve(L)
    Cinv <- tcrossprod(Linv)
    m <- Cinv %*% Y
  }

  if (mode == "basics") {
    output <- list(Cinv = Cinv, m = m)
    return(output)
  }

  if (!is.null(newX) & !is.null(newP)) {
    Xt <- newX
    Pt <- newP
    Dt <- newD
    Jt <- date2doy(Dt)

    # kXtX
    lC0 <- 0
    DD2 <- vector(mode = "list")
    for (i in 1:ndim) {
      DD2[[i]] <- abs(Xt[, i, drop = FALSE] %*% matrix(1, nrow = 1, ncol = nrow(X)) - matrix(1, nrow = nrow(Xt), ncol = 1) %*% t(X[, i, drop = FALSE]))^2
      lC0 <- lC0 - 0.5 * phi[i] * DD2[[i]]
    }
    kXtX <- tau * exp(lC0)
    Dist2 <- matrix(NA, nrow = nrow(Xt), ncol = nrow(X))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(X)) {
        Dist2[i, j] <- distMat[Pt[i], P[j]]
      }
    }
    kXtX <- kXtX * exp(-Dist2^2 * gamma1)
    tDist2 <- matrix(NA, nrow = nrow(Xt), ncol = nrow(X))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(X)) {
        tDist2[i, j] <- tdistMat[Jt[i], J[j]] #+abs(rnorm(1,0,0.01))
      }
    }
    kXtX <- kXtX * exp(-tDist2^2 * gamma2)
    kXXt <- t(kXtX)

    # kXtXt
    lC0 <- 0
    DD3 <- vector(mode = "list")
    for (i in 1:ndim) {
      DD3[[i]] <- abs(Xt[, i, drop = FALSE] %*% matrix(1, nrow = 1, ncol = nrow(Xt)) - matrix(1, nrow = nrow(Xt), ncol = 1) %*% t(Xt[, i, drop = FALSE]))^2
      lC0 <- lC0 - 0.5 * phi[i] * DD3[[i]]
    }
    kXtXt <- tau * exp(lC0)
    Dist3 <- matrix(NA, nrow = nrow(Xt), ncol = nrow(Xt))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(Xt)) {
        Dist3[i, j] <- distMat[Pt[i], Pt[j]]
      }
    }
    kXtXt <- kXtXt * exp(-Dist3^2 * gamma1)
    tDist3 <- matrix(NA, nrow = nrow(Xt), ncol = nrow(Xt))
    for (i in 1:nrow(Xt)) {
      for (j in 1:nrow(Xt)) {
        tDist3[i, j] <- tdistMat[Jt[i], Jt[j]] #+abs(rnorm(1,0,0.01))
      }
    }
    kXtXt <- kXtXt * exp(-tDist3^2 * gamma2)

    Id <- diag(1, nrow = nrow(Xt), ncol = nrow(Xt))
    mt <- kXtX %*% m
    Ct <- kXtXt - kXtX %*% Cinv %*% kXXt + ve * Id
    Ct <- (Ct + t(Ct)) / 2

    Lt <- chol(Ct)
    Ltinv <- solve(Lt)
    Ctinv <- tcrossprod(Ltinv)
  }

  if (mode == "predict") {
    output <- list(mt = matrix(mt), Ct = matrix(diag(Ct)))
    return(output)
  }

  if (mode == "optimize") {
    if (!is.null(newY)) {
      Yt <- newY
      # likelihood
      like <- -.5 * t(Yt - mt) %*% Ctinv %*% (Yt - mt) - sum(log(diag(Lt)))

      # gradient
      dl <- matrix(0, nrow = ndim + 4, ncol = 1)
      vQ <- matrix(tcrossprod(Ctinv %*% (Yt - mt)) - Ctinv, nrow = 1)
      dC <- vector(mode = "list")
      for (i in 1:ndim) {
        dC[[i]] <- -0.5 * DD3[[i]] * kXtXt + 0.5 * (DD2[[i]] * kXtX) %*% Cinv %*% kXXt + 0.5 * kXtX %*% Cinv %*% (t(DD2[[i]]) * kXXt) - 0.5 * kXtX %*% Cinv %*% (DD1[[i]] * kXX) %*% Cinv %*% kXXt
        dl[i, 1] <- .5 * vQ %*% matrix(dC[[i]], ncol = 1)
      }
      dC[[ndim + 1]] <- diag(1, nrow = nrow(Xt), ncol = nrow(Xt)) + kXtX %*% Cinv %*% Cinv %*% kXXt
      dl[ndim + 1, 1] <- .5 * vQ %*% matrix(dC[[ndim + 1]], ncol = 1)
      dC[[ndim + 2]] <- kXtXt / tau - (kXtX / tau) %*% Cinv %*% kXXt - kXtX %*% Cinv %*% (kXXt / tau) + kXtX %*% Cinv %*% (kXX / tau) %*% Cinv %*% kXXt
      dl[ndim + 2, 1] <- .5 * vQ %*% matrix(dC[[ndim + 2]], ncol = 1)
      dC[[ndim + 3]] <- kXtXt * (-Dist3^2) - (kXtX * (-Dist2^2)) %*% Cinv %*% kXXt - kXtX %*% Cinv %*% (kXXt * (-t(Dist2)^2)) + kXtX %*% Cinv %*% (kXX * (-Dist1^2)) %*% Cinv %*% kXXt
      dl[ndim + 3, 1] <- .5 * vQ %*% matrix(dC[[ndim + 3]], ncol = 1)
      dC[[ndim + 4]] <- kXtXt * (-tDist3^2) - (kXtX * (-tDist2^2)) %*% Cinv %*% kXXt - kXtX %*% Cinv %*% (kXXt * (-t(tDist2)^2)) + kXtX %*% Cinv %*% (kXX * (-tDist1^2)) %*% Cinv %*% kXXt
      dl[ndim + 4, 1] <- .5 * vQ %*% matrix(dC[[ndim + 4]], ncol = 1)
    }

    # derivative for hyperparameters wrt transformed hyperparameters -- for gradient calculation
    dpars <- matrix(c(
      (phi - phimin) * (1 - (phi - phimin) / (phimax - phimin)),
      (ve - vemin) * (1 - (ve - vemin) / (vemax - vemin)),
      (tau - taumin) * (1 - (tau - taumin) / (taumax - taumin)),
      (gamma1 - gamma1min) * (1 - (gamma1 - gamma1min) / (gamma1max - gamma1min)),
      (gamma2 - gamma2min) * (1 - (gamma2 - gamma2min) / (gamma2max - gamma2min))
    ))

    if (is.null(priors) & is.null(priors_tr)) {
      lpost <- like
      neglpost <- -lpost

      # J is gradient in parameter space - need gradient in transformed parameters
      GradLpost <- (dl) * dpars
      neglgrad <- -GradLpost
    } else if (!is.null(priors) & is.null(priors_tr)) {
      # phi
      lp_phi <- -.5 * sum(((phi - priors$E_phi)^2) / priors$V_phi)
      dlp_phi <- -(phi - priors$E_phi)^1 / priors$V_phi
      if (!is.null(priors$E_ve)) {
        lp_tau <- -.5 * sum(((tau - priors$E_tau)^2) / priors$V_tau)
        dlp_tau <- -(tau - priors$E_tau)^1 / priors$V_tau
        lp_ve <- -.5 * sum(((ve - priors$E_ve)^2) / priors$V_ve)
        dlp_ve <- -(ve - priors$E_ve)^1 / priors$V_ve
        lp_gamma1 <- -.5 * sum(((gamma1 - priors$E_gamma1)^2) / priors$V_gamma1)
        dlp_gamma1 <- -(gamma1 - priors$E_gamma1)^1 / priors$V_gamma1
        lp_gamma2 <- -.5 * sum(((gamma2 - priors$E_gamma2)^2) / priors$V_gamma2)
        dlp_gamma2 <- -(gamma2 - priors$E_gamma2)^1 / priors$V_gamma2
      }

      lp <- (lp_phi +
        lp_ve +
        lp_tau +
        lp_gamma1 +
        lp_gamma2
      )
      dlp <- matrix(c(
        dlp_phi,
        dlp_ve,
        dlp_tau,
        dlp_gamma1,
        dlp_gamma2
      ))

      lpost <- like + lp
      neglpost <- -lpost

      # J is gradient in parameter space - need gradient in transformed parameters
      GradLpost <- (dl + dlp) * dpars
      neglgrad <- -GradLpost
    }
    output <- list(neglpost = neglpost, neglgrad = neglgrad)
    return(output)
  }
}
