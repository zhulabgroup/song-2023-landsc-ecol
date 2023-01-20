path_train <- paste0(path_sub, "train/")
path_model <- paste0(path_sub, "model/")
dir.create(path_model, recursive = T)

########################
# Get new data
########################
X_all <- read_csv(paste0(path_train, "X.csv")) %>% as.matrix()
Y_all <- read_csv(paste0(path_train, "Y.csv")) %>% as.matrix()
D_all <- read_csv(paste0(path_train, "D.csv")) %>% as.matrix()
P_all <- read_csv(paste0(path_train, "P.csv")) %>% as.matrix()

missing_id <- which(rowSums(is.na(cbind(Y_all, X_all))) != 0)
if (length(missing_id) > 0) {
  X_all <- X_all[-missing_id, , drop = F]
  Y_all <- Y_all[-missing_id, , drop = F]
  P_all <- P_all[-missing_id, , drop = F]
  D_all <- D_all[-missing_id, , drop = F]
}

set.seed(42)
train_id <- sample(1:nrow(X_all), round(nrow(X_all) * 5 / 10))
valid_id <- setdiff(1:nrow(X_all), train_id)

X_train <- X_all[train_id, , drop = F]
Y_train <- Y_all[train_id, , drop = F]
P_train <- P_all[train_id, , drop = F]
D_train <- D_all[train_id, , drop = F]

X_valid <- X_all[valid_id, , drop = F]
Y_valid <- Y_all[valid_id, , drop = F]
P_valid <- P_all[valid_id, , drop = F]
D_valid <- D_all[valid_id, , drop = F]

########################
# Get basis
########################
basisnumber <- min(basisnumber, nrow(X_all))
if (basisnumber == nrow(X_all)) {
  X_basis <- X_all
  Y_basis <- Y_all
  P_basis <- P_all
  D_basis <- D_all
} else {
  J_all <- date2doy(D_all)
  set.seed(42)
  cluster <- kmeans(cbind(P_all, J_all, X_all, Y_all), basisnumber)
  cluster_id_sort <- seq(basisnumber)[order(cluster$size, decreasing = T)]
  basis_id <- rep(NA, basisnumber)
  for (i in 1:basisnumber) {
    basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
  }
  X_basis <- X_all[basis_id, , drop = F]
  Y_basis <- Y_all[basis_id, , drop = F]
  P_basis <- P_all[basis_id, , drop = F]
  D_basis <- D_all[basis_id, , drop = F]
}

write_csv(as.data.frame(X_basis), paste0(path_model, "X_basis", ".csv"))
write_csv(as.data.frame(Y_basis), paste0(path_model, "Y_basis", ".csv"))
write_csv(as.data.frame(P_basis), paste0(path_model, "P_basis", ".csv"))
write_csv(as.data.frame(D_basis), paste0(path_model, "D_basis", ".csv"))

########################
# Initialize
########################
pars_id_prev <- numeric()
pars_prev_updated <- matrix(NA, nrow = ndim + 4, ncol = 0)
pars_var_prev_updated <- vector(mode = "list", length = 0)
log_p_prev_updated <- numeric(0)

########################
# Add new sets of hyperparameters
########################
set.seed(42)
particles <- matrix(NA, nrow = ndim + 4, ncol = num_part)
for (i in 1:ndim) {
  particles[i, ] <- rtrunc(n = num_part, spec = "norm", a = phimin, b = phimax, mean = priors$E_phi[i], sd = sqrt(priors$V_phi[i]))
  particles[i, ] <- -log((phimax - phimin) / (particles[i, ] - phimin) - 1) # transform
}

particles[ndim + 1, ] <- rtrunc(n = num_part, spec = "norm", a = vemin, b = vemax, mean = priors$E_ve, sd = sqrt(priors$V_ve))
particles[ndim + 1, ] <- -log((vemax - vemin) / (particles[ndim + 1, ] - vemin) - 1) # transform

particles[ndim + 2, ] <- rtrunc(n = num_part, spec = "norm", a = taumin, b = taumax, mean = priors$E_tau, sd = sqrt(priors$V_tau))
particles[ndim + 2, ] <- -log((taumax - taumin) / (particles[ndim + 2, ] - taumin) - 1) # transform

particles[ndim + 3, ] <- rtrunc(n = num_part, spec = "norm", a = gamma1min, b = gamma1max, mean = priors$E_gamma1, sd = sqrt(priors$V_gamma1))
particles[ndim + 3, ] <- -log((gamma1max - gamma1min) / (particles[ndim + 3, ] - gamma1min) - 1) # transform

particles[ndim + 4, ] <- rtrunc(n = num_part, spec = "norm", a = gamma2min, b = gamma2max, mean = priors$E_gamma2, sd = sqrt(priors$V_gamma2))
particles[ndim + 4, ] <- -log((gamma2max - gamma2min) / (particles[ndim + 4, ] - gamma2min) - 1) # transform

pars_add <- particles
pars_id_add <- 1:5

########################
# Optimize hyperparameters
########################

Rprop_res_add <-
  foreach(
    i = 1:num_part,
    .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
  ) %dopar% {
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    set.seed(42)

    lpost <- function(p, sample_n, seed) {
      set.seed(seed)
      minibatch <- sample(1:nrow(X_train), sample_n)
      return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisD = D_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
    }
    print("start")

    fmingrad_Rprop(lpost, pars_add[, i, drop = F] + matrix(rnorm(ndim + 4, 0, 0.1)), sample_n = min(nrow(X_train), 50), maxcount = maxcount)
  }

pars_add_updated <- matrix(NA, nrow = ndim + 4, ncol = num_part)
log_p_add_updated <- rep(NA, num_part)
for (i in 1:num_part) {
  pars_add_updated[, i] <- Rprop_res_add[[i]]$xopt
  log_p_add_updated[i] <- -Rprop_res_add[[i]]$fopt
}

########################
# Estimate uncertaity of hyperparameters using log likelihood
########################
pars_var_add_updated <-
  foreach(
    i = 1:num_part,
    .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
  ) %dopar% {
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    set.seed(42)

    pars_sample <- pars_add_updated[, i, drop = F] %*% matrix(1, nrow = 1, ncol = num_sample) + t(rmvn(num_sample, mu = rep(0, ndim + 4), Sigma = diag(1, nrow = ndim + 4, ncol = ndim + 4)))
    loglik_res <- rep(NA, num_sample)

    for (k in 1:num_sample) {
      lpost <- function(p, sample_n, seed) {
        set.seed(seed)
        minibatch <- sample(1:nrow(X_train), sample_n)
        return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisD = D_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
      }

      res <- lpost(pars_sample[, k, drop = F], sample_n = min(nrow(X_train), 100), seed = k)
      loglik_res[k] <- -res$neglpost
      print(k)
    }
    loglik_res <- loglik_res - mean(loglik_res)
    loglik_res[loglik_res < -50] <- -50
    loglik_res[loglik_res > 50] <- 50
    lik_res <- exp(loglik_res)
    pars_var_add_updated_i <- cov.wt(t(pars_sample),
      wt = as.numeric(exp(loglik_res)),
      cor = FALSE, center = t(pars_add_updated[, i, drop = F])
    )$cov
    pars_var_add_updated_i
  }

pars_new <- cbind(pars_prev_updated, pars_add_updated)
pars_var_new <- c(pars_var_prev_updated, pars_var_add_updated)
log_p_new <- c(log_p_prev_updated, log_p_add_updated)
pars_id_new <- c(pars_id_prev, pars_id_add)

########################
# Choose sets of hyperparameters using validation data
########################
valid_id_subset <- vector(mode = "list")
for (n in 1:10) {
  valid_id_subset[[n]] <- sample(1:nrow(X_valid), 100)
}

log_p_valid_new <-
  foreach(
    i = 1:ncol(pars_new),
    .combine = cbind
  ) %dopar% {
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    set.seed(42)

    num_sample <- 50
    loglik_res <- rep(NA, num_sample)

    for (k in 1:num_sample) {
      minibatch <- sample(1:nrow(X_valid), min(nrow(X_valid), 100))
      res <- GPSDM(pars = pars_new[, i, drop = F], distMat = distMat, basisX = X_basis, basisP = P_basis, basisD = D_basis, basisY = Y_basis, newX = X_valid[minibatch, , drop = F], newP = P_valid[minibatch, , drop = F], newD = D_valid[minibatch, , drop = F], newY = Y_valid[minibatch, , drop = F], priors = priors, mode = c("optimize"))
      loglik_res[k] <- -res$neglpost
      print(k)
    }

    median(loglik_res)
  }

sort_id <- rev(order(log_p_new))[1:num_part]
pars <- pars_new[, sort_id, drop = F]
pars_var <- pars_var_new[sort_id]
log_p <- log_p_new[sort_id]
pars_id <- pars_id_new[sort_id]
log_p_valid <- log_p_valid_new[sort_id]
colnames(pars) <- pars_id

########################
# Save model
########################
write_csv(as.data.frame(pars), paste0(path_model, "pars", ".csv"))
for (i in 1:num_part) {
  write_csv(as.data.frame(pars_var[[i]]), paste0(path_model, "pars_var_", i, ".csv"))
}
write_csv(as.data.frame(log_p), paste0(path_model, "log_p_train", ".csv"))
write_csv(as.data.frame(log_p_valid), paste0(path_model, "log_p_valid", ".csv"))
