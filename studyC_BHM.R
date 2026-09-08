library(FDboost)

setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")

source("BHM_simulator.R")
source("stan_intro.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions_new.R")
source("EM_updating.R")
source("second_derivative_fd.R")

############## Simulation C - BHM ####################
# Parameter to vary: m = 3, 5, 8 (number of vertices)
# Fixed: sd_perc = 0.01, mu_p = rep(0, m-1), n = 500

library(compositions)

nsim  <- 50
m_params <- c(5, 8, 13)
seed0 <- 03072
n     <- 500
k     <- 23
sd0   <- 1

for (par_id in 2) {
  for (sim_id in 43:nsim) {
    
    param_c <- m_params[par_id]
    m       <- param_c
    seed    <- seed0 + sim_id * 10^5
    set.seed(seed)
    
    F_sample <- F_simulator(m = m, n = n, sd_perc = 0.01, mu_p = rep(0, m - 1))
    
    Sigma_eps <- F_sample$Sigma_eps
    Sigma_p   <- F_sample$pF$p$sigma
    mu_p      <- F_sample$pF$p$mu
    H         <- F_sample$pF$H$coefs
    f_coef    <- F_sample$nF
    D         <- F_sample$pF$H$D
    D0        <- F_sample$pF$H$D0
    get_basis <- F_sample$pF$H$basis
    
    ## Functional PCA (non-centred) ##
    Cmat <- cov(t(f_coef))
    A    <- D0 + D * 0
    M    <- D0 %*% Cmat %*% D0
    
    L       <- chol(A)
    Mat     <- solve(t(L), M)
    Mat     <- solve(L, Mat)
    Mat_fpca <- eigen(Mat)
    
    loadings  <- Mat_fpca$vectors
    pca_basis <- get_basis %*% loadings
    
    w <- diff(dom[1:2])
    for (j in 1:k) {
      norm_const      <- sqrt(sum(pca_basis[, j]^2) * w)
      pca_basis[, j]  <- pca_basis[, j] / norm_const
      loadings[, j]   <- loadings[, j] / norm_const
    }
    
    f_coef_pca <- solve(loadings) %*% f_coef
    
    D_pca  <- t(loadings) %*% D  %*% loadings
    D0_pca <- t(loadings) %*% D0 %*% loadings
    
    H0 <- diag(k)[, 1:m]
    
    EM_sample <- tryCatch({
      EM_updating(
        f_coef    = f_coef_pca,
        Sigma_eps0 = diag(sd0, k, k),
        H0        = H0,
        mu_p0     = rep(0, m - 1),
        Sigma_p0  = diag(1, m - 1, m - 1),
        B         = 50,
        D         = D_pca,
        lambdaS   = 1,
        lambdaH   = 0,
        gb        = get_basis,
        pb        = pca_basis,
        det_rate  = 0.2
      )
    },
    error = function(e) {
      message("Error at sim ", sim_id, ", m = ", m, ": ", e$message)
      NA
    })
    
    save(EM_sample, pca_basis, F_sample, get_basis,
         file = paste0("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_C/",
                       "simC_BHM_sim", sim_id, "par", param_c, "seed", seed, ".rdata"))
    print(sim_id)
  }
}
