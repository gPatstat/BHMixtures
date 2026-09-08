library(FDboost)

setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")

source("BHM_simulator.R")
source("stan_intro.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions_new.R")
source("EM_updating.R")
source("second_derivative_fd.R")

############## examples #################

# Parameters to vary:
# Simulation B
# sd=0.01,0.02,0.05,0.1 x norm of the mean density.

library(compositions)

nsim=50
sd_perc=c(0.01,0.1,0.2,0.5)
seed0=03072
m=3
n=500
k=23
sd0=1


for(par_id in 1:4){
for(sim_id in 1:nsim){
    param_a=sd_perc[par_id]
    seed=seed0+sim_id*10^5
    set.seed(seed)
    F_sample=F_simulator(m=m,n=n,sd_perc=param_a)
  
    Sigma_eps=F_sample$Sigma_eps
    Sigma_p=F_sample$pF$p$sigma
    mu_p=F_sample$pF$p$mu
    H=F_sample$pF$H$coefs
    f_coef=F_sample$nF #coefficients
    D=F_sample$pF$H$D
    D0=F_sample$pF$H$D0
    get_basis=F_sample$pF$H$basis
    
    ## Non centred PCA! ##
    # If you have scarsity in data to build the densities, use sparse-principal component
    
    D=F_sample$pF$H$D
    D0=F_sample$pF$H$D0
    
    Cmat <- cov(t(f_coef))
    A <- D0 + D*0  # try smaller penalty
    M <- D0 %*% Cmat %*% D0
    
    L <- chol(A)
    Mat <- solve(t(L), M)
    Mat <- solve(L, Mat)
    Mat_fpca <- eigen(Mat)
    
    loadings=Mat_fpca$vectors
    
    pca_basis=get_basis%*%loadings
    
    for(j in 1:k){
      norm_const=sqrt(sum(pca_basis[,j]^2)*w)
      pca_basis[,j]=pca_basis[,j]/norm_const
      loadings[,j]=loadings[,j]/norm_const
    }
    
    f_coef_pca=solve(loadings)%*%f_coef
    
    D_pca=t(loadings)%*%D%*%(loadings)
    D0_pca=t(loadings)%*%D0%*%(loadings)
    
    H0=diag(k)[,1:m]

    EM_sample=tryCatch({EM_updating(f_coef_pca, diag(sd0,k,k), H0 , rep(0,m-1), diag(1,m-1,m-1), 
                          B=50,D=D_pca,lambdaS=1,lambdaH=0,gb=get_basis,pb=pca_basis,det_rate = 0.2)
    
    },
    error = function(e) {
      message("Errors at iteration", sim_id, ": ", e$message)
      NA   
    })
    save(EM_sample,pca_basis,F_sample,get_basis, file = paste("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_A/simA_sim",sim_id,"par",param_a,"seed",seed,".rdata"))
    print(sim_id)
    }
}
