library(FDboost)

setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")
source("BHM_simulator.R")
source("Sampling_p.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions.R")
source("EM_updating.R")


############## examples #################
#set.seed(03071)

# Parameters to vary:
# Simulation B
# sd=0.01,0.02,0.05,0.1 x norm of the mean density.

library(compositions)

nsim=1
sd_perc=c(0.01,0.1,0.2,0.5)
sd0=0.1
seed=03072
m=3
n=150
k=23

set.seed(seed)
for(par_id in 2:4){
  
param_b=sd_perc[par_id]

for(sim_id in 1:nsim){
    
    F_sample=F_simulator(m=m,n=n,sd_perc=param_b)
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

    Cmat <- cov(t(f_coef))
    A <- D0 + D*1  # try smaller penalty
    M <- D0 %*% Cmat %*% D0
    
    L <- chol(A)
    Mat <- solve(t(L), M)
    Mat <- solve(L, Mat)
    Mat_fpca <- eigen(Mat)
    
    
    #D0[i,j]=integral(elem[i]*elem[j])
    loadings=Mat_fpca$vectors
    
    pca_basis=get_basis%*%loadings
    
    for(j in 1:k){
      norm_const=sqrt(sum(pca_basis[,j]^2)*w)
      pca_basis[,j]=pca_basis[,j]/norm_const
      loadings[,j]=loadings[,j]/norm_const
      #for(i in 1:n){
      #  f_coef_pca[j,i]=w*t(pca_basis[,j])%*%(get_basis)%*%f_coef[,i] #(1x100)x(100x23)x(23x1)
      #}
    }
    
    f_coef_pca=solve(loadings)%*%f_coef
    
    D_pca=t(loadings)%*%D%*%(loadings)
    D0_pca=t(loadings)%*%D0%*%(loadings)
    
    H0=diag(k)[,1:m]
    
    lamH=Mat_fpca$values[1]

    EM_sample=EM_updating(f_coef_pca, diag(sd0,k,k), H0 , rep(0,m-1), diag(1,m-1,m-1), 
                          B=20,D=D_pca,lambdaS=1,lambdaH=10^-2,gb=get_basis,pb=pca_basis)
    
    save(EM_sample,pca_basis,F_sample,get_basis, file = paste("simB_sim",sim_id,"par",param_b,"seed",seed,".rdata"))
 }
}

sd_perc=c(0.01,0.1,0.2,0.5)

seed=03072
set.seed(seed)

x11()
par(mfrow=c(4,2))
for(par_id in 1:4){
  param_b=sd_perc[par_id]
  load(paste("simB_sim",sim_id,"par",param_b,"seed",seed,".rdata"))
  verteces_display(F_sample$pF$H$coefs, get_basis, main=paste("True vertices with par_B=",sd_perc[par_id]),ylim=c(0,4))
  verteces_display(EM_sample$H,pca_basis, main=paste("Estimated vertices with par_B=",sd_perc[par_id]),ylim=c(0,4))
  #mu_true=as.matrix(t(F_sample$pF$p$mu))
  #mu_p=as.matrix(t(EM_sample$mu_p))
  #image(as.matrix(ilrInv(mu_true)))
  #image(as.matrix(ilrInv(mu_p)))
}
 
