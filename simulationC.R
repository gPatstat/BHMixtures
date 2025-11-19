# Parameters to vary:
# Simulation A
# m=2,3,5,8,13,21 - Fibonacci's series.

nsim=1
dev_params=c(0,1/5,1/3,1/2,1,2,3)
sd0=0.1
seed=03072
m=3
n=150
k=23

seed=03072
set.seed(seed)

for(par_id in 1:7){
  for(sim_id in 1:nsim){
    
    param_c=dev_params[par_id]
    mu_p=c(0,param_c)
    
    F_sample=F_simulator(m=m,n=n,sd_perc=0.01,mu_p=mu_p)
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
    A <- D0 + D*1  # try smaller penalty
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
    
    EM_sample=EM_updating(f_coef=f_coef_pca, Sigma_eps0=diag(sd0,k,k), H0=H0,
                          mu_p0=rep(0,m-1), Sigma_p0=diag(1,m-1,m-1), 
                          B=20,D=D_pca,lambdaS=1,lambdaH=0)
    
    
    save(EM_sample,pca_basis,F_sample,get_basis, file = paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
  }
}


dev_params=c(0,1/5,1/3,1/2,1,2,3)

seed=03072
set.seed(seed)

x11()
par(mfrow=c(4,2))
for(par_id in 1:4){
  param_c=dev_params[par_id]
  load(paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
  verteces_display(F_sample$pF$H$coefs, get_basis, ylim=c(0,4), main=paste("True vertices with par_C=",param_c))
  verteces_display(EM_sample$H,pca_basis, ylim=c(0,4), main=paste("Estimated vertices with par_C=",param_c))
  #mu_true=as.matrix(t(F_sample$pF$p$mu))
  #mu_p=as.matrix(t(EM_sample$mu_p))
  #image(as.matrix(ilrInv(mu_true)))
  #image(as.matrix(ilrInv(mu_p)))
}

x11()
par(mfrow=c(4,2))
for(par_id in 5:7){
  param_c=dev_params[par_id]
  load(paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
  verteces_display(F_sample$pF$H$coefs, get_basis, ylim=c(0,4), main=paste("True vertices with par_C=",param_c))
  verteces_display(EM_sample$H,pca_basis, ylim=c(0,4), main=paste("Estimated vertices with par_C=",param_c))
  #mu_true=as.matrix(t(F_sample$pF$p$mu))
  #mu_p=as.matrix(t(EM_sample$mu_p))
  #image(as.matrix(ilrInv(mu_true)))
  #image(as.matrix(ilrInv(mu_p)))
}


x11()
par(mfrow=c(4,2))
for(par_id in c(4,5,6,7)){
  param_c=dev_params[par_id]
  load(paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
  verteces_display(F_sample$pF$H$coefs, get_basis, ylim=c(0,4), main=paste("True vertices with par_C=",param_c))
  verteces_display(EM_sample$H,pca_basis, ylim=c(0,4), main=paste("Estimated vertices with par_C=",param_c))
  #mu_true=as.matrix(t(F_sample$pF$p$mu))
  #mu_p=as.matrix(t(EM_sample$mu_p))
  #image(as.matrix(ilrInv(mu_true)))
  #image(as.matrix(ilrInv(mu_p)))
}

