library(FDboost)

setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")
source("second_derivative_fd.R")
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




simul_func=function(sim_id){
  
  tryCatch({
    
    set.seed(sim_id)
    log_message <- paste("Start iteration", sim_id, "at", Sys.time(), "\n")
    write(log_message, file = "log.txt", append = TRUE)
    
    F_sample=F_simulator(m=m,n=n,sd_perc=param_b)
    Sigma_eps=F_sample$Sigma_eps
    Sigma_p=F_sample$pF$p$sigma
    mu_p=F_sample$pF$p$mu
    H=F_sample$pF$H$coefs
    f_coef=F_sample$nF #coefficients
    D=F_sample$pF$H$D
    D0=F_sample$pF$H$D0
    get_basis=F_sample$pF$H$basis
    
    ## Non centred PCA ##
    Cmat <- cov(t(f_coef))
    A <- D0 + D*0
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
    lamH=Mat_fpca$values[1]
    
    EM_sample=EM_updating(
      f_coef_pca, diag(sd0,k,k), H0, rep(0,m-1),
      diag(1,m-1,m-1),
      B=50, D=D_pca, lambdaS=1, lambdaH=10^-2,
      gb=get_basis, pb=pca_basis
    )
    
    save(EM_sample,pca_basis,F_sample,get_basis,
         file = paste("simB_sim",sim_id,"par",param_b,"seed",seed,".rdata"))
    
  },
  error=function(e){
    msg <- paste("ERROR in sim", sim_id, ":", conditionMessage(e), "at", Sys.time(), "\n")
    write(msg, file="log.txt", append=TRUE)
    
    return(NULL)   
  })
}


nsim=100
sd_perc=c(0.01,0.1,0.2,0.5)
sd0=0.1
seed=03072
m=3
n=150
k=23

numcores=6

for(par_id in 1:4){
  cl=makeCluster(numcores)
  param_b<<-sd_perc[par_id]
  #clusterExport(cl,varlist=c("param_b","nsim","m","k","n","sd0"))
  ls_env=ls()
  clusterExport(cl,varlist=ls_env)
  clusterEvalQ(cl, ({
    setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")
    library(compositions)
    library(scales)
    library(MASS)
    library(numDeriv)
    library(Matrix)
    library(FDboost)
    library(mvtnorm)
    library(base)
  }))
  
  #clusterSetRNGStream(cl, 0307)
  sim_id=1:nsim
  parLapply(cl,X=sim_id,fun=simul_func)
  
  stopCluster(cl)
}

sd_perc=c(0.01,0.1,0.2,0.5)

seed=03072
set.seed(seed)

x11()
par(mfrow=c(4,2))
for(id in 1:8){
  #param_b=sd_perc[par_id]
  load(paste("simB_sim",id,"par",0.5,"seed",seed,".rdata"))
  verteces_display(F_sample$nF, get_basis, main=paste("True vertices with par_B=",sd_perc[par_id]),ylim=c(0,4))
  #verteces_display(EM_sample$H,pca_basis, main=paste("Estimated vertices with par_B=",sd_perc[par_id]),ylim=c(0,4))
}

