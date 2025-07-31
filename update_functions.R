source("BHM_simultor.R")
source("Sampling_p.R")

# f here is the matrix of coefficients [n x k]

update_H=function(f,Sigma_eps,mu_p,Sigma_p,H_old,B=1000){
  n=1:dim(f)[2]
  W=numeric(B)
  sm=0
  ptp=0
  for(i in 1:n){
    f0=f[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H_old, mu_p, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    if(sum(W)>0){
    W=W/sum(W)
    P=sample_prop$proportions # [B x m]
    for (b in 1:B){
      sm=sm+W[b]*f0%*%t(P[b,]) # [k x 1]%*%[1 x m]=[k x m]
      ptp=ptp+W[b]*P[b,]%*%t(P[b,]) # [m x 1]%*%[1 x m]=[m x m]
    }
    }
  }
  return(sm%*%solve(ptp))
}


update_H(f_coef,Sigma_eps,mu_p,Sigma_p,matrix(0,3,3))


#create f_coef !!