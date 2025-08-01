# f_coef here is the matrix of coefficients [n x k]

update_H=function(f_coef,Sigma_eps,mu_p,Sigma_p,H_old,B=1000){
  n=dim(f_coef)[1]
  k=dim(f_coef)[2]
  m=length(mu_p)+1
  W=numeric(B)
  sm=matrix(0,k,m)
  ptp=matrix(0,m,m)
  for(i in 1:n){
    f0=f_coef[i,] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H_old, mu_p, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    if(sum(W)>0){
    W=W/sum(W)
    P=matrix(ilrInv(sample_prop$proportions),1000,3) # [B x m]
    for (b in 1:B){
      sm=sm+W[b]*f0%*%t(P[b,]) # [k x 1]%*%[1 x m]=[k x m]
      ptp=ptp+W[b]*P[b,]%*%t(P[b,]) # [m x 1]%*%[1 x m]=[m x m]
    }
    }
  }
  return(sm%*%solve(ptp))
}

