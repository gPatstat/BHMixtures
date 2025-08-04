
#####################################################
############ Coordinate-wise updating ###############
#####################################################


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



update_mu_p=function(f_coef,Sigma_eps,mu_p_old,Sigma_p,H,B=1000){
  n=dim(f_coef)[1]
  k=dim(f_coef)[2]
  m=length(mu_p)+1
  W=numeric(B)
  mu0=numeric(m-1)
  sumW=0
  for(i in 1:n){
    f0=f_coef[i,] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H, mu_p_old, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    X=sample_prop$proportions # [B x (m-1)]
    mu0=mu0+t(X)%*%W
  }
  return(m0/sumW)
}

update_Sigma_p=function(f_coef,Sigma_eps,mu_p,Sigma_p_old,H,B=1000){
  n=dim(f_coef)[1]
  k=dim(f_coef)[2]
  m=length(mu_p)+1
  W=numeric(B)
  sigma0=matrix(0,m-1,m-1)
  sumW=0
  for(i in 1:n){
    f0=f_coef[i,] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H, mu_p, Sigma_p_old,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    X=sample_prop$proportions # [B x (m-1)]
    sigma0=sigma0+(X-mu_p)%*%t(X-mu_p)%*%W
  }
  return(sigma0/sumW)
}

update_Sigma_eps=function(f_coef,Sigma_eps_old,mu_p,Sigma_p,H,B=1000){
  n=dim(f_coef)[1]
  k=dim(f_coef)[2]
  m=length(mu_p)+1
  W=numeric(B)
  sigma0=matrix(0,k,k)
  sumW=0
  for(i in 1:n){
    f0=f_coef[i,] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H, mu_p, Sigma_p_old,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    P=matrix(ilrInv(sample_prop$proportions),1000,3)  # [B x m]
    for(b in 1:B){
      sigma0=sigma0+(f0-H%*%P[b,])%*%t(f0-H%*%P[b,])*W[b]
    }
  }
  return(sigma0/sumW)
}
