
#####################################################
############ Coordinate-wise updating ###############
#####################################################


# f_coef here is the matrix of coefficients [n x k]

update_H=function(f_coef,Sigma_eps, H_old, mu_p, Sigma_p,B=100, lambda=10^-4,D){
  n=dim(f_coef)[2]
  k=dim(f_coef)[1]
  m=length(mu_p)+1
  W=numeric(B)
  C=matrix(0,k,m)
  A_r=matrix(0,m,m)
  A_l=2*lambda*(Sigma_eps%*%D)
   # print(dim(A_l))
   # print(dim(A_r))
  for(i in 1:n){
    f0=f_coef[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H_old, mu_p, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    if(sum(W)>0){
    #W=W/sum(W)
    X=sample_prop$proportions
    for (b in 1:B){
      p_x=as.numeric(ilrInv(X[b,]))
      C=C+W[b]*f0%*%t(p_x) # [k x 1]%*%[1 x m]=[k x m]
      A_r=A_r+W[b]*p_x%*%t(p_x) # [m x 1]%*%[1 x m]=[m x m]

    }
    }
  }
  vec_C=matrix(C,k*m,1,byrow=F)

  den=t(A_r)%x%diag(k)-diag(m)%x%A_l
  den=+diag(m)%x%A_l+t(A_r)%x%diag(k)
  vec_solution=solve(den)%*%vec_C
  solution=matrix(vec_solution,k,m,byrow=F)

  return(solution)
}



update_mu_p=function(f_coef,Sigma_eps, H, mu_p_old, Sigma_p,B=100){
  n=dim(f_coef)[2]
  k=dim(f_coef)[1]
  m=length(mu_p)+1
  W=numeric(B)
  mu0=matrix(0,m-1,1)
  mu_p_old=as.numeric(mu_p_old)
  sumW=0
  for(i in 1:n){
     f0=f_coef[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H, mu_p_old, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    X=sample_prop$proportions # [B x (m-1)]
    mu0=mu0+t(X)%*%W
  }
  return((mu0/sumW))
}

update_Sigma_p=function(f_coef,Sigma_eps, H, mu_p, Sigma_p_old,lambda=100,B=100){
  n=dim(f_coef)[2]
  k=dim(f_coef)[1]
  m=length(mu_p)+1
  W=numeric(B)
  sigma0=matrix(0,m-1,m-1)
  sumW=0
  for(i in 1:n){
     f0=f_coef[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H, mu_p, Sigma_p_old,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    X=sample_prop$proportions # [B x (m-1)]
    for(b in 1:B){
      sigma0=sigma0+(X[b,]-mu_p)%*%t(X[b,]-mu_p)*W[b] # colum vec %*% row vec* scalar
    }
  }##
  return((sigma0/sumW)+lambda*diag(1,m-1))
}

update_Sigma_eps=function(f_coef,Sigma_eps_old, H, mu_p, Sigma_p,B=100){
  n=dim(f_coef)[2]
  k=dim(f_coef)[1]
  m=length(mu_p)+1
  W=numeric(B)
  sigma0=matrix(0,k,k)
  sumW=0
  for(i in 1:n){
    f0=f_coef[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps_old, H, mu_p, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    sumW=sumW+sum(W)
    X=sample_prop$proportions  # [B x m]
    for(b in 1:B){
      p_x=as.numeric(ilrInv(X[b,]))
      sigma0=sigma0+(f0-H%*%p_x)%*%t(f0-H%*%p_x)*W[b] # colum vec %*% row vec
    }
  }
  return(sigma0/sumW)
}





