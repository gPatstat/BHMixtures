
##### Without penalization


update_H=function(f_coef,Sigma_eps, H_old, mu_p, Sigma_p,B=100){
  n=dim(f_coef)[2]
  k=dim(f_coef)[1]
  m=length(mu_p)+1
  W=numeric(B)
  sm=matrix(0,k,m)
  ptp=matrix(0,m,m)
  for(i in 1:n){
    f0=f_coef[,i] # [k x 1]
    sample_prop=importance_sampling_p(f0,Sigma_eps, H_old, mu_p, Sigma_p,B=B)
    W=sample_prop$wgt        # [B x 1]
    if(sum(W)>0){
      #W=W/sum(W)
      X=sample_prop$proportions
      for (b in 1:B){
        p_x=as.numeric(ilrInv(X[b,]))
        sm=sm+W[b]*f0%*%t(p_x) # [k x 1]%*%[1 x m]=[k x m]
        ptp=ptp+W[b]*p_x%*%t(p_x) # [m x 1]%*%[1 x m]=[m x m]
      }
    }
  }
  
  ptp=t(ptp)%x%diag(k)
  sm=matrix(sm,k*m,1)
  
  solution=matrix(solve(ptp)%*%sm,k,m)
  
  return(solution)
}





H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H)


x11()
par(mfrow=c(4,2))
verteces_display(H, get_basis)
verteces_display(H_new, get_basis)


