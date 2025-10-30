#####################################################
############### Importance sampling #################
#####################################################


# Sigma_eps: covariance matrix of the discretized residuals [k x k]
# H: matrix of vertices' coefficients in clr [k x m]
# f: vector of coefficients of a density in clr [k x 1]
# mu_p:  [(m-1) x 1]
# Sigma_p: [(m-1) x (m-1)]



posterior_builder=function(f0,Sigma_eps, H, mu_p, Sigma_p){
  library(compositions)
  options(digits = 22)  # aumenta cifre stampate
  inv_Sigma_eps=solve(Sigma_eps)
  inv_Sigma_p=solve(Sigma_p)
  f0=as.matrix(f0)
  func=function(x){
    p_x=as.numeric(ilrInv(x))
    val=-1/2*(t(x-mu_p)%*%inv_Sigma_p%*%(x-mu_p)+log(det(Sigma_p))+t(f0-H%*%p_x)%*%inv_Sigma_eps%*%(f0-H%*%p_x)+log(det(Sigma_eps)))
    return(-as.numeric(val))
  }
  return(func)
}


#install.packages("numDeriv")  # if not installed


importance_sampling_p=function(f0,Sigma_eps, H, mu_p, Sigma_p,B=100){
  library(scales)
  library(MASS)
  library(numDeriv)
  library(Matrix)
  mu_p=as.numeric(mu_p)
  post_prob=posterior_builder(f0,Sigma_eps, H, mu_p, Sigma_p)
  mu_opt=optim(par=mu_p,fn=post_prob)$par
  Hess=numDeriv::hessian(post_prob,mu_opt)
  if(all(eigen(Hess)$values>0)){
    Mat=Hess
  }else{
    Mat=solve(Sigma_p)
    mu_opt=mu_p
  }
  wgt=numeric(B)
  x=matrix(NA,B,length(mu_p))
  for(b in 1:B){
    x[b,]=mvrnorm(1,mu_opt,solve(Mat))
    wgt[b]=-post_prob(x[b,])+1/2*t(x[b,]-mu_opt)%*%Mat%*%(x[b,]-mu_opt)-1/2*log(det(Mat))
  }
  wgt=exp(wgt)
  if(sum(wgt)>0){
    wgt=wgt/sum(wgt)
  }
return(list(proportions=x,wgt=wgt))
}


#prova=importance_sampling_p(f_coef_pca[,10], diag(sd0,k,k), H0 , rep(0,m-1), diag(1,m-1,m-1),B=100)
#prova

