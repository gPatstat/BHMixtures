#######################################
######## BHMixtures simulator #########
#######################################


#######################################
######## Verteces simulator ###########
#######################################

source("clr2density.R")
source("second_derivative_fd.R")
library(FDboost)
library(compositions)

dom=seq(0,1,length.out=100)
w=dom[2]-dom[1]
simcoef=function(a,b,dm=dom,B=100){
  
  stopifnot(is.numeric(dm), length(dm) >= 2)
  
  w <- dm[2] - dm[1]

  #multin=rmultinom(10^5,B,1/B)
  beta=rbeta(n=B,a,b)
  histbeta=hist(beta,breaks=seq(0,1,0.01),plot=F)
  wbin=1/length(histbeta$counts)
  counts=histbeta$counts
  
  sum0=sum(counts==0)
  counts=ifelse(counts==0,0.01,counts)
  tot_counts=B+sum0*0.01
  ddens=counts/(tot_counts*wbin)
  dclrs=list(x=matrix(clr(ddens,w=wbin,inverse=F),1,length(histbeta$mids)),t=histbeta$mids)
  
  cclrs=FDboost(x~1, 
                #use bbsc() in time formula to ensure integrate-to-zero constraint 
                timeformula= ~bbsc(t,df=4, 
                                   knots=20,
                                   boundary.knots=c(0,1),
                                   degree=3,
                                   lambda=10^3),
                data=dclrs,offset=0, 
                control=boost_control(mstop=100))
  t=histbeta$mids
  basis=bbsc(t,df=4,knots=20,boundary.knots=c(0,1),degree=3)
  get_basis=extract(basis,"design",asmatrix=T) # 10 breaks x 11 elements of basis
  coefs=as.numeric(cclrs$coef(which=1)[[1]])
  new_basis <- bbsc(dm, df = 4, knots = 20, boundary.knots = c(0, 1), degree = 3,lambda=10^5)
  # Extract design matrix from the new basis
  get_basis_new <- extract(new_basis, "design", asmatrix = TRUE)
  D0=t(get_basis_new)%*%get_basis_new*w
  
  sdf=function(fx){
    second_derivative_fd(fx,dm)
  }
  
  basis_dd=na.omit(apply(get_basis_new,2,sdf))
  D <- t(basis_dd) %*% basis_dd * w  # approximate integral
  
  return(list(coefs=coefs,basis=get_basis_new, D=D,D0=D0))
}

H_simulator <- function(m,k=23){
  library(FDboost)
  library(mvtnorm)
  mat <- matrix(0, k, m)
  for (j in 1:m) {
    a=j+1
    b=m-j+2
    beta_sample=simcoef(a,b,dom)
    mat[,j]=beta_sample$coef
  }
  return(list(coefs=mat,basis=beta_sample$basis, D=beta_sample$D,D0=beta_sample$D0))
}


verteces_display <- function(H, get_basis,add=F,ylim=NULL,main="Densities",xaxt = NULL, yaxt = NULL){
  k=dim(H)[1]
  library(FDboost)
  library(MASS)
  dom=seq(0,1,length.out=100)
  w=diff(dom[1:2])
  #clr basis dim: 100 x k
  clrs=get_basis%*%H
  densities=clr2density(clrs,w=w)
  matplot(dom,densities,type="l",main=main,add=add,ylim=ylim,lty=1:k,lwd=2,
          col="black", xaxt=xaxt,yaxt=yaxt)
}


#  Example
#H_sample=H_simulator(m=4,k=23)
#verteces_display(H_sample)



#######################################
###### Proportions simulator ##########
#######################################

#mu_p
#Sigma_p

p_simulator<-function(m,n,mu_p,var=4){
  library(compositions)
  
  #A <- matrix(rnorm((m-1)^2), (m-1), (m-1))
  #Sigma_p= t(A) %*% A
  Sigma_p=diag(var,m-1,m-1)
  library(MASS)
  p_ilr=rmvnorm(n,mu_p,Sigma_p)
  p=matrix(0,m,n)
  for(i in 1:n){
    p[,i]=ilrInv(p_ilr[i,])
  }
  return(list(mu=mu_p,sigma=Sigma_p,pilr=p_ilr,p=p))
}






#######################################
###### Pure mixture simulator #########
#######################################


pF_simulator<- function(m,n,mu_p){
  H_sample=H_simulator(m) #k x m
  psim=p_simulator(m,n,mu_p)
  pp=psim$p # m x n
  
  pF=as.matrix(H_sample$coefs%*%as.matrix(pp)) # k x n
  
  return(list(pF=pF, H=H_sample,p=psim))
}


#######################################
###### Noisy mixture simulator #########
#######################################


F_simulator<- function(m,k=23,n,sd_perc,mu_p=rep(0,m-1)){
  #A <- matrix(rnorm((k)^2,sd=sd), k, k)
  #Sigma_eps= t(A) %*% A
  pF_sim=pF_simulator(m=m,n=n,mu_p)
  pF=pF_sim$pF
  get_basis=pF_sim$H$basis
  
  sd=sd_perc*sd(pF)
  Sigma_eps=diag(sd^2,k)
  coef_noise=t(rmvnorm(n,rep(0,k),Sigma_eps))
  
  nF=pF+coef_noise
  
  return(list(pF=pF_sim,nF=nF,Sigma_eps=Sigma_eps,coef_noise=coef_noise,basis=pF_sim$H$basis,D=pF_sim$H$D,D0=pF_sim$H$D0))
}




