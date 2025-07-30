#######################################
######## BHMixtures simulator #########
#######################################



#######################################
######## Verteces simulator ###########
#######################################

source("clr2density.R")



beta2coef=function(a,b,domain){
library(FDboost)
beta=rbeta(10000,a,b)
histbeta=hist(beta,breaks=10,plot=F)
wbin=1/length(histbeta$counts)
counts=ifelse(histbeta$counts==0,1,histbeta$counts)


ddens=histbeta$counts/((10000+sum(histbeta$counts==0))*wbin)
dclrs=list(x=matrix(clr(ddens,w=wbin,inverse=F),1,length(histbeta$mids)),t=histbeta$mids)

cclrs=FDboost(x~1, 
              #use bbsc() in time formula to ensure integrate-to-zero constraint 
              timeformula= ~bbsc(t,df=4, 
                                 knots=10,
                                 boundary.knots=c(0,1),
                                 degree=3,
                                 lambda=10^2),
              data=dclrs,offset=0, 
              control=boost_control(mstop=100))
t=histbeta$mids
basis=bbsc(t,df=4,knots=10,boundary.knots=c(0,1),degree=3)
get_basis=extract(basis,"design",asmatrix=T) # 10 breaks x 11 elements of basis
coefs=as.numeric(cclrs$coef(which=1)[[1]])
new_basis <- bbsc(domain, df = 4, knots = 10, boundary.knots = c(0, 1), degree = 3,lambda=10^2)
# Extract design matrix from the new basis
get_basis_new <- extract(new_basis, "design", asmatrix = TRUE)
return(list(coefs=coefs,basis=get_basis_new))
}

H_simulator <- function(m, k){
  library(FDboost)
  library(mvtnorm)
  mat <- matrix(0, k, m)
  for (j in 1:m) {
  S=runif(1,1,4)
    for(s in 1:S){
      a=runif(1,0.1,5)
      b=runif(1,0.1,5)
      beta_sample=beta2coef(a,b,domain)
      mat[,j]=mat[,j]+beta_sample$coef
    }
  }
  return(list(coefs=mat,basis=beta_sample$basis))
}


verteces_display <- function(H, get_basis){
  k=dim(H)[1]
  library(FDboost)
  library(MASS)
  domain=seq(0,1,length.out=100)
  w=diff(domain[1:2])
  #clr basis dim: 100 x k
  clrs=get_basis%*%H
  densities=clr2density(clrs,w=w)
  par(mfrow=c(1,2))
  matplot(domain,clrs,type="l",main="clrs")
  matplot(domain,densities,type="l",main="densities")
}


#  Example
#H_sample=H_simulator(m=4,k=13)
#verteces_display(H_sample)



#######################################
###### Proportions simulator ##########
#######################################

#mu_p
#Sigma_p

p_simulator<-function(m,n){
library(clusterGeneration)
library(compositions)

mu_p=rnorm((m-1),sd=1)
A <- matrix(rnorm((m-1)^2), (m-1), (m-1))
Sigma_p= t(A) %*% A

library(MASS)
p_ilr=rmvnorm(n,mu_p,Sigma_p)
p=apply(p_ilr,1,ilrInv)
return(list(mu=mu_p,sigma=Sigma_p,pilr=p_ilr,p=p))
}






#######################################
###### Pure mixture simulator #########
#######################################


pF_simulator<- function(m,k=13,n){
H_sample=H_simulator(m,k=13) #k x m
p_sample=p_simulator(m,n)$p # m x n

vert=(H_sample$basis)%*%(H_sample$coefs) # 100 x k
pF=as.matrix(vert%*%as.matrix(p_sample)) # k x n

par(mfrow=c(1,2))
matplot(domain,vert,type="l")
matplot(domain,(pF),type="l")
return(list(pF=pF,vert=vert, H=H_sample,p=p_sample))
}


pF_simulator(m=4,n=100)


#######################################
###### Noisy mixture simulator #########
#######################################


F_simulator<- function(m,k=13,n,sd){
  A <- matrix(rnorm((k)^2,sd=sd), k, k)
  Sigma_eps= t(A) %*% A
  pF_sim=pF_simulator(m=m,n=n)
  pF=pF_sim$pF
  Hs=pF_sim$H
  BH_noise=(Hs$basis)%*%t(rmvnorm(n,rep(0,k),Sigma_eps))
  
  nF=pF+BH_noise
  
  Gvert=clr2density(pF_sim$vert,w)
  pG=clr2density(pF,w)
  nG=clr2density(nF,w)
  
  par(mfrow=c(1,3))
  matplot(domain,Gvert,type="l",main="Vertices")
  matplot(domain,pG,type="l",main="Pure mixtures")
  matplot(domain,nG,type="l",main="Noisy mixtures")
  
  library(rgl)
  plot3D(t(pF_sim$p))
  
  return(list(pF=pF_sim,nF=nF,pG=pG,nG=nG))
}


F_simulator(m=3,n=100,sd=0.05)


