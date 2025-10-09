#######################################
######## BHMixtures simulator #########
#######################################



#######################################
######## Verteces simulator ###########
#######################################
#setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")
source("clr2density.R")

domain=seq(0,1,length.out=100)
w=diff(domain[1:2])

simcoef=function(domain,B=4){
library(FDboost)

a=runif(1,0.5,5)
b=runif(1,0.5,5)
beta=rbeta(1000,a,b)
histbeta=hist(beta,breaks=seq(0,1,0.1),plot=F)
wbin=1/length(histbeta$counts)
counts=histbeta$counts
  
#for(l in 1:(B-1)){
#a=runif(1,0.5,5)
#b=runif(1,0.5,5)
#beta=rbeta(1000,a,b)
#histbeta=hist(beta,breaks=seq(0,1,0.1),plot=F)
#wbin=1/length(histbeta$counts)
#counts=counts+histbeta$counts
#}

sum0=sum(counts==0)
counts=ifelse(counts==0,1,counts)
tot_counts=1000*B+sum0
ddens=counts/(tot_counts*wbin)
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
basis_dd <- extract(new_basis, "design", derivative = 2, asmatrix = TRUE)  # (200 x nbasis) matrix

D <- t(basis_dd) %*% basis_dd * w  # approximate integral

return(list(coefs=coefs,basis=get_basis_new, D=D))
}

H_simulator <- function(m,k=13){
  library(FDboost)
  library(mvtnorm)
  mat <- matrix(0, k, m)
  for (j in 1:m) {
    beta_sample=simcoef(domain, B=4)
    mat[,j]=beta_sample$coef
  }
  return(list(coefs=mat,basis=beta_sample$basis, D=beta_sample$D))
}


verteces_display <- function(H, get_basis,add=F){
  k=dim(H)[1]
  library(FDboost)
  library(MASS)
  domain=seq(0,1,length.out=100)
  w=diff(domain[1:2])
  #clr basis dim: 100 x k
  clrs=get_basis%*%H
  densities=clr2density(clrs,w=w)
  matplot(domain,densities,type="l",main="densities",add=add)
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
p=matrix(0,m,n)
for(i in 1:n){
  p[,i]=ilrInv(p_ilr[i,])
}
return(list(mu=mu_p,sigma=Sigma_p,pilr=p_ilr,p=p))
}






#######################################
###### Pure mixture simulator #########
#######################################


pF_simulator<- function(m,n){
H_sample=H_simulator(m) #k x m
psim=p_simulator(m,n)
pp=psim$p # m x n

pF=as.matrix(H_sample$coefs%*%as.matrix(pp)) # k x n

return(list(pF=pF, H=H_sample,p=psim))
}


#######################################
###### Noisy mixture simulator #########
#######################################


F_simulator<- function(m,k=13,n,sd){
  #A <- matrix(rnorm((k)^2,sd=sd), k, k)
  #Sigma_eps= t(A) %*% A
  Sigma_eps=diag(sd,k)
  pF_sim=pF_simulator(m=m,n=n)
  pF=pF_sim$pF
  coef_noise=t(rmvnorm(n,rep(0,k),Sigma_eps))

  
  nF=pF+coef_noise
    
  return(list(pF=pF_sim,nF=nF,Sigma_eps=Sigma_eps,coef_noise=coef_noise,basis=pF_sim$H$basis,D=pF_sim$H$D))
}


F_sample=F_simulator(m=4,n=100,sd=0.02)






