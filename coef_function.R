coef_f=function(f,k=13){  
  library(FDboost)
  n=dim(f)[2]
  coefs=matrix(0,n,k)
  for(id in 1:n){
    dclrs=list(x=t(f[,id]),t=domain)
    
    cclrs=FDboost(x~1, 
                  #use bbsc() in time formula to ensure integrate-to-zero constraint 
                  timeformula= ~bbsc(t,df=4, 
                                     knots=10,
                                     boundary.knots=c(0,1),
                                     degree=3,
                                     lambda=10^2),
                  data=dclrs,offset=0, 
                  control=boost_control(mstop=100))
    
    coefs[id,]=as.numeric(cclrs$coef(which=1)[[1]])
  }
  return(coefs)
}
