clr2density=function(clrs,w){
  dens=matrix(NA,dim(clrs)[1],dim(clrs)[2])
  for(i in 1:(dim(dens)[2])){
    int=sum(exp(clrs[,i]))*w
    dens[,i]=exp(clrs[,i])/int
  }
  return(dens)
}
