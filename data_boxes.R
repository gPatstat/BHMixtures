
setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")

source("vertex_error.R")

nsim=50

data_boxA=numeric()
setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_A")
seed0=03072
sd_perc=c(0.01,0.1,0.2,0.5)
for(par_id in 1:4){
  param_b=sd_perc[par_id]
  for(sim_id in 1:nsim){
  seed=seed0+10^5*sim_id
  tryCatch({
  load(paste("simB_sim",sim_id,"par",param_b,"seed",seed,".rdata"))
  if(!is.complex(pca_basis)){
  mat1=pca_basis%*%EM_sample$H
  mat2=get_basis%*%F_sample$pF$H$coefs
  #print(is.complex(pca_basis))
  data_boxA=rbind(data_boxA,c(error_vertex(mat1,mat2),param_b))
  }
  },
  error=function(e){
    msg <- paste("ERROR in sim", sim_id, ":", conditionMessage(e), "at", Sys.time(), "\n")
    write(msg, file="log.txt", append=TRUE)
    return(NULL)   
  })
  }
}

data_boxB=numeric()
setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_B")
seed0=03072
dev_params=c(1/2,1,2)
for(par_id in 1:3){
  param_c=dev_params[par_id]
  for(sim_id in 1:nsim){
    seed=seed0+10^5*sim_id
    tryCatch({
    load(paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
    mat1=pca_basis%*%EM_sample$H
    mat2=get_basis%*%F_sample$pF$H$coefs
    #print(error_vertex(mat1,mat2))
    data_boxB=rbind(data_boxB,c(error_vertex(mat1,mat2),param_c))
    },
    error=function(e){
      msg <- paste("ERROR in sim", sim_id, ":", conditionMessage(e), "at", Sys.time(), "\n")
      write(msg, file="log.txt", append=TRUE)
      
      return(NULL)   
    })
  }
}


base=data_boxA[data_boxA[,2]==0.01,1]

#data_boxA=rbind(data_boxA, cbind(base,rep(3,nsim)))
data_boxB=rbind(cbind(base,rep(0,nsim)), data_boxB)

x11()
boxplot(data_boxA[,1]~data_boxA[,2],outline=F,xlab="par_A",ylab="L2 error - best vertex permutation", main="Simulation A",ylim=c(0,1.4))
boxplot(data_boxA_lm[,1]~data_boxA_lm[,2],outline=F,xlab="par_A",ylab="L2 error - best vertex permutation", main="Simulation A",ylim=c(0,1.4))
boxplot(data_boxB[,1]~data_boxB[,2],outline=F,xlab="par_B",ylab="L2 error - best vertex permutation", main="Simulation B",ylim=c(0,1.4))
boxplot(data_boxB_lm[,1]~data_boxB_lm[,2],outline=F,xlab="par_B",ylab="L2 error - best vertex permutation", main="Simulation B",ylim=c(0,1.4))

