

setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")
data_boxA=numeric()
nsim=100
seed=03072
for(par_id in 1:3){
  m=m_params[par_id]
  param_a=m_params[par_id]
  for(sim_id in 1:nsim){
    tryCatch({
    load(paste("simA_sim",sim_id,"par",param_a,"seed",seed,".rdata"))
    if(!is.complex(pca_basis)){
        mat1=pca_basis%*%EM_sample$H
        mat2=get_basis%*%F_sample$pF$H$coefs
    data_boxA=rbind(data_boxA,c(error_vertex(mat1,mat2),param_a))
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

seed=03072
for(par_id in 1:4){
  param_b=sd_perc[par_id]
  for(sim_id in 1:nsim){
  tryCatch({
  load(paste("simB_sim",sim_id,"par",param_b,"seed",seed,".rdata"))
  if(!is.complex(pca_basis)){
  mat1=pca_basis%*%EM_sample$H
  mat2=get_basis%*%F_sample$pF$H$coefs
  #print(is.complex(pca_basis))
  data_boxB=rbind(data_boxB,c(error_vertex(mat1,mat2),param_b))}
  },
  error=function(e){
    msg <- paste("ERROR in sim", sim_id, ":", conditionMessage(e), "at", Sys.time(), "\n")
    write(msg, file="log.txt", append=TRUE)
    return(NULL)   
  })
  }
}

data_boxC=numeric()
seed=03072
for(par_id in 1:4){
  param_c=dev_params[par_id]
  for(sim_id in 1:nsim){
    tryCatch({
    load(paste("simC_sim",sim_id,"par",param_c,"seed",seed,".rdata"))
    mat1=pca_basis%*%EM_sample$H
    mat2=get_basis%*%F_sample$pF$H$coefs
    #print(error_vertex(mat1,mat2))
    data_boxC=rbind(data_boxC,c(error_vertex(mat1,mat2),param_c))
    },
    error=function(e){
      msg <- paste("ERROR in sim", sim_id, ":", conditionMessage(e), "at", Sys.time(), "\n")
      write(msg, file="log.txt", append=TRUE)
      
      return(NULL)   
    })
  }
}


x11()
par(mfrow=c(3,1))
boxplot(data_boxA[,1]~data_boxA[,2],outline=F)
boxplot(data_boxB[,1]~data_boxB[,2],outline=F)
boxplot(data_boxC[,1]~data_boxC[,2],outline=F)
