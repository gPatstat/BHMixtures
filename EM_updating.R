EM_updating = function(f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0, 
                       B = 100, D, lambdaS = 10, lambdaH = 10,gb,pb,det_rate=1) {
  ret=NULL
  V=ilrBase(D=(dim(H0)[2]))
  for (b in 1:B) {
    t1=Sys.time()
    # Save previous values
    #mu_p_prev <- mu_p0
    #Sigma_p_prev <- Sigma_p0
    
    H_prev=H0
    n_sub=ceiling(dim(f_coef)[2]*det_rate)
    
    # Update parameters
    f_coef_sub=f_coef[,sample(1:(dim(f_coef)[2]),n_sub)]
    Sigma_eps0 <- update_Sigma_eps(f_coef = f_coef_sub, Sigma_eps0, H0, mu_p0, Sigma_p0,V=V)
    print(paste("Sigma_eps done at iteration ",b))
    
    f_coef_sub=f_coef[,sample(1:(dim(f_coef)[2]),n_sub)]
    H0 <- update_H(f_coef = f_coef_sub, Sigma_eps0, H0, mu_p0, Sigma_p0, D = D, lambda = lambdaH,V=V)
    print(paste("H done at iteration ",b))
    
    f_coef_sub=f_coef[,sample(1:(dim(f_coef)[2]),n_sub)]
    mu_p0 <- update_mu_p(f_coef = f_coef_sub, Sigma_eps0, H0, mu_p0, Sigma_p0,V=V)
    print(paste("mu_p done at iteration ",b))
    
    f_coef_sub=f_coef[,sample(1:(dim(f_coef)[2]),n_sub)]
    Sigma_p0 <- update_Sigma_p(f_coef = f_coef_sub, Sigma_eps0, H0, mu_p0, Sigma_p0, lambda = lambdaS,V=V)
    print(paste("Sigma_p done at iteration ",b))
    
    # Compute relative changes
    #mu_change <- sqrt(sum((mu_p0-mu_p_prev)^2)) / sqrt(sum((mu_p_prev + 1e-8)^2))
    #Sigma_change <- norm(as.matrix(Sigma_p0-Sigma_p_prev), type = "F") / (norm(as.matrix(Sigma_p_prev), type = "F") + 1e-8)
    H_change <- norm(as.matrix(H0-H_prev), type = "F") / (norm(as.matrix(H_prev), type = "F") + 1e-8)
    
    # Print iteration number
    #cat("Iteration:", b, " | mu_change:", round(mu_change, 4), " | Sigma_change:", round(Sigma_change, 4), "\n")
    cat("Iteration:", b, " | H:", round(H_change, 4), "\n")

    
    if (H_change < 0.05) {
      cat("Converged after", b, "iterations (changes < 5%)\n")
      break
    }
    
    print(Sys.time()-t1)
    ret=list(mu_p = mu_p0, Sigma_p = Sigma_p0, H = H0, Sigma_eps = Sigma_eps0)
  }
  return(ret)
}


