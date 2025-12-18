EM_updating = function(f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0, 
                       B = 100, D, lambdaS = 10, lambdaH = 10,gb,pb) {
  for (b in 1:B) {
    #t1=Sys.time()
    # Save previous values
    mu_p_prev <- mu_p0
    Sigma_p_prev <- Sigma_p0
    
    # Update parameters
    
    Sigma_eps0 <- update_Sigma_eps(f_coef = f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0)
    
    
    
    H0 <- update_H(f_coef = f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0, D = D, lambda = lambdaH)
    
    mu_p0 <- update_mu_p(f_coef = f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0)
    
    Sigma_p0 <- update_Sigma_p(f_coef = f_coef, Sigma_eps0, H0, mu_p0, Sigma_p0, lambda = lambdaS)
    
    # Compute relative changes
    mu_change <- sqrt(sum((mu_p0-mu_p_prev)^2)) / sqrt(sum((mu_p_prev + 1e-8)^2))
    Sigma_change <- norm(as.matrix(Sigma_p0-Sigma_p_prev), type = "F") / (norm(as.matrix(Sigma_p_prev), type = "F") + 1e-8)
    
    # Print iteration number
    #cat("Iteration:", b, " | mu_change:", round(mu_change, 4), " | Sigma_change:", round(Sigma_change, 4), "\n")
    
    # Break condition: if both changes < 5%
    if (mu_change < 0.05 && Sigma_change < 0.05) {
      cat("Converged after", b, "iterations (changes < 5%)\n")
      break
    }
    
    # Optional visualization (only if needed)
    #if (b%%5==1) {
     # x11()
    #  par(mfrow = c(4, 2))
    #  image(t(loadings) %*% Sigma_eps %*% loadings)
    #  image(Sigma_eps0)
    #  image(Sigma_p)
    #  image(Sigma_p0)
    #  image(as.matrix(mu_p))
    #  image(as.matrix(mu_p0))
    #  verteces_display(H, gb)
    #  verteces_display(H0, pb)
    #}
    #print(Sys.time()-t1)
  }
  
  return(list(mu_p = mu_p0, Sigma_p = Sigma_p0, H = H0, Sigma_eps = Sigma_eps0))
}
