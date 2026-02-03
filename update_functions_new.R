
#####################################################
############ Coordinate-wise updating ###############
#####################################################

#install.packages("progress")


#####################################################
############ Coordinate-wise updating ###############
#####################################################


#install.packages("progress")
safe_stan_draws <- function(input_list, mod, m) {
  tryCatch({
    fit <- mod$sample(
      data = input_list,
      refresh = 0,
      show_messages = FALSE,
      chains = 1,
      iter_warmup = 1000,
      iter_sampling = 1000,
      thin = 10
    )
    
    # ---- Convergence check: divergent transitions ----
    
    n_div <- fit$diagnostic_summary()$num_divergent
    
    if (n_div > 0) {
      message("Divergent transitions detected: ", n_div)
      return(NULL)
    }
    # --------------------------------------------------
    sampler_df <- as.data.frame(fit$sampler_diagnostics(format = "df"))
    X_df <- suppressWarnings(as.matrix(fit$draws("x",format="df")))
    P_df <- suppressWarnings(as.matrix(fit$draws("p",format="df")))
    
    n_div_draws <- sum(sampler_df$divergent__ == 1)
    
    if (n_div_draws > 0) {
      message("Removed divergent draws: ", n_div_draws)
    }
    
    keep <- sampler_df$divergent__ == 0
    
    X_df <- X_df[keep, ]
    P_df <- P_df[keep, ]
    
    # discard log post column (first)
    X <- data.frame(X_df[,-c(dim(X_df)[2]-0:2)])
    P <- data.frame(P_df[,-c(dim(P_df)[2]-0:2)])
    
    # ---- Remove non-finite rows of X ----
    bad_rows <- rowSums(!is.finite(as.matrix(X))) > 0
    n_bad <- sum(bad_rows)
    
    if (n_bad > 0)
      message("Draws con valori non finiti rimossi: ", n_bad)
    
    X <- X[!bad_rows, , drop = FALSE]
    P <- P[!bad_rows, , drop = FALSE]
    
    # If everything is non-finite
    if (n_bad/(dim(X_df)[1]) > 0.05)
      stop("More than 5% of the samples is non-finite")
    
    return(list(X,P))
  })
}



# f_coef here is the matrix of coefficients [n x k]
update_H <- function(f_coef, Sigma_eps, H_old, mu_p, Sigma_p,
                     lambda = 0, D, V) {
  n = dim(f_coef)[2]
  k = dim(f_coef)[1]
  m = length(mu_p) + 1
  
  #pb <- progress_bar$new(
  #  format = "  updating Sigma_eps [:bar] :percent datapoints in :elapsedfull",
  #  total = n, clear = FALSE, width= 60)
  
  C   = matrix(0, k, m)
  A_r = matrix(0, m, m)
  A_l = 2 * lambda * (Sigma_eps %*% D)
  
  L_p   = t(chol(Sigma_p))
  L_eps = t(chol(Sigma_eps))
  

  input_list = list(
    m = m, k = k, n=n, f0 = t(f_coef),
    H = H_old, mu_p = as.vector(mu_p),
    L_p = L_p, L_eps = L_eps, V = V
  )
  
  draws <- safe_stan_draws(input_list, mod, m)
  X=draws[[1]]
  P=draws[[2]]
  
  if (is.null(X)) return(NULL)
  
  B=dim(X)[1]
  
  for(id in 1:n) {
    #pb$tick()
    for (b in 1:B) {
      p_x=as.numeric(P[b,(id+(n)*((1:m)-1))])
      C=C+f_coef[,id]%*%t(p_x) # [k x 1]%*%[1 x m]=[k x m] 
      A_r=A_r+p_x%*%t(p_x) # [m x 1]%*%[1 x m]=[m x m]
    }
  }
    
  # vectorize C
  vec_C=matrix(C,k*m,1,byrow=F) 
  # prepare denominator
  den=+diag(m)%x%A_l+t(A_r)%x%diag(k) 
  
  # get vectorized solution and "matricize" it
  vec_solution=solve(den)%*%vec_C 
  solution=matrix(vec_solution,k,m,byrow=F)
  return(solution)
}

#update_H(f_coef_pca,Sigma_eps0, H0, mu_p0, Sigma_p0, 0,D2)

update_mu_p <- function(f_coef, Sigma_eps, H, mu_p, Sigma_p, V) {
  
  n = dim(f_coef)[2]
  k = dim(f_coef)[1]
  m = length(mu_p) + 1
  mu0 = matrix(0, m - 1, 1)
  
  #pb <- progress_bar$new(
  #  format = "  updating Sigma_eps [:bar] :percent datapoints in :elapsedfull",
  #  total = n, clear = FALSE, width= 60)
  
  L_p   = t(chol(Sigma_p))
  L_eps = t(chol(Sigma_eps))

  input_list = list(
    m = m, k = k, n=n, f0 = t(f_coef),
    H = H, mu_p = as.vector(mu_p),
    L_p = L_p, L_eps = L_eps, V = V
  )
  
  draws <- safe_stan_draws(input_list, mod, m)
  X=draws[[1]]
  P=draws[[2]]
  
  if (is.null(X)) return(NULL)
  
  B=dim(X)[1]
  
  for(id in 1:n) {
    #pb$tick()
    mu0 = mu0 + colSums(as.matrix(X[,id+(n)*((1:(m-1))-1)]))
  }
  
  return(mu0/(B*n))
}

#update_mu_p(f_coef_pca,Sigma_eps, H, mu_p0, Sigma_p, V)

update_Sigma_p <- function(f_coef, Sigma_eps, H, mu_p,
                           Sigma_p_old, lambda = 0, V) {
  n = dim(f_coef)[2]
  k = dim(f_coef)[1]
  m = length(mu_p) + 1
  
  #pb <- progress_bar$new(
  #  format = "  updating Sigma_eps [:bar] :percent datapoints in :elapsedfull",
  #  total = n, clear = FALSE, width= 60)
  
  sigma0 = matrix(0, m - 1, m - 1)
  
  L_p   = t(chol(Sigma_p_old))
  L_eps = t(chol(Sigma_eps))
  
  input_list = list(
    m = m, k = k, n=n, f0 = t(f_coef),
    H = H, mu_p = as.vector(mu_p),
    L_p = L_p, L_eps = L_eps, V = V
  )
  
  draws <- safe_stan_draws(input_list, mod, m)
  X=draws[[1]]
  P=draws[[2]]
  
  if (is.null(X)) return(NULL)
  
  B=dim(X)[1]
  for(id in 1:n) {
    #pb$tick()
    for (b in 1:B) {
      d = as.numeric(X[b,id+(n)*((1:(m-1))-1)]) - mu_p
      sigma0 = sigma0 + d %*% t(d)
    }
  }
  
  return(sigma0/(B*n)+lambda*diag(1,(m-1)))
}



update_Sigma_eps <- function(f_coef, Sigma_eps,
                             H, mu_p, Sigma_p, V) {
  n = dim(f_coef)[2]
  k = dim(f_coef)[1]
  m = length(mu_p) + 1
  
  #pb <- progress_bar$new(
  #  format = "  updating Sigma_eps [:bar] :percent datapoints in :elapsedfull",
  #  total = n, clear = FALSE, width= 60)
  
  sigma0 = matrix(0, k, k)
  
  L_p   = t(chol(Sigma_p))
  L_eps = t(chol(Sigma_eps))
  
  input_list = list(
      m = m, k = k, n=n, f0 = t(f_coef),
      H = H, mu_p = as.vector(mu_p),
      L_p = L_p, L_eps = L_eps, V = V
    )
    
    draws <- safe_stan_draws(input_list, mod, m)
    X=draws[[1]]
    P=draws[[2]]
    
    if (is.null(X)) return(NULL)
    
    B=dim(X)[1]
    
    for (id in 1:n) {
      #pb$tick()
      for(b in 1:B){
      p_x=as.numeric(P[b,id+(n)*((1:m)-1)])
      r   = f_coef[,id] - H %*% p_x
      sigma0 = sigma0 + r %*% t(r)
      }
    }
    
    #message("update_Sigma_eps: i = ", i)
    nsamp=n*B

  return(sigma0 /nsamp)
}

