#######################################
######## BHMixtures simulator #########
#######################################


H_simulator <- function(m, k, mu, sigma) {
  library(mvtnorm)
  
  # Check if length of mu is k
  if (length(mu) != k) {
    stop(paste0("Error: length of mu must be ", k, ", but got ", length(mu), "."))
  }
  
  # Check if sigma is a k x k matrix
  if (!all(dim(sigma) == c(k, k))) {
    stop(paste0("Error: sigma must be a ", k, " x ", k, " matrix."))
  }
  
  # Check if sigma is positive definite
  if (!all(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values > 0)) {
    stop("Error: sigma must be positive definite.")
  }
  
  # Generate matrix of vertices
  mat <- matrix(NA, k, m)
  for (j in 1:m) {
    mat[, j] <- rmvnorm(1, mu, sigma)
  }
  
  return(mat)
}
