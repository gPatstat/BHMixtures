
min_perm_H1 <- function(H1, H2) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Install package 'clue': install.packages('clue')")
  }
  
  # Controllo dimensioni
  if (!all(dim(H1) == dim(H2))) {
    stop("H1 and H2 have to be of the same dimensions")
  }
  
  m <- ncol(H1)
  
  # Cost matrix: suared distance among columns
  cost <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      cost[i, j] <- sum((H1[, i] - H2[, j])^2)
    }
  }
  
  cost=t(cost)
  
  # Solve linear sum assignment problem
  perm <- clue::solve_LSAP(cost)
  
  # Best permutation
  return(as.vector(perm))
}

vertex_reorder=function(H_est,H_true){
  perm <- min_perm_H1(H_est, H_true)
  H_est_perm <- H_est[, perm]
  return(H_est_perm)
}

error_vertex=function(H_est,H_true){
  H_est=vertex_reorder(H_est,H_true)
  fr_err=sqrt(sum(w*(H_est-H_true)^2/dim(H_est)[2]))
  return(fr_err)
}

