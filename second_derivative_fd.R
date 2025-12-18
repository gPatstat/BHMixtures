second_derivative_fd <- function(fx) {
  x=domain
  dx <- diff(x)
  if (max(dx) - min(dx) > 1e-8)
    stop("x is not evenly spaced. Use the spline method instead.")
  
  dx <- dx[1]
  n  <- length(fx)
  
  d2 <- numeric(n)
  d2[2:(n-1)] <- (fx[3:n] - 2 * fx[2:(n-1)] + fx[1:(n-2)]) / dx^2
  d2[c(1, n)] <- NA  # undefined at boundaries
  return(d2)
}
