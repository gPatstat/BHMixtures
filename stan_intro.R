

# the stan file is "model.stan"
# input is given iteration by iteration:

# int<lower=2> m;                 // number of vertices
# vector[m] f0;                   // n-th clr function
# matrix[m, m] H;                 // current vertices matrix
# vector[m-1] mu_p;               // current mu_p
# matrix[m-1, m-1] Sigma_p;       // current Sigma_p
# matrix[m, m] Sigma_eps;         // current Sigma_eps
# matrix[m, m-1] V;               // base ilr

#input_list=list(m=m, f0=f0, H=H0, mu_p=mu_p0, Sigma_p=Sigma_p0,Sigma_eps=Sigma_eps0,V=V)

# usage: fit <- stan(file = 'model.stan', data = input_list)





#install.packages("cmdstanr",
#              repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
#)
library(cmdstanr)
#cmdstan_version()
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#cmdstanr::install_cmdstan()
#install.packages("C:/Users/test/Downloads/cmdstan-2.28.2.tar.gz", repos = NULL, type = "source")

source("update_functions_new.R")
mod <- cmdstan_model(stan_file = "model_stan_new.stan")
#fit <- mod$sample(data = input_list)


