


EM_updating=function(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0,B=100){
  for(b in 1:B){
    mu_p0=update_mu_p(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    H0=update_H(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    Sigma_p0=update_Sigma_p(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    Sigma_eps0=update_Sigma_p(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
  }
  return(mu_p=mu_p0,Sigma_p=Sigma_p0,H=H0,Sigma_eps=Sigma_eps0)
}