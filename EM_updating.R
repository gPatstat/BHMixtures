


EM_updating=function(f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0,B=100,D){
  for(b in 1:B){
    Sigma_eps0=update_Sigma_eps(f_coef=f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    H0=update_H(f_coef=f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0,D=D)
    mu_p0=update_mu_p(f_coef=f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    Sigma_p0=update_Sigma_p(f_coef=f_coef,Sigma_eps0, H0, mu_p0, Sigma_p0)
    print(b)
    x11()
    par(mfrow=c(4,2))
    image(Sigma_eps)
    image(Sigma_eps0)
    image(Sigma_p)
    image(Sigma_p0)
    image(as.matrix(mu_p))
    image(as.matrix(mu_p0))
    verteces_display(H, get_basis)
    verteces_display(H0, get_basis)
    
  }
  return(list(mu_p=mu_p0,Sigma_p=Sigma_p0,H=H0,Sigma_eps=Sigma_eps0))
}
