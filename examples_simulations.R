library(FDboost)

source("BHM_simultor.R")
source("Sampling_p.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions.R")


############## examples #################

F_sample=F_simulator(m=3,n=150,sd=0.01)

Sigma_eps=F_sample$Sigma_eps
Sigma_p=F_sample$pF$p$sigma
mu_p=F_sample$pF$p$mu
H=F_sample$pF$H$coefs
f=F_sample$nF


##################################################################

id=1
dclrs=list(x=t(f[,id]),t=domain)

cclrs=FDboost(x~1, 
              #use bbsc() in time formula to ensure integrate-to-zero constraint 
              timeformula= ~bbsc(t,df=4, 
                                 knots=10,
                                 boundary.knots=c(0,1),
                                 degree=3,
                                 lambda=10^2),
              data=dclrs,offset=0, 
              control=boost_control(mstop=100))

coefs=as.numeric(cclrs$coef(which=1)[[1]])

##################################################################


imps=importance_sampling_p(coefs,Sigma_eps, H, mu_p, Sigma_p,B=100)
imps$wgt


library(plotly)

# 1. All samples (converted to 3D matrix)
pts_main <- t(F_sample$pF$p$p)  # Matrix N × 3

# 2. ilrInv of imps$proportions (assume it's D-1 dimensional ilr coords)
pts_red <- t(apply(imps$proportions, 1, ilrInv))

# 3. Special point (just one row of pts_main)
pt_blue <- pts_main[id, , drop = FALSE]


fig <- plot_ly()

# Add main sample points
fig <- fig %>% add_trace(
  x = pts_main[,1], y = pts_main[,2], z = pts_main[,3],
  type = "scatter3d", mode = "markers",
  marker = list(color = "black", size = 3),
  name = "Samples"
)

# Add red points
fig <- fig %>% add_trace(
  x = pts_red[,1], y = pts_red[,2], z = pts_red[,3],
  type = "scatter3d", mode = "markers",
  marker = list(color = "red", size = 4),
  name = "Posterior samples"
)

# Add blue point
fig <- fig %>% add_trace(
  x = pt_blue[,1], y = pt_blue[,2], z = pt_blue[,3],
  type = "scatter3d", mode = "markers",
  marker = list(color = "blue", size = 6, symbol = "diamond"),
  name = "True proportion"
)

# Set axis limits
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "X", range = c(0, 1)),
    yaxis = list(title = "Y", range = c(0, 1)),
    zaxis = list(title = "Z", range = c(0, 1))
  )
)

fig


###############################################################################

f_coef=coef_f(f)

H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=t(f_coef[sample(1:(dim(f_coef)[1]),3),]))

get_basis=F_sample$pF$H$basis

x11()
par(mfrow=c(4,2))
verteces_display(H, get_basis)
verteces_display(H_old, get_basis)
verteces_display(H_new, get_basis)


H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new)
verteces_display(H_new, get_basis)

x11()
par(mfrow=c(1,2))
H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new)
verteces_display(H_new, get_basis)

H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new)
verteces_display(H_new, get_basis)

x11()
par(mfrow=c(1,2))
H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new)
verteces_display(H_new, get_basis)
