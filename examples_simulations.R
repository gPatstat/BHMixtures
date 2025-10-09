library(FDboost)

source("BHM_simultor.R")
source("Sampling_p.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions.R")
source("EM_updating.R")


############## examples #################

m=3
n=500

F_sample=F_simulator(m=m,n=n,sd=0.01)

Sigma_eps=F_sample$Sigma_eps
Sigma_p=F_sample$pF$p$sigma
mu_p=F_sample$pF$p$mu
H=F_sample$pF$H$coefs
f_coef=F_sample$nF #coefficients
D=F_sample$pF$H$D

get_basis=F_sample$pF$H$basis


x11()
par(mfrow=c(1,2))
image(Sigma_eps)
image(cov(t(F_sample$coef_noise)))


x11()
par(mfrow=c(1,2))
verteces_display(H%*%(F_sample$pF$p$p),get_basis)
verteces_display(f_coef,get_basis)


##################################################################


coefs=f_coef[,id]

imps=importance_sampling_p(coefs,Sigma_eps, H, mu_p, Sigma_p,B=100)

wgt_scaled <- imps$wgt 

# Create a color vector with transparency (using rgb)
# Assuming 2D projections for plotting
library(scales)
colors <- rgb(0, 0, 1, alpha = wgt_scaled)  # blue points with varying alpha

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
  marker = list(
    color = wgt_scaled,
    colorscale = list(c(0, 'rgb(255,200,200)'), c(1, 'rgb(255,0,0)')),
    colorbar = list(title = "Weight"),
    size = 4
  ),
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

H_old=(f_coef[,sample(1:(dim(f_coef)[1]),3)])
H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_old,D=D,lambda=100)
H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new,D=D,lambda=100)

get_basis=F_sample$pF$H$basis

x11()
par(mfrow=c(4,2))
verteces_display(H, get_basis)
verteces_display(H_old, get_basis)
verteces_display(H_new, get_basis)

H_new=update_H(f_coef,Sigma_eps,mu_p,Sigma_p,H_old=H_new,D=D,lambda=0)
verteces_display(H_new, get_basis)

#################################################################

k=13
pca=princomp(t(f_coef))
#H0=(f_coef[,sample(1:150,3)])
H0=pca$loadings[,1:m]
sd0=pca$sdev[m+1]

x11()
verteces_display(H0,get_basis)

EM_sample=EM_updating(f_coef, diag(sd0,k,k), H0 , rep(0,m-1), diag(1,m-1,m-1), B=20,D=D,lambda=100)

x11()
par(mfrow=c(1,2))
verteces_display((H), get_basis)
verteces_display(EM_sample$H, get_basis)


sample_p=mvrnorm(150,EM_sample$mu_p,EM_sample$Sigma_p)
sample_err=mvrnorm(150,rep(0,k),EM_sample$Sigma_eps)

sample_p=apply(sample_p,1,ilrInv)
EM_H=EM_sample$H


sample_from_EM=t(sample_p)%*%(t(EM_sample$H))


x11()
par(mfrow=c(1,2))
verteces_display(F_sample$pF$pF,get_basis)
verteces_display(t(sample_from_EM),get_basis)
mu_p
EM_sample$mu_p



Sigma_p
EM_sample$Sigma_p
