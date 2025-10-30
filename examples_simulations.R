library(FDboost)

source("BHM_simulator.R")
source("Sampling_p.R")
source("clr2density.R")
source("coef_function.R")
source("update_functions.R")
source("EM_updating.R")


############## examples #################
#set.seed(03071)

# Parameters to vary:
# Simulation A
# m=2,3,5,8,13,21 - Fibonacci's series.
# Simulation B
# sd=0.01,0.02,0.05,0.1 x norm of the mean density.
# Simulation C
# Deviation from the ilr_inv(0,0,0,0,0)

# Qualitative assessment: 
# Show an example.

# Quantitative assessment:
# Run 100 simulations for each scenario
# One sim. is around 30 min 
# 100*15*30= 1500 min = 750 ore

set.seed(03072)
m=2
n=150



F_sample=F_simulator(m=m,n=n,sd_perc=0.01)

Sigma_eps=F_sample$Sigma_eps
Sigma_p=F_sample$pF$p$sigma
mu_p=F_sample$pF$p$mu
H=F_sample$pF$H$coefs
f_coef=F_sample$nF #coefficients
D=F_sample$pF$H$D
D0=F_sample$pF$H$D0

get_basis=F_sample$pF$H$basis

x11()
par(mfrow=c(1,3))
verteces_display(H%*%(F_sample$pF$p$p),get_basis, main="Density mixtures")
verteces_display(f_coef,get_basis, main="Density mixtures with noise")
verteces_display(H,get_basis, main="Vertices")


##################################################################

#id=10
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

## Non centred PCA! ##
# If you have scarsity in data to build the densities, use sparse-principal component

D=F_sample$pF$H$D
D0=F_sample$pF$H$D0

m
k=23
sd0=10^-6

Cmat <- cov(t(f_coef))
A <- D0 + D * 1e4  # try smaller penalty
M <- D0 %*% Cmat %*% D0

# Solve generalized eigenproblem (A^-1 M)
Mat_fpca <- eigen(solve(A, M), symmetric = TRUE)

#D0[i,j]=integral(elem[i]*elem[j])
loadings=Mat_fpca$vectors

pca_basis=get_basis%*%loadings


x11()
matplot(1:100,pca_basis,type="l")

for(j in 1:k){
  norm_const=sqrt(sum(pca_basis[,j]^2)*w)
  pca_basis[,j]=pca_basis[,j]/norm_const
  loadings[,j]=loadings[,j]/norm_const
  #for(i in 1:n){
  #  f_coef_pca[j,i]=w*t(pca_basis[,j])%*%(get_basis)%*%f_coef[,i] #(1x100)x(100x23)x(23x1)
  #}
}

f_coef_pca=solve(loadings)%*%f_coef

x11()
matplot(1:100,pca_basis,type="l")


#y=a*1*x/a=x



D_pca=t(loadings)%*%D%*%(loadings)
D0_pca=t(loadings)%*%D0%*%(loadings)


x11()
par(mfrow=c(1,2))
verteces_display(f_coef,get_basis)
verteces_display(f_coef_pca,pca_basis)

H0=diag(k)[,1:m]

sd0=0.01
EM_sample=EM_updating(f_coef_pca, diag(sd0,k,k), H0 , rep(0,m-1), diag(1,m-1,m-1), 
                      B=20,D=D_pca,lambdaS=1,lambdaH=0)



x11()
par(mfrow=c(1,2))
verteces_display((H), get_basis,ylim=c(0,3))
verteces_display(EM_sample$H, pca_basis,ylim=c(0,3))



sample_p=mvrnorm(150,EM_sample$mu_p,EM_sample$Sigma_p)
sample_err=mvrnorm(150,rep(0,k),EM_sample$Sigma_eps)

sample_p=apply(sample_p,1,ilrInv)
EM_H=EM_sample$H


sample_from_EM=t(sample_p)%*%(t(EM_sample$H))+sample_err

x11()
par(mfrow=c(2,2))
verteces_display(H, get_basis, ylim=c(0,4), main="True vertices")
verteces_display(F_sample$pF$pF,get_basis,ylim=c(0,4),main="Dataset")
verteces_display(EM_H,pca_basis, ylim=c(0,4), main="Estimated vertices")
verteces_display(t(sample_from_EM),pca_basis,ylim=c(0,4), main="Simulated new dataset")
mu_p
EM_sample$mu_p



Sigma_p
EM_sample$Sigma_p


x11()
par(mfrow=c(1,2))
image(t(solve(loadings))%*%Sigma_eps%*%solve(loadings))
image(EM_sample$Sigma_eps)


range(EM_sample$Sigma_eps)
range(Sigma_eps)


library(ggplot2)
library(reshape2)
library(gridExtra)

# Example matrices (replace with your own)
# mat1 <- t(solve(loadings)) %*% Sigma_eps %*% solve(loadings)
# mat2 <- EM_sample$Sigma_eps

# Compute global min/max for shared color scale
vmin <- min(mat1, mat2)
vmax <- max(mat1, mat2)

# Helper function to plot heatmap
plot_heatmap <- function(mat, title) {
  df <- melt(mat)
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = heat.colors(100),
      limits = c(vmin, vmax)  # fixed scale across both plots
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(title)
}

p1 <- plot_heatmap(mat1, "Reconstructed Σₑ")
p2 <- plot_heatmap(mat2, "Estimated Σₑ")

# Combine with a single shared legend
gridExtra::grid.arrange(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  ggplotGrob(p1)$grobs[[which(sapply(ggplotGrob(p1)$grobs, function(x) x$name) == "guide-box")]],
  ncol = 3,
  widths = c(1, 1, 0.3)
)



