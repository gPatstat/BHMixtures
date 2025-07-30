#####################################################
############### Importance sampling #################
#####################################################


# Sigma_eps: covariance matrix of the discretized residuals [k x k]
# H: matrix of vertices' coefficients in clr [k x m]
# f: vector of coefficients of a density in clr [k x 1]
# mu_p:  [(m-1) x 1]
# Sigma_p: [(m-1) x (m-1)]



posterior_builder=function(f,Sigma_eps, H, mu_p, Sigma_p){
  options(digits = 22)  # aumenta cifre stampate
  inv_Sigma_eps=solve(Sigma_eps)
  inv_Sigma_p=solve(Sigma_p)
  f=as.matrix(f)
  func=function(x){
    library(compositions)
    p_x=as.numeric(ilrInv(x))
    val=-1/2*(t(x-mu_p)%*%inv_Sigma_p%*%(x-mu_p)+t(f-H%*%p_x)%*%inv_Sigma_eps%*%(f-H%*%p_x))
    return(-as.numeric(val))
  }
  return(func)
}


#install.packages("numDeriv")  # if not installed


importance_sampling_p=function(f0,Sigma_eps, H, mu_p, Sigma_p,B=100){
  library(MASS)
  library(numDeriv)
  post_prob=posterior_builder(f0,Sigma_eps, H, mu_p, Sigma_p)
  mu_opt=optim(par=(rep(0,length(mu_p))),fn=post_prob)$par
  Hess=numDeriv::hessian(post_prob,mu_opt)
  w=numeric(B)
  x=matrix(NA,B,length(mu_p))
  for(b in 1:B){
    x[b,]=mvrnorm(1,mu_opt,solve(Hess))
    w[b]=-post_prob(x[b,])+1/2*t(x[b,]-mu_opt)%*%Hess%*%(x[b,]-mu_opt)
  }
return(list(proportions=x,w=exp(w)))
}


F_sample=F_simulator(m=3,n=150,sd=0.05)

Sigma_eps=F_sample$Sigma_eps
Sigma_p=F_sample$pF$p$sigma
mu_p=F_sample$pF$p$mu
H=F_sample$pF$H$coefs
f=F_sample$nF

##################################################################


dclrs=list(x=t(f[,1]),t=domain)

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
imps$w

plot3d(t(F_sample$pF$p$p),xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))
plot3d(t(apply(imps$proportions,1,ilrInv)),add=T,col="red")
plot3d(t(F_sample$pF$p$p)[1,],add=T,col="blue")


x11()
lines(domain,f[,1])


#################################################################
library(plotly)

# 1. All samples (converted to 3D matrix)
pts_main <- t(F_sample$pF$p$p)  # Matrix N × 3

# 2. ilrInv of imps$proportions (assume it's D-1 dimensional ilr coords)
pts_red <- t(apply(imps$proportions, 1, ilrInv))

# 3. Special point (just one row of pts_main)
pt_blue <- pts_main[1, , drop = FALSE]


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
  name = "ilrInv(imps)"
)

# Add blue point
fig <- fig %>% add_trace(
  x = pt_blue[,1], y = pt_blue[,2], z = pt_blue[,3],
  type = "scatter3d", mode = "markers",
  marker = list(color = "blue", size = 6, symbol = "diamond"),
  name = "Reference Point"
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

