library(raster)
library(rasterVis)
library(xtable)
library(plyr)
library(reshape2)
library(fields)
library(spam64)

source("glb_specs.r")
source("lambda_helpers.r")

# Load in data and MCMC output =================================================

X <- readRDS("../data/lambda_covars.rds")
Boundaries <- readRDS("../data/Boundaries.rds")

Boundary <- Boundaries$Boundary
BoundaryNA <- Boundaries$BoundaryNA
Boundary.us <- Boundaries$Boundary.us
ind <- which(Boundary.us[] == 1)

# Extract parameters

MCMC.logistic.chains <- readRDS("../data/lambda_mcmc.rds")
status <- sum(!is.na(MCMC.logistic.chains[[3]]))
burn <- floor(status/2)
thin <- 1/10
ind.est <- seq(burn, status, by=1/thin)
nmcmc <- length(ind.est)

gamma.est <- MCMC.logistic.chains$gamma[ind.est]
beta.est <- MCMC.logistic.chains$beta[ind.est, ]
theta.est <- MCMC.logistic.chains$theta[ind.est]
kappa.est <- MCMC.logistic.chains$kappa[ind.est]
tau.est <- MCMC.logistic.chains$tau[ind.est]
K.est <- MCMC.logistic.chains$K[ind.est]
p.est <- MCMC.logistic.chains$p[ind.est, ]

# Compute posterior mean diffusion coefficeint

delta <- X %*% t(beta.est) %>% rowMeans()

saveRDS(delta, "../output/delta.rds")

# Loop over MCMC output and compute lambda =====================================

lambda.mean <- brick(nrows = y, 
                     ncols = x, 
                     nl = length(time.frame),
                     xmn = xmin, 
                     xmx = xmax, 
                     ymn = ymin, 
                     ymx = ymax, 
                     crs=NA)

values(lambda.mean) <- 0.0

for(i in 1:nmcmc){
  
  print(i)
  
  gamma <- gamma.est[i]
  beta <- beta.est[i, ]
  theta <- theta.est[i]
  kappa <- kappa.est[i]
  tau <- tau.est[i]
  K <- K.est[i]
  
  # Diffusion rate for pde
  delta.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  delta.r[] <- exp(X %*% beta)
  
  # Growth rate for pde
  gamma.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  gamma.r[] <- rep(gamma, q)
  
  # Diffusion rate for homogenized pde
  delta.r.inv <- 1/delta.r
  delta.r.inv <- focal(delta.r.inv, 
                       w=matrix(rep(1/smooth.fact^2, smooth.fact^2), nrow=smooth.fact), 
                       na.rm=T, pad=T, padValue=0)
  delta.bar <- aggregate(1/delta.r.inv, fact=us.fact, na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2]) 
  delta.vec <- delta.bar[ind]
  
  # Growth rate for homogenized pde
  gamma.bar <- aggregate(gamma.r, fact=us.fact, na.rm=TRUE, fun=mean)
  gamma.vec <- gamma.bar[ind]
  
  # Carrying capacity for homogenized pde
  delta.r.inv2 <- 1/delta.r^2
  delta.r.inv2 <- focal(delta.r.inv2, 
                        w=matrix(rep(1/smooth.fact^2, smooth.fact^2), nrow=smooth.fact, ncol=smooth.fact), 
                        na.rm=T, pad=T, padValue=0)
  k.bar <- aggregate(K*delta.r.inv/delta.r.inv2, fact=us.fact, na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2])
  k.vec <- k.bar[ind]
  
  # First-order neighborhood matrix
  NN <- neighborhood(Boundary.us)
  
  # Propagator matrix
  H <- propagator(NN, Boundary.us, delta.vec, dx, dy, dt)
  
  # Initial condition (D is in km)
  D <- rdist(data.frame(SpatialPoints(cell)), matrix(d, 1, 2))/1000
  lambda0.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  lambda0.r[] <- theta*exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))
  c0.r <- raster(nrows=y/us.fact, ncols=x/us.fact, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  c0.r[] <- extract(delta.r*lambda0.r, SpatialPoints(Boundary.us))
  c0.vec <- c0.r[ind]

  # Calculate c.all 
  c.all <- calcc(H, c0.vec, gamma.vec, k.vec, time.steps, dt)[, keep]
  c.m <- matrix(NA, nrow=q/us.fact^2, ncol=length(time.frame))
  c.m[ind, ] <- c.all
  c.r <- brick(nrows=y/us.fact, ncols=x/us.fact, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  c.r <- setValues(c.r, c.m)
  
  # Calculate lambda.all from c.all
  lambda.all <- disaggregate(c.r, us.fact)/delta.r
  lambda.all <- lambda.all*BoundaryNA
  
  lambda.mean <- lambda.mean + 1 / nmcmc * lambda.all
}

writeRaster(lambda.mean, "../output/lambda_mean.grd", bandorder = "BIL")
