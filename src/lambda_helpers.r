###
### C++ function for propagation
###

code <- '
arma::mat Hmat = Rcpp::as<arma::mat>(H);
arma::vec c0vec = Rcpp::as<arma::vec>(c0);
arma::vec gvec = Rcpp::as<arma::vec>(gamma);
arma::vec kvec = Rcpp::as<arma::vec>(ktilde);
double dt = Rcpp::as<double>(deltat);
int n = Rcpp::as<int>(timesteps);
int k = Hmat.n_rows;
arma::mat call(k, n);
call.col(0) = Hmat*c0vec;
for (int i = 1; i < n; ++i) {
    call.col(i) = Hmat*call.col(i - 1) + dt*gvec%call.col(i - 1) - 
                  dt*gvec%(pow(call.col(i - 1), 2)/kvec);
}
return Rcpp::wrap(call);
'
calcc <- cxxfunction(signature(H="numeric",
                               c0="numeric",
                               gamma="numeric",
                               ktilde="numeric",
                               timesteps="numeric",
                               deltat="numeric"),
                     body=code, 
                     plugin="RcppArmadillo")

###
### Calculate neighborhood matrix
###

neighborhood <- function(boundary) {
  
  ###
  ### Input:
  ### boundary: raster of homogenized boundary, 0 is land, 1 is water
  ###
  ### Output:
  ### NN: matrix containing neighborhood information of cells in water,
  ###     cell number if neighbor is in water, 0 if no-flux boundary (neighbor 
  ###     on land), NA if fatal boundary (neighbor out-of-bound)
  ###
  
  ind <- which(boundary[] == 1) # cells in water
  NN <- matrix(NA, length(ind), 4)
  for (i in 1:length(ind)) {
    cell <- ind[i]
    adj <- adjacent(boundary, cell)[, 2] # extract adjacent cells
    
    ln <- adj[which((adj + 1) == cell)] # left neighbor
    if (length(ln) > 0) {
      if (!(ln %in% ind)) {
        ln <- 0
      }
    } else {
      # ln <- NA
      ln <- 0 # change all BC to no-flux 
    }
    
    
    rn <- adj[which((adj - 1) == cell)] # right neighbor
    if (length(rn) > 0) {
      if (!(rn %in% ind)) {
        rn <- 0
      }
    } else {
      # rn <- NA
      rn <- 0 # change all BC to no-flux
    }
    
    bn <- adj[which((adj - dim(boundary)[2]) == cell)] # bottom neighbor
    if (length(bn) > 0) {
      if (!(bn %in% ind)) {
        bn <- 0
      }
    } else {
      # bn <- NA
      bn <- 0 # change all BC to no-flux
    }
    
    tn <- adj[which((adj + dim(boundary)[2]) == cell)] # top neighbor
    if (length(tn) > 0) {
      if (!(tn %in% ind)) {
        tn <- 0
      }
    } else {
      # tn <- NA
      tn <- 0 # change all BC to no-flux
    }
    
    NN[i, ] <- c(ln, rn, bn, tn)
    
  }
  
  return (NN)
}

###
### Calculate propogator matrix H
###

propagator <- function(NN, boundary, delta, dx, dy, dt) {
  
  ###
  ### Input:
  ### NN: neighborhood matrix
  ### boundary: homogenized boundary
  ### delta: motility vector
  ### dx: large scale resolution in x
  ### dy: large scale resolution in y
  ### dt: slow scale resolution in time
  ###
  ### Output:
  ### H: propagator matrix of intensity (expected abundance) lambda
  ###
  
  ind <- which(boundary[] == 1) # cells in water
  n <- nrow(NN)
  H <- matrix(0, n, n)
  for (i in 1:n) {
    is.water <- ifelse(NN[i, ] > 0, 1, 0)
    H[i, i] <- sum(c(1 - 2*delta[i]*(dt/dx^2 + dt/dy^2),
                     dt/dx^2*delta[i]*(1 - is.water[1:2]),
                     dt/dy^2*delta[i]*(1 - is.water[3:4])),
                   na.rm=TRUE)
    H[i, which(ind == NN[i, 1])] <- dt/dx^2*delta[i]*is.water[1]
    H[i, which(ind == NN[i, 2])] <- dt/dx^2*delta[i]*is.water[2]
    H[i, which(ind == NN[i, 3])] <- dt/dy^2*delta[i]*is.water[3]
    H[i, which(ind == NN[i, 4])] <- dt/dy^2*delta[i]*is.water[4]
  }
  
  return (H)
}

###
### C++ function for propagation (exponential)
###

code <- '
arma::mat Hmat = Rcpp::as<arma::mat>(H);
arma::vec c0vec = Rcpp::as<arma::vec>(c0);
arma::vec gvec = Rcpp::as<arma::vec>(gamma);
double dt = Rcpp::as<double>(deltat);
int n = Rcpp::as<int>(timesteps);
int k = Hmat.n_rows;
arma::mat call(k, n);
call.col(0) = Hmat*c0vec;
for (int i = 1; i < n; ++i) {
call.col(i) = Hmat*call.col(i - 1) + dt*gvec%call.col(i - 1);
}
return Rcpp::wrap(call);
'
calcc.exp <- cxxfunction(signature(H="numeric",
                                   c0="numeric",
                                   gamma="numeric",
                                   timesteps="numeric",
                                   deltat="numeric"),
                         body=code, 
                         plugin="RcppArmadillo")

###
### Simulate abundance
###

simAbundance <- function(parameters){
  ###
  ### Dependencies
  ###
  
  Boundary <- Boundaries$Boundary
  BoundaryNA <- Boundaries$BoundaryNA
  Boundary.us <- Boundaries$Boundary.us
  ind <- which(Boundary.us[] == 1)
  
  gamma <- parameters$gamma
  beta <- parameters$beta
  theta <- parameters$theta
  kappa <- parameters$kappa
  tau <- parameters$tau
  K <- parameters$K
  
  ###
  ### Diffusion rate for pde
  ###
  
  delta.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  delta.r[] <- exp(X%*%beta)
  # delta.r[] <- inv.logit(X%*%beta)*res^2/dt
  
  ###
  ### Growth rate for pde
  ###
  
  gamma.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  gamma.r[] <- rep(gamma, q)
  
  ###
  ### Diffusion rate for homogenized pde
  ###
  
  delta.r.inv <- 1/delta.r
  delta.r.inv <- focal(delta.r.inv, 
                       w=matrix(rep(1/smooth.fact^2, smooth.fact^2), nrow=smooth.fact), 
                       na.rm=T, pad=T, padValue=0)
  delta.bar <- aggregate(1/delta.r.inv, fact=us.fact, na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2]) 
  delta.vec <- delta.bar[ind]
  
  ###
  ### Growth rate for homogenized pde
  ###
  
  gamma.bar <- aggregate(gamma.r, fact=us.fact, na.rm=TRUE, fun=mean)
  gamma.vec <- gamma.bar[ind]
  
  ###
  ### Carrying capacity for homogenized pde
  ###
  
  delta.r.inv2 <- 1/delta.r^2
  delta.r.inv2 <- focal(delta.r.inv2, 
                        w=matrix(rep(1/smooth.fact^2, smooth.fact^2), nrow=smooth.fact, ncol=smooth.fact), 
                        na.rm=T, pad=T, padValue=0)
  k.bar <- aggregate(K*delta.r.inv/delta.r.inv2, fact=us.fact, na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2])
  k.vec <- k.bar[ind]
  
  ###
  ### First-order neighborhood matrix
  ###
  
  NN <- neighborhood(Boundary.us)
  
  ###
  ### Propagator matrix
  ###
  
  H <- propagator(NN, Boundary.us, delta.vec, dx, dy, dt)
  
  ###
  ### Initial condition (D is in km)
  ###
  
  D <- rdist(data.frame(SpatialPoints(cell)), matrix(d, 1, 2))/1000
  lambda0.r <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  lambda0.r[] <- theta*exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))
  c0.r <- raster(nrows=y/us.fact, ncols=x/us.fact, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  c0.r[] <- extract(delta.r*lambda0.r, SpatialPoints(Boundary.us))
  c0.vec <- c0.r[ind]
  
  ###
  ### Calculate c.all 
  ###
  
  c.all <- calcc(H, c0.vec, gamma.vec, k.vec, time.steps, dt)[, keep]
  c.m <- matrix(NA, nrow=q/us.fact^2, ncol=length(time.frame))
  c.m[ind, ] <- c.all
  c.r <- brick(nrows=y/us.fact, ncols=x/us.fact, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=NA)
  c.r <- setValues(c.r, c.m)
  
  ###
  ### Calculate lambda.all from c.all
  ###
  
  lambda.all <- disaggregate(c.r, us.fact)/delta.r
  lambda.all <- lambda.all*BoundaryNA
  
  ###
  ### Expected abundance
  ###
  
  EA.v <- numeric(length(keep))
  for (i in 1:length(keep)){
    EA.v[i]=sum(lambda.all[[i]][], na.rm=TRUE)
  }
  
  ###
  ### Generate abundance (N) from expected abundance (lambda)
  ###
  
  N <- cell
  N <- stack(mget(rep("N", length(keep))))
  set.seed(1234)
  for (i in 1:length(keep)) {
    N[[i]][]<- suppressWarnings(rnbinom(n=length(lambda.all[[i]][]),
                                        size=tau,
                                        mu=lambda.all[[i]][]))
  }
  n.tot.v <- N[]
  n.tot <- unname(tapply(n.tot.v, (seq_along(n.tot.v) - 1)%/%q, sum, na.rm=TRUE))
  
  ###
  ### Save results
  ###
  
  sim.out <-  list('lambda.eq'=lambda.all[[length(time.frame)]][],
                   'N.eq'=N[[length(time.frame)]][],
                   'lambda.tot'=EA.v,
                   'N.tot'=n.tot)
  
  return(sim.out)
  
}
