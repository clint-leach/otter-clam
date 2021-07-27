###
### Subroutines and Packages
###

required.packages <- c("coda",
                       "fBasics",
                       "fields",
                       "ggmap",
                       "ggplot2",
                       "gridExtra",
                       "gstat",
                       "inline",
                       "maptools",
                       "raster",
                       "rasterVis",
                       "RColorBrewer",
                       "RcppArmadillo",
                       "rgdal",
                       "rgeos",
                       "msm",
                       "doParallel",
                       "gtools")
lapply(required.packages, library, character.only=TRUE)

###
### Spatio-temporal settings
###

dt <- 1/400 
time.frame <- seq(1993, 2018) 
time.steps <- 1/dt*length(time.frame)
keep <- seq(1, time.steps, time.steps/length(time.frame))
us.fact <- 10
smooth.fact <- 15
res <- 400
dx <- res*us.fact
dy <- res*us.fact
d <- c(445000,6470000) # Perry's original value
xmin <- 400000
xmax <- 460000
ymin <- 6468000
ymax <- 6540000
x <- 150
y <- 180
q <- x*y
cell <- raster(nrows=y, ncols=x, xmn=xmin, xmx=xmax,
               ymn=ymin, ymx=ymax, crs=NA)
cell[] <- 1:q
extent <- c(xmin, xmax, ymin, ymax)
st.info <- list(dt=dt,
                time.frame=time.frame,
                us.fact=us.fact,
                smooth.fact=smooth.fact,
                res=res,
                d=d,
                extent=extent)

###
### Simulation settings
###

theta <- 500
kappa <- 60
beta0 <- 18
beta1 <- -1.5
beta2 <- 0.8
beta3 <- -0.3
beta4 <- 1
beta <- c(beta0, beta1, beta2, beta3, beta4)
gamma <- 0.25
p <- 0.75
tau <- 0.5
K <- 5

###
### Homogenization settings
###

us.fact <- 10

###
### Data collection settings
###

n.tr <- 50
n.isu <- 100
years <- c(seq(1999, 2004), 2006, 2012) # ISU data
data.years <- c(1993, seq(1996, 2006), 2009, 2010, 2012, 2017, 2018) # any data

###
### MCMC settings
###

n.iter <- 15000
checkpoint <- 1000
n.cluster <- 8


# Priors
mu.beta <- 0
var.beta <- 10^2

lo.gamma <- 0
hi.gamma <- 0.5

mu.theta <- 100
var.theta <- 200^2

mu.kappa <- 10
var.kappa <- 100^2

a.p <- 30
b.p <- 10

lo.tau <- 0
hi.tau <- 1

lo.K <- 0
hi.K <- 10^2

priors <- list(beta.prior=c(mu.beta, var.beta),
               gamma.prior=c(lo.gamma, hi.gamma),
               theta.prior=c(mu.theta, var.theta),
               kappa.prior=c(mu.kappa, var.kappa),
               p.prior=c(a.p, b.p),
               tau.prior=c(lo.tau, hi.tau),
               K.prior=c(lo.K, hi.K))

# Tuning parameters
gamma.tune <- 0.005
beta.tune <- c(0.05, 0.01, 0.01, 0.01, 0.01)
theta.tune <- 50
kappa.tune <- 0.1
tau.tune <- 0.01
K.tune <- 0.5

tuning <- list(gamma=gamma.tune,
               beta=beta.tune,
               theta=theta.tune,
               kappa=kappa.tune,
               tau=tau.tune,
               K=K.tune)

# Starting values
gamma.init <- 0.25
beta.init <- c(17, -1, 1, -1, 1)
theta.init <- 500
kappa.init <- 60
p.init <- rep(0.75, length(time.frame))
tau.init <- 0.5
K.init <- 5

inits <- list(gamma=gamma.init,
              beta=beta.init,
              theta=theta.init,
              kappa=kappa.init,
              tau=tau.init,
              p=p.init,
              K=K.init)

# Parameters to monitor
parameters <- c("gamma",
                "beta",
                "kappa",
                "theta",
                "tau",
                "p",
                "K",
                "n.tot",
                "lambda.tot",
                "score",
                "tuners",
                "N")

