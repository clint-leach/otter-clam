library(raster)
library(xtable)
library(plyr)
library(reshape2)
library(fields)
library(spam64)

# Loading in infauna sampling data =============================================

# Intertidal data
data <- read.csv("../data/intertidal_prey.csv") %>% 
  subset(!site  %in% c("FrancisE", "GeikieE", "DrakeE", "LesterE", "NetlandE", 
                       "StrawberryE", "WhidbeyE", "JohnsonE"))

# Getting which years each site was sampled (to fill in zeros)
samples <- ddply(data, .(site), summarise,
                 years = c(1993:2018),
                 sampled = years %in% Year) %>% 
  reshape2::acast(years ~ site, value.var = "sampled")

# Creating a summary table of all the sites and their locations
sites <- ddply(data, .(site), summarise, 
               Longitude = Longitude[1], 
               Latitude = Latitude[1],
               subtidal = type[1] == "subtidal")

# Making spatial points object from site list 
crdref <- CRS('+proj=longlat +datum=NAD83')
sitepts <- SpatialPoints(dplyr::select(sites, Longitude, Latitude), proj4string = crdref)

# Extracting lambda at each clam site ==========================================

Boundaries <- readRDS("../data/Boundaries.rds")
Boundary <- Boundaries$Boundary

lambda.all <- readRDS("../output/lambda_mean.rds")

# Aligning the CRS
crs(lambda.all) <- crs(Boundary)
sitepts_lambda <- spTransform(sitepts, crs(Boundary))

# Extracting lambdas
lambda <- extract(lambda.all, sitepts_lambda)
colnames(lambda) <- c(1993:2018)

lambda <- melt(lambda, varnames = c("site", "Year"), value.name = "lambda") %>% 
  mutate(site = sites$site[site]) %>% 
  subset(!is.na(lambda))

joint <- join(lambda, subset(data, spp == "SAG"))

SAG_array <- reshape2::acast(joint, Year ~ site, value.var = "AvgOfAvgOfbiomass")
otter_array <- reshape2::acast(joint, Year ~ site, value.var = "lambda")

# Filling in zeros
samples <- samples[, colnames(samples) %in% colnames(otter_array)]
SAG_array[samples] <- ifelse(is.na(SAG_array[samples]), 0.0, SAG_array[samples])

# Loading in and extracting current speed at each infauna site =================

rms_field <- raster("../data/current/w001001.adf") %>% 
  crop(Boundary)

rms_pts <- rasterToPoints(rms_field, spatial = TRUE)

# Kriging with 'fields' package with compact covar function
fit <- fastTps(rms_pts@coords, log(rms_pts@data), lon.lat = FALSE, lambda = 0, theta = 400, k = 2)

# Predicting current speed at clam sites
sitepts_rms <- spTransform(sitepts, crs(rms_field))
site_preds <- predict.fastTps(fit, xnew = sitepts_rms@coords) %>% exp()

covars <- mutate(sites, intercept = 1, rms = site_preds) %>% 
  subset(site %in% colnames(otter_array)) %>% 
  dplyr::select(intercept, Latitude, rms) %>% 
  as.matrix()

# Generating prediction locations, covariates, and distances ===================

# Subset lambda by nearshore (and not NA)
X <- readRDS("../data/lambda_covars.rds")
lambda_near <- lambda.all[X[, 2] > 0]
obs <- !is.na(rowSums(lambda_near))
lambda_near <- lambda_near[obs, ]

# Create prediction coordinates and align with RMS data
predpts <- coordinates(lambda.all, spatial = TRUE) %>% 
  magrittr::extract(X[, 2] > 0) %>% 
  magrittr::extract(obs) %>% 
  spTransform(crs(rms_field))

pred_lonlat <- spTransform(predpts, crdref)@coords

# Generate kriged predictions for each of the prediction sites
pred_rms <- predict.fastTps(fit, xnew = predpts@coords) %>% exp()

# Making full X at all prediction sites
pred_covars <- cbind(1, pred_lonlat[, 2], pred_rms)

# Computing distance matrices (using the current data CRS)
wmatch <- sites$site %in% colnames(otter_array)
site_to_site <- pointDistance(sitepts_rms, sitepts_rms, lonlat = FALSE, allpairs = TRUE)[wmatch, wmatch]
preds_to_sites <- pointDistance(predpts, sitepts_rms, lonlat = FALSE, allpairs = TRUE)[, wmatch]
preds_to_preds <- pointDistance(predpts, predpts, lonlat = FALSE, allpairs = TRUE)

# Saving all the data to load into Julia for model fitting =====================

# Some things to help with plotting
sites <- mutate(sites, x = sitepts_lambda@coords[, 1], y = sitepts_lambda@coords[, 2])
pred_sites <- coordinates(lambda.all) %>% 
  magrittr::extract(X[, 2] > 0, ) %>% 
  magrittr::extract(obs, )
  
boundary <- rasterToPoints(Boundaries$BoundaryNA) %>% as.data.frame()

# Wrapping it all up in a list
datalist <- list(lambda = otter_array,
                 y = SAG_array,
                 X = covars,
                 X_all = pred_covars,
                 lambda_all = t(lambda_near),
                 obs_sites = sites,
                 pred_sites = pred_sites,
                 Doo = site_to_site,
                 Duo = preds_to_sites,
                 Duu = preds_to_preds,
                 bay_bound = boundary)

# And saving
saveRDS(datalist, "../data/all.rds")