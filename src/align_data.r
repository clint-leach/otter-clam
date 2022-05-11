library(raster)
library(xtable)
library(plyr)
library(reshape2)
library(magrittr)

# Loading in infauna sampling data =============================================

# Raw data , downloaded from https://doi.org/10.5066/P9LODH0Z and stored in 'data' directory
counts <- read.csv("../data/seaotter_preySampling_zeroPopulated_glacierbay_2022.csv") %>% 
  subset(species_code == "LES") %>% 
  ddply(.(site_name, year, quad), summarise,
        total = sum(count))

# Sampling Sites 
sites <- read.csv("../data/seaotter_preySampling_sites_glacierbay_2022.csv") %>% 
  ddply(.(site_name), summarise, 
        latitude = latitude[1],
        longitude = longitude[1],
        site_selection = site_selection[1],
        site_zone = site_zone[1],
        nsamples = length(year))


# Making spatial points object from site list 
crdref <- CRS(SRS_string = "EPSG:4269")
sitepts <- SpatialPoints(dplyr::select(sites, longitude, latitude), proj4string = crdref)

# Extracting lambda at each clam site ==========================================

lambda.all <- brick("../output/lambda_mean.grd")

# Aligning the CRS
crs(lambda.all) <- "EPSG:26708"
sitepts_lambda <- spTransform(sitepts, wkt(lambda.all))

# Extracting lambdas
cells <- raster::extract(lambda.all, sitepts_lambda, cellnumbers = TRUE)[, 1]
lambda <- raster::extract(lambda.all, sitepts_lambda)
colnames(lambda) <- c(1993:2018)

lambda <- melt(lambda, varnames = c("site_name", "year"), value.name = "lambda") %>% 
  mutate(site_name = sites$site_name[site_name]) %>% 
  subset(!is.na(lambda))

joint <- join(lambda, counts)

SAG_array <- reshape2::acast(joint, quad ~ year ~ site_name, value.var = "total")[1:20, , ]
otter_array <- reshape2::acast(lambda, year ~ site_name, value.var = "lambda")

# Loading in and extracting current speed at each infauna site =================

rms_field <- raster("../data/current/w001001.adf") %>% 
  crop(lambda.all)

rms_pts <- rasterToPoints(rms_field, spatial = TRUE)
sitepts_rms <- spTransform(sitepts, wkt(rms_field))

# Find and extract closest current layer cell to each of the sample sites
dists <- pointDistance(rms_pts, sitepts_rms)
closest <- apply(dists, 2, which.min)
site_preds <- rms_pts@data[closest, 1]

covars <- mutate(sites, 
                 rms = site_preds) %>% 
  subset(site_name %in% colnames(otter_array)) %>% 
  dplyr::select(latitude, rms) %>% 
  as.matrix()

# Generating prediction locations, covariates, and distances ===================

# Loading in otter covariates
X <- readRDS("../data/lambda_covars.rds")

# Subset lambda by nearshore (and not NA)
obs_nearshore = which(X[, 2] > 0 & !is.na(rowSums(lambda.all[])))
lambda_near <- lambda.all[obs_nearshore]

# Create prediction coordinates and align with RMS data
predpts <- coordinates(lambda.all, spatial = TRUE) %>% 
  magrittr::extract(obs_nearshore) %>% 
  spTransform(wkt(rms_field))

pred_lonlat <- spTransform(predpts, wkt(sitepts))@coords

# Extracting current speed from nearest raster cell
dists <- pointDistance(rms_pts, predpts)
closest <- apply(dists, 2, which.min)
pred_rms <- rms_pts@data[closest, 1]

# Making full X at all prediction sites
pred_covars <- cbind(pred_lonlat[, 2], pred_rms)

# Computing distance matrices (using the current data CRS)
wmatch <- sites$site %in% colnames(otter_array)
site_to_site <- pointDistance(sitepts_rms, sitepts_rms, lonlat = FALSE, allpairs = TRUE)[wmatch, wmatch]
preds_to_sites <- pointDistance(predpts, sitepts_rms, lonlat = FALSE, allpairs = TRUE)[, wmatch]
preds_to_preds <- pointDistance(predpts, predpts, lonlat = FALSE, allpairs = TRUE)

# Saving all the data to load into Julia for model fitting =====================

# Some things to help with plotting
sites <- mutate(sites, 
                x = sitepts_lambda@coords[, 1], 
                y = sitepts_lambda@coords[, 2]) %>% 
  subset(site_name %in% colnames(otter_array))

pred_sites <- coordinates(lambda.all) %>% 
  magrittr::extract(obs_nearshore, )

Boundaries <- readRDS("../data/Boundaries.rds")
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
saveRDS(datalist, "../output/all.rds")
