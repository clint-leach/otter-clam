library(raster)
library(xtable)
library(plyr)
library(reshape2)
library(fields)
library(spam64)
library(magrittr)

# Loading in infauna sampling data =============================================

# Processed avg biomass data (which informs sampling)
biomass <- read.csv("../data/intertidal_prey.csv") %>% 
  subset(!site  %in% c("FrancisE", "GeikieE", "DrakeE", "LesterE", "NetlandE", 
                       "StrawberryE", "WhidbeyE", "JohnsonE"))

# Creating a summary table of all the sites and their locations
sites <- ddply(biomass, .(site), summarise, 
               Longitude = Longitude[1], 
               Latitude = Latitude[1],
               subtidal = type[1] == "subtidal",
               nquad = subtidal * 10 + 10)

# Years each site was sampled
years <- ddply(biomass, .(site), summarise,
               year = unique(Year))

# Filling in quadrats
quadrats <- ddply(sites, .(site), summarise,
                  quad  = 1:nquad)

# Raw count by size and quadrat data for Saxidomus giganteus
counts_wide <- read.csv("../data/SAG_count.csv")

# Processing data to get total by site, quadrat, and year
sample_info <- dplyr::select(counts_wide, c("TypeSiteYearSiteQuadTaxa", "type", "year", "site", "quad", "Taxa")) %>% 
  dplyr::rename(sample_id = TypeSiteYearSiteQuadTaxa)

counts <- counts_wide[, 2:115] %>% as.matrix() %>% 
  melt(varnames = c("sample_id", "size"), value.name = "count") %>% 
  mutate(sample_id = sample_info$sample_id[sample_id], 
         size = substring(size, 2) %>% as.numeric,
         biomass = count * 0.000214 * size ^ 2.78) %>% 
  join(sample_info) %>% 
  ddply(.(year, site, quad), summarise, 
              total = sum(count))

# Filling in zeros
counts <- join(years, quadrats) %>% 
  join(counts) %>% 
  mutate(total = ifelse(is.na(total), 0, total))

# Making spatial points object from site list 
crdref <- CRS(SRS_string = "EPSG:4269")
sitepts <- SpatialPoints(dplyr::select(sites, Longitude, Latitude), proj4string = crdref)

# Extracting lambda at each clam site ==========================================

lambda.all <- brick("../output/lambda_mean.grd")

# Aligning the CRS
crs(lambda.all) <- "EPSG:26708"
sitepts_lambda <- spTransform(sitepts, wkt(lambda.all))

# Extracting lambdas
lambda <- extract(lambda.all, sitepts_lambda)
colnames(lambda) <- c(1993:2018)

lambda <- melt(lambda, varnames = c("site", "year"), value.name = "lambda") %>% 
  mutate(site = sites$site[site]) %>% 
  subset(!is.na(lambda))

joint <- join(lambda, counts)

SAG_array <- reshape2::acast(joint, quad ~ year ~ site, value.var = "total")[1:20, , ]
otter_array <- reshape2::acast(lambda, year ~ site, value.var = "lambda")

# Loading in and extracting current speed at each infauna site =================

rms_field <- raster("../data/current/w001001.adf") %>% 
  crop(lambda.all)

rms_pts <- rasterToPoints(rms_field, spatial = TRUE)

# Kriging with 'fields' package with compact covar function
fit <- fastTps(rms_pts@coords, log(rms_pts@data), lon.lat = FALSE, lambda = 0, theta = 400, k = 2)

# Predicting current speed at clam sites
sitepts_rms <- spTransform(sitepts, wkt(rms_field))
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
  spTransform(wkt(rms_field))

pred_lonlat <- spTransform(predpts, wkt(sitepts))@coords

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
sites <- mutate(sites, 
                x = sitepts_lambda@coords[, 1], 
                y = sitepts_lambda@coords[, 2]) %>% 
  subset(site %in% colnames(otter_array))

counts <- subset(counts, site %in% colnames(otter_array))

pred_sites <- coordinates(lambda.all) %>% 
  magrittr::extract(X[, 2] > 0, ) %>% 
  magrittr::extract(obs, )

Boundaries <- readRDS("../data/Boundaries.rds")
boundary <- rasterToPoints(Boundaries$BoundaryNA) %>% as.data.frame()

# Wrapping it all up in a list
datalist <- list(lambda = otter_array,
                 y = SAG_array,
                 X = covars,
                 X_all = pred_covars,
                 lambda_all = t(lambda_near),
                 byquad = counts,
                 obs_sites = sites,
                 pred_sites = pred_sites,
                 Doo = site_to_site,
                 Duo = preds_to_sites,
                 Duu = preds_to_preds,
                 bay_bound = boundary)

# And saving
saveRDS(datalist, "../data/all.rds")
