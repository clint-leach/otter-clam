library(raster)
library(xtable)
library(plyr)
library(reshape2)
library(magrittr)
library(sf)

# Loading in infauna sampling data =============================================

# Raw data , downloaded from https://doi.org/10.5066/P9LODH0Z and stored in 'data' directory
counts <- read.csv("../data/seaotter_preySampling_zeroPopulated_glacierbay_2022.csv") %>% 
  subset(species_code == "SAG") %>% 
  ddply(.(site_name, year, quad), summarise,
        total = sum(count))

# Sampling Sites 
sites <- read.csv("../data/seaotter_preySampling_sites_glacierbay_2022.csv") %>% 
  ddply(.(site_name), summarise, 
        latitude = latitude[1],
        longitude = longitude[1],
        site_selection = site_selection[1],
        site_zone = site_zone[1],
        nsamples = length(year)) %>% 
  subset(latitude <= 58.75745) %>% 
  subset(!(latitude  > 58.72960 & longitude < -136.3606))

# Making spatial points object from site list 
crdref <- CRS(SRS_string = "EPSG:4269")
sitepts <- SpatialPoints(dplyr::select(sites, longitude, latitude), proj4string = crdref)

# Extracting lambda at each clam site ==========================================

lambda.all <- brick("../output/lambda_mean.grd")

# Aligning the CRS
crs(lambda.all) <- "EPSG:26708"
sitepts_lambda <- spTransform(sitepts, wkt(lambda.all))

# Extracting lambdas
lambda <- raster::extract(lambda.all, sitepts_lambda, buffer = 800) %>% 
  laply(colMeans, na.rm = T)

lambda <- melt(lambda, varnames = c("site_name", "year"), value.name = "lambda") %>% 
  mutate(site_name = sites$site_name[site_name], year = c(1993:2018)[year]) %>% 
  subset(year < 2013)

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
sites$current <- rms_pts@data[closest, 1]

# Shoreline complexity =========================================================

# Shoreline length
site_area <- sites %>% 
  dplyr::select(site_name, longitude, latitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(26908) %>% 
  st_buffer(dist = 1000)

# Glacier Bay coastline, Project PH6502 downloaded from https://nsde.ngs.noaa.gov/
gb <- st_read(dsn = "../data/PH6502/historicl1.shp") %>% 
  st_transform(crs = 26908)

# Coastline in vicinity of each site
coast <- st_intersection(site_area, gb) 
coast$length <- st_length(coast)

sites <- ddply(coast, .(site_name), summarise,
              shoreline = sum(length)) %>% 
  join(sites)

covars <- sites %>% 
  dplyr::select(latitude, current, shoreline) %>% 
  as.matrix()

# Generating prediction locations, covariates, and distances ===================

# Loading in otter covariates
X <- readRDS("../data/lambda_covars.rds")

# Subset lambda by nearshore (and not NA)
obs_nearshore = which(X[, 2] > 0 & !is.na(rowSums(lambda.all[])))
lambda_near <- lambda.all[obs_nearshore]

# Computing distance matrices (using the current data CRS)
wmatch <- sites$site_name %in% colnames(otter_array)
site_to_site <- pointDistance(sitepts_rms, sitepts_rms, lonlat = FALSE, allpairs = TRUE)[wmatch, wmatch]

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
                 lambda_all = t(lambda_near),
                 obs_sites = sites,
                 pred_sites = pred_sites,
                 Doo = site_to_site,
                 bay_bound = boundary)

# And saving
saveRDS(datalist, "../output/all.rds")
