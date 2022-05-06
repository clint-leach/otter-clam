library(plyr)
library(dplyr)
library(magrittr)

# Loading in infauna sampling data =============================================

# Sampling Sites (reading in 'SitesInfo' sheet of Xcel file, exported as a csv)
sites <- read.csv('../data/site_samples.csv') %>% 
  mutate(site_selection = dplyr::recode(Random.or.Selected,
                                        "Random - not resampled" = "Random",
                                        "Selected \"High Clam Densities\"" = "Selected",
                                        "Selected Clam Beds" = "Selected",
                                        "Selected Clam Beds - not resampled" = "Selected")) %>% 
  dplyr::rename(site_name = Site.Name) %>%  
  ddply(.(site_name), summarise, 
        latitude = Latitude[1],
        longitude = Longitude[1],
        site_selection = site_selection[1],
        site_zone = Zone[1],
        nsamples = length(Sample.Year)
        ) %>% 
  mutate(latitude = as.numeric(latitude), 
         longitude = as.numeric(longitude))

# Raw data (reading in 'GlacierBayPreyData2021' sheet of Xcel file, exported as a csv)
raw <- read.csv("../data/prey_data.csv") 

# Species table
spp <- ddply(raw, .(spp), summarise,
             species_name = speciesname[1],
             species_type = prey[1]) %>% 
  subset(spp != "NON")

# Summing over duplicate rows
raw <- ddply(raw, .(site_name, year, quad, spp, size), summarise,
             count = sum(count))

# List of quadrats for each sampling event
quad <- ddply(raw, .(site_name, year), summarise,
              quad = unique(quad))  

# Size range observed for each species
sizes <- ddply(raw, .(spp), summarise,
               max = max(size, na.rm = T)) %>% 
  ddply(.(spp), summarise,
        size = 14:max)

observables <- ddply(sizes, .(spp, size), summarise,
                     site_name = sites$site_name) %>% 
  join(quad)

# Checking for unsampled quadrats (to cross-reference with known missings)
observables %>% 
  ddply(.(site_name, year), summarise,
        nquad = length(unique(quad))) %>% 
  subset(!nquad %in% c(10, 20))

# Joining the nonzero counts with the observables to fill in zeros
counts <- join(observables, raw) %>% 
  mutate(count = ifelse(is.na(count), 0.0, count))

# Saving output
write.csv(counts,"../output/zero_augmented.csv", row.names = FALSE)
write.csv(sites, "../output/sites.csv", row.names = FALSE)
write.csv(spp, "../output/spp.csv", row.names = FALSE)


