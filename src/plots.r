library(ggplot2)
library(patchwork)
library(plyr)
library(magrittr)
library(reshape2)
library(ggmap)

# Reading in data
sag <- readRDS('../data/SAG_array.rds')
data <- readRDS('../data/otter_array.rds')

# Subsetting to sites with multiple obs
obs <- colSums(sag > 0, na.rm = T)
wts <- which(obs > 1)

sag <- sag[, wts] %>% melt(varnames = c("year", "site"), value.name = "biomass")
data <- data[, wts] %>% melt(varnames = c("year", "site"), value.name = "lambda")

years <- unique(data$year)
site_names <- unique(data$site)

# Reading in output
post <- readRDS("../output/chain_u0maxK.rds")

# Plotting dynamics ============================================================
keep = 25001:50000

u <- post$u[, , keep] %>% 
  melt(varnames = c("year", "site", "iter"), value.name = "biomass") %>% 
  mutate(year = years[year], site = site_names[site])

usum <- ddply(u, .(year, site), summarise,
              med = median(biomass),
              max = quantile(biomass, 0.975),
              min = quantile(biomass, 0.025))

usum %>% 
  ggplot(aes(year, med)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = sag, inherit.aes = FALSE, aes(year, biomass), color = "blue") +
  facet_wrap(~site, scales = "free_y") + 
  theme_classic()

data %>% 
  ggplot(aes(year, lambda)) + 
  geom_line() + 
  facet_wrap(~site, scales = "fixed")

u %>% 
  subset(iter %in% sample(1:25000, 100)) %>% 
  ggplot(aes(year, biomass, group = iter)) + 
  geom_line(alpha = 0.1) + 
  geom_point(data = sag, inherit.aes = FALSE, aes(year, biomass), color = "blue") +
  facet_wrap(~site, scales = "free_y") + 
  theme_classic()

# Plotting local growth rates ==================================================

r <- post$r[, keep] %>% 
  melt(varnames = c("site", "iter"), value.name = "log_r") %>% 
  mutate(r = exp(log_r), site = site_names[site])

r %>% 
  ggplot(aes(site, log_r)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("log(r)")

# Plotting local carrying capacities ===========================================

K <- post$K[, keep] %>% 
  melt(varnames = c("site", "iter"), value.name = "log_K") %>% 
  mutate(K = exp(log_K), site = site_names[site])

K %>% 
  ggplot(aes(site, exp(log_K))) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("Carrying capacity (g wet weight)")
  
# Mapping carrying capacities ==================================================

sites <- read.csv("../data/intertidal_prey.csv") %>% 
  ddply(.(site), summarise,
        Latitude = Latitude[1],
        Longitude = Longitude[1])

Ksum <- ddply(K, .(site), summarise,
              log_K = median(log_K),
              K = median(K)) %>% 
  join(sites)
  
bbox =  c(left = min(sites$Longitude) - 0.01, 
         bottom = min(sites$Latitude) - 0.01,
         right = max(sites$Longitude) + 0.01,
         top = max(sites$Latitude) + 0.01)

gb <- get_stamenmap(bbox, 
                    source = "stamen",
                    maptype = "toner-lite", force = TRUE)

ggmap(gb) + 
  geom_point(aes(Longitude, Latitude, color = K), 
             data = Ksum) + 
  scale_color_viridis_c() + 
  theme_classic()
                                                                                                                                                                                                                                                                                                                                                

