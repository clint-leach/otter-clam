library(ggplot2)
library(ggrepel)
library(patchwork)
library(plyr)
library(dplyr)
library(lubridate)
library(magrittr)
library(reshape2)
library(gganimate)
library(ggsn)
library(sf)

# Reading in data
all <- readRDS('../output/all.rds')

# Otter data  
data <-all$lambda %>% 
  melt(varnames = c("year", "site_name"), value.name = "lambda")

# S. giganteus data
sag <- all$y %>% 
  melt(varnames = c("quad", "year", "site_name"), value.name = "count") %>% 
  join(data)

# Years and sites
years <- unique(data$year)
site_names <- unique(data$site_name)

# Bay boundary raster, for plotting
boundary <- all$bay_bound

# Reading in output
post <- readRDS("../output/chain_main.rds")

# Figure1: Map of site locations ===============================================

# Sites with observed clams
nonzero <- ddply(sag, .(site_name), summarise,
                 total = sum(count, na.rm = T)) %>% 
  subset(total > 0) 

# Cleaning up labeling
left <- subset(all$obs_sites, x < median(all$obs_sites$x) & site_name %in% nonzero$site_name)
right <- subset(all$obs_sites, x >= median(all$obs_sites$x) & site_name %in% nonzero$site_name)

# Otter background
lambda <- t(all$lambda_all) %>% 
  melt(varnames = c("site", "year"), value.name = "lambda") %>% 
  mutate(year = years[year], scaled = ifelse(lambda < 0.5, 0.5, lambda)) %>% 
  cbind(all$pred_sites)

pdf(file = "../output/figures/Figure1.pdf",
    width = 5.5, height = 4.5)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.2) +
  geom_raster(aes(x, y, fill = log(scaled)), data = subset(lambda, year == 2012)) +
  scale_fill_viridis_c(name = expression(log(lambda)), option = "E") +
  geom_point(aes(x, y), data = all$obs_sites, size = 0.5, pch = 4) +
  geom_text_repel(aes(x, y, label = site_name), data = left,
                  size = 2.0, 
                  max.overlaps = Inf,
                  direction = "y", 
                  force_pull = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.75,
                  segment.size = 0.1,
                  box.padding = 0.2,
                  hjust = 0,
                  nudge_x = min(boundary$x) - 6000 - left$x) + 
  geom_text_repel(aes(x, y, label = site_name), data = right,
                  size = 2.0, 
                  max.overlaps = Inf,
                  direction = "y", 
                  force_pull = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.75,
                  segment.size = 0.1,
                  box.padding = 0.2,
                  hjust = 1,
                  nudge_x = max(boundary$x) + 15000 - right$x) + 
  scale_x_continuous(expand = expansion(mult = 0.2)) + 
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() +
  north(x.min = min(boundary$x), x.max= max(boundary$x), y.min = min(boundary$y), y.max = max(boundary$y), symbol = 10, scale = 0.05) +
  annotate("text", x = max(boundary$x) - 1400, y = max(boundary$y), label = "N", size = 2.0) + 
  scalebar(x.min = min(boundary$x), x.max= max(boundary$x), y.min = min(boundary$y), y.max = max(boundary$y), 
           dist = 5, dist_unit = "km", transform = FALSE, location = "bottomleft", st.size = 2.0, height = 0.02, border.size = 0.1) 


dev.off()

# Figure 2: timing of samples ==================================================

# Years of site samples
site_samples <- read.csv("../data/seaotter_preySampling_sites_glacierbay_2022.csv")

# Years of otter transect flights
transect_years <- c(seq(1999, 2004), 2006, 2012)
dist_years <- c(1993:1998, 2004:2005, 2009:2010)


pdf(file = "../output/figures/Figure2.pdf",
    width = 6.5, height = 5)

site_samples %>% 
  subset(site_name %in% nonzero$site_name) %>% 
  ggplot(aes(site_name, year)) + 
  geom_hline(yintercept = transect_years, size = 0.5) +
  geom_hline(yintercept = dist_years, linetype = 2, size = 0.25) +
  geom_point(size = 1) + 
  ylab("year") + 
  xlab("prey sampling site") + 
  theme_bw() + 
  scale_y_continuous(breaks = c(1993:2012), limits = c(1993, 2012), minor_breaks = NULL) + 
  facet_grid(. ~ site_zone, scales = "free_x", space = "free_x") +
  theme(strip.background = element_blank(), panel.grid.major.y = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(axis.text = element_text(size = 7))

dev.off()

# Figure S1: maps of transects and sites =======================================

# Otter sampling
otters <- read.csv("../data/seaOtter_aerialSurvey_southeast_1999_2012/GLBA sea otter strip transect counts in WGS84 1999-2012.csv") %>% 
  rename_with(tolower) %>% 
  mutate(total = adults + pups, 
         date = ymd(date),
         group_x = na_if(group_x, ".") %>% as.numeric(),
         group_y = na_if(group_y, ".") %>% as.numeric())

otter_pts <- otters %>%
  subset(total > 0) %>%
  rowwise() %>%
  mutate(point_geom = list(st_point(c(group_x, group_y)))) %>%
  mutate(point_geom = st_sfc(point_geom, crs = st_crs(4326)) %>% 
           st_transform(crs = 26908))

# Transect locations
transects <- read.csv("../data/seaOtter_aerialSurvey_southeast_1999_2012/GLBA sea otter survey transect coordinates in WGS84 1999-2012.csv") %>% 
  rename_with(tolower) %>% 
  rename(transect_end = transect) %>% 
  mutate(transect = stringr::str_sub(transect_end, 1, 3) %>% as.numeric())

# Transects that were actually flown
flown <- otters %>% 
  ddply(.(year, replicate, stratum, transect), summarise,
        date = date[1]) %>% 
  left_join(transects)

# Creating spatial transect lines object
transect_lines <- ddply(flown, .(year, replicate, stratum, transect), summarise,
                        date = date[1],
                        transect_geom = list(st_linestring(cbind(x_coordinate, y_coordinate)))) %>% 
  st_as_sf(crs = 4326) %>% 
  st_transform(crs = 26908)

# Sites
sitepts <- site_samples %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 26908)

# Glacier bay shoreline
gb <- st_read(dsn = "../data/PH6502/historicl1.shp") %>% 
  st_transform(crs = 26908)

pdf(file = "../output/figures/FigureS1.pdf",
    width = 7, height = 8)

# Maps
ggplot(gb) + 
  geom_sf(size = 0.2) + 
  geom_sf(data = transect_lines, color = "gray50", size = 0.1) + 
  geom_sf(data = sitepts, color = "blue", pch = 4, size = 0.5) + 
  facet_wrap(~year, ncol = 4) + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(axis.text = element_text(size = 7))

dev.off()

# Figure 3: posterior preds of z against lambda ================================

# Ordering sites by latitude and east-west
ordered <- all$obs_sites

ordered$region <- NA
ordered$region[ordered$site_name %in% c("71", "67", "Berg", "Berg_PCH", "PCH_230", "Fingers", "Fingers_PCH")] <- 1
ordered$region[ordered$site_name %in% c("Rush_PCH", "Johnson", "Drake", "86", "Geikie_PCH", "Geikie", "170")] <- 2
ordered$region[ordered$site_name %in% c("58", "55", "Secret", "Secret_PCH", "52", "Boulder_PCH", "Strawberry")] <- 3
ordered$region[ordered$site_name %in% c("249", "Triangle_PCH", "243", "46", "233", "43", "40")] <- 4
ordered$region[ordered$site_name %in% c("229", "Leland", "221", "30", "Puffin", "Sturgess")] <- 5

ordered <- arrange(ordered, region, desc(latitude))

z <- post$zpred %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  ddply(.(year, site_name), summarise,
        mean = mean(abundance),
        med = median(abundance),
        sd = sd(abundance),
        max = quantile(abundance, 0.1),
        min = quantile(abundance, 0.9)) %>% 
  join(data)


# Time series
z %>% 
  subset(site_name %in% nonzero$site_name) %>%
  ggplot(aes(year, mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_jitter(data = subset(sag, site_name %in% nonzero$site_name), 
              inherit.aes = FALSE, aes(year, count), size = 1.0, width = 0.0, height = 0, alpha = 0.5) +
  facet_wrap(~as.factor(site_name), scales = "free_y", ncol = 5, dir = "v") + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab(expression(italic("S. gigantea")~abundance/quadrat)) + 
  xlab("Year")

pdf(file = "../output/figures/Figure3.pdf",
    width = 6.5, height = 8)

# Phase portrait
z %>% 
  subset(site_name %in% nonzero$site_name) %>%
  ggplot(aes(lambda, mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = subset(sag, site_name %in% nonzero$site_name), inherit.aes = FALSE, aes(lambda, count), size = 0.2, alpha = 0.7, shape = 20) +
  facet_wrap(~as.factor(site_name), scales = "fixed", ncol = 5, dir = "v") + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab(expression(italic("S. gigantea")~abundance/quadrat)) + 
  xlab("otter density")

dev.off()

# Figure S2: goodness of fit ===================================================

zmean <- post$zmean %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "site_mean") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  subset(!is.na(site_mean)) %>% 
  ddply(.(year, site_name), summarise,
        mu_mean = mean(site_mean),
        mu_max = quantile(site_mean, 0.975),
        mu_min = quantile(site_mean, 0.025))

zvar <- post$zvar %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "site_var") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  subset(!is.na(site_var)) %>% 
  ddply(.(year, site_name), summarise,
        var_mean = mean(site_var),
        var_max = quantile(site_var, 0.975),
        var_min = quantile(site_var, 0.025))

sagsum <- sag %>% 
  ddply(.(site_name, year), summarise,
        mean_y = mean(count, na.rm = T),
        var_y = var(count, na.rm = T)) %>% 
  join(zmean) %>% 
  join(zvar)

sagsum %>% 
  ggplot(aes(mean_y, mu_mean)) + 
  geom_linerange(aes(mean_y, ymin = mu_min, ymax = mu_max), color = "gray20", size = 0.1) +
  geom_point(size = 0.25) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic() + 
  xlab("observed mean") + 
  ylab("posterior predicted mean") -> figa

sagsum %>% 
  ggplot(aes(var_y, var_mean)) + 
  geom_linerange(aes(var_y, ymin = var_min, ymax = var_max), color = "gray20", size = 0.1) +
  geom_point(size = 0.25) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic() + 
  xlab("observed variance") + 
  ylab("posterior predicted variance") -> figb

pdf(file = "../output/figures/FigureS2.pdf",
    width = 8, height = 4)

figa + figb

dev.off()

# Computing posterior predictive p-value for deviance ==========================

D_y <- rep(0, 20000)
D_rep <- rep(0, 20000)
for(k in 1:20000){
  u <- post$u[, , k]
  sigma <- post$sigma[, k]
  
  for(i in 1:length(site_names)){
    for(t in 1:length(years)){
      
      n <- sigma[i] / (1 - sigma[i]) * u[t, i]
      
      D_y[k] <- D_y[k] - 2 * dnbinom(all$y[, t, i], size = n, p = sigma[i], log = T) %>% sum(na.rm = T)
      
      y_rep <- rnbinom(20, size = n, p = sigma[i])
      
      D_rep[k] <- D_rep[k] - 2 * dnbinom(y_rep[!is.na(all$y[, t, i])], size = n, p = sigma[i], log = T) %>% sum(na.rm = T)
    }
  }
}

p = mean(D_y > D_rep)

# Plots of latent true biomass =================================================

u <- post$u %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  ddply(.(year, site_name), summarise,
        med = median(abundance),
        max = quantile(abundance, 0.975),
        min = quantile(abundance, 0.025)) %>% 
  join(all$obs_sites) %>% 
  join(data)

# Time series by site
u %>% 
  ggplot(aes(year, med)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = sag, inherit.aes = FALSE, aes(year, count), size = 0.5) +
  facet_wrap(~site_name, scales = "free_y", nrow = 10) + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("Abundance per quadrat")


# One-dimensional dynamics by site
u %>% 
  subset(site_name %in% nonzero$site_name) %>% 
  ggplot(aes(site_name, med, color = year)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("posterior median abundance")

# Joint phase portrait
u %>% 
  subset(site_name %in% nonzero$site_name) %>% 
  subset(site_name != "Sealers") %>% 
  ggplot(aes(lambda, med, group = site_name)) + 
  geom_line(alpha = 1.0) +
  scale_color_viridis_c("log(a)") +
  geom_point(aes(size = year), alpha = 0.8) +
  scale_radius(range = c(0.1, 3), guide = "none") + 
  theme_classic() + 
  ylab("posterior median clam abundance") + 
  xlab("otter density")

# Figure 4: plotting initial conditions regressions ============================

means <- colMeans(all$X)
sds <- apply(all$X, 2, sd)

X_pred_nat <- apply(all$X, 2, function(x) seq(min(x), max(x), length.out = 200))

X_pred <- sweep(X_pred_nat, 2, means, "-") %>%
  sweep(2, 1/sds, "*")

X_r_current <- cbind(1, 0, X_pred[, "current"], X_pred[, "current"] ^ 2, 0)
X_r_latitude <- cbind(1, X_pred[, "latitude"], 0, 0, 0)
X_r_shoreline <- cbind(1, 0, 0, 0, X_pred[, "shoreline"])

beta_r <- post$beta_r
beta_r_end <- post$beta_end

# Current 
n_current <- X_r_current %*% beta_r %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
        ) %>% 
  cbind(X_pred_nat)

n_end_current <- X_r_current %*% beta_r_end %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
  ) %>% 
  cbind(X_pred_nat)

# Latitude
n_latitude <- X_r_latitude %*% beta_r %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
  ) %>% 
  cbind(X_pred_nat)

n_end_latitude <- X_r_latitude %*% beta_r_end %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
  ) %>% 
  cbind(X_pred_nat)

# Shoreline
n_shoreline <- X_r_shoreline %*% beta_r %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
  ) %>% 
  cbind(X_pred_nat)

n_end_shoreline <- X_r_shoreline %*% beta_r_end %>% 
  melt(varnames = c("pred", "iter"), value.name = "n") %>% 
  ddply(.(pred), summarise,
        n_med = median(n),
        n_low = quantile(n, 0.1),
        n_high = quantile(n, 0.9)
  ) %>% 
  cbind(X_pred_nat)

# Plots
n_latitude %>% 
  ggplot(aes(latitude, n_med)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(latitude, ymin = n_low, ymax = n_high), alpha = 0.5, fill = "gray60") + 
  geom_line(aes(latitude, n_med), data = n_end_latitude, color = "#67a9cf", size = 1) +
  geom_ribbon(aes(latitude, ymin = n_low, ymax = n_high), alpha = 0.5, data = n_end_latitude, fill = "#67a9cf") + 
  theme_classic() + 
  ylim(c(-2.5, 4)) +
  ggtitle("(a)") +
  xlab(expression("latitude ("*degree*N*")")) + 
  ylab("log(clam abundance/quadrat)") -> figa

n_current %>% 
  ggplot(aes(current, n_med)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(current, ymin = n_low, ymax = n_high), alpha = 0.5, fill = "gray60") + 
  geom_line(aes(current, n_med), data = n_end_current, color = "#67a9cf", size = 1) +
  geom_ribbon(aes(current, ymin = n_low, ymax = n_high), alpha = 0.5, data = n_end_current, fill = "#67a9cf") + 
  theme_classic() + 
  theme(axis.title.y = element_blank()) +
  ylim(c(-2.5, 4)) +
  ggtitle("(b)") +
  xlab("current speed (m/s)") + 
  ylab("log(predicted clam abundance/quadrat)") -> figb

n_shoreline %>% 
  ggplot(aes(shoreline, n_med)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(shoreline, ymin = n_low, ymax = n_high), alpha = 0.5, fill = "gray60") + 
  geom_line(aes(shoreline, n_med), data = n_end_shoreline, color = "#67a9cf", size = 1) +
  geom_ribbon(aes(shoreline, ymin = n_low, ymax = n_high), alpha = 0.5, data = n_end_shoreline, fill = "#67a9cf") + 
  theme_classic() + 
  theme(axis.title.y = element_blank()) +
  ylim(c(-2.5, 4)) +
  ggtitle("(c)") +
  xlab("length of shoreline (m)") + 
  ylab("log(predicted clam abundance/quadrat)") -> figc

pdf(file = "../output/figures/Figure4.pdf",
    width = 7, height = 3)

figa + figb + figc + plot_layout(ncol = 3, nrow = 1)

dev.off()

# Plotting linear predictors ===================================================

# Abundance
eta_r_pred <- post$eta_r_pred %>% 
  melt(varnames = c("site_name", "iter"), value.name = "eta_r_pred") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name))

eta_r <- post$eta_r %>% 
  melt(varnames = c("site_name", "iter"), value.name = "eta_r") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  join(eta_r_pred)

eta_r %>% 
  ggplot(aes(eta_r_pred, eta_r)) + 
  geom_bin2d(bins = 100) + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_fill_viridis_c() + 
  xlab(expression(X[r] * beta[r])) + 
  ylab(expression(eta[r])) + 
  theme_classic() + 
  facet_wrap(~site_name, scales = "fixed")

# Attack rate
eta_a_pred <- post$eta_a_pred %>% 
  melt(varnames = c("site_name", "iter"), value.name = "eta_a_pred") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name))

eta_a <- post$eta_a %>% 
  melt(varnames = c("site_name", "iter"), value.name = "eta_a") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  join(eta_a_pred)

eta_a %>% 
  ggplot(aes(eta_a_pred, eta_a)) + 
  geom_bin2d(bins = 50) + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_fill_viridis_c() + 
  xlab(expression(X[a] * beta[a])) + 
  ylab(expression(eta[a])) + 
  theme_classic() +
  facet_wrap(~site_name, scales = "fixed")

eta_a %>% 
  ddply(.(site_name), summarise,
        med_main = median(eta_a),
        med_pred = median(eta_a_pred),
        upper_main = quantile(eta_a, 0.9),
        lower_main = quantile(eta_a, 0.1),
        upper_pred = quantile(eta_a_pred, 0.9),
        lower_pred = quantile(eta_a_pred, 0.1)) %>% 
  join(all$obs_sites) %>% 
  subset(site_name %in% nonzero$site_name) %>%
  ggplot(aes(med_main, med_pred, color = as.numeric(shoreline))) + 
  geom_point() +
  geom_linerange(aes(ymin = lower_pred, ymax = upper_pred), alpha = 0.5) +  
  geom_linerange(aes(xmin = lower_main, xmax = upper_main), alpha = 0.5) + 
  scale_color_viridis_c(option = "A") +
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic() + 
  xlab(expression(X[a] * beta[a])) + 
  ylab(expression(log(a)))

# Plotting site growth rates ==================================================

r <- post$eta_r %>% 
  melt(varnames = c("site_name", "iter"), value.name = "r") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  mutate(n = exp(r))

# MCMC chains
r %>% 
  ggplot(aes(iter, r)) + 
  geom_line() + 
  facet_wrap(~site_name, scales = "free_y")

# Site summaries
rsum <- ddply(r, .(site_name), summarise,
              n_med = median(n),
              n_low = quantile(n, 0.025),
              n_high = quantile(n, 0.975),
              r_med = median(r),
              r_low = quantile(r, 0.025),
              r_high = quantile(r, 0.975)) %>% 
  join(all$obs_sites)

rsum %>% 
  ggplot(aes(site_name, n_med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = n_low, ymax = n_high)) +
  scale_color_discrete() +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("r")

# Plotting site attack rates ==================================================

a <- post$eta_a %>% 
  melt(varnames = c("site_name", "iter"), value.name = "eta_a") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name), 
         a = exp(eta_a)) 

# MCMC chains
a %>% 
  ggplot(aes(iter, eta_a)) + 
  geom_line() + 
  facet_wrap(~site_name, scales = "free_y")

# Site summaries
asum <- ddply(a, .(site_name), summarise,
              eta_a_med = median(eta_a),
              eta_a_low = quantile(eta_a, 0.1),
              eta_a_high =quantile(eta_a, 0.9)) %>% 
  join(all$obs_sites)

asum %>% 
  ggplot(aes(site_name, eta_a_med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = eta_a_low, ymax = eta_a_high)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab(expression(eta[a]))

# Figure 5: attack rate vs. shoreline length ===================================

# Covariate main effects
X_a_shoreline <- cbind(1, 0, X_pred[, "shoreline"])

a_shoreline <- X_a_shoreline %*% post$beta_a %>% 
  melt(varnames = c("pred", "iter"), value.name = "eta_a") %>% 
  ddply(.(pred), summarise,
        eta_a_med = median(eta_a),
        eta_a_low = quantile(eta_a, 0.1),
        eta_a_high = quantile(eta_a, 0.9)
  ) %>% 
  cbind(X_pred_nat)

pdf(file = "../output/figures/Figure5.pdf",
    width = 6, height = 6)

a_shoreline %>% 
  ggplot(aes(shoreline, eta_a_med)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(shoreline, ymin = eta_a_low, ymax = eta_a_high), alpha = 0.25, fill = "gray60") + 
  geom_point(aes(as.numeric(shoreline), eta_a_med), data = asum, size = 1) +
  geom_linerange(aes(as.numeric(shoreline), ymin = eta_a_low, ymax = eta_a_high), data = asum, alpha = 0.5, size = 0.25) +
  theme_classic() + 
  xlab("length of shoreline (m)") + 
  ylab("log(attack rate)")

dev.off()

# Plotting nu ==================================================================

nu <- data.frame(iter = 1:length(post$nu), nu = post$nu)

mean(nu$nu)
quantile(nu$nu, c(0.05, 0.95))

nu %>% 
  ggplot(aes(iter, nu)) + 
  geom_line()

nu %>% 
  ggplot(aes(nu)) + 
  geom_histogram() +
  geom_vline(xintercept = median(nu$nu), color = "blue") + 
  theme_classic() 

# Plotting sigma ===============================================================

sigma <- post$sigma %>% 
  melt(varnames = c("site_name", "iter"), value.name = "sigma") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name))

# MCMC chains
sigma %>% 
  ggplot(aes(iter, sigma)) + 
  geom_line() + 
  facet_wrap(~site_name)


# Plotting beta_a ==============================================================

covars_a <- factor(c("intercept", "current", "shoreline"), 
                   levels = c("intercept", "current", "shoreline"))

beta_a <- post$beta_a %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta") %>% 
  mutate(coeff = covars_a[coeff])

beta_a %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(~coeff) + 
  ylab(expression(beta[a]))

ddply(beta_a, .(coeff), summarise, 
      mean = mean(beta), 
      lower = quantile(beta, 0.025),
      upper = quantile(beta, 0.975))

beta_a %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(~coeff) + 
  ylab(expression(beta[a]))

beta_a %>% 
  ggplot(aes(beta)) + 
  geom_histogram() + 
  facet_grid(.~coeff) + 
  geom_vline(xintercept = 0) +
  theme_classic() +
  xlab(expression(beta[a]))

pairs(t(post$beta_a))

# Plotting beta_r ==============================================================

covars <- factor(c("intercept", "latitude", "current", "current_sq", "shoreline"), 
                 levels = c("intercept", "latitude", "current", "current_sq", "shoreline"))

beta_r <- post$beta_r %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta") %>% 
  mutate(coeff = covars[coeff])

ddply(beta_r, .(coeff), summarise, 
      mean = mean(beta), 
      lower = quantile(beta, 0.025),
      upper = quantile(beta, 0.975))

beta_r %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(~coeff) + 
  ylab(expression(beta[r]))

beta_r %>% 
  ggplot(aes(beta)) + 
  geom_histogram() + 
  facet_grid(~coeff, scales = "free") + 
  theme_classic() + 
  geom_vline(xintercept = 0) +
  xlab(expression(beta[r]))

beta_end <- post$beta_end %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta_end") %>% 
  mutate(coeff = covars[coeff]) %>% 
  join(beta_r) %>% 
  mutate(smaller = abs(beta_end) < abs(beta),
         larger = beta_end > beta)

ddply(beta_end, .(coeff), summarise, 
      mean = mean(beta_end), 
      lower = quantile(beta_end, 0.05),
      upper = quantile(beta_end, 0.95))

beta_end %>% 
  ddply(.(coeff), summarise,
        psmaller = mean(smaller),
        plarger = mean(larger))


# Figure S3: persistence of small clams ========================================

# Summing over quadrats to get size distributions
sizedist <- read.csv("../data/seaotter_preySampling_zeroPopulated_glacierbay_2022.csv") %>% 
  subset(spp == "SAG") %>% 
  ddply(.(site, year, size, spp), summarise,
        total = sum(count),
        nquad = length(unique(quad)),
        density = total / nquad)

pdf(file = "../output/figures/FigureS3.pdf",
    width = 7, height = 9)

sizedist %>% 
  subset(site %in% nonzero$site_name) %>%
  subset(!(site %in% c("Drake", "Sturgess", "243", "249"))) %>% 
  mutate(cm = floor(size / 10)) %>% 
  ddply(.(site, year, cm), summarise,
        density = sum(density)) %>% 
  ggplot(aes(cm, density, color = year, group = year)) + 
  geom_line(size = 0.5) + 
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "E") +
  theme_classic() +
  scale_x_continuous(breaks = c(4, 8, 12)) +
  theme(axis.text = element_text(size = 7)) +
  facet_wrap(~ site, scales = "free", ncol = 5) + 
  ylab("mean count/quadrat")

dev.off()

# Figure S4: posterior of nu ===================================================

pdf(file = "../output/figures/FigureS4.pdf",
    width = 4, height = 4)

nu %>% 
  ggplot(aes(nu)) + 
  geom_histogram() +
  geom_vline(xintercept = median(nu$nu), color = "blue") + 
  theme_classic() +
  xlab(expression(nu))

dev.off()
