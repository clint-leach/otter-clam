library(ggplot2)
library(ggrepel)
library(patchwork)
library(plyr)
library(magrittr)
library(reshape2)
library(gganimate)

# Reading in data
all <- readRDS('../data/all.rds')
sag <- all$byquad
data <-all$lambda

data <- data %>% melt(varnames = c("year", "site"), value.name = "lambda")

ordered <- arrange(all$obs_sites, Latitude)

sag <- mutate(sag, site = factor(site, levels = ordered$site))

years <- unique(data$year)
site_names <- unique(data$site)

subtidal_names <- c("Berg", "BlueMouse", "Fingers", "Geikie", "Johnson", "Leland", 
                    "Puffin", "Sealers", "Secret", "Strawberry", "Wachusett")
subtidal <- site_names %in% subtidal_names

boundary <- all$bay_bound

# Reading in output
post <- readRDS("../output/chain_aprior.rds")
pred <- readRDS("../output/prediction_aprior.rds")

# Looking at temp data =========================================================

ocean <- read.csv("../data/SEAN_OC_through_2013.csv") %>% 
  mutate(date_gmt = mdy(date_gmt),
         time_gmt = hm(time_gmt),
         datetime = date_gmt + time_gmt,
         year = year(datetime),
         month = year + month(datetime))


ocean %>% 
  subset(depth == 1) %>% 
  ggplot(aes(datetime, salinity)) + 
  geom_line() + 
  facet_wrap(~station)

annual <- ocean %>% 
  ddply(.(station, depth, year), summarise,
        meant = mean(temperature),
        maxt = max(temperature),
        mint = min(temperature),
        vart = var(temperature), 
        sal = mean(salinity),
        latitude = latitude[1])

annual %>% 
  subset(depth == 0) %>% 
  ggplot(aes(latitude, meant)) + 
  geom_path(aes(color = year, group = year), size = 1, alpha = 0.7) + 
  scale_color_viridis_c() + 
  theme_classic()

annual %>% 
  subset(depth == 0) %>% 
  ggplot(aes(year, meant)) + 
  geom_line() + 
  facet_wrap(~station, scales = "free_y")

annual %>% 
  subset(depth == 0 & year == 1993) %>% 
  ggplot(aes(latitude, sal)) + 
  geom_point(aes(color = year))


# Figure1: Map of site locations ===============================================

nearshore <- as.data.frame(all$pred_sites)

pdf(file = "../output/Figure1.pdf",
    width = 5, height = 4)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_point(aes(x, y), data = nearshore, pch = 20, size = 0.1, color = "gray50", inherit.aes = FALSE) +
  geom_point(aes(x, y), data = all$obs_sites, size = 0.5, pch = 4) +
  geom_text_repel(aes(x, y, label = site), data = subset(all$obs_sites, x < median(all$obs_sites$x)),
                  size = 2.0, 
                  max.overlaps = Inf,
                  direction = "y", 
                  force_pull = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.75,
                  segment.size = 0.1,
                  box.padding = 0.2,
                  hjust = 0,
                  nudge_x = min(boundary$x) - 6000 - subset(all$obs_sites, x < median(all$obs_sites$x))$x) + 
  geom_text_repel(aes(x, y, label = site), data = subset(all$obs_sites, x > median(all$obs_sites$x)),
                  size = 2.0, 
                  max.overlaps = Inf,
                  direction = "y", 
                  force_pull = 0,
                  min.segment.length = 0,
                  segment.alpha = 0.75,
                  segment.size = 0.1,
                  box.padding = 0.2,
                  hjust = 1,
                  nudge_x = max(boundary$x) + 15000 - subset(all$obs_sites, x > median(all$obs_sites$x))$x) + 
  # scale_fill_continuous(guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.2)) + 
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() 

dev.off()

# Figure 2: posterior preds of z ===============================================

z <- post$z %>% 
  melt(varnames = c("year", "site", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site = factor(site_names[site], levels = ordered$site)) %>% 
  ddply(.(year, site), summarise,
        mean = mean(abundance),
        med = median(abundance),
        max = quantile(abundance, 0.95),
        min = quantile(abundance, 0.05))

pdf(file = "../output/Figure2.pdf",
    width = 6.5, height = 8)

z %>% 
  ggplot(aes(year, mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = sag, inherit.aes = FALSE, aes(year, total), size = 0.25) +
  facet_wrap(~as.factor(site), scales = "fixed", nrow = 10) + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab(expression(italic("S. giganteus")~abundance/quadrat)) + 
  xlab("Year")

dev.off()

# Plots of latent true biomass =================================================

u <- post$u %>% 
  melt(varnames = c("year", "site", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site = factor(site_names[site], levels = ordered$site)) %>% 
  ddply(.(year, site), summarise,
        med = median(abundance),
        max = quantile(abundance, 0.975),
        min = quantile(abundance, 0.025)) %>% 
  join(cbind(all$obs_sites, all$X[, 2:3]))

u %>% 
  ggplot(aes(year, med)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = sag, inherit.aes = FALSE, aes(year, total), size = 0.5) +
  facet_wrap(~site, scales = "free_y", nrow = 10) + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("Abundance per quadrat")


u %>% 
  join(all$obs_sites) %>% 
  ggplot(aes(site, med, color = year)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("posterior median abundance")

u %>% 
  subset(year %in% c(1993, 2018)) %>% 
  ggplot(aes(rms, med)) + 
  geom_point() + 
  geom_line(aes(group = Latitude), color = "black", alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("posterior median abundance")


# Plotting state space =========================================================

joint <- join(data, u) %>% 
  join(sag) %>% 
  join(all$obs_sites)

joint <- mutate(joint, site = factor(site, levels = ordered$site))


joint %>% 
  ggplot(aes(lambda, med, group = site)) + 
  geom_line(alpha = 0.8) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  # geom_point(aes(lambda, total), size = 2) +
  theme_classic() + 
  facet_wrap(~site, scales = "fixed", nrow = 11)  +
  ylab("Biomass (g wet weight)")

joint %>% 
  ggplot(aes(lambda, med, group = site)) + 
  scale_color_discrete(guide = F) + 
  geom_line(alpha = 0.2) +
  geom_point(aes(size = year), alpha = 0.8) + 
  scale_radius(range = c(0.1, 5)) + 
  theme_classic() + 
  ylab("Clam abundance") + 
  xlab("Otter abundance")

# Plotting site growth rates ==================================================

covars <- cbind(all$obs_sites, all$X[, 2:3]) %>% 
  join(subset(data, year == 2018))

r <- post$r %>% 
  melt(varnames = c("site", "iter"), value.name = "r") %>% 
  mutate(site = factor(site_names[site], levels = ordered$site))

r %>% 
  ggplot(aes(iter, r)) + 
  geom_line() + 
  facet_wrap(~site, scales = "free_y")

rsum <- ddply(r, .(site), summarise,
              r_med = median(r),
              r_low = quantile(r, 0.025),
              r_high =quantile(r, 0.975)) %>% 
  join(covars)

rsum %>% 
  ggplot(aes(site, r_med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = r_low, ymax = r_high)) +
  scale_color_discrete() +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("r") -> fig

rsum %>% 
  ggplot(aes(Latitude, r_med, label = site)) + 
  geom_label() + 
  theme_classic() 

rsum %>% 
  ggplot(aes(rms, r_med)) + 
  geom_point() + 
  theme_classic() 

# Plotting site attack rates ==================================================

a <- post$a %>% 
  melt(varnames = c("site", "iter"), value.name = "a") %>% 
  mutate(site = factor(site_names[site], levels = ordered$site))

a %>% 
  ggplot(aes(iter, a)) + 
  geom_line() + 
  facet_wrap(~site, scales = "free_y")

asum <- ddply(a, .(site), summarise,
              med = median(a),
              low = quantile(a, 0.1),
              high =quantile(a, 0.9)) %>% 
  join(covars)

asum %>% 
  ggplot(aes(site, med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = low, ymax = high)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("a")


asum %>% 
  ggplot(aes(rms, med, label = site)) + 
  geom_label() + 
  theme_classic() 

asum %>% 
  ggplot(aes(rms, med)) + 
  geom_point() + 
  theme_classic()


# Figure 3: maps of r, a, current, and lambda ==================================

current <- data.frame(all$pred_sites, current = all$X_all[, 3])

rpred <- pred$r %>% 
  aaply(.(1), median) %>% 
  melt(varnames = c("site"), value.name = "r") %>% 
  cbind(all$pred_sites)

apred <- pred$a %>% 
  aaply(.(1), median) %>% 
  melt(varnames = c("site"), value.name = "a") %>% 
  cbind(all$pred_sites)

lambda <- t(all$lambda_all) %>% 
  melt(varnames = c("site", "year"), value.name = "lambda") %>% 
  mutate(year = years[year], scaled = ifelse(lambda < 0.1, 0.1, lambda)) %>% 
  cbind(all$pred_sites)

pdf(file = "../output/Figure3.pdf",
    width = 6.5, height = 6.5)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = r), data = rpred) + 
  scale_fill_viridis_c(name = "r", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() + 
  ggtitle("(a)") -> figa

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill =a), data = apred) + 
  scale_fill_viridis_c(name = "a", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() +
  ggtitle("(b)") -> figb

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = current), data = current) + 
  scale_fill_viridis_c(name = "current\nspeed", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() + 
  ggtitle("(c)") -> figc

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = log(lambda)), data = subset(lambda, year == 2018)) + 
  scale_fill_viridis_c(name = expression(log(lambda)), option = "E") +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  coord_equal() + 
  ggtitle("(d)") -> figd

figa + figb + figc + figd + plot_layout(ncol = 2)

dev.off()

# Figure 4: maps of clam biomass through time ==================================

pred_sites <- as.data.frame(all$pred_sites) %>% 
  mutate(site = 1:nrow(all$pred_sites))

upred <- pred$u %>% 
  aaply(.(1, 2), median) %>% 
  melt(varnames = c("year", "site"), value.name = "u") %>% 
  mutate(year = years[year]) %>%
  join(pred_sites)

# Multifacet plot of biomass
plot_years = seq(1993, 2018, by = 5)

pdf(file = "../output/Figure4.pdf",
    width = 6, height = 9)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = u), data = subset(upred, year %in% plot_years)) + 
  scale_fill_viridis_c(name = "n", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~year, ncol = 2) + 
  coord_equal()

dev.off()

# Movie of biomass
boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = truncated), data = upred) + 
  geom_point(aes(x, y), data = all$obs_sites, pch = 4, size = 1, color = "white") + 
  scale_fill_viridis_c(name = "log(n)", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  transition_time(year) +
  ggtitle('{frame_time}')

anim_save('clams.gif')

# Mapping difference
udiff <- upred %>% 
  subset(year %in% c(1993, 2018)) %>% 
  ddply(.(site), summarise,
        prop = u[2] / u[1],
        diff = u[2] - u[1],
        x = x[1],
        y = y[1])

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = diff), data = udiff) + 
  scale_fill_viridis_c(name = "n", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank()) +
  coord_equal()

# Multifacet plot of uncertainty
boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, color = sigma), data = subset(upred, year %in% 1998)) + 
  geom_point(aes(x, y), data = all$obs_sites, pch = 4, size = 2) + 
  scale_fill_viridis_c(name = expression(sd(B[2018])), option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank()) + 
  facet_wrap(~year, nrow = 3)

# Plotting biomass as function of environment

upred <- join(upred, Xall) %>% 
  join(lambda)


upred %>% 
  subset(year %in% c(1993, 2018)) %>% 
  ggplot(aes(latitude, u, color = factor(year))) + 
  geom_point() + 
  theme_classic() + 
  facet_grid(year ~.)

upred %>% 
  subset(year %in% c(1993, 2018)) %>% 
  ddply(.(site), summarise,
        diff = diff(u), 
        latitude = latitude[1],
        current = current[1], 
        lambda = lambda[2]) %>% 
  ggplot(aes(lambda, diff)) + 
  geom_point() +
  theme_classic()

# Movie of lambda ==============================================================

lambda_years = 1993:2018

lambda <- t(all$lambda_all) %>% 
  melt(varnames = c("site", "year"), value.name = "lambda") %>% 
  mutate(year = lambda_years[year], scaled = ifelse(lambda < 0.01, 0.01, lambda)) %>% 
  cbind(all$pred_sites)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = log(scaled)), data = lambda) + 
  scale_fill_viridis_c(name = "otter\nabundance", option = "B") +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  # facet_wrap(~year)
  transition_time(year) +
  ggtitle('{frame_time}')

anim_save('otters.gif')

# Biomass flux map ===================================================

flux <- pred$flux %>% 
  aaply(.(1, 2), median) %>% 
  melt(varnames = c("year", "site"), value.name = "flux") %>% 
  mutate(year = years[year]) %>%
  join(pred_sites)

plot_years <- seq(1993, 2018, by = 5)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_raster(aes(x, y, fill = flux), data = subset(flux, year %in% plot_years)) + 
  scale_fill_viridis_c(name = "log(flux)", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  coord_equal() + 
  facet_wrap(~year)


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

# Plotting kappa ===============================================================

kappa <- data.frame(iter = 1:length(post$kappa), kappa = post$kappa)

mean(kappa$kappa)
quantile(kappa$kappa, c(0.05, 0.95))

kappa %>% 
  ggplot(aes(iter, kappa)) + 
  geom_line()

kappa %>% 
  ggplot(aes(kappa)) + 
  geom_histogram() +
  geom_vline(xintercept = median(kappa$kappa), color = "blue") + 
  theme_classic() 

# Plotting sigma ===============================================================

sigma <- post$sigma %>% 
  melt(varnames = c("iter"), value.name = "sigma")

mean(sigma$sigma)
quantile(sigma$sigma, c(0.05, 0.95))

sigma %>% 
  ggplot(aes(1:length(sigma), sigma)) + 
  geom_line()

sigma %>% 
  ggplot(aes(sigma)) + 
  geom_histogram() +
  geom_vline(xintercept = median(sigma$sigma), color = "blue") + 
  theme_classic() 

# Plotting beta_a ==============================================================

covars_a <- factor(c("intercept", "current", "current_sq"), levels = c("intercept", "current", "current_sq"))

beta_a <- post$beta_a %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta") %>% 
  mutate(coeff = covars_a[coeff]) %>% 
  join(kappa)

ddply(beta_a, .(coeff), summarise, 
      mean = mean(beta), 
      lower = quantile(beta, 0.05),
      upper = quantile(beta, 0.95))

beta_a %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(~coeff) + 
  ylab(expression(beta[a]))

beta_a %>% 
  ggplot(aes(beta)) + 
  geom_histogram() + 
  facet_grid(~coeff) + 
  theme_classic() + 
  xlab(expression(beta[a]))

beta_a %>% 
  subset(coeff == "intercept") %>% 
  ggplot(aes(kappa, beta)) + 
  geom_point()

# Plotting beta_r ==============================================================

covars <- factor(c("intercept", "latitude", "current", "current_sq"), levels = c("intercept", "latitude", "current", "current_sq"))

beta_r <- post$beta_r %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta") %>% 
  mutate(coeff = covars[coeff]) %>% 
  join(nu)

ddply(beta_r, .(coeff), summarise, 
      mean = mean(beta), 
      lower = quantile(beta, 0.05),
      upper = quantile(beta, 0.95))

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
  xlab(expression(beta[r]))

beta_r %>% 
  subset(coeff == "intercept") %>% 
  ggplot(aes(nu, beta)) + 
  geom_point()

# Other prey ===================================================================

comm <- read.csv("../data/intertidal_prey.csv")


comm %>% 
  subset(site %in% c("Puffin", "Leland", "Strawberry", "Boulder_PCH")) %>% 
  ggplot(aes(Year, AvgOfAvgOfbiomass, color = site)) + 
  geom_line() + 
  facet_wrap(~spp, scales = "free")

comm %>% 
  subset(site %in% c("Puffin", "Leland", "Strawberry", "Boulder_PCH")) %>% 
  ddply(.(site, Year), summarise,
        total = sum(AvgOfAvgOfbiomass)) %>% 
  ggplot(aes(Year, total, color = site)) + 
  geom_line()

