library(ggplot2)
library(ggrepel)
library(patchwork)
library(plyr)
library(magrittr)
library(reshape2)
library(gganimate)
library(ggsn)

# Reading in data
all <- readRDS('../output/all.rds')

# Ordering sites by latitude and east-west
ordered <- all$obs_sites

ordered$region <- NA
ordered$region[ordered$site_name %in% c("Rush_PCH", "71", "67", "Berg", "Berg_PCH", "PCH_230", "Fingers", "Fingers_PCH")] <- 1
ordered$region[ordered$site_name %in% c("86", "Geikie_PCH", "Geikie", "BlueMouse", "Strawberry", "Boulder_PCH", "Johnson", "Drake")] <- 2
ordered$region[ordered$site_name %in% c("58", "55", "Secret", "Secret_PCH", "52", "233", "Triangle_PCH", "229")] <- 3
ordered$region[ordered$site_name %in% c("249", "243", "46", "40", "43", "Sealers", "170", "211")] <- 4
ordered$region[ordered$site_name %in% c("Leland", "221", "Sturgess", "Puffin", "30")] <- 5

ordered <- arrange(ordered, region, desc(latitude))

# S. giganteus data
sag <- all$y %>% 
  melt(varnames = c("quad", "year", "site_name"), value.name = "count") %>% 
  mutate(site_name = factor(site_name, levels = ordered$site_name))

# Otter data  
data <-all$lambda %>% 
  melt(varnames = c("year", "site_name"), value.name = "lambda")

# Years and sites
years <- unique(data$year)
site_names <- unique(data$site_name)

# Bay boundary raster, for plotting
boundary <- all$bay_bound

# Reading in output
post <- readRDS("../output/chain.rds")

# Figure1: Map of site locations ===============================================

nearshore <- as.data.frame(all$pred_sites)

nonzero <- ddply(sag, .(site_name), summarise,
                 total = sum(count, na.rm = T)) %>% 
  subset(total > 0)

left <- subset(all$obs_sites, x < median(all$obs_sites$x) & site_name %in% nonzero$site)
right <- subset(all$obs_sites, x > median(all$obs_sites$x) & site_name %in% nonzero$site)

pdf(file = "../output/figures/Figure1.pdf",
    width = 5.5, height = 4.5)

boundary %>% 
  ggplot(aes(x, y)) + 
  geom_raster(fill = "steelblue3", alpha = 0.5) +
  geom_point(aes(x, y), data = nearshore, pch = 20, size = 0.1, color = "gray50", inherit.aes = FALSE) +
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
  # scale_fill_continuous(guide = "none") +
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

# Figure 2: posterior preds of z ===============================================

z <- post$z %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  ddply(.(year, site_name), summarise,
        mean = mean(abundance),
        med = median(abundance),
        max = quantile(abundance, 0.95),
        min = quantile(abundance, 0.05))


pdf(file = "../output/figures/Figure2.pdf",
    width = 6.5, height = 8)

z %>% 
  subset(site_name %in% nonzero$site_name) %>%
  ggplot(aes(year, mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  geom_point(data = subset(sag, site_name %in% nonzero$site_name), inherit.aes = FALSE, aes(year, count), size = 0.25) +
  facet_wrap(~as.factor(site_name), scales = "fixed", ncol = 5, dir = "v") + 
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab(expression(italic("S. giganteus")~abundance/quadrat)) + 
  xlab("Year")

dev.off()

# Plots of latent true biomass =================================================

u <- post$u %>% 
  melt(varnames = c("year", "site_name", "iter"), value.name = "abundance") %>% 
  mutate(year = years[year], site_name = factor(site_names[site_name], levels = ordered$site_name)) %>% 
  ddply(.(year, site_name), summarise,
        med = median(abundance),
        max = quantile(abundance, 0.975),
        min = quantile(abundance, 0.025)) %>% 
  join(cbind(all$obs_sites, all$X))


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
  ggplot(aes(site_name, med, color = year)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("posterior median abundance")

# Changes against covariates
u %>% 
  subset(year %in% c(1993, 2018)) %>% 
  ggplot(aes(latitude, med)) + 
  geom_point() + 
  geom_line(aes(group = latitude), color = "black", alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("posterior median abundance")


# Plotting state space =========================================================

joint <- join(data, u) %>% 
  join(sag) %>% 
  join(all$obs_sites)

joint <- mutate(joint, site_name = factor(site_name, levels = ordered$site_name))


joint %>% 
  ggplot(aes(lambda, med, group = site_name)) + 
  geom_line(alpha = 0.8) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + 
  theme_classic() + 
  facet_wrap(~site_name, scales = "free_y", nrow = 11)  +
  ylab("Biomass (g wet weight)")

joint %>% 
  subset(max > 0) %>% 
  ggplot(aes(lambda, med, group = site_name)) + 
  geom_line(alpha = 1.0) +
  scale_color_viridis_c() +
  geom_point(aes(size = year), alpha = 0.8) +
  scale_radius(range = c(0.1, 5)) + 
  theme_classic() + 
  ylab("Clam abundance") + 
  xlab("Otter abundance")

# Plotting site growth rates ==================================================

covars <- cbind(all$obs_sites, all$X)

r <- post$r %>% 
  melt(varnames = c("site_name", "iter"), value.name = "r") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name))

# MCMC chains
r %>% 
  ggplot(aes(iter, r)) + 
  geom_line() + 
  facet_wrap(~site_name, scales = "free_y")

# Site summaries
rsum <- ddply(r, .(site_name), summarise,
              r_med = median(r),
              r_sd = sd(r),
              r_low = quantile(r, 0.025),
              r_high =quantile(r, 0.975)) %>% 
  join(covars)

rsum %>% 
  ggplot(aes(site_name, r_med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = r_low, ymax = r_high)) +
  scale_color_discrete() +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("r")

# Growth by covariates
rsum %>% 
  ggplot(aes(latitude, r_med, label = site_name)) + 
  geom_label() + 
  theme_classic() 

rsum %>% 
  ggplot(aes(rms, r_med, label = site_name)) + 
  geom_label() + 
  theme_classic() 

# Plotting site attack rates ==================================================

a <- post$a %>% 
  melt(varnames = c("site_name", "iter"), value.name = "a") %>% 
  mutate(site_name = factor(site_names[site_name], levels = ordered$site_name))

# MCMC chains
a %>% 
  ggplot(aes(iter, a)) + 
  geom_line() + 
  facet_wrap(~site_name, scales = "free_y")

a %>% 
  ggplot() +
  geom_histogram(aes(a)) + 
  facet_wrap(~site_name)

# Site summaries
asum <- ddply(a, .(site_name), summarise,
              a_med = median(a),
              a_sd = sd(a),
              a_low = quantile(a, 0.1),
              a_high =quantile(a, 0.9),
              p_positive = mean(a > 0)) %>% 
  join(covars)

asum %>% 
  ggplot(aes(site_name, a_med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = a_low, ymax = a_high)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("a")

asum %>% 
  ggplot(aes(site_name, p_positive)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("a")


# Figure 3: maps of r, a, current, and lambda ==================================

pred <- readRDS("../output/predictions.rds")

current <- data.frame(all$pred_sites, current = all$X_all[, 2])

rpred <- pred$r %>% 
  aaply(.(1), median) %>% 
  melt(varnames = c("site"), value.name = "r") %>% 
  cbind(all$pred_sites)

apred <- pred$a %>%
  # aaply(.(1), function(x) mean(x > 0)) %>%
  aaply(.(1), sd) %>%
  melt(varnames = c("site"), value.name = "a") %>%
  cbind(all$pred_sites)


lambda <- t(all$lambda_all) %>% 
  melt(varnames = c("site", "year"), value.name = "lambda") %>% 
  mutate(year = years[year], scaled = ifelse(lambda < 0.1, 0.1, lambda)) %>% 
  cbind(all$pred_sites)

pdf(file = "../output/figures/Figure3.pdf",
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
  scale_fill_viridis_c(name = "P(a > 0)", option = "E") +
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

pdf(file = "../output/figures/Figure4.pdf",
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
  geom_raster(aes(x, y, fill = u), data = upred) + 
  scale_fill_viridis_c(name = "log(n)", option = "E") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  transition_time(year) +
  ggtitle('{frame_time}')

anim_save('clams.gif')


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

covars_a <- factor(c("intercept", "current"), levels = c("intercept", "current"))

beta_a <- post$beta_a %>% 
  melt(varnames = c("coeff", "iter"), value.name = "beta") %>% 
  mutate(coeff = covars_a[coeff]) %>% 
  join(kappa) %>% 
  join(nu)

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

beta_a %>% 
  subset(coeff == "intercept") %>% 
  ggplot(aes(nu, beta)) + 
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
