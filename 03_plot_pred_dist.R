# ET Mason
# 26 October 2023

# Generate predicted density plots based off best model fit for calculating
# standardized index of abundance

# load libraries ----------------------------------------------------------

library(INLA)
library(sdmTMB)
library(tidyverse)
library(ggplot2)
library(sf)
library(sp)
library(visreg)
library(marmap)
library(lubridate)
# library(s2)

# 01 load data ---------------------------------------------------------------
#needs to be updated
load(file="01_tidy_data/newdf4preds_m4a.2_JSOES_Chinook_salmon_yearling_Interior_Sp.RData") # this includes the model fit from m4a.2 distance and spatial varying day of year 
#load(file="01_tidy_data/newdf4preds_Yearlings.RData") # this includes the model fit from m3a

# predict on new data
set.seed(011823)

pred.spatial = predict(fit, newdata = new_df, return_tmb_object = FALSE) # return_tmb_object is FALSE (it's true when we want an idex)

# let's save it for plotting later
save(pred.spatial,file="01_tidy_data/spatialpreds_m4a2.RData")

# plot predicted density --------------------------------------------------

load(file="01_tidy_data/spatialpreds_m4a2.RData")
#load(file="01_tidy_data/spatialpreds_m3a.RData")

# plot simple maps
plot_map <- function(dat, column) {
  ggplot(dat, aes(Lon.km, Lat.km, fill = {{ column }})) +
    geom_tile() +
    coord_fixed()
}


plot_map(pred.spatial, exp(est)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible
    na.value = "yellow", limits = c(0, quantile(exp(pred.spatial$est), 0.995))
  ) +
  facet_wrap(~year) +
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = paste("maximum estimated density =", round(max(exp(pred.spatial$est))))
  ) 

plot_map(pred.spatial, epsilon_st) +
  scale_fill_gradient2() +
  facet_wrap(~year) +
  ggtitle("Spatiotemporal random effects only")


# plot pretty maps --------------------------------------------------------


library(rgeos)
library(rnaturalearth)
library(rnaturalearthhires)
# remotes::install_github("r-spatial/s2")


map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "United States of America")
# Crop the polygon for plotting and efficiency:
st_bbox(map_data) # find the rough coordinates


ggplot(na_coast_proj) + geom_sf()

sf::st_boundary(na_coast_proj)

par(mfrow =c(1,1))

ggplot() + 
  geom_tile(data = pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = exp(est))) +
  scale_fill_viridis_c( trans = "sqrt",
                        # trim extreme high values to make spatial variation more visible
                        na.value = "yellow", limits = c(0, quantile(exp(pred.spatial$est), 0.995))) +
  geom_sf(data=na_coast_proj,fill = "antique white") +
  facet_wrap(~year) +
  theme_light() +
  labs(fill = "Predicted\nabundance") +
  labs(x = "Longitude", y = "Latitude")+
  ggtitle("Prediction (fixed effects + all random effects)\nYearlings",
          subtitle = paste("maximum estimated abundance =", round(max(exp(pred.spatial$est)))))

ggsave(filename = "03_plots/Yearling_PredPlot_byYr.png",
       dpi = 600,
       width=25,
       height=25,
       units="cm")

ggplot() + 
  geom_raster(data = pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = epsilon_st)) +
  scale_fill_gradient2() +
  geom_sf(data=na_coast_proj,fill = "antique white") +
  facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude")+
  theme_light() +
  labs(title = "Spatiotemporal random effects",
       subtitle = "Yearlings")


ggsave(filename = "03_plots/Yearling_EpsPlot_byYr.png",
       dpi = 600,
       width=25,
       height=25,
       units="cm")


# Need to talk to Eric Ward and Cam Freshwater about these
# coefficient yday varying in space

ggplot()+
  geom_raster(data=pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = zeta_s_scale_yday)) +
  scale_fill_gradient2()+
  geom_sf(data=na_coast_proj,fill = "antique white") +
   facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude")+
  theme_light() +
  labs(title = "Spatially varying yday",
       subtitle = "Yearlings")

# spatially varying intercept
ggplot()+
  geom_raster(data=pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = omega_s)) +
  scale_fill_gradient2()+
  geom_sf(data=na_coast_proj,fill = "antique white") +
   facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude")+
  theme_light() +
  labs(title = "Spatially varying intercept of scale_yday",
       subtitle = "Yearlings")


