# ET Mason
# 26 October 2023

# run sdmTMB with Chinook yearling JSOES data (Rusicka file) - this is just an exercise
# in refamiliarizing myself with sdmTMB; for our purposes, we will use a JSOES
# data extract from Cheryl

# Load libraries -----------------------------------------------------------
library(INLA)
library(fmesher)
library(sdmTMB)
library(tidyverse)
library(ggplot2)
library(visreg)
library(marmap)
library(lubridate)
library(sf)
library(sp)
library(raster)
library(grDevices)



# Example model ------------


## 01 load data ---------------------------------------------------------------

load(file="01_tidy_data/yearling_dat.Rdata") # data formatted with format_raw_data.R

dat <- df %>% 
  dplyr::rename(density = 'numPkm2') %>% 
  dplyr::select(density,yday,month,week,year,year_f,depth,scale_depth,Lon.km,
                Lat.km,Station.Code)

names(dat)


## 02 fit model ---------------------------------------------------------------

# To generate a standardized index of abundance


# fit the model

# start simple using year (factor) and scaled depth as a smooth spline
# start with both the spatial and spatiotemporal random fields on
# we can later test if there is more support for a different random effects structure

m1a <- sdmTMB(
  density ~ 0 +     # zero gets rid of the intercept
    year_f + 
    s(scale_depth), # we may want to replace with, or include separately, distance to shore (the JSOES data files from Cheryl should include this information) - talk to Brian Burke and Dave Huff - they have ideas about this
  mesh = bspde,
  time = "year",
  # name of time variable
  spatial = "on",
  # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "iid",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  reml = TRUE,      # we're going to compare different random effects structures using the same fixed effects
  silent = FALSE
)


fit <- m1a # year as factor, smoothed depth, spatial on, iid on


## 03 check fit ---------------------------------------------------------------

fit
sanity(fit)
AIC(fit) # AIC = 16430.07 (reml TRUE)

# tidy table of parameters for fixed and random effects
tidy(fit,"ran_pars",conf.int = TRUE)
tidy(fit,"fixed",conf.int = TRUE)

# show model fits
visreg::visreg(fit,"scale_depth")



## 04 prediction grid --------------------------------------------------------

# Alternatively, this section can be replaced with its own script file, 
# the output of which we incorporate into section 05 or other run model scripts.

# we need to use the model to predict onto a spatial grid (geographic region of
# interest) in order to generate a standardized index of abundance

# since we don't already have a grid, we need to build one...
# it's important that we also calculate the area of each grid cell since some
# grid cells along the border of our grid will not be the same area -- we
# need to be able to account for this later on so the model knows how to 
# accurately extrapolate densities to
# come up with annual total abundance across the grid

# see https://github.com/pbs-assess/sdmTMB/discussions/151

# create a polygon based off the mesh boundaries
new_poly <- st_convex_hull(st_union(bspde$mesh_sf))
plot(new_poly)
st_crs(new_poly)

# remove land from polygon
poly <- st_difference(new_poly, na_coast_proj)
st_crs(poly) 
plot(poly)
st_bbox(poly)


# create a grid within the polygon using the raster package
resolution <- 2000 # poly is in meters; 2000 is pretty coarse: 2000m x 2000m = 4 km2
r <- raster::raster(as(poly, "Spatial"), resolution = resolution)
rr <- raster::rasterize(as(poly, "Spatial"), r, getCover = TRUE)
plot(rr)

# convert to data frame and calculate grid cell area for each row
grid <- as.data.frame(raster::rasterToPoints(rr))
grid$area <- grid$layer * resolution * resolution # meters squared
grid <- dplyr::filter(grid, area > 0) |> 
  dplyr::select(-layer)

grid$Lon.km <- grid$x / 1000
grid$Lat.km <- grid$y / 1000
grid$area.km <- grid$area / 1e06


# plot to see a display of grid cell areas coded by color
ggplot(grid, aes(Lon.km, Lat.km, colour = area.km)) +
  geom_tile(width = 10, height = 10, fill = NA) + # in meters
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()


# convert to a shapefile using the appropriate crs (we want UTMs) - this is so
# we can fill in the corresponding missing depth data using a bathymetry raster
grid_sdf <- st_as_sf(grid, coords = c("x", "y"), crs = "epsg:32610")


# Get bathymetric data, bounding box = -126,42,-123,49
JSOES_bathy = marmap::as.raster(getNOAA.bathy(-126, -123, 49, 42, res=0.25, keep=TRUE)) # resolution is in minutes
plot(JSOES_bathy)
st_crs(JSOES_bathy)

library(raster)
bathy <- projectRaster(JSOES_bathy, crs = "epsg:32610")
plot(bathy)
st_crs(bathy)

# see https://www.neonscience.org/resources/learning-hub/tutorials/extract-values-rasters-r for more information
depth_max <- raster::extract(bathy,           # raster layer
                             grid_sdf,        # SPDF with centroids for buffer
                             buffer = 1000,   # buffer size, units depend on CRS; here, extract all depths within a 1000m radius of the grid cell centroid
                             fun=max,         # what value to extract
                             df=TRUE)         # return a dataframe? 

# translate land or NAs to 0 meters depth (not underwater) - we'll filter these out later
depth_max.a <- depth_max %>% 
  mutate(depth = case_when(layer > 0 ~ 0,
                           is.na(layer) ~ 0,
                           TRUE ~ layer),
         depth = depth * -1)

# bind depth field with grid data frame
new_df.depth <- cbind(depth_max.a,grid) 

# now expand our data frame to include all covariates in our model (year, day, and
# lat/lon in kilometers) - this is a necessary step for prediction

# replication factor  
yrs <- sort(as.vector(unique(dat$year))) 
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, new_df.depth, simplify = FALSE)) 

# add year covariate
new_df$year_f <-  as.factor(rep(yrs, each = nrow(new_df.depth)))
new_df$year <-  rep(yrs, each = nrow(new_df.depth))

# add lat/lon in km
new_df$Lon.km <- new_df$x / 1000
new_df$Lat.km <- new_df$y / 1000

# create new area field in km^2
new_df$area_km2 <- new_df$area / 1e06


# let's also add a day column to predict on -- we don't currently have yday in the model,
# but later (below) we will end up testing a model with yday as a covariate 
# pick any yday within the range of existing data (day 166 is roughly June 15th)

new_df$yday <- 166 


# grab mean and stdev of depth column in dat (survey loc depths)
# use these to scale prediction grid depths
mu.depth <- mean(dat$depth)     # [1] 112.5454
sd.depth <- sd(dat$depth)       # [1] 115.0395

# do the same for yday
mu.yday <- mean(dat$yday)     # [1] 197.4743
sd.yday <- sd(dat$yday)       # [1] 49.7790

new_df <- new_df %>% 
  filter(depth > 0) %>% 
  mutate(scale_depth = (depth - mu.depth) / sd.depth,
         scale_yday = (yday - mu.yday) / sd.yday)



## 05 load prediction grid and calculate index----------------------------------------------------


# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# now make index
# create a vector of grid cell areas for accurate interpolation, otherwise, 
# get_index() assumes all grid cell areas are equal
area_vec <- new_df$area_km2
index = get_index(pred.index, bias_correct = TRUE, area = area_vec)


## 06 plot temporal trends ---------------------------

library(patchwork)                  # Load patchwork package

# this is just a custom function for publication quality plots
theme_Publication <- function(base_size=20, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour="#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}



# on the linear scale
p1 <- ggplot(index, aes(year, est)) + 
  geom_ribbon(aes(ymin=lwr,ymax=upr), fill = "grey90",alpha=0.5) +
  geom_line(linewidth = 0.5, color = "grey30")+
  ggtitle("Yearlings")+
  xlim(1998,2022)+
  xlab("Year") + ylab("Abundance")+
  theme_Publication() 
p1


# natural log
p2 <- ggplot(index, aes(year, log_est)) +
  geom_ribbon(aes(ymin=log_est-2*se,ymax=log_est+2*se), fill = "grey90",alpha=0.5) +
  geom_line()+
  xlab("Year") + ylab("Ln abundance")+
  theme_Publication()+
  xlim(1998,2022)
p2

Y <- (p1 + p2)


Y

ggsave(filename = "03_plots/SIA_Yrlng_2x2Km_DepthScaled.png",
       dpi = 600,
       width=25,
       height=15,
       units="cm")

# so, now that we've seen the step by step process, we can test different
# random effects structures. Once we've determined the best one, we can use that same
# structure to test different fixed effects models (adding in other covariates)

# Fit different models ----------------------------------------------------

# spatial field on, spatiotemporal field off
m1b <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth),
  mesh = bspde,
  time = "year",
  # name of time variable
  spatial = "on",
  # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "off",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  reml = TRUE,
  silent = FALSE
) 

# spatial field off, spatiotemporal field on
m1c <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth),
  mesh = bspde,
  time = "year",
  # name of time variable
  spatial = "off",
  # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "iid",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  reml = TRUE,
  silent = FALSE
) 

# compare re structure model fits here

fit <- m1c # substitute models

fit
sanity(fit)
AIC(fit)

# m1a, passed sanity check, AIC = 16430.07
# m1b, passed sanity check, AIC = 16581.04 (not improved over m1a)
# m1c, passed sanity check, AIC = 16546.4 (not improved over m1a)



# rerun model m1a using reml FALSE so we can compare different fixed effects models using
# the same random effects structure (makes for a more appropriate AIC comparison)

m1a.remlFALSE <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth),
  mesh = bspde,
  time = "year",
  # name of time variable
  spatial = "on",
  # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "iid",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  reml = FALSE,      # we're going to compare different fixed effects models using the same re structure
  silent = FALSE
)

fit <- m1a.remlFALSE

fit
sanity(fit)
AIC(fit)  # 16431.62...this is what we want to compare with below


# for our new model, let's also account for day of the year
m2a <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth) +
    s(yday),
  mesh = bspde,
  time = "year",
   spatial = "on",
  spatiotemporal = "iid",
  family = tweedie(),
  data = dat,
  reml = FALSE,
  silent = FALSE
)

fit <- m2a

fit
sanity(fit)
AIC(fit)  # 16430.86, doesn't improve things much (effect could be swamped by fall data)

# let's try scaling yday to see if that helps things

dat$scale_yday <- as.numeric(scale(dat$yday))
hist(dat$scale_yday)
plot(dat$scale_yday,dat$density)
plot(dat$scale_depth,dat$density)

# model structure is the same as m2a, but yday is now scaled
m2b <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth) +
    s(scale_yday),
  mesh = bspde,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  family = tweedie(),
  data = dat,
  reml = FALSE,
  silent = FALSE
)

fit <- m2b

fit
sanity(fit)
AIC(fit)  # 16430.86, same result


# we could try including scale_yday as spatially varying, recognizing that the spatial distribution of yearlings as a function of yday likely differs from year to year
m3a <- sdmTMB(
  density ~ 0 + 
    year_f + 
    s(scale_depth) +
    s(scale_yday),
  spatial_varying = ~ scale_yday,    # to account for spatial variation in where they are along their northward migration by day of year 
  mesh = bspde,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  reml = FALSE,
  silent = FALSE
)

fit <- m3a
fit
sanity(fit)
AIC(fit)  # 15986.8...that's a BIG improvement over m1a!! 

# This approach might be a good way of accounting for the influence of the timing of the northward migration on the data.
# For example, a year where a lot of fish are dying could look similar to a year
# where fish ahave quickly moved through the system...need to account for this
# Lisa early on mentioned incorporating the timing of smolt passage at Bonneville Dam 
# (directly relates to a management measure because flow in main stem can be "fixed", flow is related to passage time, and passage time is related to survival)
# potentially using the DART data (see raw data files "smoltdaily")
# if we know the approx. travel time from Bonneville Dam to the open ocean, then
# we could potentially back-calculate from the JSOES survey date what the likely passage date was (given the distribution of dates from the DART data)
# Lisa says McMichael et al. 2011 paper has estimate of travel time, "PNNL-20443 Migratory Behavior and Survival of Juvenile Salmonids in the Lower Columbia River, Estuary, and Plume in 2010 FINAL Report"

# there are likely other approaches as well...

## 01 check fit ---------------------------------------------------------------


fit <- m3a
tidy(fit,"ran_pars",conf.int = TRUE)
tidy(fit,"fixed",conf.int = TRUE)

# show model fits
visreg::visreg(fit,"scale_depth")
visreg::visreg(fit,"scale_yday")

# let's save this fit, along with the prediction grid
save(fit,new_df, file="01_tidy_data/newdf4preds_Yearlings.RData")


## 02 load prediction grid and calculate index----------------------------------------------------

load(file="01_tidy_data/newdf4preds_Yearlings.RData")

# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# now make index
# create a vector of grid cell areas for accurate interpolation, otherwise, 
# get_index() assumes all grid cell areas are equal
area_vec <- new_df$area_km2
index.m3a = get_index(pred.index, bias_correct = TRUE, area = area_vec)


# we can save these predictions for plotting temporal trends
save(pred.index,index.m3a, file="01_tidy_data/SIA_output_m3a.RData")

## 03 plot temporal trends ---------------------------

index <- index.m3a

# on the linear scale
p1 <- ggplot(index, aes(year, est)) + 
  # geom_pointrange(aes(ymin=lwr,ymax=upr))+
  geom_ribbon(aes(ymin=lwr,ymax=upr), fill = "grey90",alpha=0.5) +
  geom_line(linewidth = 0.5, color = "grey30")+
  ggtitle("Yearlings")+
  # ylim(0,3000000)+
  xlim(1998,2022)+
  xlab("Year") + ylab("Abundance")+
  theme_Publication() 
p1


# natural log
p2 <- ggplot(index, aes(year, log_est)) +
  geom_ribbon(aes(ymin=log_est-2*se,ymax=log_est+2*se), fill = "grey90",alpha=0.5) +
  # geom_pointrange(aes(ymin=log_est-2*se,ymax=log_est+2*se)) +
  geom_line()+
  # ggtitle("Yearlings")+
  xlab("Year") + ylab("Ln abundance")+
  theme_Publication()+
  xlim(1998,2022)
p2

Y <- (p1 + p2)


Y

ggsave(filename = "03_plots/SIA_Yrlng_2x2Km_m3a.png",
       dpi = 600,
       width=25,
       height=15,
       units="cm")

