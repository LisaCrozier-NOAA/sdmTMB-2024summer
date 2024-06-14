# Lisa Crozier
# Jan 16 2024

#Input:
#  myspecies=c("Chinook.salmon_yearling_Interior_Sp","California.market.squid__", "Sablefish__", "Sea.nettle__" , "Water.jelly__")
# myspecies=c("CanMagM","LarCop","EpacAd","TspinAd")#,"TspinFur","TspinJuv","EpacFur","EpacJuv")
#   spp=1
#   load(file=paste0("01_tidy_data/",myspecies[spp],"_dat.test.Rdata"),verbose=T) # data formatted with format_raw_data.R
#Output:
#write.csv(index,file=paste0("Output data for Noble/",myspecies[spp],".yeareffects.m1a.csv"))


#summary
#m1b bestmodel for market squid, but none of the other species converged
#Sablefish__ m1a failed to converge -- no models converged, might have done something wrong in the set up. Should repeat
#might need coarser mesh because of fewer datapoints?

#CanMag m1a passed all checks, almost the lowest AIC (m1b)
#EpacAd m1a did not converge 

#"CanMagM"
# df      AIC
# m1a 30 868.5553
# m1b 31 865.3884 
# m1c 30 934.0190
# best m1a, passed m1b.remlFALSE sanity check, m1b failed

# run sdmTMB with JSOES data extract from Cheryl

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
library(distancetocoast)
library(viridis)
library(patchwork)                  # Load patchwork package



# Example model ------------


## 01 load data ---------------------------------------------------------------
#spp=3;
myspecies[spp]
load(file=paste0("01_tidy_data/",myspecies[spp],"_dat.test.Rdata"),verbose=T) # data formatted with format_raw_data.R

dat <- df %>% 
  dplyr::rename(density = 'numPkm2') %>% 
  dplyr::select(density,yday,month,week,year,year_f,depth,scale_depth,
                scale_dist2shore, Lon.km, Lat.km)

names(dat)
mysum<-aggregate(dat$density,list(dat$year),sum)
names(mysum)<-c("Y",myspecies[spp]);mysum
#write.csv(mysum,file=paste0("Output data for Noble/",myspecies[spp],".yearsum.csv"))

## 02 fit model ---------------------------------------------------------------

# To generate a standardized index of abundance


# fit a model

# start simple using year (factor) and scaled depth as a smooth spline (you could also try depth as a linear variable)
# start with just the spatiotemporal random field on
# we can later test if there is more support for a different random effects structure (spatial, temporal, spatiotemporal)
#sdmTMB prefers km (bspde), plotting needs meters (bspde_1000)
#reml is false (default), that's for comparing different fixed effects with the same random effects,
#reml = true is when fixed effects are the same, but different random effects

m1a <- sdmTMB(
  density ~ 0 +     # zero gets rid of the intercept
    year_f + 
    s(scale_depth), # we may want to replace with, or include separately, distance to shore - talk to Brian Burke and Dave Huff - they have ideas about this
  #  (scale_depth), # we may want to replace with, or include separately, distance to shore - talk to Brian Burke and Dave Huff - they have ideas about this
  mesh = bspde,
  time = "year",
  # name of time variable
  spatial = "off",
  # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "iid",
  family = tweedie(),
  # could also be one of the delta_ families
  data = dat,
  # reml = TRUE,      # we're going to compare different random effects structures using the same fixed effects
  silent = FALSE
)


fit <- m1a # year as factor, smoothed depth, spatial off, iid on
#Sablefish__ m1a failed to converge


## 03 check fit ---------------------------------------------------------------
#matern effect is the distance at which data are not correlated (in km)
#if you don't get all green checks in sanity check, try a different model

fit
sanity(fit)
AIC(fit) # AIC = 4555.311 (reml FALSE)

myspecies[spp]

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
resolution <- 2000 # poly is in meters; 2000 is pretty fine: 2000m x 2000m = 4 km2
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


# If you decide to use bottom depth, you'll need to grab depth data for each grid cell
# Get bathymetric data for each grid cell, bounding box = -126,42,-123,49
JSOES_bathy = marmap::as.raster(getNOAA.bathy(-126, -123, 49, 42, res=0.25, keep=TRUE)) # resolution is in minutes
plot(JSOES_bathy)
st_crs(JSOES_bathy)

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

# Here, we're grabbing distance to shore for each grid cell (the shapefile is going to be in UTM (meters))
# library(distancetocoast)
# library(viridis)

# grab raster
dist <- crop(distance_to_coastline_10, extent(-126, -123, 43, 49)) # -126, -124, 43, 49
plot(dist,col = viridis::viridis(64))

st_crs(dist)

dist_proj <- raster::projectRaster(dist, crs = 4326)

# calculate mean distance to shore for each grid cell
dist_mean <- raster::extract(dist_proj,        # raster layer
                             grid_sdf,         # SPDF with centroids for buffer
                             fun=mean,         # what value to extract
                             df=TRUE)          # return a dataframe? 

dist_mean.b <- dist_mean %>% 
  mutate(distance = layer/1000) 

# bind distance field with grid data frame
new_df.distance <- cbind(dist_mean.b,grid)

# now bind with new_df.depth

new_df.all <- new_df.distance %>% 
  left_join(new_df.depth %>%
              dplyr::select(ID,depth))


# now expand our data frame to include all covariates in our model (year, day, and
# lat/lon in kilometers) - this is a necessary step for prediction

# replication factor  
yrs <- sort(as.vector(unique(df$year))) 
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, new_df.all, simplify = FALSE)) 

# add year covariate
new_df$year_f <-  as.factor(rep(yrs, each = nrow(new_df.all)))
new_df$year <-  rep(yrs, each = nrow(new_df.all))

# add lat/lon in km
new_df$Lon.km <- new_df$x / 1000
new_df$Lat.km <- new_df$y / 1000

# create new area field in km^2
new_df$area_km2 <- new_df$area / 1e06


# let's also add a day column to predict on -- we don't currently have yday in the model,
# but later (below) we will end up testing a model with yday as a covariate 
# pick any yday within the range of existing data (day 166 is roughly June 15th)

new_df$yday <- 166 


# grab mean and stdev of distance column in dat (survey loc distances)
# use these to scale prediction grid depths
mu.dist <- mean(df$dist)     # [1] 27.48336
sd.dist <- sd(df$dist)       # [1] 17.09894


# do the same for depth, if necessary
mu.depth <- mean(df$depth)     # [1] 114.44797
sd.depth <- sd(df$depth)       # [1] 117.99652

# do the same for yday
mu.yday <- mean(df$yday)     # [1] 164.95731
sd.yday <- sd(df$yday)       # [1] 14.42126

new_df <- new_df %>% 
  dplyr::filter(depth > 0,
                distance > 0) %>% 
  mutate(scale_dist2shore = (distance - mu.dist) / sd.dist,    # use same name as covariate used in model
         scale_depth = (depth - mu.depth) / sd.depth,
         scale_yday = (yday - mu.yday) / sd.yday)



## 05 load prediction grid and calculate index----------------------------------------------------


# make predictions: return_tmb_object = TRUE -- standardized index of abundance per year; for pretty maps, use return_tmb_object = FALSE
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# now make index
# create a vector of grid cell areas for accurate interpolation, otherwise, 
# get_index() assumes all grid cell areas are equal
area_vec <- new_df$area_km2
index = get_index(pred.index, bias_correct = TRUE, area = area_vec)
index$Species<-myspecies[spp]
index$Model<-"m1a"
index$Node<-"Condition2"
m1a.yeareffects<-index
m1a.yeareffects

#write.csv(index,file="Output data for Noble/IntSpr.yeareffects.m1a.csv")
write.csv(m1a.yeareffects,file=paste0("Output data for Noble/",myspecies[spp],".yeareffects.m1a.csv"))

## 06 plot temporal trends ---------------------------

#library(patchwork)                  # Load patchwork package

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
#   ggtitle(myspecies[spp])+
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

Y <- (p1 + p2) +
  plot_annotation(
  title = paste('m1a',myspecies[spp]))


Y

ggsave(filename = paste0('03_plots/SIA_',myspecies[spp],"2x2Km_m1a.png"),
       dpi = 600,
       width=25,
       height=15,
       units="cm")

# so, now that we've seen the step by step process, we can test different
# random effects structures. Once we've determined the best one, we can use that same
# structure to test different fixed effects models (adding in other covariates)


