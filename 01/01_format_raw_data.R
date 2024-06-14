# 26 October 2023
# Exploring Chinook salmon yearling data, 1998-2020

# Note from L. Crozier
# "hatcheries release fall yearlings, and there is also a portion of the Snake 
# River fall Chinook ESU that holds over until the following spring to migrate, 
# so they "naturally" go out as yearlings. I put naturally in quotes because we 
# think that is a response to the dams (probably Dworshak dam cooling the rearing 
# areas, rather than Lower Granite reservoir, which is what I used to think 
# caused that behavior). So it does require genetics to separate the yearlings, 
# but there are many more spring than fall. And almost all of the subs are fall"

# install INLA, fmesher
options(repos=c(
  inlabruorg = "https://inlabru-org.r-universe.dev",
  INLA = "https://inla.r-inla-download.org/R/testing",
  CRAN = "https://cran.rstudio.com"))
install.packages("fmesher")
# 
 install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
 install.packages("tidyverse")
 
# load libraries -----------------------------------------------------------
library(INLA)
library(fmesher)
library(inlabru)
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


# 01 load data ---------------------------------------------------------------

# JSOES trawl data for Chinook yearlings - Rusicka data file
rawdat <- read.csv(file="Trawl_Erica_Mason_11.27.23.csv", header=T)

str(rawdat)


# 02 format data ----------------------------------------------------------

d <- rawdat %>% 
  dplyr::select(Sample.Date,Station.Code,Transect.Code,Start.Lat..decimal.degrees.,
         Start.Long..decimal.degrees.,Depth..m.,population.density..num.km.2.) %>% 
  rename(date_ch = 'Sample.Date',
         latdd = 'Start.Lat..decimal.degrees.',
         londd = 'Start.Long..decimal.degrees.',
         depth_m = 'Depth..m.',
         numPkm2 = 'population.density..num.km.2.') %>% 
  mutate(date = lubridate::mdy(date_ch),
         week = epiweek(date),
         yday = yday(date),
         yday2 = (yday(date))^2,
         year = year(date),
         month = month(date))

# check for NAs -- sdmTMB does not allow for NAs in covariate data (however,
# year can be a smoothed covariate if missing some years)

which(is.na(d$date)) == TRUE
which(is.na(d$depth_m)) == TRUE     # 4 missing depths
which(is.na(d$numPkm2)) == TRUE
which(is.na(d$latdd)) == TRUE
which(is.na(d$londd)) == TRUE


# format GPS coordinates for sdmTMB

# project our spatial data:
# WGS84 (EPSG: 4326) is the CRS commonly used by organizations that provide 
# GIS data for the entire globe or many countries. CRS used by Google Earth. 
# So, if you have lat/long GPS data, WGS84 is the right CPS to use. 

# Go get missing depth data
# Create coordinate df and convert into an sf object and assign CRS
# coords for rows with missing depths

coords <- d %>% 
  filter(is.na(depth_m)==TRUE) %>% 
  dplyr::select(date,latdd,londd)

coords$ID <- seq(1,nrow(coords),1)

check <- d %>% 
  filter(is.na(depth_m)==FALSE)

coords_sdf <- st_as_sf(coords, coords = c("londd", "latdd"), crs = "epsg:4326")

st_crs(coords_sdf) #check that CRS was properly assigned

# Get bathymetric data
range(d$latdd)     # 43.73462 48.34717
range(d$londd)     # -125.4235 -123.9738
JSOES_bathy = marmap::as.raster(getNOAA.bathy(-125, -123, 49, 43, res=0.25, keep=TRUE))
plot(JSOES_bathy)
st_crs(JSOES_bathy)

library(raster)
bathy <- projectRaster(JSOES_bathy, crs = "epsg:4326")
st_crs(bathy)

depth_max <- raster::extract(bathy,           # raster layer
                             coords_sdf,      # SPDF with centroids for buffer
                             # buffer = ,      # buffer size, units depend on CRS
                             fun=max,         # what value to extract
                             df=TRUE)         # return a dataframe? 


depth_max.a <- depth_max %>% 
  mutate(btm_depth = layer * -1)

depth_max.b <- left_join(depth_max.a,coords) # grab associated data in d

dat.test <- d %>% 
  left_join(depth_max.b %>%
              dplyr::select(date,btm_depth,latdd,londd))

dat.test.update <- dat.test %>%              # create a new column depth w/o missing data
  mutate(depth = case_when(depth_m >0 ~ depth_m,
                           TRUE ~ btm_depth)) 

# join this table with d
d.updated <- dat.test.update %>% 
  left_join(d %>%
              dplyr::select(year,month,week,yday,date,
                            depth_m,latdd,londd,numPkm2))


# format GPS locations for sdmTMB (need to be in meters)

# first we need to convert lat lon to UTMs
# WGS84 (EPSG: 4326) is the CRS commonly used by organizations that provide 
# GIS data for the entire globe or many countries. CRS used by Google Earth. 
# So, if you have lat/long GPS data, WGS84 is the right CRS to use. 
# To convert to UTMs, then we need WGS 84 / UTM zone 10N (CRS = 32610)
utm_crs <- get_crs(d.updated, c("londd", "latdd")) # check that

d.updated <- add_utm_columns(
  d.updated,
  ll_names = c("londd", "latdd"),
  ll_crs = 4326,
  utm_names = c("Lon.km", "Lat.km"),
  utm_crs = utm_crs,
  units = c("km")
)

# add a year as factor column and scale the depth column

df <- d.updated %>%
  arrange(year,month,yday) %>% 
  mutate(
    year_f = as.factor(year),
    scale_depth = scale(as.numeric(depth))[ , 1])

head(df)

df$rownum <- seq(1,nrow(df),1)

write.csv(df,file="01_tidy_data/df_yearlings.csv")


df %>% 
  group_by(year) %>% 
  summarise(n=n(),min= min(numPkm2), max= max(numPkm2), mean = mean(numPkm2))


# 03 make a mesh -------------------------------------------------------------

# make a mesh -- which is approximating the spatial surface
# note that the data frame 'd' already contains lat and lon in UTMs -- which is what
# we want to use

mesh = make_mesh(df, xy_cols = c("Lon.km","Lat.km"), cutoff = 10, type = "kmeans")

# smaller cutoff values lead to more complexity -- and we can check the number of knots,
# or random effects being estimated with the following code
mesh$mesh$n
plot(mesh)


# code below is used to account for any overlap of mesh with land

library(rnaturalearth)

map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "United States of America")

# Crop the polygon for plotting and efficiency:
st_bbox(map_data) # find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -123, ymin = 43, xmax = -125, ymax = 49)))) # -125, -123, 49, 43,

na_coast_proj <- sf::st_transform(na_coast, crs = 32610)

# Project our survey data coordinates:
survey <- df %>% dplyr::select(londd, latdd, numPkm2) %>%
  st_as_sf(crs = 4326, coords = c("londd", "latdd")) %>% 
  st_transform(32610)

# Plot our coast and survey data:
ggplot(na_coast_proj) +
  geom_sf() +
  geom_sf(data = survey, size = 0.5) 

# Add on the barrier mesh component:
# Note that the original mesh was created using km coordinates,
# thus, the proj_scaling needs to be x1000 to convert to meters

library(sdmTMBextra)

bspde <- add_barrier_mesh(
  mesh, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

# In the above, the grey dots are the centre of triangles that are in the
# ocean. The red crosses are centres of triangles that are over land. The
# spatial range will be assumed to be 0.1 (`range_fraction`) over land compared
# to over water.

# We can make a more advanced plot if we want:

# for plotting in UTMs
df$X1000 <- df$Lon.km * 1000
df$Y1000 <- df$Lat.km * 1000

mesh_1000 <- make_mesh(df, xy_cols = c("X1000", "Y1000"), cutoff = 10 * 1000)

bspde_1000 <- add_barrier_mesh(
  mesh_1000, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1, plot = TRUE
)

mesh_df_water <- bspde_1000$mesh_sf[bspde$normal_triangles, ]
mesh_df_land <- bspde_1000$mesh_sf[bspde$barrier_triangles, ]

library(inlabru)      # to convert mesh to spatial feature

blues <- RColorBrewer::brewer.pal(5, "Blues")

ggplot(na_coast_proj) +
  geom_sf(fill = "antique white") +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green") +
  inlabru::gg(bspde_1000$mesh) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.key = element_rect(fill = "aliceblue")) +
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(breaks = seq(-126, -123, by = 1)) +
  ggtitle('sdmTMB Mesh for JSOES Trawls',subtitle = 'Yearling')

ggsave(filename = "03_plots/Yearling_mesh.png",
       dpi = 600,
       width=25,
       height=15,
       units="cm")

ggplot(na_coast_proj) +
  geom_sf(fill = "antique white") +
  inlabru::gg(bspde_1000$mesh) +
  geom_point(data=df,alpha = 0.8, aes(Lon.km * 1000, Lat.km * 1000, size = numPkm2),colour = blues[4])+
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.key = element_rect(fill = "aliceblue")) +
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(breaks = seq(-126, -123, by = 1)) +
  labs(subtitle = "Yearlings")+
  guides(size= guide_legend(title="Standardized density\n(no./km^2)"))

ggsave(filename = "03_plots/Yearling_densities.png",
       dpi = 600,
       width=25,
       height=15,
       units="cm")

ggplot(na_coast_proj) +
  geom_sf(fill = "antique white") +
  inlabru::gg(bspde_1000$mesh) +
  geom_point(data=df,alpha = 0.8, aes(Lon.km * 1000, Lat.km * 1000, size = numPkm2,group=year),colour = blues[4])+
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.key = element_rect(fill = "aliceblue")) +
  facet_wrap(~ year)+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(breaks = seq(-126, -123, by = 1)) +
  labs(subtitle = "JSOES Yearlings")+
  guides(size= guide_legend(title="Standardized density\n(no./km^2)"))

ggsave(filename = "03_plots/Yearling_densitiesByYr.png",
       dpi = 600,
       width=25,
       height=15,
       units="cm")

save(na_coast_proj,mesh_df_land,mesh_df_water,bspde,bspde_1000,df,file="01_tidy_data/yearling_dat.Rdata")


