##Code to replicate results presented in the manuscript "Sampling strategies matter to accurately estimate response curves’ parameters in species distribution models"
##Extraction and processing of altitude and climatic data


library(elevatr) #to get elevation data
library(rnaturalearth)
library(raster)
library(sf)
library(car)
library(mapview)
library(ggplot2)
library(ggpubr)
library(rcartocolor) #palettes
library(parallel) #for uesampling and Unif_sampl (https://jonlefcheck.net/2013/06/13/using-parallel-processing-in-r/)
#also install rnaturalearthhires from github
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)


##GET ALTITUDE AND CLIMATIC DATA FOR THE STUDY AREA------------------------------------------------------------------------------------------------

#Get regional boundaries of Abruzzo 
Italy <- ne_states(country = "italy", returnclass = "sf")
Abruzzo <- Italy[Italy$region == "Abruzzo", ]

#Get elevation data from: https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution
#Latitude should be approx 45°
#z = 10 -> z (zoom) and latitude determine resolution of the elevation layer
#see attr(Elevation_Abr, which = "sources") to get sources for altitude data
#the function does the job, but it also gives weird internal errors..
Elevation_Abr <- elevatr::get_elev_raster(locations = Abruzzo, z = 10, clip = "locations") #elevations in meters
#this next step is probably not needed
Elevation_Abr <- crop(Elevation_Abr, Abruzzo)

#See where we have huge holes according to the DEM (#weird -600 m elevation..)
#Holes_Abr object was deleted to save space
Holes_Abr <- reclassify(Elevation_Abr, rcl = matrix(c(-600, 0, 1, 
                                                      0, 2833, 0), byrow = T, nrow = 2))
mapview(Holes_Abr) #holes are at the seaside, and we don't need this area

#Set elevation values < 0 to 0
Elevation_Abr[Elevation_Abr < 0] <- 0

Elevation_Abr.df <- as.data.frame(Elevation_Abr, xy = T)
Elevation_Abr.df <- na.omit(Elevation_Abr.df)
colnames(Elevation_Abr.df)[3] <- "Elev"

#Get a view of elevation in Abruzzo
#copied and pasted from https://rspatialdata.github.io/elevation.html
ggplot() +
  geom_raster(data = Elevation_Abr.df, aes(x = x, y = y, fill = Elev)) +
  geom_sf(data = Abruzzo, color = "white", fill = NA) +
  coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
  scale_fill_viridis_c() +
  labs(title = "Elevation in Abruzzo", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")

#load CHELSA data
#downloaded from https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2F
#on 16/03/2022 (first time)
Chelsa_bio1 <- raster("~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio1_1981-2010_V.2.1.tif")
Chelsa_bio12 <- raster("~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio12_1981-2010_V.2.1.tif")

#the following CHELSA stack was created to avoid sampling NAs cells when performing the systematic approach
#unfortunately, there's a similar issue with the layers of the species pres./abs. So the problem is only partly solved
#by the use of this stack (see below - systematic sampling).
#For this reason, for the systematic (and for the topographic approach) we sample a (slightly) larger
#number of units than the desired sampling effort. Excess points (if any) are then removed to match the desired sampling effort. 

#add half degree of longitude to xmin and xmax (st_bbox(Abruzzo)), and half degree to lat (ymin)
Chelsa.stack.syst <- crop(stack(Chelsa_bio1, Chelsa_bio12),
                          y = extent(13.01199-.5, 14.76341+.5, 41.68691-.5, 42.90156))

#stack bioclimatic layers
Chelsa.stack <- stack(Chelsa_bio1, Chelsa_bio12)

#crop and mask them to only cover Abruzzo
Chelsa.stack <- crop(Chelsa.stack, y = extent(st_bbox(Abruzzo)))
Chelsa.stack <- mask(Chelsa.stack, mask = as(Abruzzo, "Spatial"))

#set names of Chelsa.stack.syst
names(Chelsa.stack) <- c("Bio1", "Bio12")

#no need to keep the global Bio1 and Bio12
rm(Chelsa_bio1, Chelsa_bio12)

#check correlation between Bio1 and Bio12
layerStats(Chelsa.stack, stat = "pearson", na.rm = T)

#set names of Chelsa.stack.syst
names(Chelsa.stack.syst) <- names(Chelsa.stack)


##DELINEATE AREA OF INTEREST (AOI - all lands between approx. 500 and 1800 m asl)-----------------------------------------------------------------------------

#Delineate AOI: all land between 500 and 1800 
Elev_AOI <- reclassify(Elevation_Abr, rcl = matrix(c(0, 600, NA,
                                                     600, 1800, 1,
                                                     1800, 2900, NA),
                                                   nrow = 3, byrow = T),
                       include.lowest = T, right = F)

#aggregate some cells to lower resolution so it takes less to coerce to poly
6156680/(25^2)

Elev_AOI <- raster::aggregate(Elev_AOI, fact = 25)

Elev_AOI <- rasterToPolygons(Elev_AOI, fun = function(x) {x==1}, dissolve = T)

Elev_AOI <- st_as_sf(Elev_AOI)

#get proportions of cells within 500 and 1800 m asl
quantile(extract(Elevation_Abr, as(Elev_AOI, "Spatial"), df = T)[[2]], probs = seq(0, 1, 0.025), na.rm = T)

#Plot AOI on top of elevation
ggplot() +
  geom_raster(data = Elevation_Abr.df, aes(x = x, y = y, fill = Elev)) +
  geom_sf(data = Elev_AOI, color = "white", fill = NA) +
  coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
  scale_fill_viridis_c() +
  labs(title = "Elevation in Abruzzo", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")

#Get climatic data for AOI 
Chelsa.AOI <- mask(Chelsa.stack, mask = as(Elev_AOI, "Spatial"))

#Sequences of Bio1 and Bio12 values. These will be used for generating true response curves for the species in the simulations
Bio1_seq.AOI <- seq(from = minValue(Chelsa.AOI$Bio1), to = maxValue(Chelsa.AOI$Bio1), length.out = 100)
Bio12_seq.AOI <- seq(from = minValue(Chelsa.AOI$Bio12), to = maxValue(Chelsa.AOI$Bio12), length.out = 100)