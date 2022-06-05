#Sample cells with high standard deviation of topographic vars
#compute topographic variables using raster:terrain
library(raster)
library(mapview)
library(ggpubr)
library(sf)

#Elevation_Abr: source should be SRTM, and spatial res at 45° lat ~ 100 m
Slope_Abr <- raster::terrain(x = Elevation_Abr, opt = "slope", unit = "degrees", neighbors = 8)
#aspect computed in radians to be used for northn and eastn
Aspect_Abr <- raster::terrain(x = Elevation_Abr, opt = "aspect", unit = "radians", neighbors = 8)
N_Abr <- cos(Aspect_Abr) 
E_Abr <- sin(Aspect_Abr)

Topogr_Abr <- stack(Elevation_Abr, Slope_Abr, N_Abr, E_Abr)

names(Topogr_Abr) <- c("Elev", "Slope", "N", "E")

Topogr_var_Abr <- stack(lapply(1:nlayers(Topogr_Abr), function(i) {
  z_sc_ly <- calc(x = Topogr_Abr[[i]], fun = function(.) . - cellStats(Topogr_Abr[[i]], "mean"))
  z_sc_ly <- calc(x = z_sc_ly, fun = function(.) ./cellStats(Topogr_Abr[[i]], "sd"))
  #res <- focal(z_sc_ly, w = matrix(1, nrow = 5, ncol = 5), fun = sd, na.rm = T)
  #here use aggregate instead of focal, so that the resolution matches the same of the bioclim layers
  res <- aggregate(z_sc_ly, fact = 14, fun = sd, na.rm = T)
  return(res)
  })
  )

names(Topogr_var_Abr) <- c("Elev_sd", "Slope_sd", "N_sd", "E_sd")

Topogr_het_Abr <- sum(Topogr_var_Abr)

plot(Topogr_het_Abr, col = heat.colors(200))

#Now, I need to set to NAs all cell values with sd lower than a given value
#first reduce the extent to the AOI
Topogr_het_Abr.AOI <- raster::mask(Topogr_het_Abr, mask = as(Elev_AOI, "Spatial"))

plot(Topogr_het_Abr.AOI)

#let's say, I retain all locations having value above the 2nd quartile - the distr is pretty symmetric
raster::hist(Topogr_het_Abr.AOI)
quantile(getValues(Topogr_het_Abr.AOI), na.rm = T) #2.0307561

Topogr_het_Abr.AOI[Topogr_het_Abr.AOI < 2.0307561] <- NA

#check -> ok
cellStats(Topogr_het_Abr.AOI, "min")
cellStats(Topogr_het_Abr.AOI, "countNA") #so there are 31524 - 25788: 5736 non NAs

mapview(Topogr_het_Abr.AOI)
#par(mfrow = c(1,2))
#plot(Elevation_Abr, asp = 1); plot(Topogr_het_Abr.AOI, asp = 1)

#get an idea of how good is Topogr_het_Abr.AOI to represent highly heterog places
Coords_to_use <- na.omit(data.frame(extract(Topogr_het_Abr.AOI, coordinates(Topogr_het_Abr.AOI), df = T),
                   coordinates(Topogr_het_Abr.AOI)))[c("x", "y")]

Topogr_het.df <- na.omit(data.frame(extract(Topogr_Abr, Coords_to_use, df = T), Coords_to_use))

Topogr_het.df$ID <- NULL

Topogr_het.plot <- data.frame(Val = c(Topogr_het.df$Elev, Topogr_het.df$Slope,
                                    Topogr_het.df$N, Topogr_het.df$E),
                            Var = rep(names(Topogr_Abr),
                                      each = nrow(Topogr_het.df)))

ggarrange(plotlist = lapply(unique(Topogr_het.plot$Var), function(vars) {
  ggplot(Topogr_het.plot[Topogr_het.plot$Var == vars, ], aes(x = Var, y = Val)) +
    geom_boxplot() +
    ylab(NULL) + xlab(NULL) +
    theme_minimal()
}), nrow = 2, ncol = 2)

Topogr_het.sp <- st_as_sf(Topogr_het.df, coords = c("x", "y"))

Topogr_het.sp <- st_set_crs(Topogr_het.sp, value = 4326)

mapview(Elevation_Abr) + mapview(Topogr_het.sp)

#-------------------Alt transects..this is super slow
#See links below to have an idea on how to delineate transects
#https://www.rdocumentation.org/packages/raster/versions/3.5-15/topics/adjacent
#https://stackoverflow.com/questions/69916871/find-8-neighbors-of-a-point-in-a-raster-in-r

#Altitudinal transects
#library(MBHdesign)
#library(MASS) #for data("volcano")
#library(fields) #for image.plot
#data("volcano")

#Elev_Abr_AOI <- as.matrix(Elevation_Abr)

#Elev_Abr_AOI[Elev_Abr_AOI < 600 | Elev_Abr_AOI > 1800] <- NA

#great
#range(Elev_Abr_AOI[!is.na(Elev_Abr_AOI)])

#
#num_lon <- nrow(Elev_Abr_AOI)
#num_lat <- ncol(Elev_Abr_AOI)
#image.plot(x = seq_len(num_lon), y = seq_len(num_lat), z = Elev_Abr_AOI,
#           main = "Mountain Height (m)", asp = 1)

#
#pot.sites <- expand.grid(x = 1:num_lon, y = 1:num_lat)
#pot.sites$alt <- as.vector(Elev_Abr_AOI)

#
#alt.tr.control <- list(transect.pattern = "line", transect.nPts = 10,
#                     line.length = 7, nRotate = 11, mc.cores = 6)

#alt.tr.constraints <- findDescendingTrans(
#  potential.sites = pot.sites[,c("x","y")], bathy = pot.sites$height,
#  in.area = rep(TRUE, nrow(pot.sites)), control = alt.tr.control)


