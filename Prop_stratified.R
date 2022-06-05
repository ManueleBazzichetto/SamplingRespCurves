library(ecospat)
library(classInt)
#library(rgdal) #this is to generate points according to different sampling strategies
#but I need to find an alternative as gdal would retire at the end of 2023

#Generate different sampling

Bio1.rcl <- ecospat.rcls.grd(Chelsa.stack$Bio1, no.classes = 3)
Bio12.rcl <- ecospat.rcls.grd(Chelsa.stack$Bio12, no.classes = 3)
Bio1_12.rcl <- Bio1.rcl + Bio12.rcl

#try to get a better reclassification of raster layers for stratification (based on quartiles)
#this time the stratification focuses only on the AOI
quantile(Chelsa.AOI$Bio1)

Bio1.rcl.qrt <- reclassify(Chelsa.AOI$Bio1, rcl = matrix(c(-0.36, 7.35, 1,
                                                             7.35, 9.55, 2,
                                                             9.55, 11.25, 3,
                                                             11.25, 14.26, 4), nrow = 4,
                                                           byrow = T))

quantile(Chelsa.AOI$Bio12)

Bio12.rcl.qrt <- reclassify(Chelsa.AOI$Bio12, rcl = matrix(c(631, 924.6, 1,
                                                             924.6, 1055.5, 2,
                                                             1055.5, 1169.8, 3,
                                                             1169.8, 1565, 4), nrow = 4,
                                                           byrow = T))

#multiply by 10 to have unique combo of categories
Bio12.rcl.qrt <- Bio12.rcl.qrt*10

plot(Bio1.rcl.qrt + Bio12.rcl.qrt)

par(mfrow=c(1, 2))
plot(Bio1_12.rcl); plot(Bio1.rcl.qrt + Bio12.rcl.qrt)

#create layer
Bio1_12.rcl.qrt.AOI <- (Bio1.rcl.qrt + Bio12.rcl.qrt)

freq(Bio1_12.rcl)
freq(Bio1.rcl.qrt)
freq(Bio12.rcl.qrt)
freq(Bio1_12.rcl.qrt.AOI)

#see histogram
hist(Bio1_12.rcl, breaks = 4, col = heat.colors(maxValue(Bio1_12.rcl) - minValue(Bio1_12.rcl)))
plot(Bio1_12.rcl, col = rainbow(4), asp = 1)

hist(Bio1_12.rcl.qrt.AOI, breaks = 40, col = heat.colors(maxValue(Bio1_12.rcl.qrt.AOI) - minValue(Bio1_12.rcl.qrt.AOI)))
plot(Bio1_12.rcl.qrt.AOI, col = rainbow(40), asp = 1)

#transform raster to polygon to see boundary of strata
Contour.rcl.qrt.AOI <- rasterToContour(Bio1_12.rcl.qrt.AOI)
Contour.rcl.qrt.AOI <- st_as_sf(Contour.rcl.qrt.AOI)

plot(Bio1_12.rcl.qrt.AOI)
plot(st_geometry(Contour.rcl.qrt.AOI), add = T)

#there's also an ext arg to limit the extent for generating points
#plot(raster::sampleRandom(x = Bio1_12.rcl, size = 100, xy = T, sp = T))

#plot(raster::sampleRegular(x = Bio1_12.rcl, size = 100, ext = extent(st_bbox(Abruzzo)), xy = T, sp = T))

#plot(raster::sampleStratified(x = Bio1_12.rcl, size = 10, xy = T, sp = T))

#ncells for each stratum (we will use this info to weight the n° of points for each stratum for stratifies sampling)
raster::freq(Bio1_12.rcl)

#trials for the function
Bio1_12.rcl.df <- na.omit(as.data.frame(Bio1_12.rcl, xy = T))
table(Bio1_12.rcl.df$layer)

#new strat layer
Bio1_12.rcl.qrt.AOI.df <- as.data.frame(Bio1_12.rcl.qrt.AOI, xy = T, na.rm = T)
table(Bio1_12.rcl.qrt.AOI.df$layer)

Strat_raster <- function(x, N) {
  Strata_ncell <- raster::freq(x)
  Strata_ncell <- Strata_ncell[!is.na(Strata_ncell)[, 1], ]
  Strata_ncell <- setNames(Strata_ncell[, 2], as.character(Strata_ncell[, 1]))
  Strata_ncell <- Strata_ncell/sum(Strata_ncell)
  Rast_df <- na.omit(as.data.frame(x, xy = T))
  Rast_df$layer <- as.character(Rast_df$layer)
  Points_df <- do.call(rbind, lapply(names(Strata_ncell), function(nm) {
    Subs <- Rast_df[Rast_df$layer == nm, ]
    Subs <- Subs[sample(x = nrow(Subs), size = floor(Strata_ncell[[nm]]*N), replace = F), ]
    return(Subs)
  }))
  rownames(Points_df) <- seq_len(nrow(Points_df))
  Points_df$layer <- NULL
  return(Points_df)
}

#little less points are sampled because of floor within the function - but we need to round
Example_500 <- replicate(n = 4, Strat_raster(Bio1_12.rcl.qrt.AOI, 500), simplify = F)
Example_500 <- do.call(rbind, lapply(seq_along(Example_500), function(.) {
  df <- Example_500[[.]]
  df$trial <- .
  return(df)
  }))

Example_500 <- st_as_sf(Example_500, coords = c("x", "y"))
st_crs(Example_500) <- 4326

plot(Bio1_12.rcl.qrt.AOI, asp = 1)
plot(as(Example_500[Example_500$trial == 4, ], "Spatial"), add = T, asp = 1)
