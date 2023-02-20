##Code to replicate results presented in the manuscript "Sampling strategies matter to accurately estimate response curves’ parameters in species distribution models"
##Data and (R) functions for implementing the following sampling approaches: stratified, uniform, topographic, proximity-to-road
##Random and systematic approaches are directly implemented in the simulations' script

##!!Please, load packages reported at the top of the "Bio_Elev_data" script!!

#STRATIFIED sampling--------------------------------------------------------------------------------------------------------

#The layer used for the (proportional) stratified approach is based on the quartiles of the bioclimatic variables (Bio1, Bio12)
#The layer is created for the area of interest

#get quartiles for Bio1
quantile(Chelsa.AOI$Bio1)

#reclassify layer using quartiles
Bio1.rcl.qrt <- reclassify(Chelsa.AOI$Bio1, rcl = matrix(c(-0.36, 7.35, 1,
                                                           7.35, 9.55, 2,
                                                           9.55, 11.25, 3,
                                                           11.25, 14.26, 4), nrow = 4,
                                                         byrow = T))
#get quartiles for Bio12
quantile(Chelsa.AOI$Bio12)

#reclassify layer using quartiles
Bio12.rcl.qrt <- reclassify(Chelsa.AOI$Bio12, rcl = matrix(c(631, 924.6, 1,
                                                             924.6, 1055.5, 2,
                                                             1055.5, 1169.8, 3,
                                                             1169.8, 1565, 4), nrow = 4,
                                                           byrow = T))

#multiply Bio12.rcl.qrt by 10 to have a unique combination of categories with Bio1
Bio12.rcl.qrt <- Bio12.rcl.qrt*10

#create stratified layer
Bio1_12.rcl.qrt.AOI <- (Bio1.rcl.qrt + Bio12.rcl.qrt)

#function used to perform the stratified sampling in the simulations. The function is internally called
#by the Strat_fit function for D. sperandii, and analogous function for D. tundrae
Strat_raster <- function(x, N) {
  Strata_ncell <- raster::freq(x)
  Strata_ncell <- Strata_ncell[!is.na(Strata_ncell[, 1]), ]
  Strata_names <- as.character(Strata_ncell[, 1])
  Strata_ncell <- cumsum(Strata_ncell[, 2]/sum(Strata_ncell[, 2]))
  Npts <- 0
  for (i in seq_along(Strata_ncell)){
    Npts[i] <- floor(N*Strata_ncell[i]) - sum(Npts)
  }
  Strata_ncell <- setNames(Npts, Strata_names)
  Rast_df <- as.data.frame(x, xy = T, na.rm = T)
  Rast_df$layer <- as.character(Rast_df$layer)
  Points_df <- do.call(rbind, lapply(names(Strata_ncell), function(nm) {
    Subs <- Rast_df[Rast_df$layer == nm, ]
    Subs <- Subs[sample(x = nrow(Subs), size = Strata_ncell[[nm]], replace = F), ]
    return(Subs)
  }))
  rownames(Points_df) <- seq_len(nrow(Points_df))
  Points_df$layer <- NULL
  return(Points_df)
}

#ROAD PROXIMITY sampling--------------------------------------------------------------------------------------------------------

library(osmdata) #get road data for Abruzzo
library(rgdal) #for raster projection - note that rgdal will retire by the end of 2023 (use the 'project' function of the terra package to project raster layers)

#first generate overpass query (to be able to download osmdata)
#overpass queries need a bbox
st_bbox(Abruzzo)

Over_q <- opq(bbox = c(13.01199, 41.68691, 14.76341, 42.90156))

available_features()

#Get only major roads, otherwise gets complicated to work with all kind of roads (include all but tertiary)
#this excludes links and other particular features
Abr_highway <- add_osm_features(opq = Over_q, features = c("\"highway\"=\"motorway\"",
                                                           "\"highway\"=\"trunk\"",
                                                           "\"highway\"=\"primary\"",
                                                           "\"highway\"=\"secondary\""))
Abr_highway <- osmdata_sf(Abr_highway)

names(Abr_highway)

lapply(Abr_highway[4:8], nrow) #just keep lines!

plot(st_geometry(Abr_highway$osm_lines))

Abr_highway <- Abr_highway$osm_lines

Abr_highway <- Abr_highway[c("osm_id", "name", "geometry")]

#write lines on disk and then call them back as an object to save memory
st_write(Abr_highway, dsn = "Abruzzo_roads.shp")

#then dissolve all roads in 1 feature in QGIS (as this step in R takes too long)
#and call back the dissolved object

Abr_highway <- st_read("Abruzzo_roads_diss.shp")

Abr_highway$name <- NULL
Abr_highway$osm_id <- 1

mapview(Abr_highway)

#project the shp on planar CRS (to compute distances)
Abr_highway.proj <- st_transform(Abr_highway, crs = 32633)

#project layer (D.sperandii.bin is created in the Dianthus_sperandii script) on which distances from roads will be computed
Base_layer.dist <- projectRaster(D.sperandii.bin, crs = crs("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"))

Base_layer.dist[] <- NA

Abr_highway_dist <- raster::rasterize(x = Abr_highway.proj, Base_layer.dist)

Abr_highway_dist <- raster::distance(Abr_highway_dist)

#first assign non-zero val to cells with 0 distance from roads - I set 1 (m), so cells at 0 distance are at 1
Abr_highway_dist[Abr_highway_dist == 0] <- 1

#re-project it and mask it - mask is done at this stage to avoid inflating distances at the boundaries
Abr_highway_dist.geo <- projectRaster(Abr_highway_dist, crs = crs(Chelsa.stack$Bio1))
Abr_highway_dist.geo <- mask(Abr_highway_dist.geo, mask = as(Abruzzo, "Spatial"))

#transform it in km
Abr_highway_dist.geo <- Abr_highway_dist.geo/1000

hist(Abr_highway_dist.geo)

#p decreases as negative exponential - distances close to 0 have prob close to 1
Abr_highway_dist.geo.prb <- calc(Abr_highway_dist.geo, fun = function(x) exp(-x))

mapview(Abr_highway_dist.geo.prb)

Abr_highway_dist.geo.prb.df <- as.data.frame(Abr_highway_dist.geo.prb, xy = T)

colnames(Abr_highway_dist.geo.prb.df)[3] <- "Road_distance"

ggplot() +
  geom_raster(data = Abr_highway_dist.geo.prb.df, aes(x = x, y = y, fill = Road_distance)) +
  geom_sf(data = Abruzzo, color = "white", fill = NA) +
  coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
  scale_fill_viridis_c() +
  labs(title = "Probability of sampling a cell", x = "Longitude", y = "Latitude", fill = "Prob.")

hist(Abr_highway_dist.geo.prb.df$Road_distance)

#mask Abr_highway_dist.geo.prb to match AOI
#this is the layer that will be used in the simulations for the road proximity sampling
#specifically, it will be used in the function Proximity_fit for D. sperandii, and in the analog function for D. tundrae
Abr_highway_AOI <- mask(Abr_highway_dist.geo.prb, mask = as(Elev_AOI, "Spatial"))
plot(Abr_highway_AOI)

#UNIFORM sampling--------------------------------------------------------------------------------------------------------

#Function used (in Unif_sampl) to uniformly sample points within the environmental space
#This function is also used in a parallel research project. Here, I only use 'part of it', i.e. basically the part working with the first three arguments of the function
uesampling2.0 <- function(sdf, grid.res, n.tr = 5, n.prev = NULL, sub.ts = FALSE, n.ts = 5, plot_proc = FALSE, print_proc = FALSE) {
  if(!require(sf)) install.packages('sf')
  if(!require(assertive)) install.packages('assertive')
  if(!(all(st_is(sdf, "POINT")))) {
    stop("sdf object must have a sf POINT GEOMETRY class")
  }
  if(!is_numeric(n.tr)) stop(paste(n.tr, "is not of class 'numeric'.", sep = " "))
  if(!is_logical(plot_proc)) stop("plot_proc is not of class 'logical'; it has class 'numeric'.")
  grid <- st_make_grid(sdf, n = grid.res)
  sdf$ID <- row.names(sdf)
  res <- do.call(rbind, lapply(seq_len(length(grid)), function(i) {
    if(isTRUE(print_proc)) message(paste("Processing tile", i, sep = " "))
    if(isTRUE(plot_proc)) {
      if(i == 1) {
        plot(grid, border = "black")
        plot(grid[i], col = "green", add = TRUE)
      } else {
        plot(grid[i], col = "green", add = TRUE)
      }
    }
    subs <- sdf[grid[i], ]
    if(nrow(subs) <= n.tr) {
      return(subs)
    } else {
      subs <- subs[sample(nrow(subs), n.tr, replace = FALSE), ]
      return(subs)
    }
  }))
  #dded to account for prevalence set by the user
  if(!is.null(n.prev)) {
    if(nrow(res) < n.prev) {
      n.cell <- length(grid) #removed seq_len()
      dif.abs <- n.prev - nrow(res) #n° of missing absences to reach prevalence
      #there are less missing absences than grid cells
      if(dif.abs < length(grid)) {
        Subs <- sdf[!(sdf$ID %in% res$ID), ] #this select all remaining absences 
        while(dif.abs > 0) { #as long as there are missing absences for reaching prevalence....
          ID_cell <- sample(n.cell, size = 1) #take a cell sampling its ID
          abs <- Subs[grid[ID_cell], ] #extract obs from the sdf (Subs) corresponding to the ID
          if(nrow(abs) == 0) next #if the resulting subsetted sdf is empty, try again (next)
          abs <- abs[sample(nrow(abs), 1), ] #otherwise get one obs from it (from the ones remaining that are not included in res)
          absID <- abs$ID #get the obs ID
          dif.abs <- (dif.abs - 1) #update number of missing absences
          Subs <- Subs[!(Subs$ID %in% absID), ] #exlcude the sampled absences from Subs - to avoid pseudoreplicates
          res <- rbind(res, abs) #add the new row to res
          if(nrow(Subs) == 0) { #but if there are not obs left in Subs..the exit
            message("There are not points left in the environmental space to reach prevalence")
            break #check the break
          }
        }
      } else {
        #in this case there are more absences to select than cells
        Subs <- sdf[!(sdf$ID %in% res$ID), ] #avoid pseudoreplicates (as above)
        #dif.abs <- n.prev - nrow(res) #useless
        ratio <- floor(n.prev/length(grid)) #if there are more missing absences than cells..then how many complete run over the grid can we complete?
        while(dif.abs > 0) { #as long as dif.abs is not 0 go on
          if(ratio == 0) { #but if you reach the max number of complete iter. just stop..there will be few absences left to reach the prevalence
            message("There are not points left in the environmental space to reach prevalence")
            break
          }
          new_set <- do.call(rbind, lapply(seq_len(n.cell), function(.) { #here, basically, add new run to the one done by the function above (before considering exact prevalence)
            subs <- Subs[grid[.], ]
            if(nrow(subs) == 0) {
              return(subs)
            } else {
              subs <- subs[sample(nrow(subs), 1), ]
              return(subs)
            }
          }))
          Subs <- Subs[!(Subs$ID %in% new_set$ID), ]
          res <- rbind(res, new_set)
          ratio <- (ratio - 1)
          dif.abs <- dif.abs - nrow(new.set) #update dif absences
        }
      }
    } else {
      message("Increase prevalence by randomly excluding some absence from the dataset")
    }
  }
  if(isTRUE(sub.ts)) {
    abs_val <- sdf[!(sdf$ID %in% res$ID), ]
    res_val <- do.call(rbind, lapply(rev(seq_len(length(grid))), function(i) {
      if(isTRUE(print_proc)) message(paste("Processing tile", i, sep = " "))
      if(isTRUE(plot_proc)) plot(grid[i], col = "red", add = TRUE)
      subs <- abs_val[grid[i], ]
      if(nrow(subs) <= n.ts) {
        return(subs)
      } else {
        subs <- subs[sample(nrow(subs), n.ts, replace = FALSE), ]
        return(subs)
      }
    }))
    return(list(Bkg.tr = res, Bkg.ts = res_val))
  } else {
    return(res)
  }
}

#Function internally calling "uesampling2.0". Notice that the n.tr argument is set so that an higher number
#of points are sampled within the environmental space than the number set by the user. The sample size is then adjusted
#internally in Unif_sampl to return the exact number of points to match the sampling effort.
Unif_sampl <- function(x, N, rsl) {
  Un_smp <- uesampling2.0(sdf = x, grid.res = rsl, n.tr = switch(as.character(N),
                                                                 '200' = 3, '250' = 4,
                                                                 '300' = 5, '350' = 6,
                                                                 '400' = 7, '450' = 7,
                                                                 '500' = 8))
  Un_smp.dim <- nrow(Un_smp)
  if(isTRUE(Un_smp.dim > N)) {
    Un_smp <- Un_smp[-sample(x = Un_smp.dim, size = (Un_smp.dim - N), replace = F), ]
  } 
  st_geometry(Un_smp) <- NULL
  Un_smp <- Un_smp[c("x", "y")]
  return(Un_smp)
}

#Chelsa layer used for the uniform sampling
Chelsa.AOI.df <- na.omit(as.data.frame(Chelsa.AOI, xy = T))
#add scaled Bio1 and Bio12 for environmental space -> scaled Bio1 and Bio12 will be the coordinates of 
#the raster cells in the environmental space
Chelsa.AOI.df$Bio1.sc <- with(Chelsa.AOI.df, scale(Bio1))
Chelsa.AOI.df$Bio12.sc <- with(Chelsa.AOI.df, scale(Bio12))
Chelsa.AOI.df.sp <- st_as_sf(Chelsa.AOI.df, coords = c("Bio1.sc", "Bio12.sc"))

#SYSTEMATIC sampling--------------------------------------------------------------------------------------------------------

#For the systematic sampling, the Elev_AOI.proj is created as this will be systematically sampled
#in the Systematic_fit function for D. sperandii, and analogous function for D. tundrae
Elev_AOI.proj <- st_transform(x = Elev_AOI, crs = 32633)

#TOPOGRAPHIC sampling--------------------------------------------------------------------------------------------------------

#Derive layers describing topography in Abruzzo
#Elevation_Abr: source data can be accessed through attr(Elevation_Abr, which = "sources")
Slope_Abr <- raster::terrain(x = Elevation_Abr, opt = "slope", unit = "degrees", neighbors = 8)
#aspect computed in radians to be used for northn and eastn
Aspect_Abr <- raster::terrain(x = Elevation_Abr, opt = "aspect", unit = "radians", neighbors = 8)
N_Abr <- cos(Aspect_Abr) 
E_Abr <- sin(Aspect_Abr)

Topogr_Abr <- stack(Elevation_Abr, Slope_Abr, N_Abr, E_Abr)

names(Topogr_Abr) <- c("Elev", "Slope", "N", "E")

#check this to make the two raster have same resolution and match as much as possible bioclimatic layers
#first topographic layers are reported on the same scale by standardization
#then the sd of aggregating cells is computed as a measure of cell specific heterogeneity
#(in an area of size comparable to the cell size of the CHELSA bioclimatic layers)
Topogr_var_Abr <- stack(lapply(1:nlayers(Topogr_Abr), function(i) {
  z_sc_ly <- calc(x = Topogr_Abr[[i]], fun = function(.) . - cellStats(Topogr_Abr[[i]], "mean"))
  z_sc_ly <- calc(x = z_sc_ly, fun = function(.) ./cellStats(Topogr_Abr[[i]], "sd"))
  res <- aggregate(z_sc_ly, fact = 14, fun = sd, na.rm = T)
  return(res)
  })
  )

names(Topogr_var_Abr) <- c("Elev_sd", "Slope_sd", "N_sd", "E_sd")

#The layers are summed up
Topogr_het_Abr <- sum(Topogr_var_Abr)

plot(Topogr_het_Abr, col = heat.colors(200))

##added on January 9/01/2023 - Topogr_het_Abr is resampled to spatially match Chelsa.stack
#remaining code (below) was updated accordingly
Topogr_het_Abr <- resample(Topogr_het_Abr, Chelsa.stack$Bio1, method = "bilinear") #change Rev

#Now, I need to set to NAs all cell values with sd lower than a given value
#first reduce the extent to the AOI
Topogr_het_Abr.AOI <- raster::mask(Topogr_het_Abr, mask = as(Elev_AOI, "Spatial"))
plot(Topogr_het_Abr.AOI)

#I retain all locations having the sum of the standard deviations above the 2nd quartile - notice that the distribution is pretty symmetric
raster::hist(Topogr_het_Abr.AOI)
quantile(getValues(Topogr_het_Abr.AOI), na.rm = T) #2.0257340 #change Rev

#set all values below the median as NAs, and use this layer in the Topo_fit function for D.sperandii,
#and analogous function for D. tundrae
Topogr_het_Abr.AOI[Topogr_het_Abr.AOI < 2.0257340] <- NA
