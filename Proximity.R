library(osmdata)
library(rgdal) #for raster projection

##------------------------------------------get openstreetmap data


#first generate overpass query (to be able to download osmdata)
#overpass queries need a bbox
st_bbox(Abruzzo)

Over_q <- opq(bbox = c(13.01199, 41.68691, 14.76341, 42.90156))

available_features()

#Here!! Get only major roads, otherwise gets complicated to work with all kind of roads (include all but tertiary)
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

#then dissolve all roads in 1 feature in QGIS (as on R it takes ages)
#and call back the dissolved object

Abr_highway <- st_read("Abruzzo_roads_diss.shp")

Abr_highway$name <- NULL
Abr_highway$osm_id <- 1

mapview(Abr_highway)

#project the shp on planar CRS (to compute distances)
Abr_highway.proj <- st_transform(Abr_highway, crs = 32633)

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
Abr_highway_AOI <- mask(Abr_highway_dist.geo.prb, mask = as(Elev_AOI, "Spatial"))
plot(Abr_highway_AOI)

#check how points would be sampled
plot(Abr_highway_AOI)
Prox_samling <- function(prox_lay, n = 300) {
  Prox_df <- na.omit(as.data.frame(prox_lay, xy = T))
  Prox_df <- Prox_df[sample(nrow(Prox_df), size = n, replace = F, prob = Prox_df$layer), ] #check
  Prox_df <- st_as_sf(Prox_df, coords = c("x", "y"))
  return(Prox_df)
  }

Check_prox.sampl <- do.call(rbind, replicate(n = 10, Prox_samling(prox_lay = Abr_highway_AOI, n = 50), simplify = F))
Check_prox.sampl$ID <- as.character(rep(seq_len(10), each = 50))
st_crs(Check_prox.sampl) <- 4326

#for plotting
Abr_highway_AOI.df <- as.data.frame(Abr_highway_AOI, xy = T)

ggplot() +
  geom_raster(data = Abr_highway_AOI.df, aes(x = x, y = y, fill = layer), alpha = .5) +
  geom_sf(data = Elev_AOI, color = "white", fill = NA) +
  geom_sf(data = Check_prox.sampl, aes(color = ID)) +
  coord_sf(xlim = c(13.3, 14), ylim = c(41.8, 42.4)) +
  scale_fill_viridis_c() +
  scale_color_discrete() +
  labs(title = NULL, x = "Longitude", y = "Latitude", fill = "Probability")

hist(Check_prox.sampl[Check_prox.sampl$ID == "1", ][["layer"]])

#check correlation between probability of sampling a cell and finding a species (compute joint prob)
all.equal(coordinates(D.sperandii.prob), coordinates(Abr_highway_AOI)) #different number of coordinates!
Cor_dist_prob <- na.omit(data.frame(extract(stack(D.sperandii.prob, D.sperandii.bin),
                                            coordinates(D.sperandii.prob.AOI), df = T),
                               extract(stack(Abr_highway_dist.geo, Abr_highway_AOI),
                                       coordinates(D.sperandii.prob.AOI), df = T), 
                               coordinates(D.sperandii.prob.AOI)))

colnames(Cor_dist_prob)[c(2, 3, 5, 6)] <- c("Proba", "PA", "RoadDist", "ProbSamplCell")

Cor_dist_prob$ID.1 <- NULL

with(Cor_dist_prob, plot(Proba, RoadDist))
with(Cor_dist_prob, plot(Proba, ProbSamplCell))
with(Cor_dist_prob, cor(Proba, RoadDist)) #Pearson: 0.42
with(Cor_dist_prob, cor(Proba, ProbSamplCell)) #Pearson: -0.38
hist(with(Cor_dist_prob, Proba*ProbSamplCell))
summary(with(Cor_dist_prob, Proba*ProbSamplCell))

#get coefficient relationship P.occ vs road dist
summary(glm(PA ~ RoadDist, family = binomial, data = Cor_dist_prob))
(exp(0.22689) - 1)*100

#do the same for the rare species
Cor_dist_prob.rare <- na.omit(data.frame(extract(stack(D.tundrae.prob, D.tundrae.bin),
                                                 coordinates(D.tundrae.bin.AOI), df = T),
                                    extract(stack(Abr_highway_dist.geo, Abr_highway_AOI),
                                            coordinates(D.tundrae.bin.AOI), df = T), 
                              coordinates(D.tundrae.bin.AOI)))

colnames(Cor_dist_prob.rare)[c(2, 3, 5, 6)] <- c("Proba", "PA", "RoadDist", "ProbSamplCell")

Cor_dist_prob.rare$ID.1 <- NULL

with(Cor_dist_prob.rare, plot(Proba, RoadDist))
with(Cor_dist_prob.rare, cor(Proba, RoadDist)) #0.1199693
#not sure relationship is monotonic
with(Cor_dist_prob.rare, cor(Proba, RoadDist, method = "spearman")) #0.2145394
with(Cor_dist_prob.rare, cor(Proba, ProbSamplCell)) #-0.1342271
with(Cor_dist_prob.rare, plot(Proba, ProbSamplCell))
hist(with(Cor_dist_prob.rare, Proba*ProbSamplCell))
summary(with(Cor_dist_prob.rare, Proba*ProbSamplCell))

#very low correlation 
plot(Abr_highway_AOI)

#need to create another Abr_highway_AOI for the rare species which mimics an increasing sampling effort 
#(sampling done at longer distances from the road)

#I can reduce some distance (e.g. by some percentage of decease) to simulate cells are closer to the road
#try reducing distances by 20%
#so, for example, if a cell is at 6km, it has same prob of being sample as if it was at exp(-1*6*(0.8))
Abr_highway_dist.geo.rare <- (Abr_highway_dist.geo * .8)

Abr_dist.geo.rare.prb <- calc(Abr_highway_dist.geo.rare, fun = function(x) exp(-x))

par(mfrow = c(1, 2))
plot(Abr_highway_dist.geo.prb)
plot(Abr_dist.geo.rare.prb)

#comparison 
#sequence of distance in km
#increase of 100 m
Km_seq <- seq(cellStats(Abr_highway_dist.geo, "min"), cellStats(Abr_highway_dist.geo, "max"), .1)

plot(exp(-1*Km_seq*(0.8)) ~ exp(-1*Km_seq), ylab = "Reduced", xlab = "True dist")
abline(a = 0, b = 1)

#cut Abr_dist.geo.rare.prb for AOI
Abr_highway_AOI.rare <- mask(Abr_dist.geo.rare.prb, mask = Elev_AOI)

