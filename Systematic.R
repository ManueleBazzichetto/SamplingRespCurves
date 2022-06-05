#Systematic sampling within the geographic space
library(sf)

#Example on how to use st_sample from sf
#regular works assuming CRS is planar, better to project geometries
Elev_AOI.proj <- st_transform(x = Elev_AOI, crs = 32633)
plot(st_geometry(Elev_AOI.proj))
plot(st_sample(x = Elev_AOI.proj, size = 200, type = "regular"), add = T)



