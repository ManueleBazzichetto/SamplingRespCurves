library(elevatr) #to get elevation data
library(rnaturalearth)
library(raster)
#library(terra) #for the moment, I need it for the stratified sampling with weights
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

#tinytex::reinstall_tinytex() -> tinytex updated on 21/04/2022

##GET ALTITUDE AND CLIMATIC DATA FOR THE STUDY AREA------------------------------------------------------------

#Get regional boundaries of Abruzzo 
Italy <- ne_states(country = "italy", returnclass = "sf")
Abruzzo <- Italy[Italy$region == "Abruzzo", ]

#Get elevation data from: https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution
#Latitude should be approx 45°
#z = 10 should download SRTM
#the function does the job, but it also gives weird internal errors..
Elevation_Abr <- elevatr::get_elev_raster(locations = Abruzzo, z = 10, clip = "locations") #elevations in meters
#this next step is in the tutorial, but I don't think we need it
Elevation_Abr <- crop(Elevation_Abr, Abruzzo)

#Holes_Abr was deleted to save space
#weird -600 m elevation..
Holes_Abr <- reclassify(Elevation_Abr, rcl = matrix(c(-600, 0, 1, 
                                                      0, 2833, 0), byrow = T, nrow = 2))
mapview(Holes_Abr) #holes are at the seaside..I don't need this area..

#Set elevation values < 0 to 0
Elevation_Abr[Elevation_Abr < 0] <- 0

#nice way to convert raster to dataframe (without using RStoolbox::fortify)
Elevation_Abr.df <- as.data.frame(Elevation_Abr, xy = T)
Elevation_Abr.df <- na.omit(Elevation_Abr.df)
colnames(Elevation_Abr.df)[3] <- "Elev"

#copied and pasted from https://rspatialdata.github.io/elevation.html
ggplot() +
  geom_raster(data = Elevation_Abr.df, aes(x = x, y = y, fill = Elev)) +
  geom_sf(data = Abruzzo, color = "white", fill = NA) +
  coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
  scale_fill_viridis_c() +
  labs(title = "Elevation in Abruzzo", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")

#load CHELSA data
#downloaded from https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2F
#on 16/03/2022
Chelsa_bio1 <- raster("~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio1_1981-2010_V.2.1.tif")
Chelsa_bio12 <- raster("~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio12_1981-2010_V.2.1.tif")

Chelsa.stack <- stack(Chelsa_bio1, Chelsa_bio12)

Chelsa.stack <- crop(Chelsa.stack, y = extent(st_bbox(Abruzzo)))
Chelsa.stack <- mask(Chelsa.stack, mask = as(Abruzzo, "Spatial"))

names(Chelsa.stack) <- c("Bio1", "Bio12")

#no need to keep them
rm(Chelsa_bio1, Chelsa_bio12)

layerStats(Chelsa.stack, stat = "pearson", na.rm = T)

#Original CHELSA rasters were deleted (they were ~ 1 GB) and rasters covering Abruzzo were written
#for code reproducibility
writeRaster(Chelsa.stack$Bio1, filename = "~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio1_1981-2010_V.2.1_Abr.tif")
writeRaster(Chelsa.stack$Bio12, filename = "~/Documents/SpRespCurve/CHELSA_data/CHELSA_bio12_1981-2010_V.2.1_Abr.tif")

#load WorldClim data
load("WorldC.RData")

BioClimateData.Eu <- BioClimateData.Eu[[c(1, 12)]]

BioClimate_Abr <- raster::crop(BioClimateData.Eu, y = extent(st_bbox(Abruzzo)))
BioClimate_Abr <- raster::mask(BioClimate_Abr, mask = as(Abruzzo, "Spatial"))
BioClimate_Abr$bio1 <- BioClimate_Abr$bio1/10

rm(BioClimateData.Eu)

mapview(BioClimate_Abr$bio1)
mapview(BioClimate_Abr$bio12)

plot(getValues(BioClimate_Abr$bio1) ~ getValues(BioClimate_Abr$bio12), xlab = "Prec", ylab = "Temp")
layerStats(BioClimate_Abr, stat = "pearson", na.rm = T)

##DELINEATE AREA OF INTEREST (AOI)-----------------------------------------------------------------------------

#Delineate AOI: all land between 600 and 1800 
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

#Plot AOI on top of elevation
ggplot() +
  geom_raster(data = Elevation_Abr.df, aes(x = x, y = y, fill = Elev)) +
  geom_sf(data = Elev_AOI, color = "white", fill = NA) +
  coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
  scale_fill_viridis_c() +
  labs(title = "Elevation in Abruzzo", x = "Longitude", y = "Latitude", fill = "Elevation (meters)")

#Get an idea of the bioclimatic characteristics of AOI
BioClimate_Abr.df <- na.omit(as.data.frame(mask(BioClimate_Abr, mask = as(Elev_AOI, "Spatial")))) #BioClimate_Abr was already clipped at this stage..

ggplot(BioClimate_Abr.df, aes(y = bio1)) +
  geom_boxplot() + ylab("Bio1") +
  theme(axis.text.x.bottom = element_blank())

ggplot(BioClimate_Abr.df, aes(y = bio12)) +
  geom_boxplot() + ylab("Bio12") +
  theme(axis.text.x.bottom = element_blank())

##CREATE DIANTHUS SPERANDII (WIDE DISTRIBUTION)----------------------------------------------------------------
#A 3-POINTS APPROACH IS USED TO CREATE DIANTHUS TUNDRAE (SEE RARE_SPECIES SCRIPT)

#Create species -> Dianthus sperandii
#Use CHELSA

#optimum around 6/7 C°
#stationarity point at: b1/(2*-1*b2) = 8
#at the max/min the resp is: (4*interc*b1 - b2^2)/4*b2 
plot(-2 + 0.6*getValues(Chelsa.stack$Bio1) - 0.045*(getValues(Chelsa.stack$Bio1)^2) ~ getValues(Chelsa.stack$Bio1))
plot(plogis(-2 + 0.6*getValues(Chelsa.stack$Bio1) - 0.045*(getValues(Chelsa.stack$Bio1)^2)) ~ getValues(Chelsa.stack$Bio1))

#precipitation - prob occurrence increases with precipitation
plot(plogis(-6 + 0.0044*getValues(Chelsa.stack$Bio12)) ~ getValues(Chelsa.stack$Bio12))

#both - summing the two intercepts
plot(-8 + 0.6*getValues(Chelsa.stack$Bio1) - 0.045*(getValues(Chelsa.stack$Bio1)^2) + 0.0044*getValues(Chelsa.stack$Bio12), 
     plogis(-8 + 0.6*getValues(Chelsa.stack$Bio1) - 0.045*(getValues(Chelsa.stack$Bio1)^2) + 0.0044*getValues(Chelsa.stack$Bio12)))

True_coef <- c(-8, .6, -.045, .0044)
True_coef.df <- data.frame(True_coef, Coef = c("Int", "Temp", "Temp^2", "Pr"))

#layers
D.sperandii.link <- -8 + 0.6*Chelsa.stack$Bio1 - .045*(Chelsa.stack$Bio1^2) + .0044*Chelsa.stack$Bio12
D.sperandii.prob <- calc(D.sperandii.link, plogis)
#warnings are for NAs
D.sperandii.bin <- calc(D.sperandii.prob, fun = function(x) {rbinom(1, 1, x)})

mapview(D.sperandii.bin)

#Dfs for ggplotting
D.sperandii.prob.df <- na.omit(as.data.frame(D.sperandii.prob, xy = T))
colnames(D.sperandii.prob.df)[3] <- "Prob"

D.sperandii.bin.df <- na.omit(as.data.frame(D.sperandii.bin, xy = T))
colnames(D.sperandii.bin.df)[3] <- "PA"

ggarrange(ggplot() +
            geom_raster(data = D.sperandii.prob.df, aes(x = x, y = y, fill = Prob)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_viridis_c() +
            labs(title = "Dianthus sperandii", x = "Longitude", y = "Latitude", fill = "Probability"),
          ggplot() +
            geom_raster(data = D.sperandii.bin.df, aes(x = x, y = y, fill = PA)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_viridis_c() +
            labs(title = "Dianthus sperandii", x = "Longitude", y = "Latitude", fill = "PresAbs"), nrow = 1)

##CLIMATIC DATA ONLY FOR AREA OF INTEREST----------------------------------------------------------------------

#create layers considering only the AOI (600-1800 m), to simulate a more realistic sampling design
D.sperandii.prob.AOI <- mask(D.sperandii.prob, mask = as(Elev_AOI, "Spatial"))

D.sperandii.bin.AOI <- mask(D.sperandii.bin, mask = as(Elev_AOI, "Spatial"))
Chelsa.AOI <- mask(Chelsa.stack, mask = as(Elev_AOI, "Spatial"))

#sampling effort
Sampl_effort <- seq(from = 200, to = 500, by = 50)

##We focus only on the AOI!!!!

##SIMULATIONS OF THE DIFFERENT SAMPLING STRATEGIES

####Simulate random sampling
#get coordinates for no NA cells

Climate_crds <- na.omit(as.data.frame(Chelsa.AOI, xy = T))[c("x", "y")]

Random_fit <- function(y, x, x_crds, n = 300, min_p = 30) {
  npres <- T
  while(npres) {
    Random_plots <- x_crds[sample(nrow(x_crds), n, replace = F), ]
    Fake_df <- na.omit(extract(stack(y, x), Random_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
}

#rerun
#Random_res <- replicate(n = 500, expr = Random_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI, n = 250), simplify = F)
#Random_res.mat <- do.call(rbind, lapply(Random_res, '[[', 1))
#Random_res.coln <- colnames(Random_res.mat)
#dim(Random_res.mat) <- NULL
#Random_res.df <- data.frame(Val = Random_res.mat, Coef = rep(Random_res.coln, each = 500))
#ggplot(Random_res.df, aes(y = Val, col = Coef)) +
#  geom_boxplot() +
#  facet_wrap(~ Coef, scales = "free")
#Compute Bias
#tapply(Random_res.df$Val, Random_res.df$Coef, mean) - True_coef

#What happens if we change the number of plots?
#sampling: notice that with ncell too low, the alghoritm may not converge -> also try to introduce more stochasticity
round(Sampl_effort/length(na.omit(getValues(Chelsa.AOI$Bio1))), digits = 4)*100 #1.82 2.27 2.73 3.18 3.63 4.09 4.54(% of the total area)

Random_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Random_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI, 
                                              x_crds =  Climate_crds, n = N, min_p = 30), simplify = F)
  return(res)
})

names(Random_res) <- paste("N", Sampl_effort, sep = "_")

Random_mats <- lapply(Random_res, function(.) {
  do.call(rbind, lapply(., '[[', 1))
  })

#extract cor
Random_cor <- sapply(Random_res, function(.) {
  mean(vapply(., function(i) {
    cor_val <- i[[2]]
    return(cor_val)
  }, FUN.VALUE = numeric(1)))
})

Random_mats <- do.call(rbind, lapply(names(Random_mats), function(nm) {
  df <- Random_mats[[nm]]
  df_nm <- colnames(df)
  dim(df) <- NULL
  df <- data.frame(Val = df, Coef = rep(df_nm, each = 500), N = nm)
  return(df)
  })
  )

True_coef_nm <- setNames(True_coef, nm = unique(Random_mats$Coef))

Random_mats$True_val <- unname(True_coef_nm[Random_mats$Coef])
Random_mats$Coef <- factor(Random_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Random_mats$N <- factor(Random_mats$N, levels = paste("N", Sampl_effort, sep = "_"))

ggplot(Random_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Add stratified
#here we should use some sampling proportional to the area of the stratum
Strat_fit <- function(y, x, n = 300, strata, min_p = 30) {
  npres <- T
  while(npres) {
    Strat_plots <- Strat_raster(x = strata, N = n)
    Fake_df <- na.omit(extract(stack(y, x), Strat_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
}

#this was old layer for stratif
Bio1_12.rcl.AOI <- mask(Bio1_12.rcl, mask = Elev_AOI)
#use Bio1_12.rcl.qrt.AOI

Strat_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, Strat_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI, n = N,
                                      strata = Bio1_12.rcl.qrt.AOI, min_p = 30), simplify = F)
  return(res)
  })

names(Strat_res) <- names(Random_res)

#extract cor
Strat_cor <- sapply(Strat_res, function(.) {
  mean(vapply(., '[[', 2, FUN.VALUE = numeric(1)))
})

Strat_mats <- do.call(rbind, lapply(names(Strat_res), function(nm) {
  Mat <- do.call(rbind, lapply(Strat_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
  }))

Strat_mats$Coef <- factor(Strat_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Strat_mats$True_val <- unname(True_coef_nm[as.character(Strat_mats$Coef)])
Strat_mats$N <- factor(Strat_mats$N, levels = paste("N", Sampl_effort, sep = "_"))

ggplot(Strat_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Add sampling close to paths
Proximity_fit <- function(y, x, prox_lay, n = 300, min_p = 30) {
  Prox_df.full <- na.omit(as.data.frame(prox_lay, xy = T)) #replaced Prox_df with Prox_df.full to avoid overwriting in while
  npres <- T
  while(npres) {
    Prox_df <- Prox_df.full[sample(nrow(Prox_df.full), size = n, replace = F, prob = Prox_df.full$layer), ] #check
    Prox_df$layer <- NULL
    Fake_df <- na.omit(extract(stack(y, x), Prox_df, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
}

Prox_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Proximity_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                 prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
  })

names(Prox_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor
Prox_cor <- sapply(Prox_res, function(.) {
  mean(vapply(., '[[', 2, FUN.VALUE = numeric(1)))
})

Prox_res <- do.call(rbind, lapply(names(Prox_res), function(nm) {
  Mat <- do.call(rbind, lapply(Prox_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef_nm[Mat$Coef]
  return(Mat)
  }))

Prox_res$N <- factor(Prox_res$N, levels = paste("N", Sampl_effort, sep = "_"))
Prox_res$Coef <- factor(Prox_res$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

ggplot(Prox_res, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

##Add uniform sampling
Uniform_fit <- function(y, x, x_sdf, n = 300, rsl, min_p = 30) {
  npres <- T
  while(npres) {
    Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
    Fake_df <- na.omit(extract(stack(y, x), Uniform_points, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
}

#go parallel
cr7 <- parallel::makeCluster(7)

parallel::clusterExport(cr7, c("D.sperandii.bin.AOI", "Chelsa.AOI",
                               "Chelsa.AOI.df.sp", "Uniform_fit", 
                               "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr7, list(library(sf), library(raster)))

Unif_res <- parallel::parLapply(cl = cr7, X = Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Uniform_fit(y = D.sperandii.bin.AOI,
                                               x = Chelsa.AOI,
                                               x_sdf = Chelsa.AOI.df.sp,
                                               n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

stopCluster(cr7)

names(Unif_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor -> uniformly sampling the env space reduces correlation
Unif_cor <- sapply(Unif_res, function(.) {
  mean(vapply(., '[[', 2, FUN.VALUE = numeric(1)))
})

Unif_res <- do.call(rbind, lapply(names(Unif_res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef_nm[Mat$Coef]
  return(Mat)
}))

Unif_res$N <- factor(Unif_res$N, levels = paste("N", Sampl_effort, sep = "_"))
Unif_res$Coef <- factor(Unif_res$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

ggplot(Unif_res, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

##Add systematic sampling
#systematic sampling may under-sampling (sample less than N points)
#write something to estimate P of undersampling -> no, it's actually useless -
#the regular approach will always arrange the points in the same way..no variability
#but there actually is some variability..
#sum(sapply(seq_len(1000), function(.) {
#  Pts_syst <- st_sample(x = Elev_AOI.proj, size = 500, type = "regular")
#  isTRUE(nrow(Pts_syst) < 500)
#  }))

Systematic_fit <- function(y, x, N, poly_proj, min_p = 30) {
  npres <- T
  while(npres) {
    Pts_syst <- st_sample(x = poly_proj, size = N, type = "regular")
    Pts_syst <- st_transform(Pts_syst, crs = 4326)
    Pts_syst <- st_coordinates(Pts_syst)
    if(nrow(Pts_syst) > N) {
      Pts_syst <- Pts_syst[-sample(nrow(Pts_syst), size = (nrow(Pts_syst) - N), replace = F), ]
    }
    Fake_df <- na.omit(extract(stack(y, x), Pts_syst, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
}

Syst_res <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Systematic_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                           N = n, poly_proj = Elev_AOI.proj, min_p = 30), simplify = F)
  return(res)
})

#check variability
#hist((do.call(rbind, lapply(Syst_res[[3]], '[[', 1)))[, 1])

names(Syst_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor
Syst_cor <- sapply(Syst_res, function(.) mean(vapply(., function(i) i[[2]], FUN.VALUE = numeric(1))))

Syst_res <- do.call(rbind, lapply(names(Syst_res), function(nm) {
  Mat <- do.call(rbind, lapply(Syst_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef_nm[Mat$Coef]
  return(Mat)
}))

Syst_res$N <- factor(Syst_res$N, levels = paste("N", Sampl_effort, sep = "_"))
Syst_res$Coef <- factor(Syst_res$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

ggplot(Syst_res, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Add sampling along 'topographic' transects
Topo_fit <- function(y, x, N, topo_layer, min_p = 30) {
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T)
  Topo_lyr.df <- Topo_lyr.df[!is.na(Topo_lyr.df$layer), ]
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N, replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v))
  }

Topo_res <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 500, expr = Topo_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                            N = n., topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
})

names(Topo_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor
Topo_cor <- sapply(Topo_res, function(.) mean(vapply(., function(i) i[[2]], FUN.VALUE = numeric(1))))

Topo_res <- do.call(rbind, lapply(names(Topo_res), function(nm) {
  Mat <- do.call(rbind, lapply(Topo_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef_nm[Mat$Coef]
  return(Mat)
}))

Topo_res$N <- factor(Topo_res$N, levels = paste("N", Sampl_effort, sep = "_"))
Topo_res$Coef <- factor(Topo_res$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

ggplot(Topo_res, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

##GET AND COMPARE RESULTS--------------------------------------------------------------------------------------

##compare results
Big_comp <- rbind(data.frame(Random_mats, Type = "Random"),
                  data.frame(Strat_mats, Type = "Strat"),
                  data.frame(Prox_res, Type = "Prox"),
                  data.frame(Unif_res, Type = "Uniform"),
                  data.frame(Syst_res, Type = "Systematic"),
                  data.frame(Topo_res, Type = "Topographic"))

Big_comp$Type <- factor(Big_comp$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                  "Systematic", "Topographic"))

#Uniform seems to be quite biased..
ggplot(Big_comp, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#see shape of response curve for very bad params
#Big_comp.bad <- subset(Big_comp, Type == "Random" & N == "N_200" & Coef == "(Intercept)")
#Big_comp.bad[which.min(Big_comp.bad$Val), ]

#compute estimate of mse = bias^2 + var
#3-d array
Mse_samp <- function(df, Coef., N., Type., True_vals) {
  Vals <- df[df$Coef == Coef. & df$N == N. & df$Type == Type., "Val"]
  Mse <- mean((Vals - True_vals[Coef.])^2)
  return(Mse)
}

Mse_samp(df = Big_comp, Coef. = "(Intercept)", N. = "N_200", Type. = "Random", True_vals = True_coef_nm)

Mse_list <- lapply(as.character(unique(Big_comp$Type)), function(ty) {
  mapply(function(x, y) Mse_samp(df = Big_comp, Coef. = x, N. = y, Type. = ty, True_vals = True_coef_nm),
         x = rep(names(True_coef_nm), 7), y = rep(paste("N", Sampl_effort, sep = "_"), each = 4))
  })

names(Mse_list) <- as.character(unique(Big_comp$Type))

Mse_df <- do.call(rbind, lapply(names(Mse_list), function(nm) {
  nms <- names(Mse_list[[nm]])
  df <- data.frame(Mse_val = unname(Mse_list[[nm]]), Coef = nms,
                   N = rep(paste("N", Sampl_effort, sep = "_"), each = 4), 
                   Type = nm)
  return(df)
  }))

Mse_df$Coef <- factor(Mse_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Mse_df$Type <- factor(Mse_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                              "Systematic", "Topographic"))

ggarrange(plotlist = lapply(unique(Mse_df$Type), function(ty) {
  ggplot(data = Mse_df[Mse_df$Type == ty, ], aes(x = N, y = Mse_val)) +
  geom_point(aes(col = Coef)) +
  facet_wrap(~ Coef, nrow = 1, ncol = 5, scales = "free_y") +
    xlab(NULL) + ylab("Mse") + ggtitle(ty) +
    theme_pubclean() +
    theme(legend.position = "none")
  }
  ), ncol = 1, nrow = 6)

#ok we see something here!
#at low sample sizes, performing uniform sampling within the env space is better (high precision offsets bias)
#however, as sampling effort increases, the trade-off bias variance awards other techniques
#basically the bias of the uniform approach is not anymore compensated by its precision
#and other methods are less biased and more precise as N increases..
Mse_plot <- ggplot(Mse_df, aes(x = N, y = Mse_val, col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(Random = "#F4B95AFF", Prox = "#C70E7BFF", Strat = "#007BC3FF",
                                Uniform = "#EF7C12FF", Systematic = "#FCEA1BFF", Topographic = "#009F3FFF"),
                     labels = c(Random = "Random", Prox = "Proximity", Strat = "Stratified",
                                Uniform = "Uniform", Systematic = "Systematic", Topographic = "Topographic")) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Mean squared error") + xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) 

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/MSE_sperandii.jpeg", device = "jpeg", units = "cm",
       width = 18, height = 15, dpi = 300)

#compute bias and variance, so split MSE in its 2 components
#tapply with 3 grouping indices gives you a 3D array
#Bias
Big_comp.copy <- Big_comp

Big_comp.copy$Bias <- with(Big_comp.copy, Val - True_val)

Bias_arr <- with(Big_comp.copy, tapply(Bias, INDEX = list(Type, N, Coef), mean))

Bias_df <- data.frame(do.call(rbind, lapply(1:(dim(Bias_arr)[3]), function(i) {
  df <- Bias_arr[,,i]
  dimn <- dimnames(df)
  dim(df) <- NULL
  df.res <- data.frame(Bias = df, N = rep(dimn[[2]], each = 6), Type = rep(dimn[[1]], 7))
  return(df.res)
  })), Coef = rep(levels(Big_comp.copy$Coef), each = 42))

Bias_df$Coef <- factor(Bias_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Bias_df$Type <- factor(Bias_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                         "Systematic", "Topographic"))
Bias_plot <- ggplot(Bias_df, aes(x = N, y = Bias, group = Type, col = Type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c(Random = "#F4B95AFF", Prox = "#C70E7BFF", Strat = "#007BC3FF",
                                Uniform = "#EF7C12FF", Systematic = "#FCEA1BFF", Topographic = "#009F3FFF"),
                     labels = c(Random = "Random", Prox = "Proximity", Strat = "Stratified",
                                Uniform = "Uniform", Systematic = "Systematic", Topographic = "Topographic")) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#Variance
Var_arr <- with(Big_comp, tapply(Val, INDEX = list(Type, N, Coef), FUN = var))

Var_df <- data.frame(do.call(rbind, lapply(1:(dim(Var_arr)[3]), function(i) {
  df <- Var_arr[,,i]
  dimn <- dimnames(df)
  dim(df) <- NULL
  df.res <- data.frame(Var = df, N = rep(dimn[[2]], each = 6), Type = rep(dimn[[1]], 7))
  return(df.res)
  })), Coef = rep(levels(Big_comp.copy$Coef), each = 42))

Var_df$Coef <- factor(Var_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Var_df$Type <- factor(Var_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                               "Systematic", "Topographic"))

Var_df$RMSE <- with(Var_df, sqrt(Var))

Var_plot <- ggplot(Var_df, aes(x = N, y = RMSE, group = Type, col = Type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c(Random = "#F4B95AFF", Prox = "#C70E7BFF", Strat = "#007BC3FF",
                                Uniform = "#EF7C12FF", Systematic = "#FCEA1BFF", Topographic = "#009F3FFF"),
                     labels = c(Random = "Random", Prox = "Proximity", Strat = "Stratified",
                                Uniform = "Uniform", Systematic = "Systematic", Topographic = "Topographic")) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("SQRT(Variance)") + xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

ggarrange(Bias_plot, Var_plot, ncol = 1, nrow = 2, common.legend = T)

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/B_V_sperandii.jpeg", device = "jpeg", units = "cm",
       width = 18, height = 25, dpi = 300)

#------------------------------------
#things to try:
#predictions maps..create layers with average predicted value for the species [done]
#compare map of true prob with map of prediction for AOI - also, this can be done with scatterplot of estimated vs true proba [done]
#compare response curves for a small number of replications
#maps of sampling..create maps showing how the different approaches sample the geogr space [done]

#-----------------------------maps of predictions!

#df with true proba value
#True_pred_sperandii <- data.frame(var_val = c(Plot_bio1[Plot_bio1$Species == "D. sperandii", "Temperature"],
#                                              Plot_bio12[Plot_bio12$Species == "D. sperandii", "Precipitation"]),
#                                  var = rep(c("Bio1", "Bio12"), times = c(nrow(Plot_bio1[Plot_bio1$Species == "D. sperandii", ]),
#                                                                          nrow(Plot_bio12[Plot_bio12$Species == "D. sperandii", ]))),
#                                  rbind(Plot_bio1[Plot_bio1$Species == "D. sperandii", c("Proba"), drop = F],
#                                        Plot_bio12[Plot_bio12$Species == "D. sperandii", c("Proba"), drop = F]))

#get rid of values not present in the AOI
#minValue(Chelsa.AOI$Bio1); maxValue(Chelsa.AOI$Bio1)
#minValue(Chelsa.AOI$Bio12); maxValue(Chelsa.AOI$Bio12)

#True_pred_sperandii <- do.call(rbind, lapply(unique(True_pred_sperandii$var), function(var_nm) {
#  sub_df <- True_pred_sperandii[True_pred_sperandii$var == var_nm, ]
#  if(var_nm %in% "Bio1") {
#    sub_df <- sub_df[sub_df$var_val <= 14.25 & sub_df$var_val >= -0.35, ]
#  } else {
#    sub_df <- sub_df[sub_df$var_val <= 1564 & sub_df$var_val >= 631.3, ]
#  }
#  }))

#True_pred_sperandii was recreated follwoing the approach used in Rare_species
#Bio1_seq.AOI and Bio12_seq.AOI are in the 'Rare species' script

True_pred_sperandii.2 <- data.frame(var_val = c(Bio1_seq.AOI, Bio12_seq.AOI),
                                    var = rep(c("Bio1", "Bio12"),
                                              each = 100))

True_pred_sperandii.2$Proba <- unlist(lapply(c("Bio1", "Bio12"), function(nm) {
  if(nm %in% "Bio1") {
    mat.prd <- cbind(1, True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T],
                     (True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T]^2),
                     mean(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio12", "var_val", drop = T]))
    colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef_nm))
    return(pred_val)
  } else {
    mat.prd <- cbind(1, mean(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T]),
                     median(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T]^2),
                     True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio12", "var_val", drop = T])
    colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef_nm))
    return(pred_val)
  }
}))

#random
Random_mapping <- function(y, x, x_crds, n = 300, min_p = 30) {
  require(effects)
  npres <- T
  while(npres) {
    Random_plots <- x_crds[sample(nrow(x_crds), n, replace = F), ]
    Fake_df <- na.omit(extract(stack(y, x), Random_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  #keep same names as in the raster for predict
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  #get fitted for partial effects
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
  }

Random_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Random_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI, 
                                              x_crds =  Climate_crds, n = N, min_p = 30), simplify = F)
  return(res)
})

names(Random_map) <- paste("N", Sampl_effort, sep = "_")

Random_map.lyr <- lapply(Random_map, function(i) mean(stack(lapply(i, '[[', 1))))

Random_map.fit <- do.call(rbind, lapply(names(Random_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Random_map[[nm]]), function(.) {
    inner_df <- Random_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
  }))

ggplot(Random_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#strat

Strat_mapping <- function(y, x, n = 300, strata, min_p = 30) {
  npres <- T
  while(npres) {
    Strat_plots <- Strat_raster(x = strata, N = n)
    Fake_df <- na.omit(extract(stack(y, x), Strat_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
  }
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
}


Strat_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Strat_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI, n = N,
                                                 strata = Bio1_12.rcl.qrt.AOI, min_p = 30), simplify = F)
  return(res)
})

names(Strat_map) <- names(Random_map)

Strat_map.lyr <- lapply(Strat_map, function(i) mean(stack(lapply(i, '[[', 1))))

Strat_map.fit <- do.call(rbind, lapply(names(Strat_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Strat_map[[nm]]), function(.) {
    inner_df <- Strat_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Strat_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#prox

Proximity_mapping <- function(y, x, prox_lay, n = 300, min_p = 30) {
  Prox_df.full <- na.omit(as.data.frame(prox_lay, xy = T))
  npres <- T
  while(npres) {
    Prox_df <- Prox_df.full[sample(nrow(Prox_df.full), size = n, replace = F, prob = Prox_df.full$layer), ] #check
    Prox_df$layer <- NULL
    Fake_df <- na.omit(extract(stack(y, x), Prox_df, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
  }
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
  }

Prox_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Proximity_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                 prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

names(Prox_map) <- names(Random_map)

Prox_map.lyr <- lapply(Prox_map, function(i) mean(stack(lapply(i, function(.) .[[1]]))))

Prox_map.fit <- do.call(rbind, lapply(names(Prox_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Prox_map[[nm]]), function(.) {
    inner_df <- Prox_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Prox_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#unif

Uniform_mapping <- function(y, x, x_sdf, n = 300, rsl, min_p = 30) {
  npres <- T
  while(npres) {
    Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
    Fake_df <- na.omit(extract(stack(y, x), Uniform_points, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
}

Unif_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Uniform_mapping(y = D.sperandii.bin.AOI,
                                               x = Chelsa.AOI,
                                               x_sdf = Chelsa.AOI.df.sp,
                                               n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

names(Unif_map) <- names(Random_map)

Unif_map.lyr <- lapply(Unif_map, function(i) mean(stack(lapply(i, '[[', 1))))

Unif_map.fit <- do.call(rbind, lapply(names(Unif_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Unif_map[[nm]]), function(.) {
    inner_df <- Unif_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Unif_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#Systematic

Systematic_mapping <- function(y, x, N, poly_proj, min_p = 30) {
  npres <- T
  while(npres) {
    Pts_syst <- st_sample(x = poly_proj, size = N, type = "regular")
    Pts_syst <- st_transform(Pts_syst, crs = 4326)
    Pts_syst <- st_coordinates(Pts_syst)
    if(nrow(Pts_syst) > N) {
      Pts_syst <- Pts_syst[-sample(nrow(Pts_syst), size = (nrow(Pts_syst) - N), replace = F), ]
    }
    Fake_df <- na.omit(extract(stack(y, x), Pts_syst, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
}

Syst_map <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 100, expr = Systematic_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                  N = n, poly_proj = Elev_AOI.proj, min_p = 30), simplify = F)
  return(res)
})

names(Syst_map) <- names(Random_map)

Syst_map.lyr <- lapply(Syst_map, function(i) mean(stack(lapply(i, '[[', 1))))

Syst_map.fit <- do.call(rbind, lapply(names(Syst_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Syst_map[[nm]]), function(.) {
    inner_df <- Syst_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Syst_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#Topographic
Topo_mapping <- function(y, x, N, topo_layer, min_p = 30) {
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T)
  Topo_lyr.df <- Topo_lyr.df[!is.na(Topo_lyr.df$layer), ]
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N, replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
    }
  Mod <- glm(PA ~ Bio1 + I(Bio1^2) + Bio12, family = binomial, data = Fake_df)
  Map <- raster::predict(Chelsa.AOI, model = Mod, type = "response")
  Partial_fit <- as.data.frame(effects::predictorEffects(Mod))
  Partial_fit <- do.call(rbind, lapply(names(Partial_fit), function(nm) {
    df <- Partial_fit[[nm]]
    df$var <- nm
    colnames(df)[1] <- "var_val"
    return(df)
  }))
  return(list(Map = Map, Part_fit = Partial_fit))
}

Topo_map <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 100, expr = Topo_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                            N = n., topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
})

names(Topo_map) <- names(Random_map)

Topo_map.lyr <- lapply(Topo_map, function(i) mean(stack(lapply(i, '[[', 1))))

Topo_map.fit <- do.call(rbind, lapply(names(Topo_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Topo_map[[nm]]), function(.) {
    inner_df <- Topo_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Topo_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

##maps

par(mfrow = c(6, 3))
invisible(lapply(Random_map.lyr[c(1, 2, 7)], plot))
invisible(lapply(Strat_map.lyr[c(1, 2, 7)], plot))
invisible(lapply(Prox_map.lyr[c(1, 2, 7)], plot))
invisible(lapply(Unif_map.lyr[c(1, 2, 7)], plot))
invisible(lapply(Syst_map.lyr[c(1, 2, 7)], plot))
invisible(lapply(Topo_map.lyr[c(1, 2, 7)], plot))

compareRaster(Random_map[[1]], stack(Random_map),
              stack(Strat_map),
                    stack(Prox_map),
                          stack(Unif_map),
                                stack(Syst_map), stack(Topo_map), orig = T)

Coords_pred_map <- na.omit(as.data.frame(Random_map[[1]], xy = T))
Coords_pred_map$layer <- NULL

List_pred_map <- list(RM = Random_map.lyr, StrM = Strat_map.lyr, PxM = Prox_map.lyr,
                      UM = Unif_map.lyr, SysM = Syst_map.lyr, TM = Topo_map.lyr)


Df_pred_map <- do.call(rbind, lapply(names(List_pred_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(List_pred_map[[nm]]), function(i) {
    Proba <- extract(List_pred_map[[nm]][[i]], Coords_pred_map)
    Proba.df <- data.frame(Prob = Proba, N = paste("N", Sampl_effort[i], sep = "_"))
    return(Proba.df)
    }))
  df <- data.frame(df, Type = nm)
  })
  )

#True Proba
Vec_true_prob <- extract(D.sperandii.prob.AOI, Coords_pred_map)

#recycl
Df_pred_map$True_prob <- Vec_true_prob

ggplot(Df_pred_map, aes(x = True_prob, y = Prob, col = Type)) +
  geom_point(alpha = .1) +
  geom_abline(slope = 1, intercept = 0, col = "black", lty = "dotted") +
  facet_wrap(~ N) +
  ylab("Estimated probabilities") + xlab("True probabilities") +
  theme_minimal() +
  theme(legend.position = "top")

#try enhance visualisation using mean of p for 0.1 classes
#Cut observed proba in 10 groups by .1 intervals

Unique_combo <- unique(Df_pred_map[c("N", "Type")])

Cut_df_pred <- do.call(rbind, Map(function(x, y) {
  Df_res <- Df_pred_map[Df_pred_map$Type == x & Df_pred_map$N == y, ]
  Df_res$cut_p <- with(Df_res, cut(Prob, breaks = seq(0, 1, .1), labels = seq(1, 10)))
  #here aggregate and return
  Df_res <- aggregate(Df_res[c("Prob", "True_prob")], by = list(Df_res$cut_p), mean)
  colnames(Df_res)[1] <- "Group"
  Df_res$Type <- x
  Df_res$N <- y
  return(Df_res)
}, x = Unique_combo$Type, y = Unique_combo$N))

row.names(Cut_df_pred) <- seq_len(nrow(Cut_df_pred))

#library(paletteer); library(ggthemes)
paletteer_d("LaCroixColoR::paired")
paletteer_d("beyonce::X66") #yellow of systematic from here

#range of p -> defines number of points in calibr plot
lapply(List_pred_map, function(i) sapply(i, cellStats, "max"))

Cor_pred_plot <- ggplot(Cut_df_pred[Cut_df_pred$N %in% c("N_200", "N_250","N_300"), ], aes(x = True_prob, y = Prob, col = Type)) +
  geom_abline(slope = 1, intercept = 0, col = "black", lty = "dotted", alpha = .5) +
  geom_point(alpha = .8, cex = 3) +
  geom_line(alpha = .5) +
  scale_color_manual(values = c(PxM = "#C70E7BFF", StrM = "#007BC3FF", TM = "#009F3FFF",
                                RM = "#F4B95AFF", SysM = "#FCEA1BFF",
                                UM = "#EF7C12FF"),
                     labels = c(PxM = "Proximity", StrM = "Stratified", TM = "Topographic",
                                RM = "Random", SysM = "Systematic", UM = "Uniform")) +
  facet_wrap(~ N) +
  ylab("Estimated probabilities") + xlab("True probabilities") +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 14))

#compute correlations
lapply(unique(Df_pred_map$Type), function(typ) {
  sapply(unique(Df_pred_map$N), function(sz) {
    with(Df_pred_map[Df_pred_map$Type == typ & Df_pred_map$N == sz, ], cor(Prob, True_prob,
                                                                           method = "spearman"))
  })
})

#plots of response curves
RespCurves.df <- data.frame(rbind(Random_map.fit, Strat_map.fit, Prox_map.fit,
                    Unif_map.fit, Syst_map.fit, Topo_map.fit),
                    Typ = rep(names(List_pred_map), each = nrow(Random_map.fit)))


Sampl_label <- c("RM" = "Random", "StrM" = "Stratified", "PxM" = "Proximity",
                    "UM" = "Uniform", "SysM" = "Systematic", "TM" = "Topographic")

RespCurves.df$Typ <- factor(RespCurves.df$Typ,
                                 levels = c("RM", "StrM", "PxM", "UM", "SysM", "TM"))

ggarrange(plotlist = lapply(levels(RespCurves.df$Typ), function(typ) {
  plot_rc <- ggplot(RespCurves.df[RespCurves.df$Typ == typ, ], aes(x = var_val, y = fit)) +
    geom_line(aes(group = repl), alpha = .3) +
    geom_line(data = True_pred_sperandii.2, aes(x = var_val, y = Proba), col = "red") +
    facet_grid(N ~ var, scales = "free_x",
               labeller = as_labeller(c("N_200" = "200", "N_250" = "250", "N_300" = "300",
                                        "N_350" = "350", "N_400" = "400", "N_450" = "450", "N_500" = "500",
                                        "Bio1" = "Temp", "Bio12" = "Precip"))) +
    ylab(NULL) + xlab(Sampl_label[typ]) +
    theme_classic() +
    theme(axis.title = element_text(size = 12), strip.text = element_text(size = 10),
          axis.text.x.bottom = element_text(angle = 35, vjust = .8))
  #
  if(typ %in% c("RM", "UM")) {
    plot_rc <- plot_rc + ylab("Occurrence probability") 
    }
  return(plot_rc)
}))

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/RespCurvesSperandii.jpeg", device = "jpeg",
       width = 32, height = 25, units = "cm", dpi = 300)

##MAPS OF SAMPLING EFFORT--------------------------------------------------------------------------------------

#Random sampl

Random_sampl <- function(x_crds, n = 300) {
  Random_plots <- x_crds[sample(nrow(x_crds), n, replace = F), ]
  return(Random_plots)
}

Random_xy <- do.call(rbind, replicate(n = 100, expr = Random_sampl(x_crds = Climate_crds, n = 500), simplify = F))

Random_xy.counts <- rasterize(Random_xy, D.sperandii.prob.AOI, fun = "count")

plot(Random_xy.counts, col = heat.colors(10))

#Strat sampl - re-run on May 29th with new reclassified raster

Strat_sampl <- function(strata, n = 300) {
  Strat_plots <- Strat_raster(x = strata, N = n)
  return(Strat_plots)
}

Strat_xy <- do.call(rbind, replicate(n = 100,
                                     expr = Strat_sampl(strata = Bio1_12.rcl.qrt.AOI, n = 500),
                                     simplify = F))

Strat_xy.counts <- rasterize(Strat_xy, D.sperandii.prob.AOI, fun = "count")

plot(Strat_xy.counts, col = heat.colors(10))

#Prox sampl

Proximity_sampl <- function(prox_lay, n = 300) {
  Prox_df <- na.omit(as.data.frame(prox_lay, xy = T))
  Prox_df <- Prox_df[sample(nrow(Prox_df), size = n, replace = F, prob = Prox_df$layer), ] #check
  Prox_df$layer <- NULL
  return(Prox_df)
}

Prox_xy <- do.call(rbind, replicate(n = 100,
                                     expr = Proximity_sampl(prox_lay = Abr_highway_AOI, n = 500),
                                     simplify = F))

Prox_xy.counts <- rasterize(Prox_xy, D.sperandii.prob.AOI, fun = "count")

plot(Prox_xy.counts, col = rev(heat.colors(10)))

#Unif sampl

Uniform_sampl <- function(x_sdf, rsl, n = 300) {
  Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
  return(Uniform_points)
}

Unif_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Uniform_sampl(x_sdf = Chelsa.AOI.df.sp,
                                                      n = 500, rsl = 10),
                                    simplify = F))

Unif_xy.counts <- rasterize(Unif_xy, D.sperandii.prob.AOI, fun = "count")

plot(Unif_xy.counts, col = rev(heat.colors(10)))

#Syst sampl

Systematic_sampl <- function(poly_proj, N) {
  Pts_syst <- st_sample(x = poly_proj, size = N, type = "regular")
  Pts_syst <- st_transform(Pts_syst, crs = 4326)
  Pts_syst <- st_coordinates(Pts_syst)
  if(nrow(Pts_syst) > N) {
    Pts_syst <- Pts_syst[-sample(nrow(Pts_syst), size = (nrow(Pts_syst) - N), replace = F), ]
  }
  return(Pts_syst)
}

Syst_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Systematic_sampl(poly_proj = Elev_AOI.proj, N = 500),
                                    simplify = F))

Syst_xy.counts <- rasterize(Syst_xy, D.sperandii.prob.AOI, fun = "count")

plot(Syst_xy.counts, col = rev(heat.colors(10)))

#Topogr sampl

Topo_sampl <- function(topo_layer, N) {
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T)
  Topo_lyr.df <- Topo_lyr.df[!is.na(Topo_lyr.df$layer), ]
  Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N, replace = F), c("x", "y")]
  return(Topo_coords)
}

Topo_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Topo_sampl(topo_layer = Topogr_het_Abr.AOI, N = 500),
                                    simplify = F))

Topo_xy.counts <- rasterize(Topo_xy, D.sperandii.prob.AOI, fun = "count")

plot(Topo_xy.counts, col = rev(heat.colors(10)))

#unique map of sampling effort
Sampl_eff.maps <- list(Rnd = Random_xy.counts, Str = Strat_xy.counts, Prx = Prox_xy.counts, 
                       Unif = Unif_xy.counts, Syst = Syst_xy.counts, Topo = Topo_xy.counts)

#fill all NAs with 0 val, and mask the layer again (otherwise zeros are not represented in the map)
Sampl_eff.maps <- lapply(Sampl_eff.maps, function(.) {
  .[is.na(.)] <- 0
  Map <- mask(., mask = Elev_AOI)
  return(Map)
  })

Label_sampl_eff_legend <- setNames(c("Random", "Stratified", "Proximity",
                                     "Uniform", "Systematic", "Topographic"), nm = names(Sampl_eff.maps))

Sampl_eff.maps.plots <- lapply(names(Sampl_eff.maps), function(nm) {
  Map_df <- as.data.frame(Sampl_eff.maps[[nm]], xy = T, na.rm = T)
  Map_plot <- ggplot() +
    geom_tile(data = Map_df, aes(x = x, y = y, fill = layer)) +
    geom_sf(data = Elev_AOI, color = "white", fill = NA) +
    coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
    scale_fill_carto_c(type = "quantitative", palette = "PurpOr") +
    labs(title = Label_sampl_eff_legend[nm],
         x = "Longitude", y = "Latitude", fill = "# units") +
    theme(axis.title = element_text(size = 14), title = element_text(size = 16),
          legend.text = element_text(size = 10))
  if(nm %in% "Str") {
    Map_plot <- Map_plot + geom_sf(data = Contour.rcl.qrt.AOI, color = "black", fill = NA, alpha = .4)
  }
  return(Map_plot)
  })

ggarrange(plotlist = Sampl_eff.maps.plots)

#patchwork
(Sampl_eff.maps.plots[[1]] + Sampl_eff.maps.plots[[2]] + Sampl_eff.maps.plots[[3]])/(Sampl_eff.maps.plots[[4]] + Sampl_eff.maps.plots[[5]] + Sampl_eff.maps.plots[[6]])

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/SamplEffort.png",
       device = "png", dpi = 300, width = 30, height = 18, units = "cm")

#local correlation between Bio1 and Bio12
LocalCor_Bio <- raster::corLocal(x = Chelsa.AOI$Bio1, y = Chelsa.AOI$Bio12, ngb = 3)

LocalCor_Bio

##save data for markdown -> last save: 31/05/2022
#spatial data
save(Chelsa.AOI, Chelsa.stack, D.sperandii.prob, D.sperandii.bin, Elev_AOI, 
     file = "SpatialData.RData")
#dfs
save(Elevation_Abr.df, Mse_df, D.sperandii.bin.df, D.sperandii.prob.df, Abruzzo, Abr_highway_dist.geo.prb.df,
     file = "SpatialDF.RData")

#figs
save(Mse_plot, Bias_plot, Var_plot, Cor_pred_plot, 
     file = "ResPlots.RData")

#difference random and strat
Rnd_vs_Str <- do.call(rbind, lapply(c("Rnd", "Str"), function(nm, strata = Bio1_12.rcl.qrt.AOI) {
  mappa.df <- as.data.frame(Sampl_eff.maps[[nm]], na.rm = T, xy = T)
  mappa.val <- data.frame(extract(strata, mappa.df[c("x", "y")], df = T),
                          count_str = mappa.df$layer)
  mappa.val$typ <- nm
  mappa.val$ID <- NULL
  return(mappa.val)
  }))

Rnd_vs_Str <- na.omit(Rnd_vs_Str)

ggplot(Rnd_vs_Str, aes(x = count_str, fill = typ)) +
  geom_histogram(alpha = .5, position = 'identity') +
  facet_wrap(~ layer)


