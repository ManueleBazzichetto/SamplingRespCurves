##Code to replicate results presented in the manuscript "Sampling strategies matter to accurately estimate response curves’ parameters in species distribution models"
##Figures

#Reproduce some of the figures reported in the manuscript
library(ggplot2)
#library(patchwork)
library(ggpubr)
library(colorBlindness) #colours for cells proba

#Figure of the species' (simulated) response curves
True_coef.S.nm; True_coef.T.nm

Plot_bio1 <- data.frame(Temperature = rep(Bio1_seq, 2),
                        Proba = plogis(c(cbind(1, Bio1_seq, Bio1_seq^2)%*%c((-11.3886498 + 0.0071542*cellStats(Chelsa.stack$Bio12, 'mean')), 1.3835909, -0.0920092),
                                         cbind(1, Bio1_seq, Bio1_seq^2)%*%c((-22.17261596 + 0.00490439*cellStats(Chelsa.stack$Bio12, 'mean')), 3.77916980, -0.20996874))),
                        Species = rep(c("D. sperandii", "D. tundrae"), each = length(Bio1_seq)))

Plot_bio12 <- data.frame(Precipitation = rep(Bio12_seq, 2),
                         Proba = plogis(c(cbind(1, Bio12_seq)%*%c((-11.3886498 + 1.3835909*cellStats(Chelsa.stack$Bio1, 'mean') - 0.0920092*cellStats(Chelsa.stack$Bio1^2, 'mean')), 0.0071542),
                                          cbind(1, Bio12_seq)%*%c((-22.17261596 + 3.77916980*cellStats(Chelsa.stack$Bio1, 'mean') - 0.20996874*cellStats(Chelsa.stack$Bio1^2, 'mean')), 0.00490439))),
                         Species = rep(c("D. sperandii", "D. tundrae"), each = length(Bio12_seq)))

#bio1
Plot_bio1_fig <- ggplot(Plot_bio1, aes(x = Temperature, y = Proba, col = Species)) +
  geom_line(lwd = 2) +
  scale_color_manual(values = c("#DCE319FF", "#39568CFF")) +
  ylab("Occurrence probability") + ylim(c(0, 1)) +
  theme_pubr() +
  theme(legend.position = "top", axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14, face = "italic"))

#bio12
Plot_bio12_fig <- ggplot(Plot_bio12, aes(x = Precipitation, y = Proba, col = Species)) +
  geom_line(lwd = 2) +
  scale_color_manual(values = c("#DCE319FF", "#39568CFF")) +
  ylab("Occurrence probability") + ylim(c(0, 1)) +
  theme_pubr() +
  theme(legend.position = "top", axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14, face = "italic"))

FigSimRespCurve <- ggarrange(Plot_bio1_fig, Plot_bio12_fig, nrow = 1, ncol = 2, common.legend = T)

#Figure of the species' (simulated) distributions
#Dfs for ggplotting - D. sperandii
D.sperandii.prob.df <- as.data.frame(D.sperandii.prob, xy = T, na.rm = T)
colnames(D.sperandii.prob.df)[3] <- "Prob"

D.sperandii.bin.df <- as.data.frame(D.sperandii.bin, xy = T, na.rm = T)
colnames(D.sperandii.bin.df)[3] <- "PA"

#Dfs for ggplotting - D. tundrae
D.tundrae.prob.df <- as.data.frame(D.tundrae.prob, xy = T, na.rm = T)
colnames(D.tundrae.prob.df)[3] <- "Prob"

D.tundrae.bin.df <- as.data.frame(D.tundrae.bin, xy = T, na.rm = T)
colnames(D.tundrae.bin.df)[3] <- "PA"

ggarrange(ggplot() +
            geom_raster(data = D.sperandii.prob.df, aes(x = x, y = y, fill = Prob)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_gradient(low = "#C8C8BC", high = "#DCE319FF") +
            labs(title = "Dianthus sperandii", x = "Long", y = "Lat", fill = "Prob") +
            theme_pubr() +
            theme(panel.background = element_blank(), 
                  legend.position = "right",
                  plot.title = element_text(face = 'bold.italic', hjust = .5, size = 16),
                  axis.title = element_text(size = 14)),
          ggplot() +
            geom_raster(data = D.tundrae.prob.df, aes(x = x, y = y, fill = Prob)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_gradient(low = "#CDD2DA", high = "#39568CFF") +
            labs(title = "Dianthus tundrae", x = "Long", y = "Lat", fill = "Prob") +
            theme_pubr() +
            theme(panel.background = element_blank(), 
                  legend.position = "right",
                  plot.title = element_text(face = 'bold.italic', hjust = .5, size = 16),
                  axis.title = element_text(size = 14)),
          ggplot() +
            geom_tile(data = D.sperandii.bin.df, aes(x = x, y = y, fill = factor(PA))) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_manual(values = c("0" = "#C8C8BC", "1" = "#DCE319FF")) +
            labs(title = NULL, x = "Long", y = "Lat", fill = "PresAbs") +
            theme_pubr() +
            theme(panel.background = element_blank(), 
                  legend.position = "right",
                  plot.title = element_text(face = 'bold.italic', hjust = .5, size = 16),
                  axis.title = element_text(size = 14)),
          ggplot() +
            geom_raster(data = D.tundrae.bin.df, aes(x = x, y = y, fill = factor(PA))) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_manual(values = c("0" = "#CDD2DA", "1" = "#39568CFF")) +
            labs(title = NULL, x = "Long", y = "Lat", fill = "PresAbs") +
            theme_pubr() +
            theme(panel.background = element_blank(), 
                  legend.position = "right",
                  plot.title = element_text(face = 'bold.italic', hjust = .5, size = 16),
                  axis.title = element_text(size = 14)),
          nrow = 2, ncol = 2, labels = "auto")

##Maps of sampling effort

#RANDOM

Random_sampl <- function(x_crds, n = 300) {
  Random_plots <- x_crds[sample(nrow(x_crds), n, replace = F), ]
  return(Random_plots)
}

set.seed(1789)
Random_xy <- do.call(rbind, replicate(n = 100, expr = Random_sampl(x_crds = Climate_crds, n = 500), simplify = F))

Random_xy.counts <- rasterize(Random_xy, D.sperandii.prob.AOI, fun = "count")

plot(Random_xy.counts, col = heat.colors(10))

#STRATIFIED

Strat_sampl <- function(strata, n = 300) {
  Strat_plots <- Strat_raster(x = strata, N = n)
  return(Strat_plots)
}

set.seed(1776)
Strat_xy <- do.call(rbind, replicate(n = 100,
                                     expr = Strat_sampl(strata = Bio1_12.rcl.qrt.AOI, n = 500),
                                     simplify = F))

Strat_xy.counts <- rasterize(Strat_xy, D.sperandii.prob.AOI, fun = "count")

plot(Strat_xy.counts, col = heat.colors(20))

#ROAD PROXIMITY

Proximity_sampl <- function(prox_lay, n = 300) {
  Prox_df <- na.omit(as.data.frame(prox_lay, xy = T))
  Prox_df <- Prox_df[sample(nrow(Prox_df), size = n, replace = F, prob = Prox_df$layer), ] #check
  Prox_df$layer <- NULL
  return(Prox_df)
}

set.seed(1792)
Prox_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Proximity_sampl(prox_lay = Abr_highway_AOI, n = 500),
                                    simplify = F))

Prox_xy.counts <- rasterize(Prox_xy, D.sperandii.prob.AOI, fun = "count")

plot(Prox_xy.counts, col = rev(heat.colors(10)))

#UNIFORM

Uniform_sampl <- function(x_sdf, rsl, n = 300) {
  Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
  return(Uniform_points)
}

set.seed(1861)
Unif_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Uniform_sampl(x_sdf = Chelsa.AOI.df.sp,
                                                         n = 500, rsl = 10),
                                    simplify = F))

Unif_xy.counts <- rasterize(Unif_xy, D.sperandii.prob.AOI, fun = "count")

plot(Unif_xy.counts, col = rev(heat.colors(10)))

#SYSTEMATIC

Systematic_sampl <- function(poly_proj, N) {
  Pts_syst <- st_sample(x = poly_proj, size = N, type = "regular")
  Pts_syst <- st_transform(Pts_syst, crs = 4326)
  Pts_syst <- st_coordinates(Pts_syst)
  if(nrow(Pts_syst) > N) {
    Pts_syst <- Pts_syst[-sample(nrow(Pts_syst), size = (nrow(Pts_syst) - N), replace = F), ]
  }
  return(Pts_syst)
}

set.seed(1947)
Syst_xy <- do.call(rbind, replicate(n = 100,
                                    expr = Systematic_sampl(poly_proj = Elev_AOI.proj, N = 500),
                                    simplify = F))

Syst_xy.counts <- rasterize(Syst_xy, D.sperandii.prob.AOI, fun = "count")

plot(Syst_xy.counts, col = rev(heat.colors(10)))

#TOPOGRAPHIC

Topo_sampl <- function(topo_layer, N) {
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T)
  Topo_lyr.df <- Topo_lyr.df[!is.na(Topo_lyr.df$layer), ]
  Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N, replace = F), c("x", "y")]
  return(Topo_coords)
}

set.seed(754)
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
    scale_fill_carto_c(type = "quantitative", palette = "SunsetDark") +
    labs(title = Label_sampl_eff_legend[nm],
         x = "Longitude", y = "Latitude", fill = "# units") +
    theme_pubr() +
    theme(axis.title = element_text(size = 14), title = element_text(size = 16),
          legend.text = element_text(size = 10), legend.position = "right")
  if(nm %in% "Str") {
    Map_plot <- Map_plot + geom_sf(data = Contour.rcl.qrt.AOI, color = "black", fill = NA, alpha = .6, lty = "longdash")
  }
  return(Map_plot)
})

ggarrange(plotlist = Sampl_eff.maps.plots)


