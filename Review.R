##Code to replicate results presented in the manuscript "Sampling strategies matter to accurately estimate response curves’ parameters in species distribution models"
##New analyses carried out for review

library(raster)
library(sf)
library(ggplot2)
library(ggpubr)
library(ggdensity) #for density plots in the environmental space
library(parallel)


##MISSING COVARIATES---------------------------------------------------------------

##Test the effect of missing covariates on the performance of the sampling approaches
##The distribution of two new virtual (sub)species is simulated: Dianthus sperandii subsp thermophilus and Dianthus tundrae subsp thermophilus
##Parameters for temperature and precipition are remained the same as those of D. sperandii and tundrae
##A third environmental variable (i.e., northness) is used to generate the two subspecies
##However, northness is then removed from the GLMs fitted to the data sampled by the different sampling strategies

#resample elevation layer to then compute northness using it
Elev_resample <- resample(Elevation_Abr, Chelsa.stack$Bio1, method = "bilinear")

#compute northness
Aspect_rsmp <- raster::terrain(x = Elev_resample, opt = "aspect", unit = "radians", neighbors = 8)
N_rsmp <- cos(Aspect_rsmp)

compareRaster(N_rsmp, Chelsa.stack) #TRUE

#create stack including temperature, precipitation and northness layers
T_p_n.stack <- stack(Chelsa.stack, N_rsmp)

#transform stack to data.frame
Df_t_p_n <- as.data.frame(T_p_n.stack, xy = T, na.rm = T)

#rename columns
colnames(Df_t_p_n)[c(3, 4, 5)] <- c("Temp", "Pr", "N")

##Distribution of Dianthus sperandii subsp thermophilus (wide subspecies)--------------

#regression parameter for northness set to: -1.4

#attach probabilities for wide subspecies
Df_t_p_n.pr.w <- Df_t_p_n
Df_t_p_n.pr.w$prob <- plogis(with(Df_t_p_n, True_coef.S.nm[["(Intercept)"]] + 
                                    True_coef.S.nm[["Temp"]]*Temp + True_coef.S.nm[["I(Temp^2)"]]*(Temp^2) +
                                    True_coef.S.nm[["Pr"]]*Pr - 1.4*N))

#generate layer with true probabilities
Wide_N.prob <- rasterize(x = Df_t_p_n.pr.w[c(1, 2)], y = Chelsa.stack$Bio1, field = Df_t_p_n.pr.w$prob)

#generate layer of 1/0
set.seed(969)
Wide_N.bin <- calc(Wide_N.prob, function(.) rbinom(1, 1, .))


##Distribution of Dianthus tundrae subsp thermophilus (subspecies with restricted distribution)--------------

#regression parameter for northness set to: -2

#attach probabilities for rare subspecies
Df_t_p_n.pr.r <- Df_t_p_n
Df_t_p_n.pr.r$prob <- plogis(with(Df_t_p_n, True_coef.T.nm[["(Intercept)"]] + 
                                    True_coef.T.nm[["Temp"]]*Temp + True_coef.T.nm[["I(Temp^2)"]]*(Temp^2) +
                                    True_coef.T.nm[["Pr"]]*Pr - 2*N))

#generate layer with true probabilities
Rare_N.prob <- rasterize(x = Df_t_p_n.pr.r[c(1, 2)], y = Chelsa.stack$Bio1, field = Df_t_p_n.pr.r$prob)

#generate layer of 1/0
set.seed(782)
Rare_N.bin <- calc(Rare_N.prob, function(.) rbinom(1, 1, .))

##Simulations---------------------------------------------------------------------

#params for species
Wide_N.prms <- c(True_coef.S.nm, "N" = -1.4)
Rare_N.prms <- c(True_coef.T.nm, "N" = -2)

#mask layers of subspecies distributions and predictors
T_p_n.stack.AOI <- mask(T_p_n.stack, mask = as(Elev_AOI, "Spatial"))

names(T_p_n.stack.AOI)[3] <- "N"

Wide_N.bin.AOI <- mask(Wide_N.bin, mask = as(Elev_AOI, "Spatial"))
Rare_N.bin.AOI <- mask(Rare_N.bin, mask = as(Elev_AOI, "Spatial"))

compareRaster(Wide_N.bin.AOI, T_p_n.stack.AOI)


#Simulations for Dianthus sperandii subsp thermophilus--------------------------------------------------------------

#------random

#needed to avoid sampling location where N (and Wide_N.bin.AOI) has NA
T_p_n.AOI.xy <- as.data.frame(Wide_N.bin.AOI, xy = T, na.rm = T)[c("x", "y")]

set.seed(479)
Random.wide_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Random_fit(y = Wide_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                              x_crds = T_p_n.AOI.xy, n = N, min_p = 30), 
                   simplify = F)
  return(res)
})

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Random.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Random.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Random.wide_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Random.wide_N.mat <- do.call(rbind, lapply(names(Random.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Random.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Random.wide_N.mat$Coef <- factor(Random.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Random.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Random.wide_N.mat$Coef)])
Random.wide_N.mat$N <- factor(Random.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Random.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------strat

#need to create a new Bio1_12.rcl.qrt.AOI with NA where Wide_N.bin.AOI has NA (due to northness)
compareRaster(Bio1_12.rcl.qrt.AOI, Wide_N.bin.AOI) #T

Bio1_12.strat.AOI.wide_N <- overlay(x = Wide_N.bin.AOI, y = Bio1_12.rcl.qrt.AOI, fun = function(x, y) {
  y[is.na(x)] <- NA
  return(y)
})

set.seed(835)
Strat.wide_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, Strat_fit(y = Wide_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]], n = N,
                                      strata = Bio1_12.strat.AOI.wide_N, min_p = 30), simplify = F)
  return(res)
})

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Strat.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Strat.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Strat.wide_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Strat.wide_N.mat <- do.call(rbind, lapply(names(Strat.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Strat.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Strat.wide_N.mat$Coef <- factor(Strat.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Strat.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Strat.wide_N.mat$Coef)])
Strat.wide_N.mat$N <- factor(Strat.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Strat.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------proximity

set.seed(210)
Prox.wide_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Proximity_fit(y = Wide_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                                 prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Prox.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Prox.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Prox.wide_N.res) <- paste0("N_", Sampl_effort)

Prox.wide_N.mat <- do.call(rbind, lapply(names(Prox.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Prox.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Prox.wide_N.mat$Coef <- factor(Prox.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Prox.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Prox.wide_N.mat$Coef)])
Prox.wide_N.mat$N <- factor(Prox.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Prox.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------systematic

set.seed(520)
Syst.wide_N.res <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Systematic_fit(y = Wide_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                                  N = n, perc_inc = .07, poly_proj = Elev_AOI.proj,
                                                  min_p = 30), simplify = F)
  return(res)
})

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Syst.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Syst.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Syst.wide_N.res) <- paste0("N_", Sampl_effort)

Syst.wide_N.mat <- do.call(rbind, lapply(names(Syst.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Syst.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Syst.wide_N.mat$Coef <- factor(Syst.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Syst.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Syst.wide_N.mat$Coef)])
Syst.wide_N.mat$N <- factor(Syst.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Syst.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------topographic

set.seed(669)
Topo.wide_N.res <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 500, expr = Topo_fit(y = Wide_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                            N = n., perc_incr = .07,
                                            topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
  return(res)
})

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Topo.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Topo.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Topo.wide_N.res) <- paste0("N_", Sampl_effort)

Topo.wide_N.mat <- do.call(rbind, lapply(names(Topo.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Topo.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Topo.wide_N.mat$Coef <- factor(Topo.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Topo.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Topo.wide_N.mat$Coef)])
Topo.wide_N.mat$N <- factor(Topo.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Topo.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------uniform

T_p_n.AOI.df <- as.data.frame(T_p_n.stack.AOI, xy = T, na.rm = T)
T_p_n.AOI.df$Bio1.sc <- scale(T_p_n.AOI.df$Bio1)
T_p_n.AOI.df$Bio12.sc <- scale(T_p_n.AOI.df$Bio12)
T_p_n.AOI.df.sp <- st_as_sf(T_p_n.AOI.df, coords = c("Bio1.sc", "Bio12.sc"))

#go parallel
cr.wide_N <- parallel::makeCluster(7)

parallel::clusterExport(cr.wide_N, c("Wide_N.bin.AOI", "T_p_n.stack.AOI",
                                     "T_p_n.AOI.df.sp", "Uniform_fit", 
                                     "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr.wide_N, list(library(sf), library(raster)))

set.seed(401)
Unif.wide_N.res <- parallel::parLapply(cl = cr.wide_N, X = Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Uniform_fit(y = Wide_N.bin.AOI,
                                               x = T_p_n.stack.AOI[[c(1, 2)]],
                                               x_sdf = T_p_n.AOI.df.sp,
                                               n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

stopCluster(cr.wide_N)

#check smp size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Unif.wide_N.res[[i]], '[[', 3) == Sampl_effort[i]) == 500)

#check corr
sapply(Unif.wide_N.res, function(i) mean(sapply(i, '[[', 2)))

#name list
names(Unif.wide_N.res) <- paste0("N_", Sampl_effort)

Unif.wide_N.mat <- do.call(rbind, lapply(names(Unif.wide_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif.wide_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Unif.wide_N.mat$Coef <- factor(Unif.wide_N.mat$Coef, levels = names(Wide_N.prms)[-5])
Unif.wide_N.mat$True_val <- unname(Wide_N.prms[as.character(Unif.wide_N.mat$Coef)])
Unif.wide_N.mat$N <- factor(Unif.wide_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Unif.wide_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#compute MSE, bias and variance

Wide_N.df <- rbind(data.frame(Random.wide_N.mat, Type = "random"),
                   data.frame(Strat.wide_N.mat, Type = "stratified"),
                   data.frame(Prox.wide_N.mat, Type = "proximity"),
                   data.frame(Syst.wide_N.mat, Type = "systematic"),
                   data.frame(Topo.wide_N.mat, Type = "topographic"),
                   data.frame(Unif.wide_N.mat, Type = "uniform"))

Wide_N.df$Type <- factor(Wide_N.df$Type, levels = c("random", "proximity", "stratified", "uniform",
                                                    "systematic", "topographic"))

#check results
ggplot(Wide_N.df, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#unique combos of Coef and N
Unique_coefn.wide_N <- unique(Wide_N.df[c("Coef", "N")]) 

Wide_N.MSE.list <- lapply(as.character(unique(Wide_N.df$Type)), function(ty) {
  mapply(function(x, y) Mse_samp(df = Wide_N.df, Coef. = x, N. = y, Type. = ty, True_vals = Wide_N.prms),
         x = as.character(Unique_coefn.wide_N$Coef), y = as.character(Unique_coefn.wide_N$N))
})

names(Wide_N.MSE.list) <- as.character(unique(Wide_N.df$Type))

Wide_N.MSE.df <- do.call(rbind, lapply(names(Wide_N.MSE.list), function(ty) {
  mse.df <- data.frame(Val = unname(Wide_N.MSE.list[[ty]]),
                       Coef = names(Wide_N.MSE.list[[ty]]),
                       N = rep(paste0("N_", Sampl_effort), each = 4), 
                       Type = ty)
  return(mse.df)
}))

Wide_N.MSE.df$Coef <- factor(Wide_N.MSE.df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Wide_N.MSE.df$Type <- factor(Wide_N.MSE.df$Type, levels = c("random", "proximity", "stratified", "uniform",
                                                            "systematic", "topographic"))

#add rmse
Wide_N.MSE.df$Rmse <- with(Wide_N.MSE.df, sqrt(Val))

Wide_N.RMSE_plot <- ggplot(Wide_N.MSE.df, aes(x = N, y = Rmse, col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(random = "#F4B95AFF", proximity = "#C70E7BFF", stratified = "#007BC3FF",
                                uniform = "#EF7C12FF", systematic = "#FCEA1BFF", topographic = "#009F3FFF"),
                     labels = c(random = "Random", proximity = "Proximity", stratified = "Stratified",
                                uniform = "Uniform", systematic = "Systematic", topographic = "Topographic")) +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50)))
  )) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) 


#Simulations for Dianthus tundrae subsp thermophilus--------------------------------------------------------------


#------random

set.seed(4783)
Random.rare_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Random_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                               x_crds = T_p_n.AOI.xy, n = N, min_p = 30), simplify = F)
  return(res)
})

#check separation
sapply(Random.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Random.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Random.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Random.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Random.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Random.rare_N.mat <- do.call(rbind, lapply(names(Random.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Random.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Random.rare_N.mat$Coef <- factor(Random.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Random.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Random.rare_N.mat$Coef)])
Random.rare_N.mat$N <- factor(Random.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Random.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")


#------strat

set.seed(9456)
#used Bio1_12.strat.AOI.wide_N as same layer for Rare_N would have same NAs (at same locations)
Strat.rare_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Strat_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                              strata = Bio1_12.strat.AOI.wide_N, n = N, min_p = 30), simplify = F)
  return(res)
})

#check separation
sapply(Strat.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Strat.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Strat.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Strat.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Strat.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Strat.rare_N.mat <- do.call(rbind, lapply(names(Strat.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Strat.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Strat.rare_N.mat$Coef <- factor(Strat.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Strat.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Strat.rare_N.mat$Coef)])
Strat.rare_N.mat$N <- factor(Strat.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Strat.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------proximity

set.seed(1678) 
Prox.rare_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Proximity_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                                  prox_lay = Abr_highway_AOI, n = N, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation
sapply(Prox.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Prox.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Prox.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Prox.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Prox.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Prox.rare_N.mat <- do.call(rbind, lapply(names(Prox.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Prox.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Prox.rare_N.mat$Coef <- factor(Prox.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Prox.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Prox.rare_N.mat$Coef)])
Prox.rare_N.mat$N <- factor(Prox.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Prox.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------systematic

set.seed(4758)
Syst.rare_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Systematic_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                                   poly_proj = Elev_AOI.proj, N = N, perc_incr = .07, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation
sapply(Syst.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Syst.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Syst.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Syst.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Syst.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Syst.rare_N.mat <- do.call(rbind, lapply(names(Syst.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Syst.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Syst.rare_N.mat$Coef <- factor(Syst.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Syst.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Syst.rare_N.mat$Coef)])
Syst.rare_N.mat$N <- factor(Syst.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Syst.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------topographic

set.seed(3657)
Topo.rare_N.res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Topo_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                             topo_layer = Topogr_het_Abr.AOI, N = N, perc_incr = .07, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation
sapply(Topo.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Topo.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Topo.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Topo.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Topo.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Topo.rare_N.mat <- do.call(rbind, lapply(names(Topo.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Topo.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Topo.rare_N.mat$Coef <- factor(Topo.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Topo.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Topo.rare_N.mat$Coef)])
Topo.rare_N.mat$N <- factor(Topo.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Topo.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#------uniform

#go parallel
cr.rare_N <- parallel::makeCluster(7)

parallel::clusterExport(cr.rare_N, c("Rare_N.bin.AOI", "T_p_n.stack.AOI",
                                     "T_p_n.AOI.df.sp", "Uniform_rare", 
                                     "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr.rare_N, list(library(sf), library(raster)))

set.seed(8263)
Unif.rare_N.res <- parLapply(cl = cr.rare_N, X = Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Uniform_rare(y = Rare_N.bin.AOI, x = T_p_n.stack.AOI[[c(1, 2)]],
                                                x_sdf = T_p_n.AOI.df.sp,
                                                n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

stopCluster(cr.rare_N)

#check separation
sapply(Unif.rare_N.res, function(i) sum(sapply(i, '[[', "Sep"))) #fine

#check smp sz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Unif.rare_N.res[[i]], '[[', "SampSz") != Sampl_effort[i]))

#check corr
sapply(Unif.rare_N.res, function(i) mean(sapply(i, '[[', "Cor")))

#check #p
lapply(Unif.rare_N.res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#name list
names(Unif.rare_N.res) <- paste0("N_", Sampl_effort)

#extract coefs
Unif.rare_N.mat <- do.call(rbind, lapply(names(Unif.rare_N.res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif.rare_N.res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  Mat <- as.vector(Mat)
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  return(Mat)
}))

Unif.rare_N.mat$Coef <- factor(Unif.rare_N.mat$Coef, levels = names(Rare_N.prms)[-5])
Unif.rare_N.mat$True_val <- unname(Rare_N.prms[as.character(Unif.rare_N.mat$Coef)])
Unif.rare_N.mat$N <- factor(Unif.rare_N.mat$N, levels = paste0("N_", Sampl_effort))

#check
ggplot(Unif.rare_N.mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top")

#compute MSE, bias and variance

Rare_N.df <- rbind(data.frame(Random.rare_N.mat, Type = "random"),
                   data.frame(Strat.rare_N.mat, Type = "stratified"),
                   data.frame(Prox.rare_N.mat, Type = "proximity"),
                   data.frame(Syst.rare_N.mat, Type = "systematic"),
                   data.frame(Topo.rare_N.mat, Type = "topographic"),
                   data.frame(Unif.rare_N.mat, Type = "uniform"))

Rare_N.df$Type <- factor(Rare_N.df$Type, levels = c("random", "proximity", "stratified", "uniform",
                                                    "systematic", "topographic"))

#check results
ggplot(Rare_N.df, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#unique combos of Coef and N
Unique_coefn.rare_N <- unique(Rare_N.df[c("Coef", "N")]) 

Rare_N.MSE.list <- lapply(as.character(unique(Rare_N.df$Type)), function(ty) {
  mapply(function(x, y) Mse_samp(df = Rare_N.df, Coef. = x, N. = y, Type. = ty, True_vals = Rare_N.prms),
         x = as.character(Unique_coefn.rare_N$Coef), y = as.character(Unique_coefn.rare_N$N))
})

names(Rare_N.MSE.list) <- as.character(unique(Rare_N.df$Type))

Rare_N.MSE.df <- do.call(rbind, lapply(names(Rare_N.MSE.list), function(ty) {
  mse.df <- data.frame(Val = unname(Rare_N.MSE.list[[ty]]),
                       Coef = names(Rare_N.MSE.list[[ty]]),
                       N = rep(paste0("N_", Sampl_effort), each = 4), 
                       Type = ty)
  return(mse.df)
}))

Rare_N.MSE.df$Coef <- factor(Rare_N.MSE.df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Rare_N.MSE.df$Type <- factor(Rare_N.MSE.df$Type, levels = c("random", "proximity", "stratified", "uniform",
                                                            "systematic", "topographic"))

#add rmse
Rare_N.MSE.df$Rmse <- with(Rare_N.MSE.df, sqrt(Val))

Rare_N.RMSE_plot <- ggplot(Rare_N.MSE.df, aes(x = N, y = Rmse, col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(random = "#F4B95AFF", proximity = "#C70E7BFF", stratified = "#007BC3FF",
                                uniform = "#EF7C12FF", systematic = "#FCEA1BFF", topographic = "#009F3FFF"),
                     labels = c(random = "Random", proximity = "Proximity", stratified = "Stratified",
                                uniform = "Uniform", systematic = "Systematic", topographic = "Topographic")) +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50)))
  )) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

##EFFECT OF RANDOMLY SAMPLING THE ENVIRONMENTAL SPACE-------------------------------------------

##Sampling units are randomly collected within the environmental space
##This analysis demonstrates that randomly sampling units within the environmental and geographical space
##leads to the same effect: the over-sampling of the most frequent environmental conditions

#function for randomly sampling the environmental space
Random_of_unif <- function(n, env_df, y) {
  Npres <- T
  while(Npres) {
    Rnd_subs <- env_df[sample(nrow(env_df), size = n, replace = F), ]
    Rnd_sub.xy <- cbind(x = Rnd_subs[["x"]], y = Rnd_subs[["y"]])
    Mod_df <- na.omit(data.frame(extract(x = y, y = Rnd_sub.xy, df = T),
                                 Temp = Rnd_subs[["Bio1"]],
                                 Pr = Rnd_subs[["Bio12"]]))
    Mod_df$ID <- NULL
    colnames(Mod_df)[1] <- c("PA")
    Npres <- (sum(Mod_df$PA == 1) < 30)
  }
  Mod_coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Mod_df))
  return(list(Crds_es = st_coordinates(Rnd_subs), Mod_coef = Mod_coef, Smpsz = nrow(Mod_df)))
}

#simulations ran on Dianthus sperandii
set.seed(9061)
Random_of_unif.res <- lapply(Sampl_effort, function(n) {
  Res <- replicate(n = 500, expr = Random_of_unif(n = n, env_df = Chelsa.AOI.df.sp,
                                                  y = D.sperandii.bin.AOI),
                   simplify = F)
  return(Res)
})

#plot densities
Sampl_dens.df <- do.call(rbind, lapply(seq_along(Sampl_effort), function(i) {
  Crds_envs <- do.call(rbind, lapply(Random_of_unif.res[[i]], '[[', 1))
  Crds_envs <- data.frame(Crds_envs, N = paste0("N_", Sampl_effort[i]))
  return(Crds_envs)
}))

#plot only randomly collected samples from the env_sp
set.seed(4278)
Rnd_smpl_chelsa.aoi <- Chelsa.AOI.df[sample(nrow(Chelsa.AOI.df), 2000, replace = F), ]

Rnd_smpl_chelsa.aoi <- Rnd_smpl_chelsa.aoi[c("Bio1.sc", "Bio12.sc")]

colnames(Rnd_smpl_chelsa.aoi) <- c("X", "Y")

#see: https://rpubs.com/katzkagaya/509701
Plot_env_sp <- ggplot(Chelsa.AOI.df, aes(x = Bio1.sc, y = Bio12.sc)) +
  geom_density_2d_filled() +
  #geom_point(col = "grey", alpha = .3, pch = 1) +
  ylab("Precipitation") + xlab("Temperature") +
  ylim(range(Chelsa.AOI.df$Bio12.sc)) + xlim(range(Chelsa.AOI.df$Bio1.sc)) +
  ggtitle("Environmental space") +
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

Sampl_dens.plots <- lapply(paste0("N_", Sampl_effort), function(n) {
  ggplot(Sampl_dens.df[Sampl_dens.df$N == n, ], aes(x = X, y = Y)) +
    geom_point(data = Rnd_smpl_chelsa.aoi, aes(x = X, y = Y), col = "magenta1", alpha = .1) +
    geom_hdr(method = "kde") +
    ylab(NULL) + xlab(NULL) + 
    ylim(range(Chelsa.AOI.df$Bio12.sc)) + xlim(range(Chelsa.AOI.df$Bio1.sc)) +
    ggtitle(setNames(paste("N", Sampl_effort), paste0("N_", Sampl_effort))[[n]]) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = .5))
})

Sampl_dens.plots[[8]] <- Plot_env_sp

#no need to plot all density plots - in the end, increasing sampling effort only leads to randomly sample more points..
ggarrange(plotlist = Sampl_dens.plots[c(1, 7, 8)], nrow = 1, ncol = 3, common.legend = T, legend = "bottom",
          align = "hv", labels = "auto")

#compare mean squared error associated with random sampling within the geographical and environmental space 

#check smp size - fine
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Random_of_unif.res[[i]], '[[', 3) != Sampl_effort[i]))

#extract coefs
Rnd_unif.mat <- do.call(rbind, lapply(seq_along(Sampl_effort), function(i) {
  Coef_mat <- do.call(rbind, lapply(Random_of_unif.res[[i]], '[[', 2))
  Coef_mat.nm <- colnames(Coef_mat)
  Coef_mat <- as.vector(Coef_mat)
  Coef_mat <- data.frame(Val = Coef_mat, Coef = rep(Coef_mat.nm, each = 500), N = paste0("N_", Sampl_effort[i]))
  return(Coef_mat)
}))

Rnd_unif.mat$True_val <- unname(True_coef.S.nm[Rnd_unif.mat$Coef])

#df with unique combos of Coef and N (from which we need mse)
Rnd_unif.mat.unq <- unique(Rnd_unif.mat[c("Coef", "N")])

Rnd_unif.mse <- mapply(function(x, y) {
  Df <- Rnd_unif.mat[Rnd_unif.mat$Coef == x & Rnd_unif.mat$N == y, c("Val", "True_val")]
  Mse.val <- with(Df, mean((Val - True_val)^2))
  return(Mse.val)
}, x = Rnd_unif.mat.unq$Coef, y = Rnd_unif.mat.unq$N)

Rnd_unif.mse.df <- data.frame(Val = unname(Rnd_unif.mse), Coef = names(Rnd_unif.mse),
                              N = rep(paste0("N_", Sampl_effort), each = 4))

Rnd_unif.mse.df$Type <- "Random_env"

colnames(Rnd_unif.mse.df)[1] <- "Mse_val"

#compare mse with those of random in geographic space and uniform in the env space
#get mse for random geogr and unifor env space

Rnd_g_unif_e <- Mse_df[Mse_df$Type %in% c("Random", "Uniform"), ]

Rnd_unif.mse.df <- rbind(Rnd_unif.mse.df, Rnd_g_unif_e)

Rnd_unif.mse.df$Coef <- factor(Rnd_unif.mse.df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Rnd_unif.mse.df$N <- factor(Rnd_unif.mse.df$N, levels =  unique(Rnd_unif.mse.df$N))

Rnd_unif.RMSE_plot <- ggplot(Rnd_unif.mse.df, aes(x = N, y = sqrt(Mse_val), group = Type, col = Type)) +
  geom_line(alpha = 0.4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(Random = "#F4B95AFF", Uniform = "#EF7C12FF",
                                Random_env = "black"),
                     labels = c(Random = "Random", Uniform = "Uniform",
                                Random_env = "Random in env. space")) +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y", 
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

##INCOMPLETE SAMPLING OF THE ENVIRONMENTAL SPACE-------------------------------------------

##Test the effect of incompletely sampling the environmental space on the uniform approach
##Analysis performed on Dianthus sperandii (we assume effect is the same regardless of species' ecology)
##Two sub-spaces are sampled: all observations below mean temperature, and all observations above mean temperature

#visualise sub-spaces
with(Chelsa.AOI.df, plot(Bio1.sc, Bio12.sc, col = "grey"))

#sub-portions defined by average temperature
with(Chelsa.AOI.df[Chelsa.AOI.df$Bio1.sc > 0, ], points(Bio1.sc, Bio12.sc, col = "red"))
with(Chelsa.AOI.df[Chelsa.AOI.df$Bio1.sc < 0, ], points(Bio1.sc, Bio12.sc, col = "yellow"))

#plot for review
ggplot(Chelsa.AOI.df[Chelsa.AOI.df$Bio1 < mean(Chelsa.AOI.df$Bio1), ], aes(x = Bio1.sc, y = Bio12.sc)) +
  geom_point(col = "black", pch = 6, cex = 2, alpha = .5) +
  geom_point(data = Chelsa.AOI.df[Chelsa.AOI.df$Bio1 >= mean(Chelsa.AOI.df$Bio1), ],
             aes(x = Bio1.sc, Bio12.sc), col = "green", pch = 2, cex = 2, alpha = .5) +
  geom_vline(xintercept = 0, col = "purple", lwd = 1.5) +
  ylab("Precipitation") + xlab("Temperature") +
  theme_pubr() +
  theme(text = element_text(size = 14))

##functions to perform the INcomplete sampling of the environmental space
#these are adapted from those used before for the uniform sampling

#function to sample only sub-spaces of the whole environmental space
Incompl_unif_sampl <- function(x, N, rsl, split_v) {
  #at the right of split_v
  Un_smp.right <- uesampling2.0(sdf = x[x$Bio1 >= split_v, ], grid.res = rsl, n.tr = switch(as.character(N),
                                                                                            '200' = 3, '250' = 4,
                                                                                            '300' = 4, '350' = 5,
                                                                                            '400' = 6, '450' = 6,
                                                                                            '500' = 7))
  #at the lest of split_v
  Un_smp.left <- uesampling2.0(sdf = x[x$Bio1 < split_v, ], grid.res = rsl, n.tr = switch(as.character(N),
                                                                                          '200' = 3, '250' = 4,
                                                                                          '300' = 5, '350' = 6,
                                                                                          '400' = 7, '450' = 8,
                                                                                          '500' = 9))
  #get rid of points > N
  Un_smp.r.dim <- nrow(Un_smp.right)
  Un_smp.l.dim <- nrow(Un_smp.left)
  Un_smp.right <- Un_smp.right[-sample(x = Un_smp.r.dim, size = (Un_smp.r.dim - N), replace = F), ]
  Un_smp.left <- Un_smp.left[-sample(x = Un_smp.l.dim, size = (Un_smp.l.dim - N), replace = F), ]
  st_geometry(Un_smp.right) <- NULL
  st_geometry(Un_smp.left) <- NULL
  Un_smp <- data.frame(rbind(Un_smp.right[c("x", "y")], Un_smp.left[c("x", "y")]), Portion = rep(c("R", "L"), each = N))
  return(Un_smp)
}

#function to run GLMs on data collected in sub-spaces of the whole environmental space
Incomplete_unif <- function(y, x, x_sdf, n = 300, rsl, min_p = 30, split_v) {
  npres <- T
  while(npres) {
    Uniform_points <- Incompl_unif_sampl(x = x_sdf, N = n, rsl = rsl, split_v = split_v)
    Fake_df <- na.omit(data.frame(extract(stack(y, x), Uniform_points[c("x", "y")], df = T),
                                  Portion = Uniform_points[["Portion"]]))
    Fake_df$ID <- NULL
    colnames(Fake_df)[c(1, 2, 3)] <- c("PA", "Temp", "Pr")
    npres <- !( (sum(Fake_df[Fake_df$Portion == "R", "PA"] == 1) >= min_p) && (sum(Fake_df[Fake_df$Portion == "L", "PA"] == 1) >= min_p) )
  }
  Coef.right <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, subset = (Portion == "R")))
  Coef.left <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, subset = (Portion == "L")))
  Coef <- data.frame(rbind(Coef.right, Coef.left), Portion = c("R", "L"))
  return(list(Coef = Coef, SampSz = c('R' = nrow(Fake_df[Fake_df$Portion == "R", ]),
                                      'L' = nrow(Fake_df[Fake_df$Portion == "L", ]))))
}


#launch simulations for incomplete sampling using parallel
cr_inc <- parallel::makeCluster(7)

parallel::clusterExport(cr_inc, c("D.sperandii.bin.AOI", "Chelsa.AOI",
                                  "Chelsa.AOI.df.sp", "Incomplete_unif", 
                                  "Incompl_unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr_inc, list(library(sf), library(raster)))

set.seed(6478)
Incompl_unif.res <- parLapply(cl = cr_inc, X = Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Incomplete_unif(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                   x_sdf = Chelsa.AOI.df.sp, n = n, rsl = 10, min_p = 30,
                                                   split_v = mean(Chelsa.AOI.df.sp$Bio1)), simplify = F)
  return(res)
})

stopCluster(cr_inc)

#check smp sz - fine
sapply(seq_along(Sampl_effort), function(i) {
  smp_df <- do.call(rbind, lapply(Incompl_unif.res[[i]], '[[', 2))
  smp_df <- all(colSums(smp_df != Sampl_effort[i]) == 0)
  return(smp_df)
})

Incompl_unif.mat <- do.call(rbind, lapply(seq_along(Sampl_effort), function(i) {
  Df_res <- do.call(rbind, lapply(Incompl_unif.res[[i]], '[[', 1))
  Df_res$N <- paste0("N_", Sampl_effort[i])
  return(Df_res)
}))

colnames(Incompl_unif.mat)[c(1, 3)] <- c("(Intercept)", "I(Temp^2)")


Incompl_unif.df <-  data.frame(Val = c(Incompl_unif.mat$`(Intercept)`, Incompl_unif.mat$Temp,
                                       Incompl_unif.mat$`I(Temp^2)`, Incompl_unif.mat$Pr),
                               Coef = rep(colnames(Incompl_unif.mat)[1:4], each = nrow(Incompl_unif.mat)), 
                               Portion = rep(Incompl_unif.mat$Portion, 4), 
                               N = rep(Incompl_unif.mat$N, 4))

Incompl_unif.df$True_val <- unname(True_coef.S.nm[Incompl_unif.df$Coef])

ggplot(Incompl_unif.df, aes(x = N, y = Val, col = Portion)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val)) +
  facet_wrap(~ Coef, scales = "free_y")

#compute MSE
Incompl_unif.MSE <- do.call(rbind, lapply(unique(Incompl_unif.df$Portion), function(prt) {
  df <- Incompl_unif.df[Incompl_unif.df$Portion == prt, ]
  mse.val <- mapply(function(x, y) {
    mse.val <- df[df$Coef == x & df$N == y, c("Val", "True_val")]
    mse.val <- with(mse.val, mean((Val - True_val)^2))
    return(mse.val)
  }, x = unique(Incompl_unif.df[c("Coef", "N")])[[1]],
  y = unique(Incompl_unif.df[c("Coef", "N")])[[2]])
  mse.val <- data.frame(Mse_val = unname(mse.val), Coef = names(mse.val),
                        N = unique(Incompl_unif.df[c("Coef", "N")])[[2]],
                        Portion = prt)
  return(mse.val)
}))

#add mse vals for env space (complete)
Incompl_unif.MSE <- rbind(Incompl_unif.MSE,
                          data.frame(Rnd_g_unif_e[Rnd_g_unif_e$Type %in% "Uniform", c(1, 2, 3)], Portion = "Full"))

Incompl_unif.MSE$Coef <- factor(Incompl_unif.MSE$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Incompl_unif.MSE$N <- factor(Incompl_unif.MSE$N, levels = paste0("N_", Sampl_effort))

Incompl_unif.plot <- ggplot(Incompl_unif.MSE, aes(x = N, y = sqrt(Mse_val), group = Portion, col = Portion)) +
  geom_line(alpha = 0.4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(R = "green", L = "black",
                                Full = "#EF7C12FF"),
                     labels = c(R = "Right", L = "Left",
                                Full = "Full")) +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y", 
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

##TEST EFFECT OF SAMPLING UNITS CONSISTENTLY INCLUDED IN DATASETS FOR GLMs-----------------

library(dplyr)

##A (low) amount of sampling units is consistently included in the datasets used to fit the GLMs associated with the uniform sampling
##This happens as the density of sampling units is lower at the boundary of the environmental space
##There, the number of points within a cell of the sampling grid is lower than the number of points to be collected
##across the sampling grid
##This issue tends to downwardly bias the variance estimator for the regression coefficients (estimators)
##The following analysis investigates the effect of the sampling units consistently included in the datasets
##on the simulation results for the uniform approach

#function to extract 'parametric' variance
Param_var.unif <- function(y, x, x_sdf, n = 300, rsl, min_p = 30) {
  npres <- T
  while(npres) {
    Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
    Fake_df <- na.omit(extract(stack(y, x), Uniform_points, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p)
  }
  VarCov <- diag(vcov(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df)))
  return(list(VCov = VarCov, Smpsz = nrow(Fake_df)))
}

#Simulations for Dianthus sperandii--------------------------------------------------------------

cr_prv <- parallel::makeCluster(7)

parallel::clusterExport(cr_prv, c("D.sperandii.bin.AOI", "Chelsa.AOI",
                                  "Chelsa.AOI.df.sp", "Param_var.unif", 
                                  "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr_prv, list(library(sf), library(raster)))

set.seed(4761)
Param_var.res <- parLapply(cl = cr_prv, X = Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Param_var.unif(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                  x_sdf = Chelsa.AOI.df.sp, n = n, rsl = 10, min_p = 30),
                   simplify = F)
  return(res)
})

stopCluster(cr_prv)

#check smpsz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Param_var.res[[i]], "[[", "Smpsz") == Sampl_effort[i]))

#compute average of variances
Param_var.avg <- do.call(rbind, lapply(seq_along(Sampl_effort), function(i) {
  df <- colMeans(do.call(rbind, lapply(Param_var.res[[i]], '[[', "VCov")))
  df <- data.frame(matrix(df, nrow = 1, dimnames = list("1", names(df))), N = paste0("N_", Sampl_effort[i]))
  return(df)
}))

colnames(Param_var.avg)[c(1, 3)] <- c("(Intercept)", "I(Temp^2)")

Param_var.avg <- data.frame(Var = c(Param_var.avg$`(Intercept)`, Param_var.avg$Temp, Param_var.avg$`I(Temp^2)`, Param_var.avg$Pr),
                            Coef = rep(colnames(Param_var.avg)[1:4], each = nrow(Param_var.avg)),
                            N = rep(Param_var.avg$N, 4))


#extract data for comparison
Var_unif_e <- Var_df[c("Var", "Coef", "N", "Type")]

#add missing col
Param_var.avg$Type <- "Param_var"

Param_var.avg <- rbind(Param_var.avg, Var_unif_e)

ggplot(Param_var.avg, aes(x = N, y = sqrt(Var), group = Type, col = Type)) +
  geom_line() +
  scale_color_manual(values = setNames(c("green", "grey", "grey", "grey", "red", "grey", "grey"),
                                       unique(Param_var.avg$Type))) +
  facet_wrap(~ Coef, scales = "free_y") +
  theme_pubr()

#need to add bias^2 here (from bias^2)!

Param_var.avg <- dplyr::left_join(Param_var.avg, Bias_df, by = c("Coef", "N", "Type"))

for(cf in unique(Param_var.avg$Coef)) {
  for(n in unique(Param_var.avg$N)) {
    Param_var.avg[Param_var.avg$Type == "Param_var" & Param_var.avg$Coef == cf & Param_var.avg$N == n, "Bias"] <- Param_var.avg[Param_var.avg$Type == "Uniform" & Param_var.avg$Coef == cf & Param_var.avg$N == n, "Bias"]
  }
}

#square bias!
Param_var.avg$Bias_sq <- with(Param_var.avg, Bias^2)

#now compute mse
Param_var.avg$Mse <- with(Param_var.avg, Var + Bias_sq)

#add column for plotting
Param_var.avg$Type_simple <- Param_var.avg$Type

Param_var.avg$Type_simple <- with(Param_var.avg, ifelse(Type_simple == "Uniform", "Uniform_unc",
                                                        ifelse(Type_simple == "Param_var", "Uniform_cor",
                                                               "Other")))

Param_var.avg$Type_simple <- factor(Param_var.avg$Type_simple, levels = c("Uniform_unc", "Uniform_cor", "Other"))
Param_var.avg$Coef <- factor(Param_var.avg$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

Param_var_S <- ggplot(Param_var.avg, aes(x = N, y = sqrt(Mse), group = Type, col = Type_simple)) +
  geom_line(alpha = 0.4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = setNames(c("magenta", "blue", "grey"),
                                       levels(Param_var.avg$Type_simple)), 
                     labels = setNames(c("Uniform uncorr.", "Uniform corr.", "Other"),
                                       levels(Param_var.avg$Type_simple)), name = "Type") +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y", 
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#Simulations for Dianthus tundrae--------------------------------------------------------------

cr_prv.r <- parallel::makeCluster(7)

parallel::clusterExport(cr_prv.r, c("D.tundrae.bin.AOI", "Chelsa.AOI",
                                    "Chelsa.AOI.df.sp", "Param_var.unif", 
                                    "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr_prv.r, list(library(sf), library(raster)))

set.seed(8475) 
Param_var.res.rare <- parLapply(cl = cr_prv.r, X = Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Param_var.unif(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                  x_sdf = Chelsa.AOI.df.sp, n = n, rsl = 10, min_p = 30),
                   simplify = F)
  return(res)
})

stopCluster(cr_prv.r)

#check smpsz
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Param_var.res.rare[[i]], "[[", "Smpsz") == Sampl_effort[i]))

#compute average of variances
Param_var.avg.rare <- do.call(rbind, lapply(seq_along(Sampl_effort), function(i) {
  df <- colMeans(do.call(rbind, lapply(Param_var.res.rare[[i]], '[[', "VCov")))
  df <- data.frame(matrix(df, nrow = 1, dimnames = list("1", names(df))), N = paste0("N_", Sampl_effort[i]))
  return(df)
}))

colnames(Param_var.avg.rare)[c(1, 3)] <- c("(Intercept)", "I(Temp^2)")

Param_var.avg.rare <- data.frame(Var = c(Param_var.avg.rare$`(Intercept)`, Param_var.avg.rare$Temp, Param_var.avg.rare$`I(Temp^2)`, Param_var.avg.rare$Pr),
                                 Coef = rep(colnames(Param_var.avg.rare)[1:4], each = nrow(Param_var.avg.rare)),
                                 N = rep(Param_var.avg.rare$N, 4))

#extract data for comparison
Var_unif_e.rare <- Var_rare_df[c("Var", "Coef", "N", "Type")]

#add missing col
Param_var.avg.rare$Type <- "Param_var"

Param_var.avg.rare <- rbind(Param_var.avg.rare, Var_unif_e.rare)

#need to add bias^2 here (from bias^2)!
Param_var.avg.rare <- dplyr::left_join(Param_var.avg.rare, Bias_rare_df, by = c("Coef", "N", "Type"))

for(cf in unique(Param_var.avg.rare$Coef)) {
  for(n in unique(Param_var.avg.rare$N)) {
    Param_var.avg.rare[Param_var.avg.rare$Type == "Param_var" & Param_var.avg.rare$Coef == cf & Param_var.avg.rare$N == n, "Bias"] <- Param_var.avg.rare[Param_var.avg.rare$Type == "Uniform" & Param_var.avg.rare$Coef == cf & Param_var.avg.rare$N == n, "Bias"]
  }
}

#square bias!
Param_var.avg.rare$Bias_sq <- with(Param_var.avg.rare, Bias^2)

#now compute mse
Param_var.avg.rare$Mse <- with(Param_var.avg.rare, Var + Bias_sq)

#add column for plotting
Param_var.avg.rare$Type_simple <- Param_var.avg.rare$Type

Param_var.avg.rare$Type_simple <- with(Param_var.avg.rare, ifelse(Type_simple == "Uniform", "Uniform_unc",
                                                                  ifelse(Type_simple == "Param_var", "Uniform_cor",
                                                                         "Other")))

Param_var.avg.rare$Type_simple <- factor(Param_var.avg.rare$Type_simple, levels = c("Uniform_unc", "Uniform_cor", "Other"))
Param_var.avg.rare$Coef <- factor(Param_var.avg.rare$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))

Param_var_T <- ggplot(Param_var.avg.rare, aes(x = N, y = sqrt(Mse), group = Type, col = Type_simple)) +
  geom_line(alpha = 0.4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = setNames(c("magenta", "blue", "grey"),
                                       levels(Param_var.avg.rare$Type_simple)), 
                     labels = setNames(c("Uniform uncorr.", "Uniform corr.", "Other"),
                                       levels(Param_var.avg.rare$Type_simple)), name = "Type") +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y", 
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

