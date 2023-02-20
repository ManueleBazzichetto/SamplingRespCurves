##Code to replicate results presented in the manuscript "Sampling strategies matter to accurately estimate response curves’ parameters in species distribution models"
##Simulations for Dianthus sperandii (virtual species with wide distribution)

##!!Please, load packages reported at the top of the "Bio_Elev_data" script!!

#Set value of regression parameters - named vector of parameters created after random sampling (see below)
True_coef.S <- c(-11.3886498, 1.3835909, -0.0920092, 0.0071542)

#Generate layers of D. sperandii distribution
D.sperandii.link <- -11.3886498 + 1.3835909*Chelsa.stack$Bio1 -0.0920092*(Chelsa.stack$Bio1^2) + 0.0071542*Chelsa.stack$Bio12
D.sperandii.prob <- calc(D.sperandii.link, plogis)
#Get layer of 1/0
#warnings are for NAs
set.seed(1984)
D.sperandii.bin <- calc(D.sperandii.prob, fun = function(x) {rbinom(1, 1, x)})

#See distribution
mapview(D.sperandii.bin)

#Mask layers of species distribution to only consider AOI
D.sperandii.prob.AOI <- mask(D.sperandii.prob, mask = as(Elev_AOI, "Spatial"))
D.sperandii.bin.AOI <- mask(D.sperandii.bin, mask = as(Elev_AOI, "Spatial"))


#SIMULATIONS----------------------------------------------------------------------------------------------------------------
Sampl_effort <- seq(from = 200, to = 500, by = 50)

#RANDOM sampling

#get coords of cells that are not NAs
Climate_crds <- as.data.frame(Chelsa.AOI, xy = T, na.rm = T)[c("x", "y")]

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
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

set.seed(1930)
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

#check sample size - ok
sapply(seq_along(Random_res), function(.) {
  sum(vapply(Random_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
})

Random_mats <- do.call(rbind, lapply(names(Random_mats), function(nm) {
  df <- Random_mats[[nm]]
  df_nm <- colnames(df)
  dim(df) <- NULL
  df <- data.frame(Val = df, Coef = rep(df_nm, each = 500), N = nm)
  return(df)
})
)

True_coef.S.nm <- setNames(True_coef.S, nm = unique(Random_mats$Coef))

Random_mats$True_val <- unname(True_coef.S.nm[Random_mats$Coef])
Random_mats$Coef <- factor(Random_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Random_mats$N <- factor(Random_mats$N, levels = paste("N", Sampl_effort, sep = "_"))

#Get an idea of the distribution of coefficients' values
ggplot(Random_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 


#STRATIFIED sampling
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
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

set.seed(1956)
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

#check sample size - ok
sapply(seq_along(Strat_res), function(.) {
  sum(vapply(Strat_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
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
Strat_mats$True_val <- unname(True_coef.S.nm[as.character(Strat_mats$Coef)])
Strat_mats$N <- factor(Strat_mats$N, levels = paste("N", Sampl_effort, sep = "_"))

ggplot(Strat_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#ROAD PROXIMITY sampling

Proximity_fit <- function(y, x, prox_lay, n = 300, min_p = 30) {
  Prox_df.full <- na.omit(as.data.frame(prox_lay, xy = T)) #replaced Prox_df with Prox_df.full to avoid overwriting in while
  npres <- T
  while(npres) {
    Prox_df <- Prox_df.full[sample(nrow(Prox_df.full), size = n, replace = F, prob = Prox_df.full$layer), ] #check
    Prox_df$layer <- NULL
    Fake_df <- na.omit(extract(stack(y, x), Prox_df, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == n)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

set.seed(1986)
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

#check sample size - ok
sapply(seq_along(Prox_res), function(.) {
  sum(vapply(Prox_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
})

Prox_res <- do.call(rbind, lapply(names(Prox_res), function(nm) {
  Mat <- do.call(rbind, lapply(Prox_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef.S.nm[Mat$Coef]
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

#UNIFORM sampling

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
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

#go parallel
cr7 <- parallel::makeCluster(7)

parallel::clusterExport(cr7, c("D.sperandii.bin.AOI", "Chelsa.AOI",
                               "Chelsa.AOI.df.sp", "Uniform_fit", 
                               "Unif_sampl", "uesampling2.0", "Sampl_effort"))

parallel::clusterEvalQ(cr7, list(library(sf), library(raster)))

set.seed(1953)
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

#check sample size - ok
sapply(seq_along(Unif_res), function(.) {
  sum(vapply(Unif_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
})

Unif_res <- do.call(rbind, lapply(names(Unif_res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef.S.nm[Mat$Coef]
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

#SYSTEMATIC sampling

Systematic_fit <- function(y, x, N, perc_inc, poly_proj, min_p = 30) {
  npres <- T
  while(npres) {
    Pts_syst <- st_sample(x = poly_proj, size = N + floor(N*perc_inc), type = "regular")
    Pts_syst <- st_transform(Pts_syst, crs = 4326)
    Pts_syst <- st_coordinates(Pts_syst)
    Fake_df <- na.omit(data.frame(extract(y, Pts_syst, df = T), extract(x, Pts_syst)))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == N) 
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

#I used D.sperandii.bin and Chelsa.stack.syst to solve issue with points regularly sampled within Elev_AOI.proj
#ending up in NA cells for D.sperandii.bin.AOI and Chelsa.AOI
#This solves the thing for Chelsa.stack.syst, but not for D.sperandii.bin
#I increase a bit (7%) sampling effort (e.g. X points for each sampling effort). Excess is removed in the function

set.seed(2022)
Syst_res <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 500, expr = Systematic_fit(y = D.sperandii.bin, x = Chelsa.stack.syst,
                                                  N = n, perc_inc = .07, poly_proj = Elev_AOI.proj, min_p = 30), simplify = F)
  return(res)
})

#check variability
#hist((do.call(rbind, lapply(Syst_res[[3]], '[[', 1)))[, 1])

names(Syst_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor
Syst_cor <- sapply(Syst_res, function(.) mean(vapply(., function(i) i[[2]], FUN.VALUE = numeric(1))))

#check sample size - ok
sapply(seq_along(Syst_res), function(.) {
  sum(vapply(Syst_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
})

Syst_res <- do.call(rbind, lapply(names(Syst_res), function(nm) {
  Mat <- do.call(rbind, lapply(Syst_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef.S.nm[Mat$Coef]
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

#TOPOGRAPHIC sampling

Topo_fit <- function(y, x, N, perc_incr, topo_layer, min_p = 30) { 
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T, na.rm = T)
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N + floor(perc_incr*N), replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == N)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, SampSz = nrow(Fake_df)))
}

set.seed(1919)
Topo_res <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 500, expr = Topo_fit(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                            N = n., perc_incr = .07,
                                            topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
})

names(Topo_res) <- paste("N", Sampl_effort, sep = "_")

#extract cor
Topo_cor <- sapply(Topo_res, function(.) mean(vapply(., function(i) i[[2]], FUN.VALUE = numeric(1))))

#check sample size - ok
sapply(seq_along(Topo_res), function(.) {
  sum(vapply(Topo_res[[.]], '[[', 3, FUN.VALUE = numeric(1)) != Sampl_effort[.])
})

Topo_res <- do.call(rbind, lapply(names(Topo_res), function(nm) {
  Mat <- do.call(rbind, lapply(Topo_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500), N = nm)
  Mat$True_val <- True_coef.S.nm[Mat$Coef]
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


#GET AND COMPARE RESULTS----------------------------------------------------------------------------------------------

#Big data.frame with estimators of regression parameters for the different sampling efforts and strategies
Big_comp <- rbind(data.frame(Random_mats, Type = "Random"),
                  data.frame(Strat_mats, Type = "Strat"),
                  data.frame(Prox_res, Type = "Prox"),
                  data.frame(Unif_res, Type = "Uniform"),
                  data.frame(Syst_res, Type = "Systematic"),
                  data.frame(Topo_res, Type = "Topographic"))

Big_comp$Type <- factor(Big_comp$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                  "Systematic", "Topographic"))

#Comparison among sampling strategies
ggplot(Big_comp, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#function for computing mean squared error (MSE)
Mse_samp <- function(df, Coef., N., Type., True_vals) {
  Vals <- df[df$Coef == Coef. & df$N == N. & df$Type == Type., "Val"]
  Mse <- mean((Vals - True_vals[Coef.])^2)
  return(Mse)
}

#get MSE
Mse_list <- lapply(as.character(unique(Big_comp$Type)), function(ty) {
  mapply(function(x, y) Mse_samp(df = Big_comp, Coef. = x, N. = y, Type. = ty, True_vals = True_coef.S.nm),
         x = rep(names(True_coef.S.nm), 7), y = rep(paste("N", Sampl_effort, sep = "_"), each = 4))
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

Mse_plot <- ggplot(Mse_df, aes(x = N, y = Mse_val, col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(Random = "#F4B95AFF", Prox = "#C70E7BFF", Strat = "#007BC3FF",
                                Uniform = "#EF7C12FF", Systematic = "#FCEA1BFF", Topographic = "#009F3FFF"),
                     labels = c(Random = "Random", Prox = "Proximity", Strat = "Stratified",
                                Uniform = "Uniform", Systematic = "Systematic", Topographic = "Topographic")) +
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50)))
  )) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#change for rev
RootMSE_plot <- ggplot(Mse_df, aes(x = N, y = sqrt(Mse_val), col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
  scale_color_manual(values = c(Random = "#F4B95AFF", Prox = "#C70E7BFF", Strat = "#007BC3FF",
                                Uniform = "#EF7C12FF", Systematic = "#FCEA1BFF", Topographic = "#009F3FFF"),
                     labels = c(Random = "Random", Prox = "Proximity", Strat = "Stratified",
                                Uniform = "Uniform", Systematic = "Systematic", Topographic = "Topographic")) +
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

#compute bias
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
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#compute variance
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
  scale_x_discrete(labels = setNames(as.character(seq(200, 500, 50)),
                                     paste0('N_', as.character(seq(200, 500, 50))))) +
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("SQRT(Variance)") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))


#PREDICTED VS TRUE SPECIES RESPONSE CURVES--------------------------------------------------------------------------------------

#Compute true probabilities for D. sperandii
True_pred_sperandii.2 <- data.frame(var_val = c(Bio1_seq.AOI, Bio12_seq.AOI),
                                    var = rep(c("Bio1", "Bio12"),
                                              each = 100))

#This creates a vector of true probabilities - computations carried out in matrix form
True_pred_sperandii.2$Proba <- unlist(lapply(c("Bio1", "Bio12"), function(nm) {
  if(nm %in% "Bio1") {
    mat.prd <- cbind(1, True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T],
                     (True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T]^2),
                     mean(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio12", "var_val", drop = T]))
    #colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef.S.nm))
    return(pred_val)
  } else {
    mat.prd <- cbind(1, mean(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T]),
                     median(True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio1", "var_val", drop = T])^2,
                     True_pred_sperandii.2[True_pred_sperandii.2$var == "Bio12", "var_val", drop = T])
    #colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef.S.nm))
    return(pred_val)
  }
}))

#RANDOM
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
  return(list(Map = Map, Part_fit = Partial_fit, SampSz = nrow(Fake_df)))
}

set.seed(2890)
Random_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Random_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI, 
                                                  x_crds =  Climate_crds, n = N, min_p = 30), simplify = F)
  return(res)
})

#check sample size
vapply(seq_along(Random_map), function(.) sum(sapply(Random_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

names(Random_map) <- paste("N", Sampl_effort, sep = "_")

#Prediction maps for the different sampling efforts
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

#STRATIFIED

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
  return(list(Map = Map, Part_fit = Partial_fit, SampSz = nrow(Fake_df)))
}

set.seed(2931)
Strat_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Strat_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI, n = N,
                                                 strata = Bio1_12.rcl.qrt.AOI, min_p = 30), simplify = F)
  return(res)
})

#check sample size
vapply(seq_along(Strat_map), function(.) sum(sapply(Strat_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

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

#ROAD PROXIMITY

Proximity_mapping <- function(y, x, prox_lay, n = 300, min_p = 30) {
  Prox_df.full <- na.omit(as.data.frame(prox_lay, xy = T))
  npres <- T
  while(npres) {
    Prox_df <- Prox_df.full[sample(nrow(Prox_df.full), size = n, replace = F, prob = Prox_df.full$layer), ] #check
    Prox_df$layer <- NULL
    Fake_df <- na.omit(extract(stack(y, x), Prox_df, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == n)
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
  return(list(Map = Map, Part_fit = Partial_fit, SampSz = nrow(Fake_df)))
}

set.seed(4094)
Prox_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Proximity_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                     prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

#check sample size
vapply(seq_along(Prox_map), function(.) sum(sapply(Prox_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

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

#UNIFORM

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
  return(list(Map = Map, Part_fit = Partial_fit, SamplSz = nrow(Fake_df)))
}

set.seed(2879)
Unif_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Uniform_mapping(y = D.sperandii.bin.AOI,
                                                   x = Chelsa.AOI,
                                                   x_sdf = Chelsa.AOI.df.sp,
                                                   n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

#check sample size
vapply(seq_along(Unif_map), function(.) sum(sapply(Unif_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

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

#SYSTEMATIC

Systematic_mapping <- function(y, x, N, perc_inc, poly_proj, min_p = 30) {
  npres <- T
  while(npres) {
    Pts_syst <- st_sample(x = poly_proj, size = N + floor(N*perc_inc), type = "regular")
    Pts_syst <- st_transform(Pts_syst, crs = 4326)
    Pts_syst <- st_coordinates(Pts_syst)
    Fake_df <- na.omit(data.frame(extract(y, Pts_syst, df = T), extract(x, Pts_syst)))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == N)
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
  return(list(Map = Map, Part_fit = Partial_fit, SamplSz = nrow(Fake_df)))
}

set.seed(2948)
Syst_map <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 100, expr = Systematic_mapping(y = D.sperandii.bin.AOI, x = Chelsa.stack.syst,
                                                      N = n, perc_inc = 0.07, poly_proj = Elev_AOI.proj, min_p = 30), simplify = F)
  return(res)
})

#check sample size
vapply(seq_along(Syst_map), function(.) sum(sapply(Syst_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

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

#TOPOGRAPHIC
Topo_mapping <- function(y, x, N, perc_inc, topo_layer, min_p = 30) {
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T, na.rm = T)
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N + floor(N*perc_inc), replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df)[1] <- "PA"
    npres <- !(sum(Fake_df$PA == 1) >= min_p && nrow(Fake_df) == N)
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
  return(list(Map = Map, Part_fit = Partial_fit, SampSz = nrow(Fake_df)))
}

set.seed(2441)
Topo_map <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 100, expr = Topo_mapping(y = D.sperandii.bin.AOI, x = Chelsa.AOI,
                                                N = n., perc_inc = 0.07, topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
})

#check sample size
vapply(seq_along(Topo_map), function(.) sum(sapply(Topo_map[[.]], '[[', 3) != Sampl_effort[.]), FUN.VALUE = numeric(1))

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

#Figure showing predicted vs true response curves
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