##Code to replicate results presented in the manuscript "Sampling strategy matters to accurately estimate response curves’ parameters in species distribution models"
##Simulations for Dianthus tundrae (virtual species with restricted distribution)

##!!Please, load packages reported at the top of the "Bio_Elev_data" script!!

#Set value of regression parameters
True_coef.T.nm <- setNames(c(-22.17261596, 3.77916980, -0.20996874, 0.00490439), nm = names(True_coef.S.nm))

#Generate layers of D. tundrae distribution
D.tundrae.link <- -22.17261596 + 3.77916980*Chelsa.stack$Bio1 -0.20996874*(Chelsa.stack$Bio1^2) + 0.00490439*Chelsa.stack$Bio12
D.tundrae.prob <- calc(D.tundrae.link, plogis)

#Get layer of 1/0
#warnings are for NAs
set.seed(100)
D.tundrae.bin <- calc(D.tundrae.prob, fun = function(x) {rbinom(1, 1, x)})

#See distribution
mapview(D.tundrae.bin)

#Mask layers of species distribution to only consider AOI
D.tundrae.bin.AOI <- mask(D.tundrae.bin, mask = Elev_AOI)


#SIMULATIONS----------------------------------------------------------------------------------------------------------------

#RANDOM sampling

Random_rare <- function(y, x, x_crds, n = 300, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  npres <- T
  while(npres) {
    Random_plots <- x_crds[sample(nrow(x_crds), n, replace = F), ]
    Fake_df <- na.omit(extract(stack(y, x), Random_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

set.seed(76)
Rnd_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Random_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                               x_crds = Climate_crds, n = N, min_p = 30), simplify = F)
  return(res)
})

#check separation - fine
lapply(Rnd_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Rnd_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Rnd_rare_res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#check sample size - fine
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Rnd_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Rnd_rare_res) <- paste0("N_", Sampl_effort)

Rnd_res_mats <- do.call(rbind, lapply(names(Rnd_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Rnd_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Rnd_res_mats$Coef <- factor(Rnd_res_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Rnd_res_mats$True_val <- unname(True_coef.T.nm[as.character(Rnd_res_mats$Coef)])
Rnd_res_mats$N <- factor(Rnd_res_mats$N, levels = paste0("N_", Sampl_effort))

#Get an idea of distribution of coefficients' value
ggplot(Rnd_res_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#STRATIFIED sampling

Strat_rare <- function(y, x, n = 300, strata, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  npres <- T
  while(npres) {
    Strat_plots <- Strat_raster(x = strata, N = n)
    Fake_df <- na.omit(extract(stack(y, x), Strat_plots, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

set.seed(12)
Str_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Strat_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                              strata = Bio1_12.rcl.qrt.AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

#check separation - fine
lapply(Str_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Str_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Str_rare_res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#check sample size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Str_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Str_rare_res) <- paste0("N_", Sampl_effort)

Str_rare_mats <- do.call(rbind, lapply(names(Str_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Str_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Str_rare_mats$Coef <- factor(Str_rare_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Str_rare_mats$True_val <- unname(True_coef.T.nm[as.character(Str_rare_mats$Coef)])
Str_rare_mats$N <- factor(Str_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Str_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#ROAD PROXIMITY sampling

Proximity_rare <- function(y, x, prox_lay, n = 300, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  Prox_df.full <- na.omit(as.data.frame(prox_lay, xy = T)) #here changed Prox_df to Prox_df.full
  npres <- T
  while(npres) {
    Prox_df <- Prox_df.full[sample(nrow(Prox_df.full), size = n, replace = F, prob = Prox_df.full$layer), ] #check
    Prox_df$layer <- NULL
    Fake_df <- na.omit(extract(stack(y, x), Prox_df, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p && nrow(Fake_df) == n)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

set.seed(37)
Prx_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Proximity_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                  prox_lay = Abr_highway_AOI, n = N, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation - fine
lapply(Prx_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Prx_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Prx_rare_res, function(.) colMeans(do.call(rbind, lapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  c(N_pr = Tab[Tab$Var1 == "1", "Freq"], N_ab = Tab[Tab$Var1 == "0", "Freq"])
}))))

#check sample size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Prx_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Prx_rare_res) <- paste0("N_", Sampl_effort)

Prx_rare_mats <- do.call(rbind, lapply(names(Prx_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Prx_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Prx_rare_mats$Coef <- factor(Prx_rare_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Prx_rare_mats$True_val <- unname(True_coef.T.nm[as.character(Prx_rare_mats$Coef)])
Prx_rare_mats$N <- factor(Prx_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Prx_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#UNIFORM sampling

Uniform_rare <- function(y, x, x_sdf, n = 300, rsl, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  npres <- T
  while(npres) {
    Uniform_points <- Unif_sampl(x = x_sdf, N = n, rsl = rsl)
    Fake_df <- na.omit(extract(stack(y, x), Uniform_points, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

cr7 <- parallel::makeCluster(7)

parallel::clusterExport(cr7, c("D.tundrae.bin.AOI", "Chelsa.AOI", "Chelsa.AOI.df.sp",
                               "Unif_sampl", "Uniform_rare", "uesampling2.0",
                               "Sampl_effort"))

parallel::clusterEvalQ(cr7, list(library(brglm2), library(detectseparation), library(raster), library(sf)))

set.seed(579)
Unif_rare_res <- parLapply(cl = cr7, X = Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Uniform_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                x_sdf = Chelsa.AOI.df.sp,
                                                n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

stopCluster(cr7)

#check separation - fine
lapply(Unif_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Unif_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Unif_rare_res, function(.) {
  colMeans(do.call(rbind, lapply(., function(i) {
    Tab <- as.data.frame(i[[3]])
    c(N_pr = Tab[Tab$Var1 == "1", "Freq"], N_abs = Tab[Tab$Var1 == "0", "Freq"])
  })))
})

#check sample size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Unif_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Unif_rare_res) <- paste0("N_", Sampl_effort)

Unif_rare_mat <- do.call(rbind, lapply(names(Unif_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Unif_rare_mat$Coef <- factor(Unif_rare_mat$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Unif_rare_mat$True_val <- unname(True_coef.T.nm[as.character(Unif_rare_mat$Coef)])
Unif_rare_mat$N <- factor(Unif_rare_mat$N, levels = paste0("N_", Sampl_effort))

ggplot(Unif_rare_mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#SYSTEMATIC sampling

Systematic_rare <- function(y, x, N, perc_incr, poly_proj, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  npres <- T
  while(npres) {
    Pts_syst <- st_sample(x = poly_proj, size = N + floor(N*perc_incr), type = "regular")
    Pts_syst <- st_transform(Pts_syst, crs = 4326)
    Pts_syst <- st_coordinates(Pts_syst)
    Fake_df <- na.omit(data.frame(extract(y, Pts_syst, df = T), extract(x, Pts_syst)))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p && nrow(Fake_df) == N)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

#use D.tundrae.bin and Chelsa.stack.syst as done for D. sperandii to alleviate problem with NAs
set.seed(248)
Syst_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Systematic_rare(y = D.tundrae.bin, x = Chelsa.stack.syst,
                                                   poly_proj = Elev_AOI.proj, N = N, perc_incr = .07, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation - fine
lapply(Syst_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Syst_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Syst_rare_res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#check sample size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Syst_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Syst_rare_res) <- paste0("N_", Sampl_effort)

Syst_rare_mats <- do.call(rbind, lapply(names(Syst_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Syst_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Syst_rare_mats$Coef <- factor(Syst_rare_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Syst_rare_mats$True_val <- unname(True_coef.T.nm[as.character(Syst_rare_mats$Coef)])
Syst_rare_mats$N <- factor(Syst_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Syst_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#TOPOGRAPHIC sampling

Topo_rare <- function(y, x, N, topo_layer, perc_incr, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T, na.rm = T)
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N + floor(N*perc_incr), replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    if(nrow(Fake_df) > N) {
      Fake_df <- Fake_df[-sample(nrow(Fake_df), size = (nrow(Fake_df) - N), replace = F), ]
    }
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p && nrow(Fake_df) == N)
  }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation, SampSz = nrow(Fake_df)))
}

set.seed(1808)
Topo_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Topo_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                             topo_layer = Topogr_het_Abr.AOI, N = N, perc_incr = .07, min_p = 30),
                   simplify = F)
  return(res)
})

#check separation - fine
lapply(Topo_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Topo_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Topo_rare_res, function(.) mean(sapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  N_pr <- Tab[Tab$Var1 == "1", "Freq"]
})))

#check sample size
sapply(seq_along(Sampl_effort), function(i) sum(sapply(Topo_rare_res[[i]], '[[', 5) != Sampl_effort[i]))

names(Topo_rare_res) <- paste0("N_", Sampl_effort)

Topo_rare_mats <- do.call(rbind, lapply(names(Topo_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Topo_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Topo_rare_mats$Coef <- factor(Topo_rare_mats$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Topo_rare_mats$True_val <- unname(True_coef.T.nm[as.character(Topo_rare_mats$Coef)])
Topo_rare_mats$N <- factor(Topo_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Topo_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 


#GET AND COMPARE RESULTS----------------------------------------------------------------------------------------------

#Big data.frame with estimators of regression parameters for the different sampling efforts and strategies
Big_comp.rare <- rbind(data.frame(Rnd_res_mats, Type = "Random"),
                       data.frame(Str_rare_mats, Type = "Strat"),
                       data.frame(Prx_rare_mats, Type = "Prox"),
                       data.frame(Syst_rare_mats, Type = "Systematic"),
                       data.frame(Unif_rare_mat, Type = "Uniform"),
                       data.frame(Topo_rare_mats, Type = "Topographic"))

Big_comp.rare$Type <- factor(Big_comp.rare$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                            "Systematic", "Topographic"))

#Comparison among sampling strategies
ggplot(Big_comp.rare, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#get MSE
Mse_rare <- lapply(as.character(unique(Big_comp.rare$Type)), function(ty) {
  mapply(function(x, y) Mse_samp(df = Big_comp.rare, Coef. = x, N. = y, Type. = ty, True_vals = True_coef.T.nm),
         x = rep(names(True_coef.T.nm), 7), y = rep(paste("N", Sampl_effort, sep = "_"), each = 4))
})

names(Mse_rare) <- as.character(unique(Big_comp.rare$Type))

Mse_rare_df <- do.call(rbind, lapply(names(Mse_rare), function(nm) {
  nms <- names(Mse_rare[[nm]])
  df <- data.frame(Mse_val = unname(Mse_rare[[nm]]), Coef = nms,
                   N = rep(paste("N", Sampl_effort, sep = "_"), each = 4), 
                   Type = nm)
  return(df)
}))

Mse_rare_df$Coef <- factor(Mse_rare_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Mse_rare_df$Type <- factor(Mse_rare_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                        "Systematic", "Topographic"))

Mse_rare_plot <- ggplot(Mse_rare_df, aes(x = N, y = Mse_val, col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
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
  ylab("Mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#change for rev
RootMSE_rare_plot <- ggplot(Mse_rare_df, aes(x = N, y = sqrt(Mse_val), col = Type, group = Type)) +
  geom_line(alpha = .4) +
  geom_point(position = position_dodge2(width = .1), cex = 2.5, alpha = .7) +
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
  ylab("Root mean squared error") + xlab("Sampling effort") +
  theme_pubr() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#compute bias
Big_rare.copy <- Big_comp.rare

Big_rare.copy$Bias <- with(Big_rare.copy, Val - True_val)

Bias_rare_arr <- with(Big_rare.copy, tapply(Bias, INDEX = list(Type, N, Coef), mean)) 

Bias_rare_df <- data.frame(do.call(rbind, lapply(1:(dim(Bias_rare_arr)[3]), function(i) {
  df <- Bias_rare_arr[,,i]
  dimn <- dimnames(df)
  dim(df) <- NULL
  df.res <- data.frame(Bias = df, N = rep(dimn[[2]], each = 6), Type = rep(dimn[[1]], 7))
  return(df.res)
})), Coef = rep(levels(Big_rare.copy$Coef), each = 42))

Bias_rare_df$Coef <- factor(Bias_rare_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Bias_rare_df$Type <- factor(Bias_rare_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                          "Systematic", "Topographic"))

Bias_rare_plot <- ggplot(Bias_rare_df, aes(x = N, y = Bias, group = Type, col = Type)) +
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
Var_rare_arr <- with(Big_comp.rare, tapply(Val, INDEX = list(Type, N, Coef), FUN = var))

Var_rare_df <- data.frame(do.call(rbind, lapply(1:(dim(Var_rare_arr)[3]), function(i) {
  df <- Var_rare_arr[,,i]
  dimn <- dimnames(df)
  dim(df) <- NULL
  df.res <- data.frame(Var = df, N = rep(dimn[[2]], each = 6), Type = rep(dimn[[1]], 7))
  return(df.res)
})), Coef = rep(levels(Big_rare.copy$Coef), each = 42))

Var_rare_df$Coef <- factor(Var_rare_df$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Var_rare_df$Type <- factor(Var_rare_df$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                        "Systematic", "Topographic"))

Var_rare_df$RMSE <- with(Var_rare_df, sqrt(Var))

Var_rare_plot <- ggplot(Var_rare_df, aes(x = N, y = RMSE, group = Type, col = Type)) +
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

#Compute true probabilities for D. tundrae
True_pred_tundrae.2 <- data.frame(var_val = c(Bio1_seq.AOI, Bio12_seq.AOI),
                                  var = rep(c("Bio1", "Bio12"),
                                            each = 100))

True_pred_tundrae.2$Proba <- unlist(lapply(c("Bio1", "Bio12"), function(nm) {
  if(nm %in% "Bio1") {
    mat.prd <- cbind(1, True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T],
                     (True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T]^2),
                     MeanPrec.AOI)
    #colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef.T.nm))
    return(pred_val)
  } else {
    mat.prd <- cbind(1, MeanTemp.AOI,
                     MeanTemp.AOI^2,
                     True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio12", "var_val", drop = T])
    #colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef.T.nm))
    return(pred_val)
  }
}))

#RANDOM

set.seed(1896)
Rnd_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Random_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI, 
                                                  x_crds =  Climate_crds,
                                                  Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                  n = N, min_p = 30), simplify = F)
  return(res)
})

#check sample size
sapply(seq_along(Rnd_rare_map), function(i) sum(sapply(Rnd_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Rnd_rare_map) <- names(Random_map)

#Prediction maps for the different sampling efforts
Rnd_rare_map.lyr <- lapply(Rnd_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Rnd_rare_map.fit <- do.call(rbind, lapply(names(Rnd_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Rnd_rare_map[[nm]]), function(.) {
    inner_df <- Rnd_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Rnd_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#STRATIFIED

set.seed(1922)
Str_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Strat_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                 Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                 n = N,
                                                 strata = Bio1_12.rcl.qrt.AOI, min_p = 30), simplify = F)
  return(res)
})

#check sample size
sapply(seq_along(Str_rare_map), function(i) sum(sapply(Str_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Str_rare_map) <- names(Rnd_rare_map)

Str_rare_map.lyr <- lapply(Str_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Str_rare_map.fit <- do.call(rbind, lapply(names(Str_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Str_rare_map[[nm]]), function(.) {
    inner_df <- Str_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Str_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#ROAD PROXIMITY

set.seed(1942)
Prox_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Proximity_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                     Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                     prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

#check sample size
sapply(seq_along(Prox_rare_map), function(i) sum(sapply(Prox_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Prox_rare_map) <- names(Rnd_rare_map)

Prox_rare_map.lyr <- lapply(Prox_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Prox_rare_map.fit <- do.call(rbind, lapply(names(Prox_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Prox_rare_map[[nm]]), function(.) {
    inner_df <- Prox_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Prox_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#UNIFORM

set.seed(1936)
Unif_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Uniform_mapping(y = D.tundrae.bin.AOI,
                                                   x = Chelsa.AOI,
                                                   x_sdf = Chelsa.AOI.df.sp,
                                                   Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                   n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

#check sample size
sapply(seq_along(Unif_rare_map), function(i) sum(sapply(Unif_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Unif_rare_map) <- names(Rnd_rare_map)

Unif_rare_map.lyr <- lapply(Unif_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Unif_rare_map.fit <- do.call(rbind, lapply(names(Unif_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Unif_rare_map[[nm]]), function(.) {
    inner_df <- Unif_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Unif_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#SYSTEMATIC

set.seed(1964)
Syst_rare_map <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 100, expr = Systematic_mapping(y = D.tundrae.bin.AOI, x = Chelsa.stack.syst,
                                                      N = n, perc_inc = 0.07, poly_proj = Elev_AOI.proj,
                                                      Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                      min_p = 30), simplify = F)
  return(res)
})

#check sample size
sapply(seq_along(Syst_rare_map), function(i) sum(sapply(Syst_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Syst_rare_map) <- names(Rnd_rare_map)

Syst_rare_map.lyr <- lapply(Syst_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Syst_rare_map.fit <- do.call(rbind, lapply(names(Syst_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Syst_rare_map[[nm]]), function(.) {
    inner_df <- Syst_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Syst_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#TOPOGRAPHIC

set.seed(1952)
Topo_rare_map <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 100, expr = Topo_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                N = n., perc_inc = 0.07, topo_layer = Topogr_het_Abr.AOI,
                                                Con_val = c(Bio1 = MeanTemp.AOI, Bio12 = MeanPrec.AOI),
                                                min_p = 30),
                   simplify = F)
})

#check sample size
sapply(seq_along(Topo_rare_map), function(i) sum(sapply(Topo_rare_map[[i]], '[[', 3) != Sampl_effort[i]))

names(Topo_rare_map) <- names(Rnd_rare_map)

Topo_rare_map.lyr <- lapply(Topo_rare_map, function(i) mean(stack(lapply(i, '[[', 1))))

Topo_rare_map.fit <- do.call(rbind, lapply(names(Topo_rare_map), function(nm) {
  df <- do.call(rbind, lapply(seq_along(Topo_rare_map[[nm]]), function(.) {
    inner_df <- Topo_rare_map[[nm]][[c(., 2)]]
    inner_df$repl <- as.character(.)
    return(inner_df)
  }))
  df$N <- nm
  return(df)
}))

ggplot(Topo_rare_map.fit, aes(x = var_val, y = fit)) +
  geom_line(aes(group = repl), alpha = .3) +
  geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#Figure showing predicted vs true response curves
RespCurves.df.rare <- data.frame(rbind(Rnd_rare_map.fit, Str_rare_map.fit, Prox_rare_map.fit,
                                       Syst_rare_map.fit, Topo_rare_map.fit, Unif_rare_map.fit),
                                 Typ = rep(c("Rnd", "Str", "Prx", "Syst", "Topo", "Unif"),
                                           times = sapply(list(Rnd_rare_map.fit, Str_rare_map.fit, Prox_rare_map.fit,
                                                               Syst_rare_map.fit, Topo_rare_map.fit, Unif_rare_map.fit), nrow)))

RespCurves.df.rare$Typ <- factor(RespCurves.df.rare$Typ,
                                 levels = c("Rnd", "Str", "Prx", "Unif", "Syst", "Topo"))

Sampling_label <- c("Rnd" = "Random", "Str" = "Stratified", "Prx" = "Proximity",
                    "Unif" = "Uniform", "Syst" = "Systematic", "Topo" = "Topographic")

ggarrange(plotlist = lapply(levels(RespCurves.df.rare$Typ), function(typ) {
  plot_rc <- ggplot(RespCurves.df.rare[RespCurves.df.rare$Typ == typ, ], aes(x = var_val, y = fit)) +
    geom_line(aes(group = repl), alpha = .3) +
    geom_line(data = True_pred_tundrae.2, aes(x = var_val, y = Proba), col = "red") +
    facet_grid(N ~ var, scales = "free_x",
               labeller = as_labeller(c("N_200" = "200", "N_250" = "250", "N_300" = "300",
                                        "N_350" = "350", "N_400" = "400", "N_450" = "450", "N_500" = "500",
                                        "Bio1" = "Temp", "Bio12" = "Precip"))) +
    xlab(Sampling_label[typ]) + ylab(NULL) +
    theme_classic() +
    theme(axis.title = element_text(size = 12), strip.text = element_text(size = 10),
          axis.text.x.bottom = element_text(angle = 35, vjust = .8))
  if(typ %in% c("Rnd", "Unif")) plot_rc <- plot_rc + ylab("Occurrence probability")
  return(plot_rc)
}))