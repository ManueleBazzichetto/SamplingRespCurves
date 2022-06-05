#Rare species
lapply(c("raster", "sf"), require, character.only = T)

#try with the virtualspecies package
#RareSpParam <- formatFunctions(Bio1 = c(fun = 'quadraticFun', a = -1, b = 2, c = 0),
#                                 Bio12 = c(fun = 'quadraticFun', a = -1, b = 2, c = 0))

#RareSpTrueP <- virtualspecies::generateSpFromFun(raster.stack = Chelsa.stack, parameters = RareSpParam,
#                                  species.type = "additive", plot = T)

#so, maybe what it does is to generate the quadratic (from, eg, dnorm) and rescale between 0-1 using: (x - min)/(max - min)
#functions (e.g. quadraticFun) are applied to data, and args are passed as lists - see body of generateSpFromFun)

Bio1_seq <- seq(cellStats(Chelsa.stack$Bio1, "min"), cellStats(Chelsa.stack$Bio1, "max"), .1)
Bio12_seq <- seq(cellStats(Chelsa.stack$Bio12, "min"), cellStats(Chelsa.stack$Bio12, "max"), .1)

##CREATE DIANTHUS TUNDRAE--------------------------------------------------------------------------------------
#THREE-POINTS APPROACH (SET THREE POINTS WITH COORDINATES C(Y, X) AND FIT A QUADRATIC FUNCTION WITH LEAST SQUARES)
#THEN, GET COEFFICIENTS, AND USE THEM TO SIMULATE SPECIES RESPONSE CURVE

#remember -> low p -> low prevalence
#knowing three points, we could use LS estimation to fit a parabola
#Temp -> unimodal function (quadratic)
#k (vertex over y) is 0.5; h (vertex over x) is 7.5
#pt1 (intercept) is -4; 6
#pt2 (right end tail) is -4; 9

Coords_quadr_LS <- data.frame(Y = c(0.5, -4, -4), X = c(7.5, 6, 9))
with(Coords_quadr_LS, plot(Y ~ X))
#fit quadratic
lm(Y ~ X + I(X^2), data = Coords_quadr_LS)
#saturated model -> obs == fitted
#so the fit is:
ggplot(Coords_quadr_LS, aes(x = X, y = Y)) +
  geom_smooth(method = "lm", formula = (y ~ x + I(x^2))) +
  geom_point()

#coords vertex quadratic
Coords_quadratic <- function(a, b, c) {
  h <- -1*(b/(2*c))
  k <- (a - (b^2/(4*c)))
  return(c(h = h, k = k))
}

coefficients(lm(Y ~ X + I(X^2), data = Coords_quadr_LS))

Coords_quadratic(a = -112, b = 30, c = -2)

#a gives value of the function when vary is 0
par(mfrow = c(1, 2))
(function(a, b, c, vary) {
  plot((a + b*vary + c*(vary^2)) ~ vary, ylab = "logit", xlab = "x")
  plot(plogis(a + b*vary + c*(vary^2)) ~ vary, ylab = "p", xlab = "x")
})(a = -112, b = 30, c = -2, vary = Bio1_seq)

#Prec -> linear relationship
plot(plogis(-4.6 + 0.0044*Bio12_seq) ~ Bio12_seq)

#params for bio1: a = -112, b = 30, c = -2
#params for bio12: a = -4.6, b = 0.0044

D.tundrae.link <- (-112 - 4.6) + 30*Chelsa.stack$Bio1 - 2*(Chelsa.stack$Bio1^2) + 0.0044*Chelsa.stack$Bio12
D.tundrae.prob <- calc(D.tundrae.link, plogis)

plot(D.tundrae.prob)

#warnings are for NAs
D.tundrae.bin <- calc(D.tundrae.prob, fun = function(x) {rbinom(1, 1, x)})

#Dfs for ggplotting
D.tundrae.prob.df <- na.omit(as.data.frame(D.tundrae.prob, xy = T))
colnames(D.tundrae.prob.df)[3] <- "Prob"

D.tundrae.bin.df <- na.omit(as.data.frame(D.tundrae.bin, xy = T))
colnames(D.tundrae.bin.df)[3] <- "PA"

ggarrange(ggplot() +
            geom_raster(data = D.tundrae.prob.df, aes(x = x, y = y, fill = Prob)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_viridis_c() +
            labs(title = "Dianthus tundrae", x = "Longitude", y = "Latitude", fill = "Probability"),
          ggplot() +
            geom_raster(data = D.tundrae.bin.df, aes(x = x, y = y, fill = PA)) +
            geom_sf(data = Elev_AOI, color = "white", fill = NA) +
            coord_sf(xlim = c(13.01193, 14.76349), ylim = c(41.68703, 42.90137)) +
            scale_fill_viridis_c() +
            labs(title = "Dianthus tundrae", x = "Longitude", y = "Latitude", fill = "PresAbs"), nrow = 1)

#comparison between D. tundrae and sperandii
par(mfrow = c(1, 2))
plot(D.sperandii.bin)
plot(D.tundrae.bin)

#mask D.tundrae.bin for AOI
D.tundrae.bin.AOI <- mask(D.tundrae.bin, mask = Elev_AOI)

raster::freq(D.tundrae.bin.AOI)

#True coef
True_coef_rare <- setNames(c(-116.6, 30, -2, 0.0044), nm = unique(Rnd_res_mats$Coef))

##SIMULATE THE DIFFERENT SAMPLING APPROACHES-------------------------------------------------------------------

#Sampling! We need to assure at least 30 presences in each training set

#Random

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
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

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
Rnd_res_mats$True_val <- unname(True_coef_rare[as.character(Rnd_res_mats$Coef)])
Rnd_res_mats$N <- factor(Rnd_res_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Rnd_res_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Strat

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
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

#note that Strat_raster repeat at each iteration some computations can be done just once....
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
Str_rare_mats$True_val <- unname(True_coef_rare[as.character(Str_rare_mats$Coef)])
Str_rare_mats$N <- factor(Str_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Str_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Proximity

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
    npres <- !(sum(Fake_df$PA) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

#check using Abr_highway_AOI instead of Abr_highway_AOI.rare
Proximity_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
               prox_lay = Abr_highway_AOI, n = 200, min_p = 30)

replicate(n = 10, expr = Proximity_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                         prox_lay = Abr_highway_AOI, n = 200, min_p = 30),
          simplify = F)

#go parallel (if needed)
#cr5 <- parallel::makeCluster(5)

#parallel::clusterExport(cr5, c("D.tundrae.bin.AOI", "Chelsa.AOI",
#                               "Abr_highway_AOI.rare", "Proximity_rare", "Sampl_effort"))

#parallel::clusterEvalQ(cr5, list(library(brglm2), library(detectseparation), library(raster)))

#note that Proximity_fun repeats at each iteration a computation that can be done just once (see line 209 in Proximity_rare....
Prx_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Proximity_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                  prox_lay = Abr_highway_AOI, n = N, min_p = 30),
                   simplify = F)
  return(res)
})

#stopCluster(cr5)

#check separation - fine
lapply(Prx_rare_res, function(.) sum(sapply(., function(i) i[[4]])))

#check correlation
lapply(Prx_rare_res, function(.) mean(sapply(., function(i) i[[2]])))

#check n
lapply(Prx_rare_res, function(.) colMeans(do.call(rbind, lapply(., function(i) {
  Tab <- as.data.frame(i[[3]])
  c(N_pr = Tab[Tab$Var1 == "1", "Freq"], N_ab = Tab[Tab$Var1 == "0", "Freq"])
  }))))

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
Prx_rare_mats$True_val <- unname(True_coef_rare[as.character(Prx_rare_mats$Coef)])
Prx_rare_mats$N <- factor(Prx_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Prx_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Unif

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
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

#Try from 250
Uniform_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
             x_sdf = Chelsa.AOI.df.sp,
             n = 250, rsl = 10, min_p = 30)

#go parallel (if needed)
cr7 <- parallel::makeCluster(7)

parallel::clusterExport(cr7, c("D.tundrae.bin.AOI", "Chelsa.AOI", "Chelsa.AOI.df.sp",
                               "Unif_sampl", "Uniform_rare", "uesampling2.0",
                               "Sampl_effort"))

parallel::clusterEvalQ(cr7, list(library(brglm2), library(detectseparation), library(raster), library(sf)))


Unif_rare_res <- parLapply(cl = cr7, X = Sampl_effort[-c(1)], function(N) {
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

names(Unif_rare_res) <- paste0("N_", Sampl_effort[-c(1)])

Unif_rare_mat <- do.call(rbind, lapply(names(Unif_rare_res), function(nm) {
  Mat <- do.call(rbind, lapply(Unif_rare_res[[nm]], '[[', 1))
  Mat_nm <- colnames(Mat)
  dim(Mat) <- NULL
  Mat <- data.frame(Val = Mat, Coef = rep(Mat_nm, each = 500))
  Mat$N <- nm
  return(Mat)
}))

Unif_rare_mat$Coef <- factor(Unif_rare_mat$Coef, levels = c("(Intercept)", "Pr", "Temp", "I(Temp^2)"))
Unif_rare_mat$True_val <- unname(True_coef_rare[as.character(Unif_rare_mat$Coef)])
Unif_rare_mat$N <- factor(Unif_rare_mat$N, levels = paste0("N_", Sampl_effort[-c(1)]))

ggplot(Unif_rare_mat, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Systematic

Systematic_rare <- function(y, x, N, poly_proj, min_p = 30) {
  require(brglm2)
  require(detectseparation)
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
    npres <- !(sum(Fake_df$PA) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

Syst_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Systematic_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                   poly_proj = Elev_AOI.proj, N = N, min_p = 30),
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
Syst_rare_mats$True_val <- unname(True_coef_rare[as.character(Syst_rare_mats$Coef)])
Syst_rare_mats$N <- factor(Syst_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Syst_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#Topographic

Topo_rare <- function(y, x, N, topo_layer, min_p = 30) {
  require(brglm2)
  require(detectseparation)
  Topo_lyr.df <- as.data.frame(topo_layer, xy = T)
  Topo_lyr.df <- Topo_lyr.df[!is.na(Topo_lyr.df$layer), ]
  npres <- T
  while(npres) {
    Topo_coords <- Topo_lyr.df[sample(nrow(Topo_lyr.df), size = N, replace = F), c("x", "y")]
    Fake_df <- na.omit(extract(stack(y, x), Topo_coords, df = T))
    Fake_df$ID <- NULL
    colnames(Fake_df) <- c("PA", "Temp", "Pr")
    npres <- !(sum(Fake_df$PA) >= min_p)
    }
  Coef <- coef(glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df))
  Sep_check <- glm(PA ~ Temp + I(Temp^2) + Pr, family = binomial, data = Fake_df, method = "detect_separation")
  Cor_v <- with(Fake_df, cor(Temp, Pr))
  return(list(Coef = Coef, Cor = Cor_v, N_pa = table(Fake_df$PA), Sep = Sep_check$separation))
}

Topo_rare_res <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 500, expr = Topo_rare(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                             topo_layer = Topogr_het_Abr.AOI, N = N, min_p = 30),
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
Topo_rare_mats$True_val <- unname(True_coef_rare[as.character(Topo_rare_mats$Coef)])
Topo_rare_mats$N <- factor(Topo_rare_mats$N, levels = paste0("N_", Sampl_effort))

ggplot(Topo_rare_mats, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot() +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

##GET AND COMPARE RESULTS--------------------------------------------------------------------------------------

##compare results
Big_comp.rare <- rbind(data.frame(Rnd_res_mats, Type = "Random"),
                       data.frame(Str_rare_mats, Type = "Strat"),
                      data.frame(Prx_rare_mats, Type = "Prox"),
                      data.frame(Syst_rare_mats, Type = "Systematic"),
                      data.frame(Unif_rare_mat, Type = "Uniform"),
                      data.frame(Topo_rare_mats, Type = "Topographic"))

Big_comp.rare$Type <- factor(Big_comp.rare$Type, levels = c("Random", "Prox", "Strat", "Uniform",
                                                  "Systematic", "Topographic"))

#Uniform seems to be quite biased..
ggplot(Big_comp.rare, aes(y = Val, x = N)) +
  geom_hline(yintercept = 0, col = "black", alpha = .4) +
  geom_boxplot(aes(fill = Type)) +
  geom_hline(aes(yintercept = True_val, col = Coef)) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ Coef, scales = "free") +
  theme_classic() +
  theme(legend.position = "top") 

#MSE, variance and bias

Mse_rare <- lapply(as.character(unique(Big_comp.rare$Type))[-5], function(ty) {
  mapply(function(x, y) Mse_samp(df = Big_comp.rare, Coef. = x, N. = y, Type. = ty, True_vals = True_coef_rare),
         x = rep(names(True_coef_nm), 7), y = rep(paste("N", Sampl_effort, sep = "_"), each = 4))
})

names(Mse_rare) <- as.character(unique(Big_comp.rare$Type))[-5]

Mse_rare$Uniform <- mapply(function(x, y) Mse_samp(df = Big_comp.rare,
                                                   Coef. = x, N. = y, Type. = "Uniform",
                                                   True_vals = True_coef_rare),
                           x = rep(names(True_coef_nm), 6),
                           y = rep(paste("N", Sampl_effort[-1], sep = "_"), each = 4))

Mse_rare_df <- do.call(rbind, lapply(names(Mse_rare), function(nm) {
  nms <- names(Mse_rare[[nm]])
  if(nm %in% c("Uniform")) {
    df <- data.frame(Mse_val = unname(Mse_rare[[nm]]), Coef = nms,
                     N = rep(paste("N", Sampl_effort[-1], sep = "_"), each = 4), 
                     Type = nm)
  } else {
    df <- data.frame(Mse_val = unname(Mse_rare[[nm]]), Coef = nms,
                   N = rep(paste("N", Sampl_effort, sep = "_"), each = 4), 
                   Type = nm)
    }
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
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("Mean squared error") + xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) 

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/MSE_tundrae.jpeg", device = "jpeg", units = "cm",
       width = 18, height = 15, dpi = 300)

#Bias

Big_rare.copy <- Big_comp.rare

Big_rare.copy$Bias <- with(Big_rare.copy, Val - True_val)

Bias_rare_arr <- with(Big_rare.copy, tapply(Bias, INDEX = list(Type, N, Coef), mean)) #NA for N_200 Uniform

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
  facet_wrap(~ Coef, scales = "free_y",
             labeller = as_labeller(c("(Intercept)" = "Intercept",
                                      "Pr" = "Precipitation",
                                      "Temp" = "Temperature",
                                      "I(Temp^2)" = "Quadr. Temp."))) +
  ylab("SQRT(Variance)") + xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "top", strip.background = element_blank(),
        text = element_text(size = 14), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

ggarrange(Bias_rare_plot, Var_rare_plot, ncol = 1, nrow = 2, common.legend = T)

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/B_V_tundrae.jpeg", device = "jpeg", units = "cm",
       width = 18, height = 25, dpi = 300)

#MAPS OF PREDICTIONS AND PLOT OF CORRELATION ESTIMTED VS TRUE PROBABILITIES OF SP. OCCURRENCE -------------------------

#df with true proba value - here I should used the model that generate the species distribution
#and not the compute the fitted values separately for temp and prec
#True_pred_tundrae <- data.frame(var_val = c(Plot_bio1[Plot_bio1$Species == "D. tundrae", "Temperature"],
#                                              Plot_bio12[Plot_bio12$Species == "D. tundrae", "Precipitation"]),
#                                  var = rep(c("Bio1", "Bio12"), times = c(nrow(Plot_bio1[Plot_bio1$Species == "D. tundrae", ]),
#                                                                          nrow(Plot_bio12[Plot_bio12$Species == "D. tundrae", ]))),
#                                  rbind(Plot_bio1[Plot_bio1$Species == "D. tundrae", c("Proba"), drop = F],
#                                        Plot_bio12[Plot_bio12$Species == "D. tundrae", c("Proba"), drop = F]))

#get rid of values not present in the AOI
#minValue(Chelsa.AOI$Bio1); maxValue(Chelsa.AOI$Bio1)
#minValue(Chelsa.AOI$Bio12); maxValue(Chelsa.AOI$Bio12)

#True_pred_tundrae <- do.call(rbind, lapply(unique(True_pred_tundrae$var), function(var_nm) {
#  sub_df <- True_pred_tundrae[True_pred_tundrae$var == var_nm, ]
#  if(var_nm %in% "Bio1") {
#    sub_df <- sub_df[sub_df$var_val <= 14.25 & sub_df$var_val >= -0.35, ]
#  } else {
#    sub_df <- sub_df[sub_df$var_val <= 1564 & sub_df$var_val >= 631.3, ]
#  }
#}))

#alternative formulation of predicted value -> if it works, delete what's above (deleted)
#regressor for quadratic term of temp is kept at its median value (as done by predictorEffects, see Results_interpr)

Bio1_seq.AOI <- seq(from = minValue(Chelsa.AOI$Bio1), to = maxValue(Chelsa.AOI$Bio1), length.out = 100)
Bio12_seq.AOI <- seq(from = minValue(Chelsa.AOI$Bio12), to = maxValue(Chelsa.AOI$Bio12), length.out = 100)

True_pred_tundrae.2 <- data.frame(var_val = c(Bio1_seq.AOI, Bio12_seq.AOI),
                                  var = rep(c("Bio1", "Bio12"),
                                            each = 100))

True_pred_tundrae.2$Proba <- unlist(lapply(c("Bio1", "Bio12"), function(nm) {
  if(nm %in% "Bio1") {
    mat.prd <- cbind(1, True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T],
                     (True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T]^2),
                    mean(True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio12", "var_val", drop = T]))
    colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef_rare))
    return(pred_val)
  } else {
    mat.prd <- cbind(1, mean(True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T]),
                     median(True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio1", "var_val", drop = T]^2),
                     True_pred_tundrae.2[True_pred_tundrae.2$var == "Bio12", "var_val", drop = T])
    colnames(mat.prd) <- NULL
    pred_val <- plogis(mat.prd%*%unname(True_coef_rare))
    return(pred_val)
  }
}))


#Random
Rnd_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Random_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI, 
                                                  x_crds =  Climate_crds, n = N, min_p = 30), simplify = F)
  return(res)
})

names(Rnd_rare_map) <- names(Random_map)

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

#Strat
Str_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Strat_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI, n = N,
                                                 strata = Bio1_12.rcl.qrt.AOI, min_p = 30), simplify = F)
  return(res)
})

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

#Prox
Prox_rare_map <- lapply(Sampl_effort, function(N) {
  res <- replicate(n = 100, expr = Proximity_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                     prox_lay = Abr_highway_AOI, n = N, min_p = 30), simplify = F)
  return(res)
})

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

#Unif
Unif_rare_map <- lapply(Sampl_effort[-1], function(N) {
  res <- replicate(n = 100, expr = Uniform_mapping(y = D.tundrae.bin.AOI,
                                                   x = Chelsa.AOI,
                                                   x_sdf = Chelsa.AOI.df.sp,
                                                   n = N, rsl = 10, min_p = 30), simplify = F)
  return(res)
})

names(Unif_rare_map) <- names(Rnd_rare_map)[-1]

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
  geom_line(data = True_pred_tundrae, aes(x = var_val, y = Proba), col = "red") +
  facet_grid(N ~ var, scales = "free_x") +
  ylab("Occurrence probability") + xlab(NULL) +
  theme_classic()

#Systematic
Syst_rare_map <- lapply(Sampl_effort, function(n) {
  res <- replicate(n = 100, expr = Systematic_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                      N = n, poly_proj = Elev_AOI.proj, min_p = 30), simplify = F)
  return(res)
})

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

#Topo
Topo_rare_map <- lapply(Sampl_effort, function(n.) {
  res <- replicate(n = 100, expr = Topo_mapping(y = D.tundrae.bin.AOI, x = Chelsa.AOI,
                                                N = n., topo_layer = Topogr_het_Abr.AOI, min_p = 30),
                   simplify = F)
})

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

#plots of response curves
RespCurves.df.rare <- data.frame(rbind(Rnd_rare_map.fit, Str_rare_map.fit, Prox_rare_map.fit,
                                  Syst_rare_map.fit, Topo_rare_map.fit),
                            Typ = rep(c("Rnd", "Str", "Prx", "Syst", "Topo"),
                                      each = nrow(Rnd_rare_map.fit)))

RespCurves.df.rare <- rbind(RespCurves.df.rare, data.frame(Unif_rare_map.fit, Typ = "Unif"))

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

ggsave(filename = "~/Documents/SpRespCurve/Fig_rmarkdown/RespCurvesTundrae.jpeg", device = "jpeg",
       width = 32, height = 25, units = "cm", dpi = 300)


#Maps of predictions and correlations ------------------------------------------

compareRaster(Rnd_rare_map.lyr[[1]], stack(Rnd_rare_map.lyr[-1]),
              stack(Str_rare_map.lyr), stack(Prox_rare_map.lyr), stack(Unif_rare_map.lyr),
              stack(Syst_rare_map.lyr), stack(Topo_rare_map.lyr), orig = T)


Coords_rare_pred_map <- as.data.frame(Rnd_rare_map.lyr[[1]], xy = T, na.rm = T)
Coords_rare_pred_map <- Coords_rare_pred_map[c("x", "y")]

List_rare_pred_map <- list(Rnd = Rnd_rare_map.lyr, Str = Str_rare_map.lyr, Prx = Prox_rare_map.lyr, 
                           Unif = Unif_rare_map.lyr, Syst = Syst_rare_map.lyr, Topo = Topo_rare_map.lyr)


Df_rare_pred_map <- do.call(rbind, lapply(names(List_rare_pred_map), function(nm) {
  Df_pred <- do.call(rbind, lapply(seq_along(List_rare_pred_map[[nm]]), function(i) {
    Pred <- extract(List_rare_pred_map[[nm]][[i]], Coords_rare_pred_map)
    Pred_data <- data.frame(Prob = Pred, N = paste0("N_", if(nm %in% "Unif") Sampl_effort[i+1] else Sampl_effort[i]))
  }))
  Df_pred$Type <- nm
  return(Df_pred)
  }))

#True Proba
Vec_true_rare_prob <- extract(D.tundrae.prob, Coords_rare_pred_map)

#recycling
Df_rare_pred_map$True_prob <- Vec_true_rare_prob

ggplot(Df_rare_pred_map, aes(x = True_prob, y = Prob, col = Type)) +
  geom_point(alpha = .1) +
  geom_abline(slope = 1, intercept = 0, col = "black", lty = "dotted") +
  facet_wrap(~ N) +
  ylab("Estimated probabilities") + xlab("True probabilities") +
  theme_minimal() +
  theme(legend.position = "top")

#try enhance visualisation using mean of p for 0.1 classes
#Cut observed proba in 10 groups by .1 intervals

Unique_rare_combo <- unique(Df_rare_pred_map[c("N", "Type")])

Cut_df_rare_pred <- do.call(rbind, Map(function(x, y) {
  Df_res <- Df_rare_pred_map[Df_rare_pred_map$Type == x & Df_rare_pred_map$N == y, ]
  Df_res$cut_p <- with(Df_res, cut(Prob, breaks = seq(0, 1, .1), labels = seq(1, 10)))
  #here aggregate and return
  Df_res <- aggregate(Df_res[c("Prob", "True_prob")], by = list(Df_res$cut_p), mean)
  colnames(Df_res)[1] <- "Group"
  Df_res$Type <- x
  Df_res$N <- y
  return(Df_res)
  }, x = Unique_rare_combo$Type, y = Unique_rare_combo$N))

row.names(Cut_df_rare_pred) <- seq_len(nrow(Cut_df_rare_pred))

#library(paletteer); library(ggthemes)
paletteer_d("LaCroixColoR::paired")
paletteer_d("beyonce::X66") #yellow of systematic from here

#range of p
lapply(List_rare_pred_map, function(i) sapply(i, cellStats, "max"))

Cor_rare_pred_plot <- ggplot(Cut_df_rare_pred[Cut_df_rare_pred$N %in% c("N_200", "N_250","N_300"), ], aes(x = True_prob, y = Prob, col = Type)) +
  geom_abline(slope = 1, intercept = 0, col = "black", lty = "dotted", alpha = .5) +
  geom_point(alpha = .8, cex = 3) +
  geom_line(alpha = .5) +
  scale_color_manual(values = c(Prx = "#C70E7BFF", Str = "#007BC3FF", Topo = "#009F3FFF",
                       Rnd = "#F4B95AFF", Syst = "#FCEA1BFF",
                       Unif = "#EF7C12FF"),
                     labels = c(Prx = "Proximity", Str = "Stratified", Topo = "Topographic",
                                Rnd = "Random", Syst = "Systematic", Unif = "Uniform")) +
  facet_wrap(~ N) +
  ylab("Estimated probabilities") + xlab("True probabilities") +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 14))

##save data for markdown -> last save: 31/05/2022
#spatial data
save(D.tundrae.prob, D.tundrae.bin, 
     file = "SpatialDataRare.RData")
#dfs
save(D.tundrae.bin.df, D.tundrae.prob.df,
     file = "SpatialDFRare.RData")

#resplots
save(Mse_rare_plot, Var_rare_plot, Bias_rare_plot, Cor_rare_pred_plot,
     file = "ResPlotsRare.RData")
