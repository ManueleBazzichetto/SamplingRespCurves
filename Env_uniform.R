#UEsampling!
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

Chelsa.AOI.df <- na.omit(as.data.frame(Chelsa.AOI, xy = T))
#add scaled Bio1 and Bio12 for environmental space
Chelsa.AOI.df$Bio1.sc <- with(Chelsa.AOI.df, scale(Bio1))
Chelsa.AOI.df$Bio12.sc <- with(Chelsa.AOI.df, scale(Bio12))
Chelsa.AOI.df.sp <- st_as_sf(Chelsa.AOI.df, coords = c("Bio1.sc", "Bio12.sc"))

plot(st_geometry(Chelsa.AOI.df.sp))

#how many points for cell (considering a 10X10 grid)
Grid_uniform <- st_make_grid(x = Chelsa.AOI.df.sp, n = 10)

sapply(st_geometry(Grid_uniform), function(i) {
  nrow(Chelsa.AOI.df.sp[i, ])
})

plot(st_geometry(Chelsa.AOI.df.sp))
plot(Grid_uniform, add = T)

#see how D. sperandii occupies the env. space
Coords.D.sperandii <- as.data.frame(D.sperandii.bin.AOI, xy = T, na.rm = T)
Coords.D.sperandii.p <- Coords.D.sperandii[Coords.D.sperandii$layer == 1, ]
Chelsa.D.sperandii <- raster::extract(x = Chelsa.AOI, y = Coords.D.sperandii.p[c("x", "y")], df = T)

#standardize Bio1 and Bio12 using scaling factor used for Chelsa.AOI.df
#Bio1: center 9.170865; scale 2.555482
#Bio12: center 1049.6; scale 169.4228

Chelsa.D.sperandii$Bio1 <- with(Chelsa.D.sperandii, (Bio1 - 9.170865)/2.555482)
Chelsa.D.sperandii$Bio12 <- with(Chelsa.D.sperandii, (Bio12 - 1049.6)/169.4228)

Chelsa.D.sperandii.sp <- st_as_sf(Chelsa.D.sperandii, coords = c("Bio1", "Bio12"))

plot(st_geometry(Chelsa.AOI.df.sp))
plot(Grid_uniform, add = T)
plot(Chelsa.D.sperandii.sp, add = T, col = "green")

#see how the rare species occupies the env space
Coords.D.tundrae <- as.data.frame(D.tundrae.bin.AOI, xy = T, na.rm = T)
Coords.D.tundrae.p <- Coords.D.tundrae[Coords.D.tundrae$layer == 1, ]
Chelsa.D.tundrae <- raster::extract(x = Chelsa.AOI, y = Coords.D.tundrae.p[c("x", "y")], df = T)

#standardize Bio1 and Bio12 using scaling factor used for Chelsa.AOI.df
#Bio1: center 9.170865; scale 2.555482
#Bio12: center 1049.6; scale 169.4228

Chelsa.D.tundrae$Bio1 <- with(Chelsa.D.tundrae, (Bio1 - 9.170865)/2.555482)
Chelsa.D.tundrae$Bio12 <- with(Chelsa.D.tundrae, (Bio12 - 1049.6)/169.4228)

Chelsa.D.tundrae.sp <- st_as_sf(Chelsa.D.tundrae, coords = c("Bio1", "Bio12"))

plot(st_geometry(Chelsa.AOI.df.sp))
plot(Grid_uniform, add = T)
plot(Chelsa.D.tundrae.sp, add = T, col = "green")

#comparison
par(mfrow = c(1, 2))
plot(st_geometry(Chelsa.AOI.df.sp), main = "D. sperandii")
plot(Grid_uniform, add = T)
plot(Chelsa.D.sperandii.sp, add = T, col = "#DCE319FF")
plot(st_geometry(Chelsa.AOI.df.sp), main = "D. tundrae")
plot(Grid_uniform, add = T)
plot(Chelsa.D.tundrae.sp, add = T, col = "#39568CFF")

#for c(200, 300, 400, 500)
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 3) #203
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 4) #203
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 5) #332
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 6) #396
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 7) #458
uesampling2.0(Chelsa.AOI.df.sp, grid.res = 10, n.tr = 8) #519


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

#cool
sapply(Sampl_effort, function(n) nrow(Unif_sampl(x = Chelsa.AOI.df.sp, N = n, rsl = 10)))

