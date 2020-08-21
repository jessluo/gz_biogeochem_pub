# FUNCTIONS


library("plyr")
library("reshape2")
library("oce")
library("ggplot2")
library("stringr")
library("RColorBrewer")

`%ni%` <- Negate(`%in%`)

## {------------------- GEO - FUNCTIONS -----------------------------
as.radians <- function(theta = 0)
{
  pi <- 3.14159265358979323846
  return(theta * pi/180)
}

area_earth_1deg <- function(lat, lon){
  # calculates the surface area of the earth in a 1 degree grid from the center of the grid
  # inputs in degrees
  require(oce)
  pi <- 3.14159265358979323846
  R <- 6371.0 # Volumetric mean radius of the earth, in km, from nssdc.gsfc.nasa.gov
  lat1 <- as.radians(lat-0.5)
  lat2 <- as.radians(lat+0.5)
  lon1 <- lon-0.5 #longitudes do not need to be in radians  
  lon2 <- lon+0.5
  A <- pi / 180 * R^2 * abs(sin(lat1)-sin(lat2)) * abs(lon1-lon2)
  return(A)
}

area_earth_4pts <- function(lat1, lon1, lat2, lon2){
  # calculates the surface area of the earth within 4 points
  # inputs in degrees
  require(oce)
  pi <- 3.14159265358979323846
  lat1 <- as.radians(lat1)
  lat2 <- as.radians(lat2)
  R <- 6371.0 # Volumetric mean radius of the earth, in km, from nssdc.gsfc.nasa.gov
  A <- 2 * pi * R^2 * abs(sin(lat1)-sin(lat2)) * (abs(lon1-lon2) / 360)
  return(A)
}

gm_mean = function(x, na.rm=TRUE, na.exclude=FALSE){
  if(na.exclude){
    x <- x[which(!is.na(x))]
    x <- x[x>0]
    res <- exp(mean(log(x), na.rm=na.rm))
    return(res)
  } else {
    x <- x[x>0]
    res <- exp(mean(log(x), na.rm=na.rm))
    return(res)
  }
}

gm_sd <- function(x, na.rm = TRUE)
{
  require("EnvStats")
  vals=x[which(x > 0)]
  res <- geoSD(vals, na.rm=na.rm)
}

weighted.geomean <- function(x, w, ...)
{
  w=w[which(x>0)]
  x=x[which(x>0)]
  return(prod(x^w, ...)^(1/sum(w)))
}

#The weighted geometric sample deviation is the square root of the weighted geometric sample variance. The geometric sample variance is
#var(logdata,weights) = exp(var(logdata)/b)
#where the effective base b is b(weights) = (sumi wi)2/sumi(wi2)
weighted.geosd <- function(x, w, ...){
  exp(var(log(x))/(sum(w)^2/sum(w^2)))
}
# }

## ------------------- PROCESSING FUNCTIONS ------------------------

add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]), 
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)), 
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}


remove.jelly.NAs <- function(d){
  suppressWarnings(temp <- which(is.na(d$Chordata) & is.na(d$Cnidaria) & is.na(d$Ctenophora)))
  if (length(temp) != 0){
    d <- d[(!is.na(d$Chordata) | !is.na(d$Cnidaria) | !is.na(d$Ctenophora)),]
  }
  return(d)
}


which_closest_value <- function(x, y){
  idx <- which(abs(y-x)==min(abs(y-x)))
  return(idx)
}

# collapse values into a 1 degree grid - with the point at the center between the two grids
collapse.into.1deg <- function(d, na.rm=FALSE){
  d$lon <- floor(d$lon) + 0.5
  d$lat <- floor(d$lat) + 0.5
  d <- ddply(d, ~lat+lon, colMeans, na.rm=na.rm, .progress ="text")
  return(d)
}

process.nc.biomass.data <- function(d, times=timebnds_months){
  returndf <- NULL
  nt = length(times)
  
  if (nt==1){
    temp <- as.data.frame(cbind(lon, d))
    names(temp) <- c("lon", lat) # lat values from 
    temp <- melt(temp, id.vars = "lon", variable.name="lat", value.name=times)
    latlon <- temp[,1:2]
    returndf <- as.data.frame(temp[,3])
  }
  if (nt > 1){
    for (i in 1:length(times)){
      temp <- d[,,i]
      temp <- as.data.frame(cbind(lon, temp))
      names(temp) <- c("lon", lat) # lat values from 
      temp <- melt(temp, id.vars = "lon", variable.name="lat", value.name=times[i])
      latlon <- temp[,1:2]
      returndf <- as.data.frame(cbind(returndf, temp[,3]))
    }
  }

  names(returndf) <- times
  returndf <- cbind(latlon, returndf)
  returndf$lat <- as.numeric(as.character(returndf$lat))
  # correct the longitude to be -180 to 180
  returndf$lon <- ifelse(returndf$lon < -180, returndf$lon+360, returndf$lon)
  return(returndf)
}

summarize.into.seasons <- function(d){
  # notes that these dates start on the 18th of the previous month
  winter_months <- c("Jan", "Feb", "Mar")
  spring_months <- c("Apr", "May", "Jun")
  summer_months <- c("Jul", "Aug", "Sep")
  fall_months <- c("Oct", "Nov", "Dec")
  year <- c(winter_months, spring_months, summer_months, fall_months)
  
  if(any(winter_months %in% names(d))){
    win <- rowMeans(d[,which(names(d) %in% winter_months)], na.rm=T)
    # row sums of NAs is 0, so use rowMeans to find the NAs and put it into the array
    wintemp <- rowMeans(d[,which(names(d) %in% winter_months)], na.rm=T)
    win[which(is.na(wintemp))] <- NA
  } else {print("Error - no winter months in dataset")}
  
  if(any(spring_months %in% names(d))){
    spr <- rowMeans(d[,which(names(d) %in% spring_months)], na.rm=T)
    sprtemp <- rowMeans(d[,which(names(d) %in% spring_months)], na.rm=T)
    spr[which(is.na(sprtemp))] <- NA        
  } else {print("Error - no spring months in dataset")}
  
  if(any(summer_months %in% names(d))){
    sum <- rowMeans(d[,which(names(d) %in% summer_months)], na.rm=T)
    sumtemp <- rowMeans(d[,which(names(d) %in% summer_months)], na.rm=T)
    sum[which(is.na(sumtemp))] <- NA
  } else {print("Error - no summer months in dataset")}
  
  if(any(fall_months %in% names(d))){
    fall <- rowMeans(d[,which(names(d) %in% fall_months)], na.rm=T)
    falltemp <- rowMeans(d[,which(names(d) %in% fall_months)], na.rm=T)
    fall[which(is.na(falltemp))] <- NA
  } else {print("Error - no fall months in dataset")}  
  
  if(any(year %in% names(d))){
    ann <- rowMeans(d[,which(names(d) %in% year)], na.rm=T)
    anntemp <- rowMeans(d[,which(names(d) %in% year)], na.rm=T)
    ann[which(is.na(anntemp))] <- NA
  } else {print("Error - no months in dataset")}  
  
  output <- data.frame(d, win, spr, sum, fall, ann)
  
  # convert to average C biomass over 200 m
  output[3:ncol(output)] <- output[3:ncol(output)]/200
  
  return(output)
}

summarize.into.seasons.sum <- function(d){
  # notes that these dates start on the 18th of the previous month
  winter_months <- c("Jan", "Feb", "Mar")
  spring_months <- c("Apr", "May", "Jun")
  summer_months <- c("Jul", "Aug", "Sep")
  fall_months <- c("Oct", "Nov", "Dec")
  year <- c(winter_months, spring_months, summer_months, fall_months)
  
  if(any(winter_months %in% names(d))){
    win <- rowSums(d[,which(names(d) %in% winter_months)], na.rm=T)
    # row sums of NAs is 0, so use rowSums to find the NAs and put it into the array
    wintemp <- rowSums(d[,which(names(d) %in% winter_months)], na.rm=T)
    win[which(is.na(wintemp))] <- NA
  } else {print("Error - no winter months in dataset")}
  
  if(any(spring_months %in% names(d))){
    spr <- rowSums(d[,which(names(d) %in% spring_months)], na.rm=T)
    sprtemp <- rowSums(d[,which(names(d) %in% spring_months)], na.rm=T)
    spr[which(is.na(sprtemp))] <- NA        
  } else {print("Error - no spring months in dataset")}
  
  if(any(summer_months %in% names(d))){
    sum <- rowSums(d[,which(names(d) %in% summer_months)], na.rm=T)
    sumtemp <- rowSums(d[,which(names(d) %in% summer_months)], na.rm=T)
    sum[which(is.na(sumtemp))] <- NA
  } else {print("Error - no summer months in dataset")}
  
  if(any(fall_months %in% names(d))){
    fall <- rowSums(d[,which(names(d) %in% fall_months)], na.rm=T)
    falltemp <- rowSums(d[,which(names(d) %in% fall_months)], na.rm=T)
    fall[which(is.na(falltemp))] <- NA
  } else {print("Error - no fall months in dataset")}  
  
  if(any(year %in% names(d))){
    ann <- rowSums(d[,which(names(d) %in% year)], na.rm=T)
    anntemp <- rowSums(d[,which(names(d) %in% year)], na.rm=T)
    ann[which(is.na(anntemp))] <- NA
  } else {print("Error - no months in dataset")}  
  
  output <- data.frame(d, win, spr, sum, fall, ann)
  
  # convert to average C biomass over 200 m
  output[3:ncol(output)] <- output[3:ncol(output)]/200
  
  return(output)
}

plotting <- function(var, data, plot.title=var, plot.limits=NA, na.value="grey80", size=3, trans=NA){
  # import coastline
  # data(coastlineWorld) # oce package
  # coastline.world=data.frame(lon=coastlineWorld[["longitude"]],lat=coastlineWorld[["latitude"]]) 
  # coastline.world <- read.csv("/glade/u/home/jluo/p/tools/gshhg_world_c.csv")
  coastline.world <- read.csv("data/gshhg_world_c.csv")
  if(is.na(plot.limits)){plot.limits=NULL}
  
  require(ggplot2)
  require(RColorBrewer)

  # print(head(data)) # check
  
  if(is.na(trans)){
    p <- ggplot(mapping=aes(x=lon, y=lat)) + 
      geom_point(aes_string(colour=var), size=size, data=data, shape=15) + 
      geom_path(data=coastline.world) + 
      theme_bw() + 
      scale_color_gradientn(colours=brewer.pal(11,'Spectral'), na.value=na.value, limits=plot.limits) + 
      labs(title=plot.title) 
  }
  
  if(!is.na(trans)){
    p <- ggplot(mapping=aes(x=lon, y=lat)) + 
      geom_point(aes_string(colour=var), size=size, data=data, shape=15) + 
      geom_path(data=coastline.world) + 
      theme_bw() + 
      scale_color_gradientn(colours=brewer.pal(11,'Spectral'), na.value=na.value, limits=plot.limits, trans="log10") + 
      labs(title=plot.title) 
  }
  
  return(p)
}


# grid into 5 deg squares, calculate for the variables given
regrid_geomean <- function(d, vars){
  
  dnames <- vars
  N <- length(vars)
  
  worldlat5deg <- seq(from = -90, to = 90, by = 5)
  worldlon5deg <- seq(from = -180, to = 180, by =5)
  
  # create a 3-D matrix with N slices
  datamatrix <- array(dim=c(length(worldlat5deg), length(worldlon5deg), N))
  
  
  matlat <- matrix(nrow = length(worldlat5deg), ncol = length(worldlon5deg))
  matlon <- matrix(nrow = length(worldlat5deg), ncol = length(worldlon5deg))
  
  vec5deg <- matrix(nrow=(length(worldlat5deg)*length(worldlon5deg)), ncol = N+2)
  colnames(vec5deg) <- c("lat","lon", dnames)
  
  counter <- 1
  
  # gridding loop
  for (latloop in 1:length(worldlat5deg)) {
    for (lonloop in 1:length(worldlon5deg)) {
      
      # set a temporary dataframe
      tempd <- d[which(d$lat>=worldlat5deg[latloop] & d$lat<=worldlat5deg[latloop+1] & 
                         d$lon>=worldlon5deg[lonloop] & d$lon<=worldlon5deg[lonloop+1]),]
      
      # set the lat and lon matrix
      matlat[latloop,lonloop] <- worldlat5deg[latloop]
      matlon[latloop,lonloop] <- worldlon5deg[lonloop]
      
      vec5deg[counter,1] <- (worldlat5deg[latloop])+2.5   # latitudes
      vec5deg[counter,2] <- (worldlon5deg[lonloop])+2.5   # longitudes
      
      for (i in 1:N){
        var <- dnames[i]
        temp <- tempd[,which(names(tempd)==var)]

        datamatrix[latloop, lonloop, i] <- exp(mean(log(temp), na.rm=TRUE))
        
        vec5deg[counter, i+2] <- exp(mean(log(temp), na.rm=TRUE))
        
        remove(temp)
      }
      
      counter <- counter + 1     
    }
  }
  
  dc <- as.data.frame(vec5deg)  
  return(dc)
}


# grid into 5 deg squares, calculate mean for variables given
regrid <- function(d, vars){
  
  dnames <- vars
  N <- length(dnames)
  
  worldlat5deg <- seq(from = -90, to = 90, by = 5)
  worldlon5deg <- seq(from = -180, to = 180, by =5)
  
  # create a 3-D matrix with N slices
  datamatrix <- array(dim=c(length(worldlat5deg), length(worldlon5deg), N))
  
  
  matlat <- matrix(nrow = length(worldlat5deg), ncol = length(worldlon5deg))
  matlon <- matrix(nrow = length(worldlat5deg), ncol = length(worldlon5deg))
  
  vec5deg <- matrix(nrow=(length(worldlat5deg)*length(worldlon5deg)), ncol = N+2)
  colnames(vec5deg) <- c("lat","lon", dnames)
  
  counter <- 1
  
  # gridding loop
  for (latloop in 1:length(worldlat5deg)) {
    for (lonloop in 1:length(worldlon5deg)) {
      
      # set a temporary dataframe
      tempd <- d[which(d$lat>=worldlat5deg[latloop] & d$lat<=worldlat5deg[latloop+1] & 
                         d$lon>=worldlon5deg[lonloop] & d$lon<=worldlon5deg[lonloop+1]),]
      
      # set the lat and lon matrix
      matlat[latloop,lonloop] <- worldlat5deg[latloop]
      matlon[latloop,lonloop] <- worldlon5deg[lonloop]
      
      vec5deg[counter,1] <- (worldlat5deg[latloop])+2.5   # latitudes
      vec5deg[counter,2] <- (worldlon5deg[lonloop])+2.5   # longitudes
      
      for (i in 1:N){
        var <- dnames[i]
        temp <- tempd[,which(names(tempd)==var)]
        datamatrix[latloop, lonloop, i] <- mean(temp, na.rm = TRUE)
        
        vec5deg[counter, i+2] <- mean(temp, na.rm = TRUE)
        
        remove(temp)
      }
      
      counter <- counter + 1     
    }
  }
  
  dc <- as.data.frame(vec5deg)  
  return(dc)
}


biomass_by_biome <- function(d){
  # load biomes
  load('data/biomes_chl_mixedlayer.Rdata')
  
  # import area/volume
  av <- read.csv("data/surfaceArea_volume_1-deg.csv")
  av <- join(av,biomes_df, by=c("lat","lon"))
  av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(sum(x$vol))})
  names(av_biomesum) <- c("biome", "vol")
  
  biomass = join(d, biomes_df, by=c("lat","lon"))
  
  biomass = ddply(biomass, ~rank_phylum+biome, function(x){
    return(data.frame(Biomass.mgCm3=mean(x$Biomass.mgCm3), n=nrow(x)))
  })
  
  biomass = join(biomass, av_biomesum, by="biome")
  biomass$Biomass.gC = biomass$Biomass.mgCm3 * biomass$vol / 1000
  
  dtot = ddply(biomass, ~rank_phylum, function(x){
    dtot = data.frame(biome="TOTAL", t(colSums(x[,c("Biomass.mgCm3","vol","Biomass.gC","n")])))
  })
  
  biomass = arrange(rbind(biomass, dtot), rank_phylum, biome)
  
  return(biomass)
}


# calculate_biome_mean_UO <- function(data, taxon, biomes=biomes_df, n_biome=N_BIOME, n_ens=N_ENS, saveFiles=TRUE){
#   
#   mc_biomeres_mat <- array(NA, dim=c(n_biome, 8, n_ens))
#   
#   #print('Extract and compute averages for each biome')
#   for (i in 1:length(data)){
#     temp <- data[[i]]
#     temp$taxon <- as.character(temp$taxon)
#     temp <- temp[which(temp$taxon==taxon),] # swap out cnid, cten, chor if you want those specific fluxes
#     temp <- temp[,which(names(temp) != "taxon")]
#     temp$biome = as.numeric(temp$biome)
#     # biome means
#     temp_b <- ddply(temp, ~biome, function(x){
#       return(apply(x[,which(names(temp) %ni% c("lat","lon","biome"))], MARGIN = 2, mean))
#     })
#     mc_biomeres_mat[,,i] <- as.matrix(temp_b)
#   }
#   
#   # dimensions of mc_biomeres_mat are: [biome, rates, ensemble]
#   
#   # Compute min, max, and mean for each biome --
#   mins <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
#   maxs <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
#   uo_mean <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
#   bcis_upper <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
#   bcis_lower <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
#   
#   for (i in 1:dim(mc_biomeres_mat)[1]){ # loop over every biome
#     temp <- mc_biomeres_mat[i,,]
#     temp <- t(temp) # transpose
#     
#     rates <- temp[,2:8]
#     
#     min <- apply(rates, 2, min, na.omit=TRUE)
#     min <- c(i, min)
#     mins[i,] <- min
#     
#     max <- apply(rates, 2, max, na.omit=TRUE)
#     max <- c(i, max)
#     maxs[i,] <- max
#     
#     avg <- apply(rates, 2, mean, na.rm=TRUE) # should we be taking geometric means here??
#     avg <- c(i, avg)
#     uo_mean[i,] <- avg
#     
#     bci <- apply(rates, 2, return_quant)
#     bci_lower <- c(i, bci[1,])
#     bci_upper <- c(i, bci[2,])
#     bcis_lower[i,] <- bci_lower
#     bcis_upper[i,] <- bci_upper
#     
#   }
#   
#   if(saveFiles){
#     save(bcis_lower, file=str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_QUANT_lower.Rdata"))
#     save(bcis_upper, file=str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_QUANT_upper.Rdata"))
#   }
#   
#   # Put results together --
#   res <- cbind(melt(as.data.frame(uo_mean), id.vars = "V1", value.name = "mean"),
#                melt(as.data.frame(bcis_lower), id.vars = "V1", value.name = "lower")[,"lower"],
#                melt(as.data.frame(bcis_upper), id.vars = "V1", value.name = "upper")[,"upper"])
#   names(res) <- c("biome","var","mean","lower","upper")
#   res$var <- factor(res$var, labels=c("ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion"))
#   res$biome <- cut(res$biome, breaks=n_biome, labels=c("COAST","LC","HCSS","HCPS"))
#   
#   res <- join(res, av_biomesum, by="biome")
#   res$mean_tot <- res$mean * res$vol
#   res$lower_tot <- res$lower * res$vol
#   res$upper_tot <- res$upper * res$vol
#   
#   if(saveFiles){
#     write.csv(res, str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_MeanResults.csv"), row.names=FALSE)
#     }
#   
#   return(res)
#   
# }


return_bci <- function(data, na.rm=TRUE){
  require(simpleboot)
  if(all(data == 0)){
    return(c(0,0))
  }
  # get rid of NAs
  if (na.rm==TRUE){
    data = data[! is.na(data)]
  }
  
  set.seed(100)
  
  boot <- one.boot(data, mean, R=2000, trim=0)
  
  if(boot$t0==0){return(c(0,0))}
  
  bci <- boot::boot.ci(boot, conf=0.95, type=c("bca"))
  
  if(is.null(bci)){
    return(rep(median(boot$data), 2))
  } 
  
  limits <- bci$bca[4:5]
  return(limits)
}


return_quant <- function(data, probs=c(0.1,0.9)){
  res <- quantile(data, probs=probs)
  return(res)
}


