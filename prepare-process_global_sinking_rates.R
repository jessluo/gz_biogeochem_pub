#!/usr/bin/env Rscript
#
#     Process global sinking rates. Uses World Ocean Atlas temperature by depth
#       - Uses export decomposition model from Lebrato et al. 2011
#       - Uses sinking speeds from Lebrato et al. 2013. (conservative estimates)
#
#     
#
#----------------------------------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
library("oce")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- TRUE

`%ni%` <- Negate(`%in%`)

source("functions.R")

# }

## { Initial processing of World Ocean Atlas temperature at depth data -------------------
d <- read.csv("data/woa_annual_temp_by_depth.csv")

# fill in missing values
rd <- read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv")
nodata <- join(rd[,c("lat", "lon")], d, by=c("lat", "lon"))
nodata <- nodata[which(is.na(nodata$X0)),]
nodata <- nodata[!duplicated(nodata),] # get rid of duplicates

# save processing time by getting rid of the points that fall on land
d <- d[which(!is.na(d$X0)),]

# because all the missing values are in coastal areas that are shallow, pick the MINIMUM of the nearby values
for (i in 1:nrow(nodata)){
  tempdepth <- d[which(d$lat <= (nodata$lat[i] + 1.5) & d$lat >= (nodata$lat[i] - 1.5) & 
                                    d$lon <= (nodata$lon[i] + 1.5) & d$lon >= (nodata$lon[i] - 1.5)),]
  if(nrow(tempdepth)==0){
    tempdepth <- d[which(geodDist(d$lon, d$lat, nodata$lon[i], nodata$lat[i]) == 
                                      min(geodDist(d$lon, d$lat, nodata$lon[i], nodata$lat[i]))),]
  }
  # by doing straight colMeans without na.rm=TRUE we are pulling out the shallowest site
  # temperature variation isn't so much that it will make a huge difference. this is easiest by far
  if(is.na(mean(tempdepth$X5)) & length(which(is.na(tempdepth$X5))) / nrow(tempdepth) < 0.5){
    tempdepth <- tempdepth[!is.na(tempdepth$X5),]
  }
  nodata[i, 3:ncol(nodata)] <- colMeans(tempdepth[,3:ncol(tempdepth)]) 
}

# put together
d <- rbind(d, nodata)

# melt depths
dm <- melt(d, value.name = "temp", id.vars = c("lat", "lon"), variable.name = "depth")
dm$depth <- as.character(dm$depth)
dm$depth <- as.numeric(str_sub(dm$depth, 2, -1))

# remove all the NA temps, for faster processing through the rest of the code
dm <- dm[!is.na(dm$temp),]


#calculate the change in temperature by depth
dm <- ddply(dm, ~lat+lon, function(x){
  difftemp=diff(x$temp)
  diffdepth=diff(x$depth)
  difftemp[which(difftemp==0)] <- 1e-10
  kt <- difftemp/diffdepth
  kt <- c(NA, kt)
  x <- cbind(x, kt)  
  return(x)
}, .progress="text")

# }

## { Global gelatinous zooplankton biomass export ----------------------------
# sinking velocities array (m d-1)
C <- c(100, 800, 900, 1000, 1100, 1200)

# create empty columns for the C values
dm <- cbind(dm, data.frame(matrix(NA, nrow=nrow(dm), ncol=length(C))))
# rename names
names(dm)[(ncol(dm)-(length(C)-1)):ncol(dm)] <- C

# keep melting
dm <- melt(dm, id.vars=c("lat", "lon", "depth", "temp", "kt"), variable.name="C")
dm$C <- as.numeric(as.character(dm$C))
dm <- dm[,which(names(dm) != "value")]

# calculate the amount of time it takes to sink through each depth bin
# time[i] is calculated as the amount of time it takes to sink from depth[i-1] to depth[i]
dm <- ddply(dm, ~lat+lon, function(x){
  x$time <- NA
  uniqueC <- unique(x$C)
  for (i in 1:length(uniqueC)){
    time <- diff(x[x$C == uniqueC[i],"depth"]) / uniqueC[i]
    time <- c(NA, time)
    x$time[x$C == uniqueC[i]] <- time
  }
  return(x)
}, .progress="text")

dm <- arrange(dm, lat, lon, C, depth)

# prep for the calculation - do it on an array instead of a loop
temp <- c(NA, dm$temp[1:(nrow(dm)-1)])
dm$temp_1 <- temp
# dm$temp_1[which(is.na(dm$kt))] <- NA

# calculation of the m_ratio for each depth bin
# M_zi : M_zi.1 = exp(((-0.140*exp(0.145*T_zi.1)) / (0.145*K_Ti*C))(exp^0.145*K_Ti*deltaT*C - 1))
dm$m_ratio_db <- exp(((-0.140*exp(0.145*dm$temp_1)) / (0.145*(-dm$kt*dm$C))) * (exp(0.145*-dm$kt*dm$time*dm$C) - 1))


# split up the dataframe by latitude and longitude and sinking speed
# then calculate the total amount of biomass exported to each depth bin
dm <- ddply(dm, ~lat+lon+C, function(x){
  if(nrow(x)==1) {
    x$m_ratio=1
  } else {
    x$m_ratio = c(1,cumprod(x$m_ratio_db[2:nrow(x)]))
  }
  return(x)
}, .progress="text")


# clean up
dm <- dm[,which(names(dm) %ni% c("kt", "temp_1", "m_ratio_db"))] 
dm <- dm[complete.cases(dm) | dm$depth==0,]

# round to 4 significant figures
dm$temp <- signif(dm$temp, 4)
dm$time <- signif(dm$time, 4)
dm$m_ratio <- signif(dm$m_ratio, 4)

dm$latlon <- paste0(dm$lat, ", ", dm$lon)
save(dm, file="data/global_particle_sinking.Rdata")

# }
