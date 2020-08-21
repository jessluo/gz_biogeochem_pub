#!/usr/bin/env Rscript
#
#           Global input data processing
#           1) Load and process ESRL OI SST values
#           2) Calculate the surface area and volume under a 1-deg grid
#
#           Jessica Luo
#------------------------------------------------------------------------------------------

# load libraries
library("plyr")
library("reshape2")
library("ncdf4")
library("matlab")
library("oce")
library("ggplot2")
library("stringr")

# load functions
source("functions.R")

##{ Process SST ----------------------------------------------------------------------------
# long-term monthly mean SST from 1971-2000, in 1 degree grids, from NOAA ESRL
# http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
nc_mtemp <- nc_open("raw_data/sst.ltm.1971-2000.nc")

sst <- ncvar_get(nc_mtemp, nc_mtemp$var[[1]])
str(sst) # num[1:360, 1:180, 1:12]

lon <- nc_mtemp$dim$lon$vals
lat <- nc_mtemp$dim$lat$vals
times <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

sstdf <- process.nc.biomass.data(sst, times=times)
# change longitude to -180 to 180
sstdf$lon <- ifelse(sstdf$lon > 180, sstdf$lon-360, sstdf$lon)


write.csv(sstdf, "data/sst.ltm.1971-2000.1-deg.csv", row.names=F)

rsst <- data.frame(sstdf[,1:2],ann=rowMeans(sstdf[,3:14]))
rsst <- rsst[complete.cases(rsst),]

save(rsst, "data/sst.ltm.1-deg.annual.Rdata")
# }

##{ Calculate Area -----------------------------------------------------------------------

# pull the lat and lon values from the SST grid
# could potentially pull from other grids but this one is the most readily available and easiest to use
lat_lon <- sstdf[complete.cases(sstdf),c("lat", "lon")]

# add in 'missing'/coastal values so data points don't get tossed
rd <- read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv", as.is=TRUE)
rd$latlon <- str_c(rd$lat, ',', rd$lon)
add <- which(rd$latlon %ni% str_c(lat_lon$lat, ',', lat_lon$lon))

lat_lon <- rbind(lat_lon, rd[add, c("lat", "lon")])

# grab seafloor data
seafloor <- read.csv("data/woa_seafloor_depth.csv")
seafloor <- join(lat_lon, seafloor, by=c("lat", "lon"), type="left")

# changes the seafloor.depth values that are either 0 or NA to be 2
seafloor[which(is.na(seafloor$depth)),"depth"] <- 5.5
seafloor[seafloor$depth==0,"depth"] <- 5.5

seafloor$intDepth <- ifelse(seafloor$depth>=200, 198, seafloor$depth-3)

area.km2 <- area_earth_1deg(lat_lon$lat, lat_lon$lon) # in squared km
area.m2 <- area.km2*1E6 # squared m
vol <- area.m2 * seafloor$intDepth # cubic meters

av <- cbind(lat_lon, area.km2, area.m2, vol)

# scale back area and vol in areas where there is a coastline - assuming this file is ocean area/vol only.
coastline.world <- read.csv("data/gshhg_world_c.csv")
coastline.world$lon_x1 <- floor(coastline.world$lon) + 0.5
coastline.world$lat_x1 <- floor(coastline.world$lat) + 0.5

coastline.x1 <- ddply(coastline.world, ~lat_x1+lon_x1, function(x){return(nrow(x))})
coastline.x1 <- coastline.x1[complete.cases(coastline.x1),]
coastline.x1$multiplier=0.54
names(coastline.x1)[1:2]<-c("lat","lon")

add <- which(rd$latlon[add] %ni% str_c(coastline.x1$lat_x1,',',coastline.x1$lon_x1))
coastline.x1 <- rbind(coastline.x1, data.frame(c(rd[add, c("lat","lon")], V1=1, multiplier=0.17)))

av <- join(av,coastline.x1[,c("lat","lon","multiplier")])
av[is.na(av$multiplier),"multiplier"] <- 1


av$vol <- av$vol * av$multiplier
av$area.km2 <- av$area.km2 * av$multiplier
av$area.m2 <- av$area.m2 * av$multiplier

sum(av$area.km2)
# oceans surface should be 361,900,000 km2

write.csv(av, "data/surfaceArea_volume_1-deg.csv", row.names=F)

# }

# }
