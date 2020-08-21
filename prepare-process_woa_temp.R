#./prepare-process_sst.R
#
#     Imports and processing World Ocean Atlas hydrological data
#       - Temperature and salinity by depth
#       - Calculates sequestration depth and seafloor depth
#
#     Jessica Luo 2015-2016
#
#------------------------------------------------------------------------------------------------

library(ncdf4)
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)

NEW_DOWNLOAD=TRUE
PROCESS_RAW=TRUE

##{ Functions -----------------------------------------------------------------------------------
`%ni%` <- Negate(`%in%`)

process_nc_woa_phys_data <- function(d, depths=depths){
  returndf <- NULL
  for (i in 1:length(depths)){
    temp <- d[,,i]
    temp <- as.data.frame(cbind(lon, temp))
    names(temp) <- c("lon", lat) # rows are longitudes, columns are latitudes
    temp <- melt(temp, id.vars = "lon", variable.name="lat", value.name=as.character(depths[i]))
    latlon <- temp[,1:2]
    returndf <- as.data.frame(cbind(returndf, temp[,3]))
  }
  names(returndf) <- as.character(depths)
  returndf <- cbind(latlon, returndf)
  returndf$lat <- as.numeric(as.character(returndf$lat))
  # correct the longitude to be -180 to 180
  # returndf$lon <- ifelse(returndf$lon < -180, returndf$lon+360, returndf$lon)
  return(returndf)
}

# }

## { Temperature and salinity file download and process ---------------------------------------------
	
# # Temperature file download and save
# # download the file directly
if(NEW_DOWNLOAD){
	filename = "https://data.nodc.noaa.gov/thredds/fileServer/nodc/archive/data/0114815/public/temperature/netcdf/decav/1.00/woa13_decav_t00_01.nc"
	loc_temp = "raw_data/woa13_decav_t00_01.nc"
	download.file(filename, loc_temp)
	
	## Salinity -- commented out here because Salinity is not needed for our purposes here.
	# filename = "https://data.nodc.noaa.gov/thredds/fileServer/nodc/archive/data/0114815/public/salinity/netcdf/decav/1.00/woa13_decav_s00_01.nc"
	# loc_sal = "raw_data/woa13_decav_s00_01.nc"
	# download.file(filename, loc_sal)
}

# Process files
if(PROCESS_RAW){
	dir.create("data")
	
	# TEMPERATURE
	nc <- nc_open(loc_temp)
	# print(nc) # view the variables
	lat <- ncvar_get(nc, 'lat')
	lon <- ncvar_get(nc, 'lon')
	depths <- ncvar_get(nc, 'depth')

	data <- ncvar_get(nc, "t_an") # units are lon, lat, depth
	data <- process_nc_woa_phys_data(data, depths)

	write.csv(data, "data/woa_annual_temp_by_depth.csv", row.names=F)
	
	
	# # SALINITY
	# nc <- nc_open(loc_sal)
	# # print(nc) # view the variables
	# lat <- ncvar_get(nc, 'lat')
	# lon <- ncvar_get(nc, 'lon')
	# depths <- ncvar_get(nc, 'depth')
	#
	# data <- ncvar_get(nc, "s_an") # units are lon, lat, depth
	# data <- process_nc_woa_phys_data(data, depths)
	#
	# write.csv(data, "data/woa_annual_sal_by_depth.csv", row.names=F)
}

# }

# { Initial processing of World Ocean Atlas temperature/salinity at depth data -------------------
t <- read.csv("data/woa_annual_temp_by_depth.csv")
# s <- read.csv("data/woa_annual_sal_by_depth.csv")

# save processing time by getting rid of the points that fall on land
t <- t[which(!is.na(t$X0)),]
# s <- s[which(!is.na(s$X0)),]

# fill in missing temperature values
rd <- read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv")

notemp <- join(rd[,c("lat", "lon")], t, by=c("lat", "lon"))
notemp <- notemp[which(is.na(notemp$X0)),]
notemp <- notemp[!duplicated(notemp),]

# nosal <- join(rd[,c("lat", "lon")], s, by=c("lat", "lon"))
# nosal <- nosal[which(is.na(nosal$X0)),]
# nosal <- nosal[!duplicated(nosal),]

# because all the missing values are in coastal areas that are shallow, pick the MINIMUM of the nearby values
for (i in 1:nrow(notemp)){
  tempdepth <- t[which(t$lat <= (notemp$lat[i] + 1.5) & t$lat >= (notemp$lat[i] - 1.5) & 
                         t$lon <= (notemp$lon[i] + 1.5) & t$lon >= (notemp$lon[i] - 1.5)),]
  # saldepth <- s[which(s$lat <= (notemp$lat[i] + 1.5) & s$lat >= (notemp$lat[i] - 1.5) &
                         # s$lon <= (notemp$lon[i] + 1.5) & s$lon >= (notemp$lon[i] - 1.5)),]
  if(nrow(tempdepth)==0){
    tempdepth <- t[which(geodDist(t$lon, t$lat, notemp$lon[i], notemp$lat[i]) == 
                           min(geodDist(t$lon, t$lat, notemp$lon[i], notemp$lat[i]))),]
    # saldepth <- s[which(geodDist(s$lon, s$lat, notemp$lon[i], notemp$lat[i]) ==
                           # min(geodDist(s$lon, s$lat, notemp$lon[i], notemp$lat[i]))),]
  }
  # by doing straight colMeans without na.rm=TRUE we are pulling out the shallowest site
  # temperature variation isn't so much that it will make a huge difference. this is easiest by far
  if(is.na(mean(tempdepth$X5)) & length(which(is.na(tempdepth$X5))) / nrow(tempdepth) < 0.5){
    tempdepth <- tempdepth[!is.na(tempdepth$X5),]
    # saldepth <- saldepth[!is.na(saldepth$X5),]
  }
  notemp[i, 3:ncol(notemp)] <- colMeans(tempdepth[,3:ncol(tempdepth)]) 
  # nosal[i, 3:ncol(nosal)] <- colMeans(saldepth[,3:ncol(saldepth)])
}

# put together
t <- rbind(t, notemp)
# s <- rbind(s, nosal)

# melt depths
tm <- melt(t, value.name = "temp", id.vars = c("lon", "lat"), variable.name = "depth")
tm$depth <- as.character(tm$depth)
tm$depth <- as.numeric(str_sub(tm$depth, 2, -1))

# sm <- melt(s, value.name = "sal", id.vars = c("lon", "lat"), variable.name = "depth")
# sm$depth <- as.character(sm$depth)
# sm$depth <- as.numeric(str_sub(sm$depth, 2, -1))

# remove all the NA values, for faster processing through the rest of the code
tm <- tm[!is.na(tm$temp),]
# sm <- sm[!is.na(sm$sal),]


# export depth - 100 m
export.depth <- ddply(tm, ~lon+lat, function(x){
  if(x[nrow(x), "depth"] < 100){return(x[nrow(x), "depth"])}
  else {return(100)}
}, .progress="text")

export.depth <- rename(export.depth, replace=c("V1" = "depth"))
write.csv(export.depth, "data/woa_export_depth.csv", row.names=FALSE)

# pull out the sequestration depth
# use 1000 m (Carlson and Passow 2012)
seq.depth <- ddply(tm, ~lon+lat, function(x){
  if(x[nrow(x), "depth"] < 1000){return(x[nrow(x), "depth"])}
  else {return(1000)}
}, .progress="text")

seq.depth <- rename(seq.depth, replace=c("V1" = "depth"))
write.csv(seq.depth, "data/woa_sequestration_depth.csv", row.names=FALSE)

# find the bottom depth
seafloor.depth <- ddply(tm, ~lon+lat, function(x){
  return(depth=x[nrow(x),"depth"])
}, .progress="text")

seafloor.depth <- rename(seafloor.depth, replace=c("V1" = "depth"))
write.csv(seafloor.depth, "data/woa_seafloor_depth.csv", row.names=FALSE)



# }